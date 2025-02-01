/**
 * @file WPGen.cpp
 * @brief Implementation of the WPGen class for generating gradient waveforms and k-space trajectories.
 * 
 * This file contains the implementation of the WPGen class, which is responsible for generating
 * gradient waveforms and k-space trajectories for MRI sequences. The class supports both single Tau
 * and dual Tau configurations and can handle different spiral directions (in, out, in-out).
 * 
 * @date 2024-06-01
 * @version 1.0
 * @author Guruprasad Krishnamoorty
 */

#include "WPGen.hpp"
#include <iostream>
namespace wp{
// Constructor that takes Tau and computes number of arms inturn
WPGen::WPGen(double fov_, double res_, double reqTau_, double sMax_,
  double gMax_, double omegaMax_, double gradDwell_,
  double adcDwell_, double gamma_, bool addOuterRing_, double _alphaVD) : mFov(fov_),
      mRes(res_),
      mReqTau(reqTau_),
      slewMaxIn(sMax_),
      gradMaxIn(gMax_),
      omegaMax(omegaMax_),
      mGradRast(gradDwell_),
      mDwell(adcDwell_),
      mGamma(gamma_),
      mAddOuterRing(addOuterRing_),
      mAlphaVD(_alphaVD)
{
  double kMax = (2.0 / sqrt(M_PI)) * (0.5 / res_); /* True Resolution! */

  /* First determine number of arms */
  double freqTerm = (gamma_ * sMax_) / (3.0 * pow(omegaMax_, 3.0));
  double slewTerm = (gamma_ * pow(gMax_, 3.0)) / (6.0 * pow(sMax_, 2.0));
  double gMaxTerm = (pow(kMax, 2.0)) / (2.0 * gamma_ * gMax_);
  double wpTerm = freqTerm + slewTerm + gMaxTerm;
  nyqNumberArms = 2.0 * M_PI * fov_ * wpTerm / mReqTau;
  // Determine the actual readout duration first
  _initialize();

  // To match the required Tau with actual Tau
  double constTerm = readoutDur * nyqNumberArms;
  double currReqArms = constTerm / mReqTau;
  double currReqTau = (nyqNumberArms / currReqArms) * mReqTau;
  nyqNumberArms = 2.0 * M_PI * fov_ * wpTerm / currReqTau;

  // Initilaize again with the updated nyqNumberArms
  _initialize();
}

// Constructor that takes number of arms
WPGen::WPGen(double fov_, double res_, int nArms_, double sMax_,
  double gMax_, double omegaMax_, double gradDwell_,
  double adcDwell_, double gamma_, bool addOuterRing_, double _alphaVD) : mFov(fov_),
      mRes(res_),
      nyqNumberArms(nArms_),
      slewMaxIn(sMax_),
      gradMaxIn(gMax_),
      omegaMax(omegaMax_),
      mGradRast(gradDwell_),
      mDwell(adcDwell_),
      mGamma(gamma_),
      mAddOuterRing(addOuterRing_),
      mAlphaVD(_alphaVD)
{
  _initialize();
}

void WPGen::_initialize()
{
  _initializeParameters();
  // Find compatible constraints, so each segment does not have "negative" duration
  _findCompatibleConstraints();
  // Find timings
  _computeTimings();
}

void WPGen::_initializeParameters()
{
  mDelta = nyqNumberArms / (2.0 * M_PI * mFov);

  /* We define a spiral circle to have a diameter (2/sqrt(Pi))*(1/res) in k-space to match a square
      sampling region of width (1/res), to keep resolution comparable - we call this True Resolution

      Van Gelderen P. Comparing true resolution in square versus circular k-space sampling.
      Proceedings of the 6th Annual Meeting of ISMRM, Sydney, Australia, 1998. p 424.

      We map this to a k-space grid that is 25% larger than (FOV/RES), corresponding to 80% requested resolution
      This all works out nicely :-) but look down 4 lines if you don't want to do this
  */
  mTrueRes = sqrt(M_PI) / 2.0;
  mGridRes = 0.8;

  /* If you want to cheat, you can ignore this by uncommenting the next 2 lines
      trures = 1.0;
      gridres = 1.0;
  */

  int origGridMtx = static_cast<int>(floor(mFov / (mGridRes * mRes)));
  gridMtxSize = _favFFTLength(origGridMtx);

  coordsScaleFactor = static_cast<double>(origGridMtx) / 
                static_cast<double>(gridMtxSize);

  mRampDownComputed = false;
  mBaseSpiralComputed = false;
  mKSPComputed = false;
  mM0M1Computed = false;

  m0m1Pts = 0; mx0 = 0; mx1 = 0; my0 = 0; my1 = 0;

  if(mAddOuterRing){
    // Yoyo covers a little extra k-space with the ring, we hit Rc+delta
    // We translate that here into a new krad_max for WHIRLED PEAS to hit

    double klim = 0.5/(mTrueRes*mRes);
    rcLim = klim-mDelta;
    kradMax = sqrt(rcLim * rcLim + mDelta * mDelta);
  }
  else{
    kradMax = 0.5 / (mTrueRes * mRes);
  }
}


void WPGen::_findCompatibleConstraints()
{
  // Actual gmax is limited by the min dwell time supported by adc
  double gMax_dwell_limited = 1.0 / (mGamma * mDwell * mFov);

  if(gMax_dwell_limited > gradMaxIn)
    gMax_dwell_limited = gradMaxIn;

  double omega1 = sqrt(2.0 * mGamma * slewMaxIn / (3.0 * mDelta));
  double omega2 = 2.0 * mGamma * gMax_dwell_limited / (3.0 * mDelta);
  omegaMax = (omegaMax > 0) ? std::min(std::min(omegaMax, omega1), omega2) : std::min(omega1, omega2);

  double slew1 = gMax_dwell_limited * omegaMax;
  double slew2 = sqrt(pow(omegaMax, 4.0) * (kradMax * kradMax - mDelta * mDelta) / (mGamma * mGamma));
  slewMax = std::min(std::min(slewMaxIn, slew1), slew2);

  double grad1 = pow((slewMax * slewMax) * (kradMax * kradMax - mDelta * mDelta) / (mGamma * mGamma), 0.25);
  gradMax = std::min(gMax_dwell_limited, grad1);
  gradMaxSpiral = gradMax;
}

void WPGen::_computeTimings()
{
  // segments start at tx0, end at tx1, with total time t_X

  // Arc
  arcDur0 = 0.0;
  arcDur1 = (5.0 * M_PI + 1.0) / (6.0 * omegaMax);
  arcDur = arcDur1 - arcDur0;

  // Omega Constrained
  omegaDur0 = 1.0 / omegaMax;
  omegaDur1 = (mGamma * slewMax) / (mDelta * pow(omegaMax, 3.0));
  omegaDur = omegaDur1 - omegaDur0;

  // Slew Constrained
  slewDur0 = (2.0 * mGamma * slewMax) / (3.0 * mDelta * pow(omegaMax, 3.0));
  slewDur1 = (2.0 * mGamma * pow(gradMax, 3.0)) / (3.0 * mDelta * slewMax * slewMax);
  slewDur = slewDur1 - slewDur0;

  // Gradient Constrained
  gradDur0 = (mGamma * pow(gradMax, 3.0)) / (2.0 * mDelta * slewMax * slewMax);
  gradDur1 = (kradMax * kradMax - mDelta * mDelta) / (2.0 * mGamma * mDelta * gradMax);
  gradDur = gradDur1 - gradDur0;

  // Outer YoYo Ring (90 degrees)
  outRingDur = 0.0;
  if (mAddOuterRing) {
    thetaPrimeRing = mGamma * gradMax / sqrt(kradMax * kradMax - mDelta * mDelta);
    outRingDur = 0.5 * M_PI / thetaPrimeRing;
    outRingDur0 = gradDur1;
    outRingDur1 = gradDur1 + outRingDur;
  }
  else {
    thetaPrimeRing = 0.0;
    outRingDur = 0.0;
    outRingDur0 = 0.0;
    outRingDur1 = 0.0;
  }

  // gradient rampdown is a Hanning Window
  // This may help a little with spiral-in (??)
  rampDur = _roundGrad(0.5 * M_PI * gradMax / slewMaxIn);

  readoutDur = _roundGrad(arcDur + omegaDur + slewDur + gradDur + outRingDur);

  totalDur = readoutDur + rampDur;

  spiralReadPts = static_cast<int>(readoutDur / mGradRast);
  spiralRampPts = static_cast<int>(rampDur / mGradRast);
  totalSpiralPts = spiralReadPts + spiralRampPts;

  kspacePts = static_cast<int>(readoutDur / mDwell);

  // Constant Gradient crd variables (fast arc followed by constant gradient)
  mDwellCG = 1.0 / (mGamma * gradMaxIn * mFov);
  mArcDurCG = nyqNumberArms / (2.0 * mGamma * gradMaxIn * mFov);
  mG0CG = mDelta / (2.0 * mGamma * gradMaxIn);
  mCtaCG = 2.0 * M_PI * (mGamma * mGamma * gradMaxIn * gradMaxIn * mFov * mFov) / (nyqNumberArms * nyqNumberArms);
  mTotalCGDur = mArcDurCG + (gradDur1 - mG0CG) + outRingDur;
}

void WPGen::ComputeBaseSpiral(dVector& gradX, dVector& gradY)
{
  // Define some constants
  const double arc_ta = M_PI / (3.0 * omegaMax);
  const double arc_tb = (1.0 + 2.0 * M_PI) / (6.0 * omegaMax);

  const double cga = mDelta * omegaMax / (3.0 * mGamma);
  const double cgw = (mDelta * omegaMax * omegaMax) / mGamma;
  const double cgs = pow((3.0 * mDelta * slewMax * slewMax) / (2.0 * mGamma), (1.0 / 3.0));
  const double cgg = gradMax;

  const double cta = omegaMax / 3.0;
  const double ctw = omegaMax;
  const double cts = pow((9.0 * mGamma * slewMax) / (4.0 * mDelta), 1.0 / 3.0);
  const double ctg = sqrt(2.0 * mGamma * gradMax / mDelta);

  double momT = 0.0;
  double gmag = 0.0, t, tt, theta = 0.0;
  double t0 = readoutDur - (spiralReadPts * mGradRast) + 0.5 * mGradRast;

  gradX.resize(spiralReadPts, 0.0);
  gradY.resize(spiralReadPts, 0.0);

  // Variable density parameters
  double tNorm, densityScale;

  // Compute GRADIENT
  for (int i = 0; i < spiralReadPts; i++) {
    t = t0 + static_cast<double>(i) * mGradRast;
    tNorm = t / readoutDur; // Normalized time for variable density scaling 
    densityScale = pow(tNorm, mAlphaVD); // Smooth density modulation using a cubic polynomial

    // ARC
    if (t < arcDur1) {
      if (t < arc_ta) {
        gmag = cga * (1 - cos(3 * omegaMax * t));
        theta = cta * (t - (sin(3 * omegaMax * t) / (3 * omegaMax))) + 1 - (0.5 * M_PI);
      } else if (t < (arc_ta + arc_tb)) {
        tt = t - arc_ta;
        gmag = 2 * cga;
        theta = cta * (arc_ta + 2 * tt) + 1 - (0.5 * M_PI);
      } else {
        tt = t - arc_ta - arc_tb;
        gmag = cga * (3 - cos(3 * omegaMax * tt));
        theta = cta * (arc_ta + (2 * arc_tb) + 3 * tt - (sin(3 * omegaMax * tt) / (3 * omegaMax))) + 1 - (0.5 * M_PI);
      }
    } else {
      t = t - arcDur1 + omegaDur0;

      // FRQ
      if (t < omegaDur1) {
        gmag = cgw * t;
        theta = ctw * t;
      } else {
        t = t - omegaDur1 + slewDur0;

        // SLEW
        if (t < slewDur1) {
          gmag = cgs * pow(t, (1.0 / 3.0));
          theta = cts * pow(t, (2.0 / 3.0));
        } else {
          t = t - slewDur1 + gradDur0;

          // GRAD
          if (mAddOuterRing) {
            if (t < gradDur1) {
              gmag = cgg;
              theta = ctg * sqrt(t);
            } else {
              gmag = cgg;
              theta = ctg * sqrt(gradDur1) + thetaPrimeRing * (t - gradDur1);
            }
          } else {
            gmag = cgg;
            theta = ctg * sqrt(t);
          }
        }
      }
    }

    gmag *= densityScale; // Scale gradient magnitude with density factor

    gradX[i] = gmag * cos(theta);
    gradY[i] = gmag * sin(theta);

    // Calculate moments
    momT = mGradRast * static_cast<double>(i - (spiralReadPts + spiralRampPts));
    mx0 += mGradRast * gradX[i];
    my0 += mGradRast * gradY[i];
    mx1 += mGradRast * gradX[i] * momT;
    my1 += mGradRast * gradY[i] * momT;
  }

  m0Angle = atan2(my0, mx0);
  mBaseSpiralComputed = true;
}

void WPGen::ComputeRampDown(dVector& rampX, dVector& rampY)
{
  rampX.resize(spiralRampPts, 0.0);
  rampY.resize(spiralRampPts, 0.0);

  double t, theta = 0.0, gmag, momT;

  const double ctg = sqrt(2.0 * mGamma * gradMax / mDelta);
  const double cgg = gradMax;

  // Compute GRADIENT RAMPDOWN
  for (int i = 0; i < spiralRampPts; i++) {
    t = static_cast<double>(i) * mGradRast;

    if(!mAddOuterRing)
      theta = ctg * sqrt(t + gradDur1);
    else
      theta = ctg * sqrt(gradDur1) + thetaPrimeRing * (outRingDur1 + t - gradDur1);

    gmag = cgg * cos(0.5 * M_PI * t / rampDur);
    rampX[i] = gmag * cos(theta);
    rampY[i] = gmag * sin(theta);

    // calculate moments
    momT = mGradRast * static_cast<double>(i - spiralRampPts);
    mx0 = mx0 + mGradRast * rampX[i];
    my0 = my0 + mGradRast * rampY[i];
    mx1 = mx1 + mGradRast * rampX[i] * momT;
    my1 = my1 + mGradRast * rampY[i] * momT;
  }

  mRampDownComputed = true;
}

bool WPGen::ComputeM0M1Trapezoids(dVector& m0m1X,  dVector& m0m1Y)
{
  if(!(mRampDownComputed && mBaseSpiralComputed))
  {
    printf("Error! ComputeBaseSpiral() and ComputeRampDown() must to be called before calling this routine\n");
    return mM0M1Computed;
  }

  /* Guess at refocussing trapezoids
   * trap A followed by trapb
   * ramps are 2*tr
   * plateaus are ta and tb, include half of ramp either side
   * this way areas are ga*ta and gb*tb, respectively
   *
   *         _________
   *        /         \
   *       /           \
   *      /             \
   * _____/               \___________
   *     | |<-- ta --->| |
   *   ->| |<-tr     ->| |<-tr
   */

  // Print mx0, mx1, my0, my1
  // Initialize variables
  bool searching = true;
  double scale_a = 1.0;
  double scale_b = 1.0;
  double tra, trb, ta, tb;
  double den, aa, bb, cc, dd;
  double gxa, gxb, gya, gyb, ga_norm, gb_norm;

  // The variable smooth describes what fraction of the ramp will be turned into parabola
  // smooth goes from 0 (all linear) to 1 (all parabolic)
  // total tr increases by a factor of (1+smooth)
  double smooth = 0.5;

  while (searching) {
    searching = false;

    tra = (1+smooth) * 0.5 * scale_a * gradMaxIn / slewMaxIn;
    trb = (1+smooth) * 0.5 * scale_b * gradMaxIn / slewMaxIn;
    ta = 2.*tra;
    tb = trb;

    for (int i = 0; i < 10; i++) {
      // Calculate new gradients
      den = ta * tb * (tra + trb + 0.5 * (ta + tb));
      aa = tb * (2 * tra + trb + ta + 0.5 * tb) / den;
      bb = -tb / den;
      cc = -ta * (tra + 0.5 * ta) / den;
      dd = ta / den;

      // Calculate new gradients
      gxa = -aa * mx0 - bb * mx1;
      gxb = -cc * mx0 - dd * mx1;
      gya = -aa * my0 - bb * my1;
      gyb = -cc * my0 - dd * my1;

      // Normalize new gradients
      ga_norm = scale_a * gradMaxIn / sqrt(gxa * gxa + gya * gya);
      gb_norm = scale_b * gradMaxIn / sqrt(gxb * gxb + gyb * gyb);

      gxa = gxa * ga_norm;
      gxb = gxb * gb_norm;
      gya = gya * ga_norm;
      gyb = gyb * gb_norm;

      // Calculate trapezoid durations
      den = (gxa * gyb - gxb * gya);
      ta = (-gyb * mx0 + gxb * my0) / den;
      tb = (gya * mx0 - gxa * my0) / den;

    }

    if (ta < 2.0 * tra) {
      searching = true;
      scale_a -= 0.1;
    }

    if (tb < 2.0 * trb) {
      searching = true;
      scale_b -= 0.1;
    }

    if (scale_a < 0.15) {
      searching = false;
      printf("ERROR!!!\n");
      return false;
    }

    if (scale_b < 0.15) {
      searching = false;
      printf("ERROR!!!\n");
      return false;
    }
  }

  double tra2 = 2 * tra;
  double trb2 = 2 * trb;
  double ta2 = ta - 2.0 * tra;
  double tb2 = tb - 2.0 * trb;

  // Round durations to grast while maintaining the area 
  gxa = (gxa * (tra2 + ta2)) / (_roundGrad(tra2) + _roundGrad(ta2));
  gya = (gya * (tra2 + ta2)) / (_roundGrad(tra2) + _roundGrad(ta2));

  gxb = (gxb * (trb2 + tb2)) / (_roundGrad(trb2) + _roundGrad(tb2));
  gyb = (gyb * (trb2 + tb2)) / (_roundGrad(trb2) + _roundGrad(tb2));

  tra2 = _roundGrad(tra2);
  ta2 = _roundGrad(ta2);
  trb2 = _roundGrad(trb2);
  tb2 = _roundGrad(tb2);

    // Calculate the number of points for m0m1 waveform
  m0m1Pts = static_cast<int>((2 * tra2 + ta2 + 2 * trb2 + tb2) / mGradRast) + 1;

  m0m1X.resize(m0m1Pts,0.0);
  m0m1Y.resize(m0m1Pts,0.0);

  double tlina = (1.0 - smooth) * tra2 / (1.0 + smooth);
  double tpara = smooth * tra2 / (1.0 + smooth);
  double tlinb = (1.0 - smooth) * trb2 / (1.0 + smooth);
  double tparb = smooth * trb2 / (1.0 + smooth);

  // Variables for smooth trapezoids
  double gxa_1 = gxa * smooth / 2.0;
  double gxa_2 = gxa * (1.0 - smooth);
  double gya_1 = gya * smooth / 2.0;
  double gya_2 = gya * (1.0 - smooth);
  double gxb_1 = gxb * smooth / 2.0;
  double gxb_2 = gxb * (1.0 - smooth);
  double gyb_1 = gyb * smooth / 2.0;
  double gyb_2 = gyb * (1.0 - smooth);

  double t_trapa = tra2 + ta2 + tra2;

  for (int i = 0; i < m0m1Pts; i++) {
    double t = static_cast<double>(i) * mGradRast;

    // TRAP A
    if (t < tpara) { // parabolic start
      double tt = t / tpara;
      m0m1X[i] = tt * tt * gxa_1;
      m0m1Y[i] = tt * tt * gya_1;
    } else if (t < tpara + tlina) { // linear ramp up
      double tt = (t - tpara) / tlina;
      m0m1X[i] = gxa_1 + tt * gxa_2;
      m0m1Y[i] = gya_1 + tt * gya_2;
    } else if (t < tra2) { // parabolic end
      double tt = (tra2 - t) / tpara;
      m0m1X[i] = gxa - tt * tt * gxa_1;
      m0m1Y[i] = gya - tt * tt * gya_1;
    } else if (t < tra2 + ta2) { // plateau
      m0m1X[i] = gxa;
      m0m1Y[i] = gya;
    } else if (t < tra2 + ta2 + tpara) { // parabolic start
      double tt = (t - tra2 - ta2) / tpara;
      m0m1X[i] = gxa - tt * tt * gxa_1;
      m0m1Y[i] = gya - tt * tt * gya_1;
    } else if (t < tra2 + ta2 + tpara + tlina) { // linear ramp down
      double tt = (tra2 + ta2 + tpara + tlina - t) / tlina;
      m0m1X[i] = gxa_1 + tt * gxa_2;
      m0m1Y[i] = gya_1 + tt * gya_2;
    } else if (t < 2 * tra2 + ta2) { // parabolic end
      double tt = (2 * tra2 + ta2 - t) / tpara;
      m0m1X[i] = tt * tt * gxa_1;
      m0m1Y[i] = tt * tt * gya_1;
    }

    // TRAP B
    else if (t < t_trapa + tparb) { // parabolic start
      double tt = (t - t_trapa) / tparb;
      m0m1X[i] = tt * tt * gxb_1;
      m0m1Y[i] = tt * tt * gyb_1;
    } else if (t < t_trapa + tparb + tlinb) { // linear ramp up
      double tt = (t - (t_trapa + tparb)) / tlinb;
      m0m1X[i] = gxb_1 + tt * gxb_2;
      m0m1Y[i] = gyb_1 + tt * gyb_2;
    } else if (t < t_trapa + trb2) { // parabolic end
      double tt = (t_trapa + trb2 - t) / tparb;
      m0m1X[i] = gxb - tt * tt * gxb_1;
      m0m1Y[i] = gyb - tt * tt * gyb_1;
    } else if (t < t_trapa + trb2 + tb2) { // plateau
      m0m1X[i] = gxb;
      m0m1Y[i] = gyb;
    } else if (t < t_trapa + trb2 + tb2 + tparb) { // parabolic start
      double tt = (t - t_trapa - trb2 - tb2) / tparb;
      m0m1X[i] = gxb - tt * tt * gxb_1;
      m0m1Y[i] = gyb - tt * tt * gyb_1;
    } else if (t < t_trapa + trb2 + tb2 + tparb + tlinb) { // linear ramp down
      double tt = (t_trapa + trb2 + tb2 + tparb + tlinb - t) / tlinb;
      m0m1X[i] = gxb_1 + tt * gxb_2;
      m0m1Y[i] = gyb_1 + tt * gyb_2;
    } else if (t < t_trapa + 2 * trb2 + tb2) { // parabolic end
      double tt = (t_trapa + 2 * trb2 + tb2 - t) / tparb;
      m0m1X[i] = tt * tt * gxb_1;
      m0m1Y[i] = tt * tt * gyb_1;
    }
  }

  double gMaxTrap = 0.0; 
  for (int i = 0; i < m0m1Pts; i++) {
    double ssq = sqrt(m0m1X[i] * m0m1X[i] + m0m1Y[i]*m0m1Y[i]);
    if(ssq > gMaxTrap)
      gMaxTrap = ssq;
  }
  gradMax = std::max(gMaxTrap, gradMax) + 0.01; // 0.1 to avoid rounding errors
  
  m0m1Dur = m0m1Pts * mGradRast;
  totalDur += m0m1Dur;

  mM0M1Computed = true;
  return mM0M1Computed;
}

void WPGen::ComputeKSPnSDC(dVector& kspX,
dVector& kspY,
dVector& sdc)
{
  double t, tt, theta, alpha, krad = 0.0;
  double phi = 0.0, arm2arm, sdcPt = 0.0;

  // Define some constants
  const double arcTa = M_PI / (3.0 * omegaMax);
  const double arcTb = (1.0 + 2.0 * M_PI) / (6.0 * omegaMax);

  const double cga = mDelta * omegaMax / (3.0 * mGamma);
  const double cgw = (mDelta * omegaMax * omegaMax) / mGamma;
  const double cgs = pow((3.0 * mDelta * slewMax * slewMax) / (2.0 * mGamma), (1.0 / 3.0));
  const double cgg = gradMaxSpiral;

  const double cta = omegaMax / 3.0;
  const double ctw = omegaMax;
  const double cts = pow((9.0 * mGamma * slewMax) / (4.0 * mDelta), 1.0 / 3.0);
  const double ctg = sqrt(2.0 * mGamma * gradMaxSpiral / mDelta);

  const double csa = cga / gradMaxSpiral;
  const double csw = cgw / gradMaxSpiral;
  const double css = cgs / gradMaxSpiral;
  const double csg = cgg / gradMaxSpiral;

  double meanSDC = 0.0;
  double meanSDC2 = 0.0;

  kspX.resize(kspacePts, 0.0);
  kspY.resize(kspacePts, 0.0);
  sdc.resize(kspacePts, 0.0);

  kspTrace.resize(kspacePts, 0.0);

  double t0 = readoutDur - (kspacePts * mDwell) + 0.5 * mDwell;

  double phi0 = ctg * sqrt(gradDur1);

  // Power-law Variable Density Parameters
  double tNorm, densityScale;

  for (int i = 0; i < kspacePts; i++) {
    t = t0 + static_cast<double>(i) * mDwell;

    // Normalized time for power-law scaling
    tNorm = t / readoutDur;  // normalize time relative to the total duration
    densityScale = pow(tNorm, mAlphaVD);

    // ARC
    if (t < arcDur1) {
      if (t < arcTa) {
        theta = cta * (t - (sin(3.0 * omegaMax * t) / (3.0 * omegaMax)));
        krad = mDelta * sqrt(2.0 * (1.0 - cos(theta)));
        phi = atan2(1.0 - cos(theta), sin(theta)) + 1.0 - 0.5 * M_PI;
        arm2arm = sin(theta);
        sdcPt = csa * (1.0 - cos(3.0 * omegaMax * t)) * arm2arm;
      } else if (t < (arcTa + arcTb)) {
        tt = t - arcTa;
        theta = cta * (arcTa + 2.0 * tt);
        krad = mDelta * sqrt(2.0 * (1.0 - cos(theta)));
        phi = atan2(1.0 - cos(theta), sin(theta)) + 1.0 - 0.5 * M_PI;
        arm2arm = sin(theta);
        sdcPt = 2.0 * csa * arm2arm;
      } else {
        tt = t - arcTa - arcTb;
        theta = cta * (arcTa + (2.0 * arcTb) + 3.0 * tt - (sin(3.0 * omegaMax * tt) / (3.0 * omegaMax)));
        krad = mDelta * sqrt(2.0 * (1.0 - cos(theta)));
        phi = atan2(1.0 - cos(theta), sin(theta)) + 1.0 - 0.5 * M_PI;
        arm2arm = sin(theta);
        sdcPt = csa * (3.0 - cos(3.0 * omegaMax * tt)) * arm2arm;
      }
    } else {
      t = t - arcDur1 + omegaDur0;

      // FRQ
      if (t < omegaDur1) {
        sdcPt = csw * t;
        theta = ctw * t;
        krad = mDelta * sqrt(theta * theta + 1.0);
        phi = theta - acos(mDelta / krad);
      } else {
        t = t - omegaDur1 + slewDur0;

        // SLEW
        if (t < slewDur1) {
          sdcPt = css * pow(t, 1.0 / 3.0);
          theta = cts * pow(t, 2.0 / 3.0);
          krad = mDelta * sqrt(theta * theta + 1.0);
          phi = theta - acos(mDelta / krad);
        } else {
          t = t - slewDur1 + gradDur0;

          // GRAD
          if (mAddOuterRing) {
            if (t < gradDur1) {
              sdcPt = csg;
              theta = ctg * sqrt(t);
              krad = mDelta * sqrt(theta * theta + 1.0);
              phi = theta - acos(mDelta / krad);
            } else {
              // RING
              alpha = 0.5 * M_PI * (outRingDur1 - t) / outRingDur;
              sdcPt = csg * sin(alpha); // Sin() is ad-hoc for inter-arm distance
              theta = thetaPrimeRing * (t - gradDur1);
              krad = sqrt(rcLim * rcLim + mDelta * mDelta + 2.0 * rcLim * mDelta * sin(theta));
              phi = phi0 - acos((rcLim * sin(theta) + mDelta) / krad);
            }
          } else {
            sdcPt = csg;
            theta = ctg * sqrt(t);
            krad = mDelta * sqrt(theta * theta + 1.0);
            phi = theta - acos(mDelta / krad);
          }
        }
      }
    }

    // Apply power-law density scaling to the k-space and SDC
    krad *= densityScale;  // Adjust k-space radial distance with power-law scaling

    // Normalize K-space with mGridres * mRes
    // output k-space between -0.5 and +0.5 (our standard, but you scale as desired)
    kspX[i] = krad * cos(phi) * mGridRes * mRes * coordsScaleFactor;
    kspY[i] = krad * sin(phi) * mGridRes * mRes * coordsScaleFactor;
    sdc[i] = sdcPt;

    meanSDC += sdcPt;
    meanSDC2 += (sdcPt * sdcPt);
    kspTrace[i] = sqrt(kspX[i] * kspX[i] + kspY[i] * kspY[i]);
  }

  meanSDC /= kspacePts;
  meanSDC2 /= kspacePts;

  // Calculate SNR factor numerically
  // This is the SNR loss by not weighting all data uniformly
  snrFactor = meanSDC / sqrt(meanSDC2);
  snrFactor = meanSDC / sqrt(meanSDC2);

  mKSPComputed = true;
}

void WPGen::ComputeKSPnSDCForCG(dVector& kspX, dVector& kspY, dVector& sdc)
{
  const double cgg = gradMaxSpiral;
  const double csg = cgg / gradMaxSpiral;
  int kpts_cg = static_cast<int>(mTotalCGDur / mDwellCG);

  kspX.resize(kpts_cg, 0.0);
  kspY.resize(kpts_cg, 0.0);
  sdc.resize(kpts_cg, 0.0);

  double ctg = sqrt(2.0 * mGamma * gradMax / mDelta);
  double theta, phi = 0., arm2arm, krad = 0., phi0, alpha;

  for (int i = 0; i < kpts_cg; i++) {
    double t = static_cast<double>(i) * mDwellCG;

    // ARC
    if (t < mArcDurCG) {
      theta = mCtaCG * t * t;
      krad = mDelta * sqrt(2.0 * (1.0 - cos(theta)));
      phi = atan2(1.0 - cos(theta), sin(theta)) + 1.0 - 0.5 * M_PI;
      arm2arm = sin(theta);
      sdc[i] = (t / mArcDurCG) * arm2arm;
    } else {
      t = t - mArcDurCG + mG0CG;

      // GRAD
      if (mAddOuterRing) {
        if (t < mTotalCGDur) {
          sdc[i] = 1.0;
          theta = ctg * sqrt(t);
          krad = mDelta * sqrt(theta * theta + 1.0);
          phi = theta - acos(mDelta / krad);
        } else {
          // RING
          phi0 = ctg * sqrt(outRingDur1);
          if (t < outRingDur1) {
            alpha = 0.5 * M_PI * (outRingDur1 - t) / outRingDur;
            sdc[i] = csg * sin(alpha); // Sin() is ad-hoc for inter-arm distance
            theta = thetaPrimeRing * (t - gradDur1);
            krad = sqrt(rcLim * rcLim + mDelta * mDelta + 2.0 * rcLim * mDelta * sin(theta));
            phi = phi0 - acos((rcLim * sin(theta) + mDelta) / krad);
          }
        }
      } else {
        sdc[i] = 1.0;
        theta = ctg * sqrt(t);
        krad = mDelta * sqrt(theta * theta + 1.0);
        phi = theta - acos(mDelta / krad);
      }
    }

    kspX[i] = krad * cos(phi) * mGridRes * mRes * coordsScaleFactor;
    kspY[i] = krad * sin(phi) * mGridRes * mRes * coordsScaleFactor;
  }
}

void WPGen::ComputeTimeMap(std::vector<dVector > & timeMap)
{
  //////////////////////////////////////
  // Create a Time map for blur kernels
  // We use these for deblurring
  //////////////////////////////////////
  double mtxEdge = 0.5 / (mGridRes * mRes);
  int mtx2 = gridMtxSize / 2;
  double rcMax = sqrt(kradMax * kradMax - mDelta * mDelta);
  double rcEdge = sqrt(mtxEdge * mtxEdge - mDelta * mDelta);

  // Temporal offsets to make smooth transitions
  double tbaseFrq = ((1 + 5 * M_PI) / (6 * omegaMax)) - (1.0 / omegaMax);
  double tbaseSlw = tbaseFrq + ((mGamma * slewMax) / (3.0 * mDelta * pow(omegaMax, 3)));
  double tbaseGrd = tbaseSlw + ((mGamma * pow(gradMax, 3)) / (6.0 * mDelta * pow(slewMax, 2)));
  double tbaseEdg = tbaseGrd + (rcMax * rcEdge / (2.0 * mGamma * mDelta * gradMax));
  double c2Edg = rcMax / ((rcMax - rcEdge) * (2.0 * mGamma * mDelta * gradMax));
  double tauTop = tbaseEdg;
  double tacq, rad, rc = 0.0, rcEdgeTemp = 0.0;

  // Make sure we catch outer edges of time_out
  timeMap.resize(gridMtxSize, dVector(gridMtxSize, tauTop));

  // Make sure we catch outer edges of time_out
  for (int i = 0; i < mtx2 - 1; i++) {
    for (int j = 0; j < mtx2 - 1; j++) {

      rad = sqrt((double)(i * i + j * j)) / mFov;
      if (rad > mDelta) {
        rc = sqrt(rad * rad - mDelta * mDelta);
      }
      //-----
      // Arc |
      //-----
      if (rad < sqrt(2.0) * mDelta) {
        tacq = (M_PI / (6.0 * omegaMax)) + ((4 * M_PI + 1.0) * rad / (6.0 * sqrt(2.0) * mDelta * omegaMax));
      }
      //-----
      // FRQ |
      //-----
      else if (rc < (mGamma * slewMax) / (omegaMax * omegaMax)) {
        tacq = tbaseFrq + (rc / (mDelta * omegaMax));
      }
      //------
      // SLEW |
      //------
      else if (rc < (mGamma * gradMax * gradMax) / (slewMax)) {
        tacq = tbaseSlw + (2.0 / (3.0 * mDelta)) * sqrt(rc * rc * rc / (mGamma * slewMax));
      }
      //------
      // GRAD |
      //------
      else if (rc < rcMax) {
        tacq = tbaseGrd + (rc * rc / (2.0 * mGamma * mDelta * gradMax));
      }
      //------
      // EDGE |
      //------
      // smooth kernel phase calcs
      // quadratic that fits end of GRAD stage
      else if (rc < rcEdge) {
        rcEdgeTemp = rc - rcEdge;
        tacq = tbaseEdg + c2Edg * rcEdgeTemp * rcEdgeTemp;
      }
      else {
        tacq = tauTop;
      }
      timeMap[mtx2+i][mtx2+j] = tacq;
      timeMap[mtx2-i][mtx2+j] = tacq;
      timeMap[mtx2+i][mtx2-j] = tacq;
      timeMap[mtx2-i][mtx2-j] = tacq;
    }
  }
}

// Destructor
WPGen::~WPGen()
{
  // Add any necessary cleanup code here
}

double WPGen::_roundGrad(double dur)
{
  double roundedDur = ceil(dur / mGradRast) * mGradRast;
  return roundedDur;
}

// Concatenates two input vectors
dVector WPGen::_concatenate(const dVector& input1, const dVector& input2)
{
    dVector combined;

    combined.reserve(input1.size() + input2.size());
    combined.insert(combined.end(), input1.begin(), input1.end());
    combined.insert(combined.end(), input2.begin(), input2.end());

    return combined;
}

int WPGen::_favFFTLength(int inGridSize)
{
    std::vector<int> fftLengths = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 24, 25, 27,
            28, 30, 32, 35, 36, 40, 42, 45, 48, 49, 50, 54, 56, 60, 63, 64, 70, 72,
            75, 80, 81, 84, 90, 96, 98, 100, 105, 108, 112, 120, 125, 126, 128, 135,
            140, 144, 150, 160, 162, 168, 175, 180, 189, 192, 196, 200, 210, 216,
            224, 225, 240, 243, 245, 250, 252, 256, 270, 280, 288, 300, 315, 320,
            324, 336, 343, 350, 360, 375, 384, 392, 400, 405, 420, 432, 448, 450,
            480, 486, 490, 500, 504, 512, 525, 540, 560, 567, 576, 588, 600, 625,
            630, 640, 648, 672, 675, 700, 720, 729, 735, 750, 768, 784, 800, 810,
            840, 864, 875, 896, 900, 945, 960, 972, 980, 1000, 1008, 1024, 1050,
            1080, 1120, 1134, 1152, 1176, 1200, 1215, 1225, 1260, 1280, 1296, 1344};

    size_t length=0 ;               
    for (length = 0 ; length < fftLengths.size(); length++) {
      if (fftLengths[length] > inGridSize) {
        return fftLengths[length];
      }
    }

    return inGridSize; 
}
}