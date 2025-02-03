/**
 * @file WPGen.hpp
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


#pragma once
#ifndef WPGEN_H
#define WPGEN_H

#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include "visibility.h"

using dVector = std::vector<double>;
using dVector2D = std::vector<std::vector<double> >;
using iVector = std::vector<int>;

namespace wp{

class DLL_PUBLIC WPGen {
public:



    // Constructor that takes number of arms
    WPGen(double fov_, double res_, int nArms, double sMax_, 
        double gMax_, double omegaMax_, double gradRast_, double dwell_, double gamma_,
        bool addOuterRing_, double alphaVD = 0.0);

    // Constructor that Tau and computes number of arms inturn
    WPGen(double fov_, double res_, double reqTau, double sMax_, 
        double gMax_, double omegaMax_, double gradRast_, double dwell_, double gamma_,
        bool addOuterRing_, double alphaVD = 0.0);

    ~WPGen(); // Destructor

    void ComputeBaseSpiral(dVector& gradX, dVector& gradY);
    void ComputeRampDown(dVector& rampX, dVector& rampY);
    bool ComputeM0M1Trapezoids(dVector& m0m1X, dVector& m0m1Y);
    void ComputeKSPnSDC(dVector& kspX,   dVector& kspY, dVector& sdc);
    void ComputeTimeMap(std::vector<dVector >& timeMap);

    // Getter methods
    int GetSpiralReadPts() const { return spiralReadPts; }
    int GetSpiralRampPts() const { return spiralRampPts; }
    int GetTotalSpiralPts() const { return totalSpiralPts; }
    int GetM0M1Pts() const { return m0m1Pts; }
    int GetKspacePts() const { return kspacePts; }

    double GetOmegaMax() const { return omegaMax; }
    double GetSlewMax() const { return slewMax; }
    double GetGradMax() const { return gradMax; }

    double GetArcDur() const { return arcDur; } 
    double GetFreqDur() const { return omegaDur; }
    double GetSlewDur() const { return slewDur; }
    double GetGradDur() const { return gradDur; }
    double GetRampDur() const { return rampDur; }
    double GetActualTau() const { return readoutDur; }
    int GetNyqNumberArms() const { return static_cast<int>(ceil(nyqNumberArms)); }
    double GetSnrFactor() const { return snrFactor; }
    int GetGridMtxSize() const { return gridMtxSize; }
    double GetOutRingDur() const { return outRingDur; }
    double GetM0Angle() const { return m0Angle; }
    double GetGradRast() const { return (mGradRast); } //msec
    double GetDwell() const { return mDwell; } //msec

    dVector GetKspTrace() const {
        if(mKSPComputed) {
            return kspTrace;
        }
        else{
            dVector dummy;
            printf("Error! ComputeKSPnSDC() must to be called before calling this routine\n");
            return dummy;
        }
         
    }

    double GetTotalDur() const { 
        
        if(mM0M1Computed) {
            return totalDur;
        }
        else{
            printf("Error! ComputeM0M1Trapezoids() must to be called before calling this routine\n");
            return 0.0;
        }
         
    }

    double GetM0M1Dur() const { 
        
        if(mM0M1Computed) {
            return m0m1Dur;
        }
        else{
            printf("Error! ComputeM0M1Trapezoids() must to be called before calling this routine\n");
            return 0.0;
        }
         
    }

private:
   // Private data members
    // Member variables related to imaging parameters
    double mFov; // Field of view
    double mRes; // Resolution
    double mReqTau;
    double nyqNumberArms; // Number of Nyquist ghosting arms
    double slewMaxIn; // Maximum slew allowed by the system
    double gradMaxIn; // Maximum gradient strength allowed by the system
    double omegaMax; // Maximum omega value
    double mGradRast; // Gradient dwell time
    double mDwell; // ADC dwell time
    double mGamma; // Gyromagnetic ratio
    

    // Getter methods
    double kradMax; // Maximum krad value
    double mDelta; // Delta value

    // Member variables related to durations
    double arcDur0; // Arc duration 0
    double arcDur1; // Arc duration 1
    double arcDur; // Arc duration in msec
    bool mAddOuterRing;
    double mAlphaVD;

    double omegaDur0; // Omega constrained duration 0
    double omegaDur1; // Omega constrained duration 1
    double omegaDur; // Omega constrained duration in msec

    double slewDur0; // Slew constrained duration 0
    double slewDur1; // Slew constrained duration 1
    double slewDur; // Slew constrained duration in msec

    double gradDur0; // Grad max constrained duration 0
    double gradDur1; // Grad max constrained duration 1
    double gradDur; // Grad max constrained duration in msec

    double rampDur; // Ramp down duration in msec
    double readoutDur; // Actual readout duration in msec

    double m0m1Dur;

    double totalDur; // Total waveform duration includin spiral, ramp, trapezoids

    double outRingDur;
    double outRingDur0;
    double outRingDur1;
    double thetaPrimeRing;

    // Member variables related to spiral waveform
    int spiralReadPts; // Number of readout points in the spiral waveform
    int spiralRampPts; // Number of ramp points in the spiral waveform
    int totalSpiralPts; // Total number of points in the spiral waveform
    int m0m1Pts; // Number of points for M0 and M1 trapezoids
    int kspacePts; // Number of points in k-space

    double slewMax; // Maximum slew rate of the spiral waveform
    double gradMax; // Maximum gradient strength of the wavform including spiral and trapezoids
    double gradMaxSpiral; // Maximum gradient strength of the spiral waveform

    
    double snrFactor; // Signal-to-Noise Ratio (SNR) factor
    double mGridRes; // Grid resolution
    double mTrueRes; // True resolution
    double mx0, mx1, my0, my1; // M0 and M1 trapezoid areas
    double m0Angle;
    int gridMtxSize; // Grid Matrix size

    bool mRampDownComputed;
    bool mBaseSpiralComputed;
    bool mKSPComputed;
    bool mM0M1Computed;

    double rcLim;

    dVector kspTrace;

    

    double mTotalCGDur;
    double mArcDurCG;
    double mCtaCG;
    double mG0CG;
    double mDwellCG;

    double coordsScaleFactor;

    //Private functions
    double _roundGrad(double dur);

    void _findCompatibleConstraints();
    void _computeTimings();
    void _initializeParameters();
    void _initialize();
    dVector _concatenate(const dVector& input1, const dVector& input2);

      double smoothStep(double x) {
   return x * x * (3 - 2 * x); // Smoother transition between 0 and 1
  }

  int _favFFTLength(int inGridSize);
};
}
#endif // WPGEN_H
