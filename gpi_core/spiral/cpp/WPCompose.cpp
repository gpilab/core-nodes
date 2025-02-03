
/**
 * @file WPCompose.cpp
 * @brief Implementation of the WPCompose class for generating gradient waveforms and k-space trajectories.
 * 
 * This file contains the implementation of the WPCompose class, which is responsible for generating
 * gradient waveforms and k-space trajectories for MRI sequences. The class supports both single Tau
 * and dual Tau configurations and can handle different spiral directions (in, out, in-out).
 * 
 * @date 2024-06-01
 * @version 1.0
 * @author Guruprasad Krishnamoorty
 */
#include "WPCompose.hpp"
namespace wp{
// Constructor for single Tau

WPCompose::WPCompose(WPGen* wp1_, int numEchoes_, SpiralDirect spiralType_): 
    mWPOut(wp1_), mNEchoes(numEchoes_),
      mSpDirection(spiralType_)
{   
     // Calculate base gradient waveforms
    mWPOut->ComputeBaseSpiral(mSpiralOutX, mSpiralOutY);
    mWPOut->ComputeRampDown(mRampDownX, mRampDownY);
    mWPOut->ComputeM0M1Trapezoids(mTrapOutX, mTrapOutY);

    mGradPts = mSpiralOutX.size() + mRampDownX.size() + mTrapOutX.size();

    mSpiralInX = _flip(mSpiralOutX);
    mSpiralInY = _flip(mSpiralOutY);
    mRampUpX = _flip(mRampDownX);
    mRampUpY = _flip(mRampDownY);
    mTrapInX = _flip(mTrapOutX);
    mTrapInY = _flip(mTrapOutY);

    
    readOutDur = mWPOut->GetActualTau();
    mTimeEcho1 = 0.0;
    mTrailingDur = mWPOut->GetRampDur() + mWPOut->GetM0M1Dur();
    mLeadingDur = 0.0;
    mNyqArms = mWPOut->GetNyqNumberArms();

    if (mSpDirection == SpiralDirect::SPIRAL_INOUT) {
        mTimeEcho1 = mWPOut->GetTotalDur();
        mGradPts *= 2;
        readOutDur *= 2;

        if (mNEchoes > 1) {
            mDeltaTE = 2.0 * mWPOut->GetActualTau();
            readOutDur *= mNEchoes;
        }
        mLeadingDur = mTrailingDur;
    }

    if (mSpDirection == SpiralDirect::SPIRAL_IN) {
        mTrailingDur = 0.0;
        mLeadingDur = mWPOut->GetRampDur() + mWPOut->GetM0M1Dur();
        mTimeEcho1 = mWPOut->GetTotalDur();
    }

    if (mSpDirection == SpiralDirect::SPIRAL_OUTIN) {
         mTrailingDur = 2.0 * mWPOut->GetActualTau();
         mLeadingDur = 0.0;
         mTimeEcho1 = 0.0;

        mGradPts *= 2;
        readOutDur *= 2;

        if (mNEchoes > 1) {
            mDeltaTE = 2.0 * mWPOut->GetActualTau();
            readOutDur *= mNEchoes;
        }
    }

    mGradMax = mWPOut->GetGradMax();
    mSlewMax = mWPOut->GetSlewMax();

    if (mNyqArms % 2 != 0) 
        mNyqArms++;

    mGradRast = mWPOut->GetGradRast();

    mDualTau = false;

    mCentrNullPts = 0;
}

// Constructor for Dual Tau
WPCompose::WPCompose(WPGen* wp1_, WPGen* wp2_, int centrNullPts_): 
    mWPIn(wp1_), mWPOut(wp2_), mCentrNullPts(centrNullPts_)
{   
    mGradRast = mWPOut->GetGradRast();
    mDwell = mWPOut->GetDwell();

    
    int halfNullPts = mCentrNullPts / 2;
    dVector zeroV(halfNullPts, 0);
    // Calculate base gradient waveforms
    mWPIn->ComputeBaseSpiral(mSpiralInX, mSpiralInY);
    mSpiralInX = _flip(mSpiralInX);
    mSpiralInY = _flip(mSpiralInY);

    if(halfNullPts > 0){
        mSpiralInX = _concatenate(mSpiralInX, zeroV);
        mSpiralInY = _concatenate(mSpiralInY, zeroV);
    }

    mWPIn->ComputeRampDown(mRampUpX, mRampUpY);
    mRampUpX = _flip(mRampUpX);
    mRampUpY = _flip(mRampUpY);

    mWPIn->ComputeM0M1Trapezoids(mTrapInX, mTrapInY);
    mTrapInX = _flip(mTrapInX);
    mTrapInY = _flip(mTrapInY);

    mWPOut->ComputeBaseSpiral(mSpiralOutX, mSpiralOutY);
    mWPOut->ComputeRampDown(mRampDownX, mRampDownY);
    mWPOut->ComputeM0M1Trapezoids(mTrapOutX, mTrapOutY);

    if(halfNullPts > 0){
        mSpiralOutX = _concatenate(zeroV, mSpiralOutX);
        mSpiralOutY = _concatenate(zeroV, mSpiralOutY);
    }

    mGradPts = mTrapInX.size() + mRampUpX.size() + mSpiralInX.size() +
                mSpiralOutX.size() + mRampDownX.size() + mTrapOutX.size();

    mSpDirection = SpiralDirect::SPIRAL_INOUT; // Only In-Out support with dual Tau for now

    mNEchoes = 1; // Only 1 echo supported for now

    mTimeEcho1 = mWPIn->GetTotalDur();

    if(halfNullPts > 0){
        mTimeEcho1 += (halfNullPts * mWPIn->GetGradRast());
    }

    mDeltaTE = 0.0;

    double gMaxIn = mWPIn->GetGradMax();
    double gMaxOut = mWPOut->GetGradMax();
    double sMaxIn = mWPIn->GetSlewMax();
    double sMaxOut = mWPOut->GetSlewMax();

    mGradMax = std::max(gMaxIn, gMaxOut);
    mSlewMax = std::max(sMaxIn, sMaxOut);

    int outArms = mWPOut->GetNyqNumberArms();
    int inArms = mWPIn->GetNyqNumberArms();

    readOutDur = mWPIn->GetActualTau() 
                    + mWPOut->GetActualTau();

    if(halfNullPts > 0)
        readOutDur += (2 * halfNullPts * mWPIn->GetGradRast());

    mNyqArms = std::min(inArms, outArms);

    if (mNyqArms % 2 != 0) 
        mNyqArms++;

    mTrailingDur = mWPOut->GetRampDur() + mWPOut->GetM0M1Dur();
    mLeadingDur = mWPIn->GetRampDur() + mWPIn->GetM0M1Dur();

    mDualTau = true;
}

void WPCompose::composeTimeMap(dVector2D& timeMapIn, dVector2D& timeMapOut) 
{
    timeMapIn.clear();
    timeMapOut.clear();

    if (mDualTau) {
        mWPIn->ComputeTimeMap(timeMapIn);
    } else {
        mWPOut->ComputeTimeMap(timeMapIn);
    }
    mWPOut->ComputeTimeMap(timeMapOut);

    double maxTimeIn = 0.0;
    for (const auto& row : timeMapIn) {
        double rowMax = *std::max_element(row.begin(), row.end());
        if (rowMax > maxTimeIn) {
            maxTimeIn = rowMax;
        }
    }

    timeMapIn = _subtract(maxTimeIn, timeMapIn);

    if (mSpDirection == SpiralDirect::SPIRAL_OUT) {
        timeMapIn.clear();
    } else if (mSpDirection == SpiralDirect::SPIRAL_IN) {
        timeMapOut.clear();
    }
}

void WPCompose::composeKSPnSDC(dVector& kspX, dVector& kspY,
 dVector& sdc, bool cg_crds_sdc,
 iVector& kspLengths) 
{
    kspX.clear();
    kspY.clear();
    sdc.clear();
    kspLengths.clear();

    // 0th index containing ksp-in pts, 
    // 1st index containing central null pts, 
    // 2nd index contaings spiral out pts
    kspLengths.resize(3, 0);


    // Calculate base ksp and sdc
    mWPOut->ComputeKSPnSDC(mKspOutX, mKspOutY, mSdcOut);

    kspLengths[2] = mKspOutX.size();

    if (mDualTau) {
        mWPIn->ComputeKSPnSDC(mKspInX, mKspInY, mSdcIn);

        kspLengths[0] = mKspInX.size();

        mKspInX = _negative(_flip(mKspInX));
        mKspInY = _negative(_flip(mKspInY));
        mSdcIn = _flip(mSdcIn);

        int halfNullPts = mCentrNullPts / 2;
        int halfNullPtsKsp = halfNullPts * (mGradRast / mDwell);
        dVector zeroVKsp(halfNullPtsKsp, 0);

        // Zeros points are added in center only for standard trajectroy to be used
        if ((halfNullPtsKsp > 0) && (!cg_crds_sdc)) {

            kspLengths[1] = mCentrNullPts * (mGradRast / mDwell);

            mKspInX = _concatenate(mKspInX, zeroVKsp);
            mKspInY = _concatenate(mKspInY, zeroVKsp);
            mSdcIn = _concatenate(mSdcIn, zeroVKsp);

            mKspOutX = _concatenate(zeroVKsp, mKspOutX);
            mKspOutY = _concatenate(zeroVKsp, mKspOutY);
            mSdcOut = _concatenate(zeroVKsp, mSdcOut);
        }
    } else {
        mKspInPts = mKspOutPts;

        mKspInX = _negative(_flip(mKspOutX));
        mKspInY = _negative(_flip(mKspOutY));
        mSdcIn = _flip(mSdcOut);
    }

    if (mSpDirection == SpiralDirect::SPIRAL_OUT) {
        kspX = _repeat(mKspOutX, mNEchoes);
        kspY = _repeat(mKspOutY, mNEchoes);
        sdc = mSdcOut;
    } else if (mSpDirection == SpiralDirect::SPIRAL_IN) {
        kspX = mKspInX;
        kspY = mKspInY;
        sdc = mSdcIn;
    } else if (mSpDirection == SpiralDirect::SPIRAL_INOUT) {
        // Setup building blocks
        dVector pcKspOutX, pcKspOutY;
        _phasorConjucate(mKspOutX, mKspOutY, pcKspOutX, pcKspOutY);

        dVector pcKspInX, pcKspInY;
        _phasorConjucate(mKspInX, mKspInY, pcKspInX, pcKspInY);

        dVector kspInOutX = _concatenate(mKspInX, mKspOutX);
        dVector kspInOutY = _concatenate(mKspInY, mKspOutY);

        dVector pcKspInOutX = _concatenate(pcKspInX, pcKspOutX);
        dVector pcKspInOutY = _concatenate(pcKspInY, pcKspOutY);

        dVector sdcInOut = _concatenate(mSdcIn, mSdcOut);

        sdc = sdcInOut;

        kspX = kspInOutX; // 1st Echo
        kspY = kspInOutY; // 1st Echo

        for (int echo = 2; echo <= mNEchoes; echo++) {
            if (echo % 2 == 0) {
                // Echo is even so append a rotated conjugate
                kspX = _concatenate(kspX, pcKspInOutX);
                kspY = _concatenate(kspY, pcKspInOutY);
            } else {
                // Echo is odd so append the base waveform
                kspX = _concatenate(kspX, kspInOutX);
                kspY = _concatenate(kspY, kspInOutY);
            }
        }
    }
    else if (mSpDirection == SpiralDirect::SPIRAL_OUTIN){
    
        // Setup building blocks
        dVector pcKspInX, pcKspInY;
        _phasorConjucate(mKspInX, mKspInY, pcKspInX, pcKspInY);

        dVector kspX1 = _concatenate(mKspOutX, pcKspInX);
        dVector kspY1 = _concatenate(mKspOutY, pcKspInY);
        dVector sdc1 = _concatenate(mSdcOut, mSdcIn);

        for(int echo = 0; echo < mNEchoes; echo++)
        {
            kspX = _concatenate(kspX, kspX1);
            kspY = _concatenate(kspY, kspY1);
            sdc = _concatenate(sdc, sdc1);
        }
	}

}

void WPCompose::composeGradients(dVector& gradX, dVector& gradY) 
{
    gradX.clear();
    gradY.clear();

    if (mSpDirection == SpiralDirect::SPIRAL_OUT) {

        gradX = _concatenate(mSpiralOutX, mRampDownX, mTrapOutX);
        gradX = _repeat(gradX, mNEchoes);

        gradY = _concatenate(mSpiralOutY, mRampDownY, mTrapOutY);
        gradY = _repeat(gradY, mNEchoes);

        if(mNEchoes > 1)
        {
            mDeltaTE = mWPOut->GetTotalDur();
        }    
    } 
    else if (mSpDirection == SpiralDirect::SPIRAL_IN) {

        gradX = _concatenate(mTrapInX, mRampUpX, mSpiralInX);
        gradY = _concatenate(mTrapInY, mRampUpY, mSpiralInY);
    }
    else if (mSpDirection == SpiralDirect::SPIRAL_INOUT) {
        
        // Setup building blocks
        dVector pcSpiralOutX, pcSpiralOutY;
        _phasorConjucate(mSpiralOutX, mSpiralOutY, pcSpiralOutX, pcSpiralOutY);

        dVector pcSpiralInX, pcSpiralInY;
        _phasorConjucate(mSpiralInX, mSpiralInY, pcSpiralInX, pcSpiralInY);
        
        dVector pcRampOutX, pcRampOutY;
        _phasorConjucate(mRampDownX, mRampDownY, pcRampOutX, pcRampOutY);

        dVector pcTrapOutX, pcTrapOutY;
        _phasorConjucate(mTrapOutX, mTrapOutY, pcTrapOutX, pcTrapOutY);

        dVector spiralInOutX = _concatenate(mSpiralInX, mSpiralOutX);
        dVector spiralInOutY = _concatenate(mSpiralInY, mSpiralOutY);

        dVector pcSpiralInOutX = _concatenate(pcSpiralInX, pcSpiralOutX);
        dVector pcSpiralInOutY = _concatenate(pcSpiralInY, pcSpiralOutY);

        dVector rampDownTrapX = _concatenate(mRampDownX, mTrapOutX);
        dVector rampDownTrapY = _concatenate(mRampDownY, mTrapOutY);

        dVector pcRampDownTrapX = _concatenate(pcRampOutX, pcTrapOutX);
        dVector pcRampDownTrapY = _concatenate(pcRampOutY, pcTrapOutY);

        // Composition starts
        // Start = TrapezoidIn + RampUp
        gradX = _concatenate(mTrapInX, mRampUpX);
        gradY = _concatenate(mTrapInY, mRampUpY);

        for(int echo = 1; echo <= mNEchoes; echo++)
        {
            if(echo % 2 == 0)
            {
                // Echo is even so append a rotated conjugate
                gradX = _concatenate(gradX, pcSpiralInOutX);
                gradY = _concatenate(gradY, pcSpiralInOutY);
            }
            else{
                // Echo is odd so append the base waveform
                gradX = _concatenate(gradX, spiralInOutX);
                gradY = _concatenate(gradY, spiralInOutY);
            }
        }
        
        // Finish = Ramp down + Trapezoid
        if (mNEchoes % 2 == 0) {
            // Echoes are even so the rampdown and trapezoids needs to a rotated conjugate
            gradX = _concatenate(gradX, pcRampDownTrapX);
            gradY = _concatenate(gradY, pcRampDownTrapY);
        } 
        else {
            // Echoes are odd
            gradX = _concatenate(gradX, rampDownTrapX);
            gradY = _concatenate(gradY, rampDownTrapY);
        }
    }
    else if(mSpDirection == SpiralDirect::SPIRAL_OUTIN) {
        // Setup building blocks
        dVector pcSpiralInX, pcSpiralInY;
        _phasorConjucate(mSpiralInX, mSpiralInY, pcSpiralInX, pcSpiralInY);

        dVector gradX1 = _concatenate(mSpiralOutX, pcSpiralInX);
        dVector gradY1 = _concatenate(mSpiralOutY, pcSpiralInY);

        for(int echo = 0; echo < mNEchoes; echo++)
        {
            gradX = _concatenate(gradX, gradX1);
            gradY = _concatenate(gradY, gradY1);
        }
    }

    mTotalPts = gradX.size();
    mTotalDur = gradX.size() * mWPOut->GetGradRast();

    mEchoToEndDur = mTotalDur - mTimeEcho1;
}

dVector WPCompose::_repeat(const dVector& input, int nTimes)
{
    dVector concatenated;
    concatenated.reserve(input.size() * nTimes);
    for (int i = 0; i < nTimes; i++) {
        concatenated.insert(concatenated.end(), input.begin(), input.end());
    }
    return concatenated;
}

dVector WPCompose::_flip(const dVector& input)
{
    dVector flipped(input.rbegin(), input.rend());
    return flipped;
}

dVector WPCompose::_negative(const dVector& input)
{
    dVector output = input;
    for (size_t i = 0; i < output.size(); i++) {
        output[i] *= -1;
    }
    return output;
}

// Concatenates two input vectors
dVector WPCompose::_concatenate(const dVector& input1, const dVector& input2)
{
    dVector combined;

    combined.reserve(input1.size() + input2.size());
    combined.insert(combined.end(), input1.begin(), input1.end());
    combined.insert(combined.end(), input2.begin(), input2.end());

    return combined;
}

// Concatenates three input vectors
dVector WPCompose::_concatenate(const dVector& input1, 
const dVector& input2, const dVector& input3)
{
    dVector combined = _concatenate(input1, input2);
    combined.insert(combined.end(), input3.begin(), input3.end());

    return combined;
}

void WPCompose::_phasorConjucate(const dVector& x, const dVector& y,
dVector& phCOnjX, dVector& phCOnjY)
{
    dVector negY = _negative(y);

    double m0Angle = mWPOut->GetM0Angle();
    double phAng = M_PI + 2.0 * m0Angle;

    double phReal = cos(phAng);
    double phImag = sin(phAng);

    phCOnjY = _linearComb(negY, phReal, '+', x, phImag);
    phCOnjX = _linearComb(x, phReal, '-', negY, phImag);
}


/*
 * Performs a linear combination of two input vectors and returns the result.
 * vect3 = vect1 * a (op) vect2 * b;
 * If op is '+', performs addition; if op is '-', performs subtraction.
 */
dVector WPCompose::_linearComb(const dVector& vect1, double a, const char op, const dVector& vect2, double b)
{
    dVector vect3;
    vect3.reserve(vect1.size());
    double result;
    for (size_t i = 0; i < vect1.size(); i++) {
        switch (op) {
            case '+':
                result = vect1[i] * a + vect2[i] * b;
                break;
            case '-':
                result = vect1[i] * a - vect2[i] * b;
                break;
            default:
                // Handle invalid operator
                printf("Invalid operator: %c\n", op);
                return vect3;
        }
        vect3.push_back(result);
    }
    return vect3;
}

double WPCompose::_roundUpGrad(double dur)
{
  double roundedDur = ceil(dur / mGradRast) * mGradRast;
  return roundedDur;
}

// Destructor
WPCompose::~WPCompose() {
    // Cleanup any resources here
    // You can add your cleanup code here, such as releasing memory or closing files
    // For example:
    // delete mWPIn;
    // delete mWPOut;
}

dVector2D WPCompose::_subtract(double value, const dVector2D& matrix)
{
    dVector2D result = matrix;
    for (auto& row : result) {
        for (auto& elem : row) {
            elem = value - elem;
        }
    }
    return result;
}

}