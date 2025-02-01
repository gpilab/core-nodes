
#pragma once
#ifndef WPCOMPOSE_H
#define WPCOMPOSE_H

#include "WPGen.hpp"
#include <iostream>
#include<cmath>


namespace wp{

#define EPSILON 1e-6

class WPCompose {
public:
    enum SpiralDirect {
        SPIRAL_OUT = 1,
        SPIRAL_IN = 2,
        SPIRAL_INOUT = 3,
        SPIRAL_OUTIN = 4
    };

    // Constructor for same Tau spiral
    WPCompose(WPGen* wp1_, int numEchoes_, SpiralDirect spiralType_);
    // Constructor for dual Tau spiral
    WPCompose(WPGen* wp1_, WPGen* wp2_, int centrNullPts_);

    // Destructor
    ~WPCompose();

    void composeGradients(dVector& gradX, dVector& gradY);

    void composeKSPnSDC(dVector& kspX, dVector& kspY, dVector& sdc, bool cg_crds_sdc, iVector& kspLengths);

    void composeTimeMap(dVector2D& timeMapIn, dVector2D& timeMapOut);

    int getNEchoes() const { return mNEchoes; }

    int getKspInPts() const { return mKspInPts; }
    int getKspOutPts() const { return mKspOutPts; }

    int getGradPts() const { return mGradPts; }

    int getGridMtxSize() const { return mWPOut->GetGridMtxSize(); }

    double getTimeEcho1() const { return mTimeEcho1; }
    double getDeltaTE() const { return mDeltaTE; }

    double getTotalDur() const { return mTotalDur; }

    int getTotalPts() const { return mTotalPts; }

    double getGradMax () const { return mGradMax; }

    double getSlewMax () const { return mSlewMax; }

    double getTrailingDur () const{ return mTrailingDur;}

    double getReadOutDur () const{ return readOutDur; }

    double getLeadingDur () const{ return mLeadingDur;}

    double getEchoToEndDur () const{ return mEchoToEndDur;}

    int getNyqArms() const { return mNyqArms; }

private:
    WPGen* mWPIn;
    WPGen* mWPOut;

    int mCentrNullPts;

    int mNEchoes;
    

    SpiralDirect mSpDirection;
    dVector mSpiralOutX, mSpiralOutY;
    dVector mRampDownX, mRampDownY;
    dVector mTrapOutX, mTrapOutY;
    dVector mKspOutX, mKspOutY, mSdcOut;
    int mGradPts;
    double mTimeEcho1;
    double mDeltaTE;
    double mTotalDur;
    double readOutDur;
    int mTotalPts;
    int mNyqArms;
    double mGradMax, mSlewMax;
    double mTrailingDur; // Time from end of readout to end of waveform
    double mLeadingDur; // Time from start of waveform to the start of readout
    double mEchoToEndDur; //Time from echo to the end of waveform
    double mGradRast;
    double mDwell;
    bool mDualTau;
    int mKspInPts, mKspOutPts;

    dVector mKspInX, mKspInY, mSdcIn;

    dVector     mSpiralInX, mSpiralInY, mRampUpX, mRampUpY, mTrapInX, mTrapInY;

    dVector _concatenate(const dVector& input1, const dVector& input2);

    dVector _concatenate(const dVector& input1, const dVector& input2, const dVector& input3);

    dVector _concatenate(const dVector& input1, const dVector& input2, const dVector& input3, const dVector& input4);

    dVector _repeat(const dVector& input, int nTimes);

    dVector _flip(const dVector& input);

    dVector _negative(const dVector& input);

    double _roundUpGrad(double dur);

    void _phasorConjucate(const dVector& x, const dVector& y, dVector& phCOnjX, dVector& phCOnjY);

    dVector _linearComb(const dVector& vect1, double a, const char op, const dVector& vect2, double b);

    dVector2D _subtract(double value, const dVector2D& matrix);
};

struct SpiralGenParams {

    SpiralGenParams() = default;

    double fov = 0.0; //metre
    double res = 0.0; //metre
    double tauin = 0.0; //msec
    double tauout = 0.0; //msec
    double smax = 0.0; //mT/m/s
    double gmax = 0.0; //mT/m
    double omegmax = 0.0; //Omege
    double rast = 0.0; // msec
    double dwell = 0.0; //msec
    double gamma = 0.0; // KHz/T
    wp::WPCompose::SpiralDirect spiraldirect = wp::WPCompose::SpiralDirect::SPIRAL_OUT;
    int nechoes = 0;
    bool outring = false;
    bool addgre = false;
    double vdfactor = 0.0;

    // Overload the != operator
    bool operator!=(const SpiralGenParams& other) const {
        return std::fabs(fov - other.fov) > EPSILON || std::fabs(res - other.res) > EPSILON ||
               std::fabs(tauout - other.tauout) > EPSILON || std::fabs(smax - other.smax) > EPSILON ||
               std::fabs(gmax - other.gmax) > EPSILON || std::fabs(omegmax - other.omegmax) > EPSILON ||
               std::fabs(rast - other.rast) > EPSILON || std::fabs(dwell - other.dwell) > EPSILON ||
               std::fabs(gamma - other.gamma) > EPSILON || spiraldirect != other.spiraldirect ||
               nechoes != other.nechoes || std::fabs(tauin - other.tauin) > EPSILON ||
               outring != other.outring || addgre != other.addgre || std::fabs(vdfactor - other.vdfactor) > EPSILON;
    }

    // Overload the assignment operator
    SpiralGenParams& operator=(const SpiralGenParams& other) {
        if (this != &other) {
            fov = other.fov;
            res = other.res;
            tauin = other.tauin;
            tauout = other.tauout;
            smax = other.smax;
            gmax = other.gmax;
            omegmax = other.omegmax;
            rast = other.rast;
            dwell = other.dwell;
            gamma = other.gamma;
            spiraldirect = other.spiraldirect;
            nechoes = other.nechoes;
            outring = other.outring;
            addgre = other.addgre;
            vdfactor = other.vdfactor;
        }
        return *this;
    }

    // Equality operator
    bool operator==(const SpiralGenParams& other) const {
        return !(*this != other);
    }

    // Validation function
    bool validate() const {
        return fov > 0 && res > 0 && smax > 0 && gmax > 0 && gamma > 0;
    }    
} ;

}
#endif // WPCOMPOSE_H
