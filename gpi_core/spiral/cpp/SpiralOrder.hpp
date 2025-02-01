#pragma once
#ifndef ARMORDER_H
#define ARMORDER_H


#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <vector>

#define GOLDENAGNLE  137.5077640500378546463487


using dVector = std::vector<double>;
using iVector = std::vector<int>;

namespace wp{

class SpiralArmOrder {
public:
    enum OrderType {
        LINEAR = 1,
        SKIPP = 2,
        TWO_WAY = 3,
        MIXED = 4,
        GOLDEN = 5
    };

    SpiralArmOrder(int numberArms, int numberEcho, float angularCoverage, OrderType orderType, bool inCoherentMech, bool debug);

    void composeAngles(std::vector<dVector> &spiralArmAngles);

    virtual ~SpiralArmOrder();

private:
    int mNArms;
    int mNEcho;
    float mAngCoverage;
    OrderType mOrderType;

    int mTotalAngles;
    
    dVector armAngles;
    iVector armOrder;
    bool mInCohMech;
    bool mDebug;

   

    void _linearOrder(iVector& orderArray);
    void _skipOrder(iVector& orderArray);
    void _twowayOrder(iVector& orderArray);
    void _mixedOrder(iVector& orderArray);

    dVector _fillAngles(const iVector& orderArray);

};

}

#endif // ARMORDER_H