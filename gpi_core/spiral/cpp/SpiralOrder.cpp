

/**
 * @file SpiralOrder.cpp
 * @brief Implementation of the SpiralArmOrder class for generating spiral arm orders and angles.
 * 
 * This file contains the implementation of the SpiralArmOrder class, which is responsible for generating
 * different types of spiral arm orders and their corresponding angles for MRI sequences. The class supports
 * various order types including linear, skip, two-way, mixed, and golden angle orders. It also handles
 * coherent and incoherent mechanisms for angle generation.
 * @method _linearOrder Generates a linear order of angles.
 * @method _skipOrder Generates a skip order of angles.
 * @method _twowayOrder Generates a two-way order of angles.
 * @method _mixedOrder Generates a mixed order of angles.
 * @method _fillAngles Fills the angle array based on the order array.
 * @method composeAngles Composes the spiral arm angles based on the specified order type.
 * 
 * @date 2024-06-01
 * @version 1.0
 * @author Guruprasad Krishnamoorty

 * 
 * @note Ensure that the angular coverage is within the range [0, 2π]. If not, it will be set to 2π.
 * 
 * @example
 * @code
 * wp::SpiralArmOrder spiralOrder(4, 128, M_PI, wp::OrderType::LINEAR, true, false);
 * std::vector<dVector> spiralArmAngles;
 * spiralOrder.composeAngles(spiralArmAngles);
 * @endcode
 * 
 * @see SpiralOrder.hpp
 */


#include "SpiralOrder.hpp"
namespace wp{
SpiralArmOrder::SpiralArmOrder(int numberArms, int numberEcho, float angularCoverage, OrderType orderType, 
bool inCoherentMech, bool debug)
    : mNArms(numberArms),
      mNEcho(numberEcho),
      mAngCoverage(angularCoverage),
      mOrderType(orderType),
      mInCohMech(inCoherentMech),
      mDebug(debug)
{
    if (mAngCoverage < 0 || mAngCoverage > 2 * M_PI)
    {
        mAngCoverage = 2 * M_PI;
    }

    mTotalAngles = mNArms;

    if(inCoherentMech)
        mTotalAngles *= mNEcho;

    armAngles.resize(mTotalAngles);
    armOrder.resize(mTotalAngles);
}

SpiralArmOrder::~SpiralArmOrder(){}

void SpiralArmOrder::_linearOrder(iVector &orderArray)
{
    for (int i = 0; i < mTotalAngles; i++)
    {
        orderArray[i] = i;
        if (mDebug){
            printf("Order: %d, index: %d\n", orderArray[i], i);
        }
    }
}

void SpiralArmOrder::_skipOrder(iVector &orderArray)
{
    int currOrder = 0;
    for (int i = 0; i < mTotalAngles; i++)
    {
        orderArray[i] = currOrder;

        currOrder += 2;
        if (currOrder >= mTotalAngles)
        {
            currOrder = 1;
        }

        if (mDebug){
            printf("Order: %d, index: %d\n", orderArray[i], i);
        }
    }
}

void SpiralArmOrder::_twowayOrder(iVector &orderArray)
{
    int currOrder = 0;
    int evenCounter = 0;
    int oddCounter = mTotalAngles - 1;
    for (int i = 0; i < mTotalAngles; i++)
    {
        if (i % 2)
        {
            currOrder = oddCounter;
            oddCounter--;
        }
        else
        {
            currOrder = evenCounter;
            evenCounter++;
        }
        orderArray[i] = currOrder;
        if (mDebug){
            printf("Order: %d, index: %d\n", orderArray[i], i);
        }
    }
}

void SpiralArmOrder::_mixedOrder(iVector &orderArray)
{
    int currOrder = 0;
    int evenCounterA = 1;
    int oddCounterA = mTotalAngles - 2;
    int evenCounterB = mTotalAngles - 1;
    int oddCounterB = 0;
    for (int i = 0; i < mTotalAngles; i++)
    {
        if (i < mTotalAngles / 2)
        {
            currOrder = i % 2 ? oddCounterA : evenCounterA;
            oddCounterA -= 2;
            evenCounterA += 2;
        }
        else
        {
            currOrder = i % 2 ? oddCounterB : evenCounterB;
            oddCounterB += 2;
            evenCounterB -= 2;
        }

        orderArray[i] = currOrder;
        if (mDebug){
            printf("Order: %d, index: %d\n", currOrder, i);
        }
    }
}

dVector SpiralArmOrder::_fillAngles(const iVector &orderArray)
{
    dVector angleArray(mTotalAngles, 0.0);

    double angleInc =  mAngCoverage / static_cast<double>(mTotalAngles);
    if (mOrderType == OrderType::GOLDEN)    
        angleInc = GOLDENAGNLE * M_PI / 180.0;

    for (int i = 0; i < mTotalAngles; i++)
    {
        angleArray[i] = fmod(angleInc * orderArray[i], 2 * M_PI);
    }

    return angleArray;
}

void SpiralArmOrder::composeAngles(std::vector<dVector> &spiralArmAngles)
{
    spiralArmAngles.clear();

    iVector orderArray(mTotalAngles);

    switch (mOrderType)
    {
    case OrderType::LINEAR:
        _linearOrder(orderArray);
        break;
    case OrderType::SKIPP:
        _skipOrder(orderArray);
        break;
    case OrderType::TWO_WAY:
        _twowayOrder(orderArray);
        break;
    case OrderType::MIXED:
        _mixedOrder(orderArray);
        break;
    case OrderType::GOLDEN:
        _linearOrder(orderArray);
        break;
    default:
        // Handle invalid order type
        break;
    }

    dVector angleArray = _fillAngles(orderArray);

    spiralArmAngles.resize(mNArms, dVector(mNEcho));

    for (int i = 0; i < mNArms; i++){
        for (int j = 0; j < mNEcho; j++){

            int index = i;

            if(mInCohMech)
                index = i + j * mNArms;

            spiralArmAngles[i][j] = angleArray[index];


            if (mDebug){
                printf("Arm: %d, Echo: %d, Angle: %.2f, index: %d\n", 
                i, j, 180/M_PI * angleArray[index], index);
            }
        }
    }
}

}