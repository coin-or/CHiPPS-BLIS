/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * ALPS is distributed under the Eclipse Public License as part of the       *
 * COIN-OR repository (http://www.coin-or.org).                              *
 *                                                                           *
 * Authors:                                                                  *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Conceptual Design:                                                        *
 *                                                                           *
 *          Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *          Laszlo Ladanyi, IBM T.J. Watson Research Center                  *
 *          Matthew Saltzman, Clemson University                             *
 *                                                                           * 
 *                                                                           *
 * Copyright (C) 2001-2023, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include <cmath>
#include <cassert>

#include "Alps.h"
#include "Blis.h"

#include "BlisPseudo.h"

//#############################################################################

void 
BlisPseudocost::update(const int dir,
                       const double parentObjValue,
                       const double objValue,
                       const double solValue)
{
    
    double objDiff = objValue - parentObjValue;    

#ifdef BLIS_DEBUG
    assert(objDiff/(1.0+objValue) >= -1.0e-5);
#endif

    update(dir, objDiff, solValue);    
}

//#############################################################################

void
BlisPseudocost::update(int dir,
                       double objDiff,
                       double solValue)
{
    double fraction;
    double cost;
    
#ifdef BLIS_DEBUG
    if (objDiff < -1.0e-1) {
        std::cout << "objDiff=" << objDiff
                  << ", dir="<< dir << ", x=" << solValue <<std::endl;
        assert(0);
    }
#endif

    if (objDiff < 0.0) {
        return;
    }

    if (dir == 1) {
        fraction = ceil(solValue) - solValue;
        if (fraction >= 1.0e-5) {
            cost = objDiff / (fraction + 1.0e-9);
            upCost_ = (upCost_ * upCount_ + cost) / (1 + upCount_);
            ++upCount_;
        }
        else {
#ifdef BLIS_DEBUG
            ALPS_PRINTF("WARNING: small fraction %.10g, shouldn't happen.\n",
                        fraction);
            assert(0);
#endif
        }
    }
    else if (dir == -1) {
        fraction = solValue - floor(solValue);
        if (fraction >= 1.0e-5) {
            cost = objDiff / (fraction + 1.0e-9);
            downCost_ = (downCost_ * downCount_ + cost) / (1 + downCount_);
            ++downCount_;
        }
        else {
#ifdef BLIS_DEBUG
            ALPS_PRINTF("WARNING: small fraction %.10g, shouldn't happen.\n",
                        fraction);
            assert(0);
#endif
        }
    }
    else {
        ALPS_PRINTF("ERROR: wrong direction %d.\n", dir);
        assert(0);
    }
    
    score_ = weight_* ALPS_MIN(upCost_, downCost_) + 
        (1.0 - weight_) * ALPS_MAX(upCost_, downCost_);
    
}

//#############################################################################

void 
BlisPseudocost::update(double upCost, int upCount, 
                       double downCost, int downCount)
{
    if (upCount) {
        upCount_ += upCount;
        upCost_ = (upCost_ + upCost)/upCount;
    }

    if (downCount) {
        downCount_ += downCount;
        downCost_ = (downCost_ + downCost)/downCount_;
    }

    score_ = weight_* ALPS_MIN(upCost_, downCost_) + 
        (1.0 - weight_) * ALPS_MAX(upCost_, downCost_);
}

//#############################################################################

/** Pack pseudocost to the given object. */
AlpsReturnStatus 
BlisPseudocost::encodeTo(AlpsEncoded *encoded) const 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    
    encoded->writeRep(weight_);
    encoded->writeRep(upCost_);
    encoded->writeRep(upCount_);
    encoded->writeRep(downCost_);
    encoded->writeRep(downCount_);
    encoded->writeRep(score_);

#if 0
    std::cout << "encodeTo: weight_=" << weight_ 
              << ", upCost_=" << upCost_
              << ", upCount_=" << upCount_
              << ", downCost_=" << downCost_ 
              << ", downCount_=" << downCount_
              << ", score_=" << score_
              << std::endl;
#endif

    return status;
}

//#############################################################################

/** Unpack pseudocost from the given encode object. */
AlpsReturnStatus 
BlisPseudocost::decodeFrom(AlpsEncoded &encoded)
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    double weight, upCost, downCost, score;
    int upCount, downCount;   

    encoded.readRep(weight);
    encoded.readRep(upCost);
    encoded.readRep(upCount);
    encoded.readRep(downCost);
    encoded.readRep(downCount);
    encoded.readRep(score);

    update(upCost, upCount, downCost, downCount);

#if 0
    std::cout << "decodeFrom: weight=" << weight 
              << ", upCost=" << upCost
              << ", upCount=" << upCount
              << ", downCost=" << downCost
              << ", downCount=" << downCount 
              << ", score=" << score
              << std::endl;
#endif

    return status;
}

//#############################################################################

AlpsEncoded*
BlisPseudocost::encode() const 
{
    AlpsEncoded* encoded = new AlpsEncoded(BLIS_PSEUDO);
    
    encoded->writeRep(weight_);
    encoded->writeRep(upCost_);
    encoded->writeRep(upCount_);
    encoded->writeRep(downCost_);
    encoded->writeRep(downCount_);
    encoded->writeRep(score_);
    
    return encoded;
}

//#############################################################################

AlpsKnowledge* 
BlisPseudocost::decode(AlpsEncoded& encoded) const 
{
    double weight;
    int upCount;
    double upCost;
    int downCount;
    double downCost;
    double score;

    encoded.readRep(weight);
    encoded.readRep(upCost);
    encoded.readRep(upCount);
    encoded.readRep(downCost);
    encoded.readRep(downCount);
    encoded.readRep(score);

    BlisPseudocost *pcost = new BlisPseudocost(upCost,
                                               upCount,
                                               downCost,
                                               downCount,
                                               score);
    pcost->setWeight(weight);

    return pcost;
}

//#############################################################################
