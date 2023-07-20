/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * BLIS is distributed under the Eclipse Public License as part of the       *
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

#include "CoinHelperFunctions.hpp"

#include "BcpsBranchObject.h"
#include "BlisObjectInt.h"
#include "BlisBranchObjectInt.h"
#include "BlisModel.h"


//#############################################################################

/** Useful constructor
    Loads actual upper & lower bounds for the specified variable. */
BlisObjectInt::BlisObjectInt(int objectIndex,
                             int iColumn, 
                             double lb,
                             double ub,
                             double breakEven)
    : 
    BcpsObject()
{
    
    assert (breakEven > 0.0 && breakEven < 1.0);
    objectIndex_ = objectIndex;
    columnIndex_ = iColumn;
    originalLower_ = lb;
    originalUpper_ = ub;
    breakEven_ = breakEven;
}

//#############################################################################

// Copy constructor 
BlisObjectInt::BlisObjectInt( const BlisObjectInt & rhs)
    :
    BcpsObject(rhs)

{
    objectIndex_ = rhs.objectIndex_;
    columnIndex_ = rhs.columnIndex_;
    originalLower_ = rhs.originalLower_;
    originalUpper_ = rhs.originalUpper_;
    breakEven_ = rhs.breakEven_;
}

//#############################################################################

// Assignment operator 
BlisObjectInt & 
BlisObjectInt::operator = (const BlisObjectInt& rhs)
{
    if (this != &rhs) {
	BcpsObject::operator=(rhs);
	columnIndex_ = rhs.columnIndex_;
	breakEven_ = rhs.breakEven_;
    }
    return *this;
}

//#############################################################################

// Compute the infeasibility based on currently relax solution.
double 
BlisObjectInt::infeasibility(BcpsModel *m, int & preferredWay) const
{
    BlisModel * model = dynamic_cast<BlisModel *>(m);
    OsiSolverInterface * solver = model->solver();
    
    const double * solution =  solver->getColSolution();

    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();

    double value = solution[columnIndex_];
    
    value = CoinMax(value, lower[columnIndex_]);
    value = CoinMin(value, upper[columnIndex_]);
    
    //printf("%d %g %g %g %g\n",columnIndex_,value,lower[columnIndex_],
    //   solution[columnIndex_],upper[columnIndex_]);

    double nearest = floor(value + (1.0 - breakEven_));
    double integerTolerance = model->BlisPar()->entry(BlisParams::integerTol);
    if (nearest > value) {
	preferredWay = 1;
    }
    else {
	preferredWay = -1;
    }
    
    double weight = fabs(value - nearest);

    // normalize so weight is 0.5 at break even
    if (nearest < value) {
	weight = (0.5/breakEven_) * weight;
    }
    else {
	weight = (0.5/(1.0 - breakEven_)) * weight;
    }
    
    if (fabs(value - nearest) <= integerTolerance) {
	return 0.0;
    }
    else {
	return weight;
    }
}

//#############################################################################

// Force this object within exiting bounds, then fix the bounds at the
// the nearest integer value. Assume solution value is within tolerance of
// the nearest integer.
void 
BlisObjectInt::feasibleRegion(BcpsModel *m)
{

    BlisModel *model = dynamic_cast<BlisModel* >(m);
    OsiSolverInterface * solver = model->solver();
  
    const double * solution =  solver->getColSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();

    double value = solution[columnIndex_];

    // 1) Force value to be in bounds.
    value = CoinMax(value, lower[columnIndex_]);
    value = CoinMin(value, upper[columnIndex_]);
    
    double nearest = floor(value + 0.5);

    // 2) Fix variable at the nearest integer
    assert (fabs(value - nearest) <= 0.01);
    solver->setColLower(columnIndex_, nearest);
    solver->setColUpper(columnIndex_, nearest);
}

//#############################################################################

// Creates a branching object from this infeasible object.
BcpsBranchObject * 
BlisObjectInt::createBranchObject(BcpsModel *m, int direction) const
{
    BlisModel *model = dynamic_cast<BlisModel* >(m);
    OsiSolverInterface * solver = model->solver();
    
    //double integerTolerance = model->BlisPar()->entry(BlisParams::integerTol);
    
    const double * solution = solver->getColSolution();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    
    double value = solution[columnIndex_];
    //std::cout << "COL"<< columnIndex_ << ": x = " << value << std::endl;
    
    // Force value in bounds.
    value = CoinMax(value, lower[columnIndex_]);
    value = CoinMin(value, upper[columnIndex_]);
    
    return new BlisBranchObjectInt(model, objectIndex_, direction, value);
}

//#############################################################################

// Feasible branching. Given a integeral solution value, return a 
// branching object which would give a new feasible point in direction
// reduced cost says would be better. If no better feasible point,
// return NULL.

BcpsBranchObject * 
BlisObjectInt::preferredNewFeasible(BcpsModel *m) const
{
    BlisModel *model = dynamic_cast<BlisModel* >(m);
    OsiSolverInterface * solver = model->solver();
    
    double value = solver->getColSolution()[columnIndex_];
    
    double nearest = floor(value + 0.5);

    assert(fabs(value - nearest) <= model->BlisPar()->entry(BlisParams::integerTol));

    double dj = solver->getObjSense()*solver->getReducedCost()[columnIndex_];

    BlisBranchObjectInt * object = NULL;

    if (dj >= 0.0) {
	// Better go down
	if (nearest > originalLower_ + 0.5) {
	    // Has room to go down
	    object = new BlisBranchObjectInt(model,
                                             objectIndex_,
                                             -1,
                                             nearest - 1.0,
                                             nearest - 1.0);
	}
    } 
    else {
	// Better go up
	if (nearest < originalUpper_ - 0.5) {
	    // Has room to go up
	    object = new BlisBranchObjectInt(model, 
                                             objectIndex_, 
                                             -1,
                                             nearest + 1.0,
                                             nearest + 1.0);
	}
    }

    return object;
}
  
//#############################################################################

// Feasible branching. Given a integeral solution value, return a 
// branching object which would give a new feasible point in direction
// reduced cost says would be worse. If no better feasible point,
// return NULL.

BcpsBranchObject * 
BlisObjectInt::notPreferredNewFeasible(BcpsModel *m) const 
{
    BlisModel *model = dynamic_cast<BlisModel* >(m);
    OsiSolverInterface * solver = model->solver();
    double value = solver->getColSolution()[columnIndex_];
    
    double nearest = floor(value + 0.5);

    assert (fabs(value-nearest) <= model->BlisPar()->entry(BlisParams::integerTol));
    double dj = solver->getObjSense()*solver->getReducedCost()[columnIndex_];
    BlisBranchObjectInt * object = NULL;

    if (dj <= 0.0) {
	// Worse if go down.
	if (nearest > originalLower_ + 0.5) {
	    // Has room to go down.
	    object = new BlisBranchObjectInt(model,
                                             objectIndex_,
                                             -1,
                                             nearest - 1.0,
                                             nearest - 1.0);
	}
    } 
    else {
	// Worse if go up.
	if (nearest < originalUpper_-0.5) {
	    // Has room to go up.
	    object = new BlisBranchObjectInt(model,
                                             objectIndex_,
                                             -1,
                                             nearest + 1.0,
                                             nearest + 1.0);
	}
    }

    return object;
}

//#############################################################################

// Reset bounds to the ones in LP solver.
void 
BlisObjectInt::resetBounds(BcpsModel *m)
{
    originalLower_ = dynamic_cast<BlisModel *>
	(m)->solver()->getColLower()[columnIndex_];
    originalUpper_ = dynamic_cast<BlisModel *>
	(m)->solver()->getColUpper()[columnIndex_];
}

//#############################################################################
