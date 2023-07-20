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


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#include "CoinTime.hpp"
#include "OsiSolverInterface.hpp"
#include "CglProbing.hpp"

#include "BcpsObjectPool.h"

#include "BlisHelp.h"
#include "BlisModel.h"
#include "BlisMessage.h"
#include "BlisConGenerator.h"
#include "BlisConstraint.h"

//#############################################################################

// Normal constructor
BlisConGenerator::BlisConGenerator(BlisModel * model,
				   CglCutGenerator * generator,
				   const char * name,
				   BlisCutStrategy strategy,
				   int cutGenerationFrequency,
				   bool normal, 
				   bool atSolution, 
				   bool infeasible)
{
    model_ = model;
    generator_ = generator;
    generator_->refreshSolver(model_->solver());
    
    if (name) {
        name_ = name;
    }
    else {
        name_ = "UNKNOWN";
    }
    
    strategy_ = strategy;
    cutGenerationFrequency_ = cutGenerationFrequency;
    normal_ = normal;
    atSolution_ = atSolution;
    whenInfeasible_ = infeasible;
    numConsGenerated_ = 0;
    numConsUsed_ = 0;
    time_ = 0.0;
    calls_ = 0;
    noConsCalls_ = 0;
}

//#############################################################################

// Copy constructor 
BlisConGenerator::BlisConGenerator(const BlisConGenerator & rhs)
{
    model_ = rhs.model_;
    generator_ = rhs.generator_;
    generator_->refreshSolver(model_->solver());
    strategy_ = rhs.strategy_;
    cutGenerationFrequency_ = rhs.cutGenerationFrequency_;
    name_ = rhs.name_;
    normal_ = rhs.normal_;
    atSolution_ = rhs.atSolution_;
    whenInfeasible_ = rhs.whenInfeasible_;
    numConsGenerated_ = 0;
    numConsUsed_ = 0;
    time_ = 0.0;
    calls_ = 0;
    noConsCalls_ = 0;
}

//#############################################################################

// Assignment operator 
BlisConGenerator & 
BlisConGenerator::operator=( const BlisConGenerator& rhs)
{
    if (this != &rhs) {
        delete generator_;
        model_ = rhs.model_;
        generator_ = rhs.generator_;
        generator_->refreshSolver(model_->solver());
        strategy_ = rhs.strategy_;
	cutGenerationFrequency_ = rhs.cutGenerationFrequency_;
        name_ = rhs.name_;
        normal_ = rhs.normal_;
        atSolution_ = rhs.atSolution_;
        whenInfeasible_ = rhs.whenInfeasible_;
        numConsGenerated_ = 0;
        numConsUsed_ = 0;
        time_ = 0.0;
        calls_ = 0;
        noConsCalls_ = 0;
    }
    
    return *this;
}

//#############################################################################

// This is used to refresh any inforamtion.
// It also refreshes the solver in the con generator
// in case generator wants to do some work 

void 
BlisConGenerator::refreshModel(BlisModel * model)
{
  model_ = model;
  generator_->refreshSolver(model_->solver());
}

//#############################################################################

// Generate constraints for the model data contained in si.
// The generated constraints are inserted into and returned in the 
// constraint pool.
// Default implementation use Cgl cut generators.
bool
BlisConGenerator::generateConstraints(BcpsConstraintPool &conPool)
{
    bool status = false;
    
    OsiSolverInterface * solver = model_->solver();
    
#if defined(BLIS_DEBUG_MORE)
    std::cout << "model_->getNodeCount() = " << model_->getNodeCount()
              << std::endl;
#endif
    
    //--------------------------------------------------
    // Start to generate constraints...
    //--------------------------------------------------

    assert(generator_ != NULL);
    
    int j;
    OsiCuts newOsiCuts;
    CglProbing* generator = dynamic_cast<CglProbing *>(generator_);
    
    if (generator) {
	// It is CglProbing - return tight column bounds
        CglTreeInfo info;
	generator->generateCutsAndModify(*solver, newOsiCuts, &info);
	const double * tightLower = generator->tightLower();
	const double * lower = solver->getColLower();
	const double * tightUpper = generator->tightUpper();
	const double * upper = solver->getColUpper();
	const double * solution = solver->getColSolution();
	
	int numberColumns = solver->getNumCols();
	double primalTolerance = 1.0e-8;
	for (j = 0; j < numberColumns; ++j) {
	    if ( (tightUpper[j] == tightLower[j]) &&
		 (upper[j] > lower[j]) ) {
		// Fix column j
		solver->setColLower(j, tightLower[j]);
		solver->setColUpper(j, tightUpper[j]);
		if ( (tightLower[j] > solution[j] + primalTolerance) ||
		     (tightUpper[j] < solution[j] - primalTolerance) ) {
		    status = true;
		}
	    }
	}
    }
    else {
	// Other Cgl cut generators
	generator_->generateCuts(*solver, newOsiCuts);
    }
    
    //--------------------------------------------------
    // Create blis constraints and remove zero length row cuts.
    //--------------------------------------------------
    
    int numNewConstraints = newOsiCuts.sizeRowCuts();
    for (j = 0; j < numNewConstraints; ++j) {
	OsiRowCut & rCut = newOsiCuts.rowCut(j);
	int len = rCut.row().getNumElements();

#ifdef BLIS_DEBUG_MORE
	std::cout << "Cut " << j<<": length = " << len << std::endl;
#endif
	if (len > 0) {
	    // Create BlisConstraints from OsiCuts.
	    BlisConstraint *blisCon = BlisOsiCutToConstraint(&rCut);
	    conPool.addConstraint(blisCon);
	}
	else if (len == 0) {
	    // Empty cuts
#ifdef BLIS_DEBUG
	    std::cout << "WARNING: Empty cut from " << name_ << std::endl;
#endif
	}
	else {
#ifdef BLIS_DEBUG
	    std::cout << "ERROR: Cut length = " << len << std::endl;
#endif
	    // Error
	    assert(0);
	}
    }

    // Adjust cut strategy.
    if ( (strategy_ == BlisCutStrategyAuto) &&
	 (noConsCalls_ > BLIS_CUT_DISABLE) ) {
	strategy_ = BlisCutStrategyNone;
    }

    return status;
}

//#############################################################################
