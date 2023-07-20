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

#include <cassert>
#include <iostream>
#include <utility>
#include <cmath>
#include <vector>

#include "CoinUtility.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiCuts.hpp"

#include "AlpsKnowledge.h"
#include "AlpsEnumProcessT.h"
#include "AlpsKnowledgeBroker.h"
#include "AlpsTreeNode.h"

#include "BcpsBranchStrategy.h"

#include "Blis.h"
#include "BlisBranchObjectInt.h"
#include "BlisBranchObjectBilevel.h"
#include "BlisConstraint.h"
#include "BlisHelp.h"
#include "BlisTreeNode.h"
#include "BlisModel.h"
#include "BlisNodeDesc.h"
#include "BlisModel.h"
#include "BlisObjectInt.h"
#include "BlisParams.h"
#include "BlisSolution.h"
//#include "BlisVariable.h"

#define REMOVE_SLACK 1
#define BLIS_SLACK_MAX 4

//#############################################################################

AlpsTreeNode*
BlisTreeNode::createNewTreeNode(AlpsNodeDesc *&desc) const
{
    BlisModel* model = dynamic_cast<BlisModel*>(desc_->getModel());
    BcpsBranchStrategy *strategy = model->branchStrategy();

    BlisBranchingStrategy type = 
	static_cast<BlisBranchingStrategy>(strategy->getType());

    switch(type){
	
     case BlisBranchingStrategyMaxInfeasibility:
     case BlisBranchingStrategyPseudoCost:
     case BlisBranchingStrategyReliability:
     case BlisBranchingStrategyStrong:
	 {
	     double estimate = solEstimate_;

	     // Set solution estimate for this nodes.
	     // double solEstimate = quality_ + sum_i{min{up_i, down_i}}
	     int branchDir=dynamic_cast<BlisNodeDesc *>(desc)->getBranchedDir();
	     int branchInd=dynamic_cast<BlisNodeDesc *>(desc)->getBranchedInd();
	     double lpX = dynamic_cast<BlisNodeDesc *>(desc)->getBranchedVal();
	     double f = lpX - floor(lpX);
	     assert(f > 0.0);
	     
	     int objInd = model->getIntObjIndices()[branchInd];
	     BlisObjectInt *obj = dynamic_cast<BlisObjectInt *>(model->objects(objInd));
	     
	     if (branchDir == -1) {
		 estimate -= (1.0-f) * obj->pseudocost().getUpCost();
	     }
	     else {
		 estimate -= f * obj->pseudocost().getDownCost();
	     }
	     
#ifdef BLIS_DEBUG_MORE
	     printf("BLIS:createNewTreeNode: quality=%g, solEstimate=%g\n",
		    quality_, solEstimate_);
#endif
	     break;
	 }

     case BlisBranchingStrategyBilevel:

	break;

    }

    // Create a new tree node
    BlisTreeNode *node = new BlisTreeNode(desc);    
    desc = NULL;
	
    return node;
}

//#############################################################################

// NOTE: if rampup,
// - parent must be explicit if not NULL,
// - this node is explicit.

int
BlisTreeNode::process(bool isRoot, bool rampUp)
{
    BlisReturnStatus returnStatus = BlisReturnStatusUnknown;
    BlisLpStatus lpStatus = BlisLpStatusUnknown;
    int j, k = -1;
    int numCols, numRows, numCoreCols, numCoreRows;
    int numStartRows, origNumStartRows;
    int maxNumCons, tempNumCons = 0;

    int numIntInfs = 0;
    int numObjInfs = 0;

    int voilatedNumCons = 0;
    int origNumOldCons = 0;
    int currNumOldCons = 0;
    int newNumCons = 0;
    int maxNewNumCons = 0;

    int numAppliedCons = 0;
    int cutStrategy;
    
    // Only autmatic stategy has depth limit.
    int maxConstraintDepth = 20;

    int numPassesLeft = 0;      
    int bStatus = -1;
 
    double cutoff = ALPS_INC_MAX;
    double parentObjValue = getQuality();
    double preObjValue = -ALPS_OBJ_MAX;
    double improvement = 100.0;

    double *currLpSolution = NULL;

    bool keepOn = true;
    bool needBranch = false;
    bool lpFeasible = false;
    bool foundSolution = false;
    bool genConsHere = false;
    bool shareCon = false;
    bool shareVar = false;

    CoinWarmStartBasis::Status rowStatus;
    BlisConstraint *aCon = NULL;
    BcpsObject **newConstraints = NULL;

    BlisSolution *ipSol = NULL;
    
    int numDelRows = 0;
    int *delRow = NULL;
    int *oldConsPos = NULL;
    
    std::vector<int> delIndices;
    
    BlisModel* model = dynamic_cast<BlisModel*>(desc_->getModel());

    AlpsPhase phase = knowledgeBroker_->getPhase();

    int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);
    int hubMsgLevel = model->AlpsPar()->entry(AlpsParams::hubMsgLevel);
    int workerMsgLevel = model->AlpsPar()->entry(AlpsParams::workerMsgLevel);

    BlisParams * BlisPar = model->BlisPar();

    int maxPass = BlisPar->entry(BlisParams::cutPass);
    int quickCutPass = BlisPar->entry(BlisParams::quickCutPass);

    double tailOffTol = BlisPar->entry(BlisParams::tailOff);

    if (msgLevel >= 100){
	std::cout << std::endl;
	std::cout << "**********************************" << std::endl;
	std::cout << "* NOW PROCESSING NODE " << index_ << std::endl;
	std::cout << "**********************************" << std::endl;
	printf("Initial Bound:  %.3f\n", quality_);
    }

    //if ( (quickCutPass > 0 ) && (phase == AlpsPhaseRampup) ) {
    if (quickCutPass > 0) {
        std::cout << "+++ Quick branching; pass is "<<quickCutPass<<std::endl;
        maxPass = quickCutPass + 1;
    }
    else if (maxPass < ALPS_INT_MAX) {
        // Add one to solve initial subprolem.
        ++maxPass;
    }
    
    shareCon = BlisPar->entry(BlisParams::shareConstraints);
    shareVar = BlisPar->entry(BlisParams::shareVariables);    

    cutoff = model->getCutoff();

    numPassesLeft = model->getNumBranchResolve();
    
    // std::cout << "numPassesLeft = " << numPassesLeft << std::endl;
    
    //------------------------------------------------------
    // Check if this can be fathomed by objective cutoff.
    //------------------------------------------------------
    
#if 0
    std::cout << "parentObjValue = " << parentObjValue 
              << "cutoff = " << cutoff << std::endl;
#endif

    if (parentObjValue > cutoff && !isRoot) {
	setStatus(AlpsNodeStatusFathomed);
        //std::cout << "fathom!" <<std::endl;
	goto TERM_PROCESS;
    }
    else {
        model->setActiveNode(this);
        model->addNumNodes();
    }
    
    //------------------------------------------------------
    // Get model information and parameters.
    //------------------------------------------------------

    numCols = model->solver()->getNumCols();
    numRows = model->solver()->getNumRows();
    numCoreCols = model->getNumCoreVariables();
    numCoreRows = model->getNumCoreConstraints();
    numAppliedCons =  numRows - numCoreRows;
    
    maxNumCons = model->getMaxNumCons();

    currLpSolution = new double [numCols];

    //------------------------------------------------------
    // Decides if can generate constraints.
    //------------------------------------------------------

    // Mark if this node is root or not.
    model->isRoot_ = isRoot;

    genConsHere = false;

    cutStrategy = model->getCutStrategy();

    assert(cutStrategy != BlisCutStrategyNotSet);
    
    if (cutStrategy == BlisCutStrategyNone) {
	genConsHere = false;
    }
    else if (cutStrategy == BlisCutStrategyRoot) {
	// The original root only
	if (isRoot && (index_ == 0)) genConsHere = true;
    }
    else if (cutStrategy == BlisCutStrategyAuto) {
	if (depth_ < maxConstraintDepth) {
            genConsHere = true;
	}
    }
    else if (cutStrategy == BlisCutStrategyPeriodic) {
	genConsHere = true;
    }
    else {
	genConsHere = true;
    }

    if (genConsHere && (phase == AlpsPhaseRampup)) {
	if (!(BlisPar->entry(BlisParams::cutRampUp))) {
	    genConsHere = false;
	}
    }

    //--------------------------------------------------
    // Call HeurisBounding before solving for the first node.
    //--------------------------------------------------
    
    if (model->getNumNodes() == 1) {
	int heurStatus = callHeuristics(model, true); // before root
	if (heurStatus == 1) {
	    cutoff = model->getCutoff();
	}
	else if (heurStatus == 2) {
	    // Fathom this node
	    goto TERM_PROCESS;
	}
    }

    //======================================================
    // Restore, load and solve the subproblem.
    // (1) LP infeasible
    //     a. set status to be fathom.
    // (2) LP feasible
    //     a. MILP feasible. Check whether need update incumbent.
    //     b. LP feasible but not MIP feasible. Check whether can be 
    //        fathomed, if not, choose a branch variable.
    //======================================================

    //------------------------------------------------------
    // Extract info from this node and load subproblem into lp solver.
    //------------------------------------------------------
    
    installSubProblem(model);

    //------------------------------------------------------
    // Bounding, heuristic searching, constraint generating.
    //------------------------------------------------------

    numStartRows = model->solver()->getNumRows();
    origNumStartRows = numStartRows;

    origNumOldCons = numStartRows - numCoreRows;
    currNumOldCons = origNumOldCons;
    
#if 0
    std::cout << "PROCESS: genConsHere =" << genConsHere
	      << ", cut strategy =" << model->getCutStrategy()
	      << ", numCoreRows =" << numCoreRows
	      << ", numStartRows =" << numStartRows
	      << ", currNumOldCons =" << currNumOldCons 
              << ", index_ = " << index_ << std::endl;
#endif
    
    if (currNumOldCons > 0) {
	oldConsPos = new int [currNumOldCons];
	for (k = 0; k < currNumOldCons; ++k) {
	    oldConsPos[k] = k;
	}
    }

    if (genConsHere) {
	if (maxNumCons > ALPS_INT_MAX - 100) {
	    maxNewNumCons = 10000;
	}
	else {
	    maxNewNumCons = maxNumCons;
	}
        newConstraints = new BcpsObject* [maxNewNumCons];    
    }

    model->boundingPass_ = 0;
    while (keepOn && (model->boundingPass_ < maxPass)) {
        ++(model->boundingPass_);
        keepOn = false;
        
        //--------------------------------------------------
        // Bounding to get the quality of this node.
        //--------------------------------------------------
        
        if ( ((knowledgeBroker_->getProcType() == AlpsProcessTypeMaster) || 
              (knowledgeBroker_->getProcType() == AlpsProcessTypeSerial)) && 
             isRoot && (model->boundingPass_ == 1) ) {
            if (msgLevel >= 50) {
                model->blisMessageHandler()->message(BLIS_ROOT_PROCESS, 
                                                     model->blisMessages())
                    << model->getNumRows()
                    << model->getNumCols()
                    << CoinMessageEol;
                model->solver()->messageHandler()->setLogLevel(1);
            }
            else {
                model->solver()->messageHandler()->setLogLevel(-1);
            }
            getKnowledgeBroker()->tempTimer().start();
        }
        
	// Get lowe bound.
        lpStatus = static_cast<BlisLpStatus> (bound(model));

	if (model->boundingPass_ == 1) {
	    int iter = model->solver()->getIterationCount();
	    model->addNumIterations(iter);
            if (isRoot) {
                getKnowledgeBroker()->tempTimer().stop();
                if (((knowledgeBroker_->getProcType()==AlpsProcessTypeMaster)||
                     (knowledgeBroker_->getProcType()==AlpsProcessTypeSerial))
		    && (msgLevel >= 50)) {
                    model->solver()->messageHandler()->setLogLevel(0);
                    model->blisMessageHandler()->message(BLIS_ROOT_TIME, 
							 model->blisMessages())
                        << getKnowledgeBroker()->tempTimer().getCpuTime() 
                        << CoinMessageEol;
                }
            }
	}
        
        switch(lpStatus) {
        case BlisLpStatusOptimal:

            // Check if IP feasible 
	    ipSol = model->feasibleSolution(numIntInfs, numObjInfs);

            if (ipSol) {         
                // IP feasible
                model->storeSolution(BlisSolutionTypeBounding, ipSol);
                // Update cutoff
                cutoff = model->getCutoff();
                setStatus(AlpsNodeStatusFathomed);
		goto TERM_PROCESS;
            }
            else {
                if (quality_ > cutoff) {
                    setStatus(AlpsNodeStatusFathomed);
                    goto TERM_PROCESS;
                }
                needBranch = true;
                reducedCostFix(model);
                
                //------------------------------------------
                // Check if tailoff or have to keep on.
                //------------------------------------------

                if (model->boundingPass_ > 1) {
                    improvement = quality_ - preObjValue;
                    if (improvement > tailOffTol) {
                        // NOTE: still need remove slacks, although
                        //       tailoff.
                        keepOn = true;
                    }
                    
#if 0
                    std::cout << "PROCESS: boundingPass_["
                              << model->boundingPass_ << "], improvement=" 
                              << improvement << ", tailOffTol=" << tailOffTol
                              << ", preObj=" << preObjValue 
                              << ", newObj=" << quality_
                              << std::endl;
#endif
                }
                else {
                    keepOn = true;
                }

                if (model->boundingPass_ == maxPass) {
                    if (!numIntInfs && !numObjInfs) {
                        /* No fractional objects, have to process one more pass.
                           For applications like VRP. */
                        ++maxPass;
                        keepOn = true;
                    }
                }

                //------------------------------------------
                // Update previous objective value.
                //------------------------------------------

                preObjValue = quality_;
                
                //------------------------------------------
                // Remove non-core slack constraints. 
                //------------------------------------------

                numRows = model->getNumRows();
                
                if ( genConsHere &&
                     //(improvement > tailOffTol) && 
                     //(numRows > numStartRows) ) {
                     (numRows > numCoreRows) ) {   
                 
#if 1                  
                    if ( (numStartRows + newNumCons != numRows) ||
                         (numCoreRows+currNumOldCons +newNumCons != numRows) ){
                        
                        std::cout << "ERROR: numRows=" << numRows
                                  << "; numCoreRows=" << numCoreRows
                                  << "; numStartRows=" << numStartRows
                                  << "; newNumCons=" << newNumCons
                                  << "; currNumOldCons=" << currNumOldCons
                                  << std::endl;
                        
                        assert(numRows - numStartRows == newNumCons);
                    }
#endif
                    
                    int *oldDelMark = NULL;
                    if (currNumOldCons > 0) {
                        oldDelMark = new int [currNumOldCons];
                        CoinZeroN(oldDelMark, currNumOldCons);
                    }
                    int *newDelMark = NULL;
                    if (newNumCons > 0) {
                        newDelMark = new int [newNumCons];
                        CoinZeroN(newDelMark, newNumCons);
                    }
                    
                    const CoinWarmStartBasis* ws= 
                        dynamic_cast<CoinWarmStartBasis*>
                        (model->solver()->getWarmStart());
                    
                    // Make sure delIndices is empty.
                    assert(delIndices.size()==0);

#if REMOVE_SLACK
                    for (k = numCoreRows; k < numRows; ++k) {
                        rowStatus = ws->getArtifStatus(k);

                        if (rowStatus == CoinWarmStartBasis::basic) {
			    int count;
                            if (k < numStartRows) {
				BlisConstraint *tmpCon = 
				    model->oldConstraints()[(k-numCoreRows)];
				count = tmpCon->getNumInactive() + 1;
				tmpCon->setNumInactive(count);
				if (tmpCon->getNumInactive() > BLIS_SLACK_MAX){
				    oldDelMark[(k-numCoreRows)] = 1;
				    delIndices.push_back(k);
				}
                            }
                            else {
				BcpsObject* tmpCon = 
				    newConstraints[(k-numStartRows)];
				count = tmpCon->getNumInactive() + 1;
				tmpCon->setNumInactive(count);
				if (tmpCon->getNumInactive() > BLIS_SLACK_MAX){
				    newDelMark[(k-numStartRows)] = 1;
				    delIndices.push_back(k);
				}
                            }
                        }
                    }
#endif
                    numDelRows = static_cast<int> (delIndices.size());
		    
		    if (msgLevel > 100) {
			std::cout << "PROCESS: new cuts=" << newNumCons 
				  << ", delete slack cuts=" << numDelRows 
				  << ", numRows=" << numRows 
				  << ", numStartRows=" <<numStartRows 
				  << std::endl;
		    }
    
                    
                    if (numDelRows > 0) {
                        delRow = new int [numDelRows];
                        for (k = 0; k < numDelRows; ++k) {
                            delRow[k] = delIndices[k];
#ifdef BLIS_DEBUG
			    std::cout << "REMOVE: slack row " << delRow[k] 
				      << std::endl;
#endif
                        }
                        
                        //----------------------------------
                        // Delete from lp solver.
                        //----------------------------------
                        
                        model->solver()->deleteRows(numDelRows, delRow);
                        
                        delete [] delRow;
                        delRow = NULL;
                        delIndices.clear();

                        //----------------------------------
                        // Delete from the old cut position array.
                        //----------------------------------

                        tempNumCons = 0;
                        for (k = 0; k < currNumOldCons; ++k) {
                            if (oldDelMark[k] != 1) {
                                // Survived
                                oldConsPos[tempNumCons++] = oldConsPos[k];    
                            }
                        }
                        currNumOldCons = tempNumCons;
                        numStartRows = numCoreRows + currNumOldCons;
                        
                        //----------------------------------
                        // Delete from new cut vector.
                        //----------------------------------

                        //std::cout << std::endl;
                        tempNumCons = 0;
                        for (k = 0; k < newNumCons; ++k) {
                            if (newDelMark[k] == 1) {
                                // Deleted
#ifdef BLIS_DEBUG_MORE
                                std::cout << "delete cut " << k 
                                          << ", size=" 
                                          << dynamic_cast<BlisConstraint*>(newConstraints[k])->getSize()
                                          << std::endl;
#endif
                                
                                delete newConstraints[k];
                                newConstraints[k] = NULL;
                            }
                            else {
                                // Survived
                                newConstraints[tempNumCons++] = newConstraints[k];
                            }
                        }
                        //assert(tempNumCons + numDelRows == newNumCons);
                        numAppliedCons -= newNumCons;
                        numAppliedCons += tempNumCons;
                        newNumCons = tempNumCons;                        
                        
                        //----------------------------------
                        // Resolve to update solution info in lp solver.
                        //----------------------------------
                        
                        int easy = 2;
                        model->solver()->setHintParam(OsiDoInBranchAndCut,
                                                      true, OsiHintDo, &easy);
                        model->solver()->resolve();
                        model->solver()->setHintParam(OsiDoInBranchAndCut,
                                                      true, OsiHintDo, NULL) ;
#ifdef BLIS_DEBUG
                        if (model->solver()->getIterationCount() != 0) {
                            // TODO: maybe some cuts become slack again
#ifdef BLIS_DEBUG
                            std::cout << "SLACK: resolve changed solution!"
                                      << ", iter=" 
				      << model->solver()->getIterationCount()
				      << std::endl;
#endif
                        }
			else {
#ifdef BLIS_DEBUG
			    std::cout<<"SLACK: resolve don't changed solution!"
                                     << std::endl;
#endif
			}
			
                        double tempOV = model->solver()->getObjValue();
                        double ovDiff = fabs(quality_ - tempOV);
                        
                        if (ovDiff /(1.0 + tempOV) > 1.0e-3) {
                            std::cout << "ERROR: SLACK: quality_("<<quality_ 
                                      << ") != tempOV(" << tempOV
                                      << ")" << std::endl;
                            assert(0);
                        }
                        else {
                            std::cout << "AFTER SLACK: quality_("<<quality_ 
                                      << ") == tempOV(" << tempOV
                                      << ")" << std::endl;
                        }
#endif
                    }
                    
                    delete ws;
                    delete [] newDelMark;
                    delete [] oldDelMark;
                }
            }
            
            break;
        case BlisLpStatusAbandoned:
            assert(0);
	    returnStatus = BlisReturnStatusErrLp;
            goto TERM_PROCESS;
        case BlisLpStatusDualInfeasible:
            // FIXME: maybe also primal infeasible
#ifdef BLIS_DEBUG
	    assert(0);
#endif
	    returnStatus = BlisReturnStatusUnbounded;
            goto TERM_PROCESS;
        case BlisLpStatusPrimalInfeasible:
            setStatus(AlpsNodeStatusFathomed);
            quality_ = -ALPS_OBJ_MAX;       // Remove it as soon as possilbe
	    returnStatus = BlisReturnStatusInfeasible;
            goto TERM_PROCESS;
        case BlisLpStatusDualObjLim:
            setStatus(AlpsNodeStatusFathomed);
            quality_ = -ALPS_OBJ_MAX;       // Remove it as soon as possilbe
	    returnStatus = BlisReturnStatusOverObjLim;
            goto TERM_PROCESS;
        case BlisLpStatusPrimalObjLim:
        case BlisLpStatusIterLim:
            /* Can't say much, need branch */
            needBranch = true;
#ifdef BLIS_DEBUG
            assert(0);
#endif
	    returnStatus = BlisReturnStatusBranch;
            goto TERM_BRANCH;
            break;
        default:
#ifdef BLIS_DEBUG
            std::cout << "PROCESS: unknown status "  <<  lpStatus << std::endl;
            assert(0);
#endif
            break;
        }

        //--------------------------------------------------
        // Call heuristics.
        //--------------------------------------------------

#if 0
        std::cout << "Process " << getKnowledgeBroker()->getProcRank()
                  << " : model->boundingPass_=" << model->boundingPass_
                  << " ; maxPass = " << maxPass
                  << " ; keepOn = " << keepOn << std::endl;
#endif

        if (keepOn) {
	    int heurStatus = callHeuristics(model, false);
	    if (heurStatus == 1) {
		cutoff = model->getCutoff();
	    }
	    else if (heurStatus == 2) {
		// Fathom this node
		goto TERM_PROCESS;
	    }
        }
        
        //--------------------------------------------------
        // Generate constraints.
        //--------------------------------------------------

#if 0
        std::cout << "keepOn = " << keepOn
                  << ", geneConsHere = " << genConsHere
                  << ", numAppliedCons = " << numAppliedCons
                  << ", maxNumCons = " << maxNumCons
                  << std::endl;
#endif

        if ( keepOn && genConsHere && 
             (numAppliedCons < maxNumCons) && 
             (model->boundingPass_ < maxPass) ) {
            
            //OsiCuts newOsiCuts;
            BcpsConstraintPool newConPool;

            memcpy(currLpSolution, 
                   model->getLpSolution(),
                   numCols * sizeof(double));
            
	    // Get violated constraints that are from other processes.
            tempNumCons = newConPool.getNumConstraints();
	    getViolatedConstraints(model, currLpSolution, 
				   *(model->constraintPoolReceive()));
	    voilatedNumCons = newConPool.getNumConstraints() - tempNumCons;
	    
	    // Generate constraints (only if no violated).
	    if (voilatedNumCons == 0) {
		lpStatus = static_cast<BlisLpStatus> 
		    (generateConstraints(model, newConPool));
            
		if (lpStatus != BlisLpStatusOptimal) {
		    setStatus(AlpsNodeStatusFathomed);
		    quality_ = -ALPS_OBJ_MAX; // Remove it as soon as possilbe
		    goto TERM_PROCESS;
		}
	    }
            
	    tempNumCons = newConPool.getNumConstraints();
            
            if (tempNumCons > 0) {
		// Select and install new constraints
		applyConstraints(model, currLpSolution, newConPool);
		
		// Some weak/parallel/dense constraints might be discarded.
		tempNumCons = newConPool.getNumConstraints();
		if (tempNumCons > 0) {
		    keepOn = true;
		}
		else {
		    keepOn = false;
		}

                // Move cuts from pool to array newConstraints.
                for (k = 0; k < tempNumCons; ++k) {
                    aCon = dynamic_cast<BlisConstraint *>
			(newConPool.getConstraint(k));
                    newConstraints[newNumCons++] = aCon;
                    if (newNumCons >= maxNewNumCons) {
                        // No space, need resize
#ifdef BLIS_DEBUG
                        std::cout << "NEWCUT: resize, maxNewNumCons = " 
                                  << maxNewNumCons << std::endl;
#endif
                        maxNewNumCons *= 2;
                        BcpsObject **tempNews = new BcpsObject* [maxNewNumCons];
                        memcpy(tempNews, 
                               newConstraints,
                               newNumCons * sizeof(BcpsObject *));
                        delete [] newConstraints;
                        newConstraints = tempNews;
                    }

		    // Make a copy to send pool if share
		    if (shareCon && (voilatedNumCons == 0)) {
			if (aCon->getValidRegion() == BcpsValidGlobal) {
			    model->constraintPoolSend()->
				addConstraint(new BlisConstraint(*aCon));
			}
#if 0
			std::cout << "+++ Num of send new constraint = " 
				  << model->constraintPoolSend()->getNumConstraints()
				  << std::endl;
#endif
		    }
                }
		
		newConPool.clear();
                numAppliedCons += tempNumCons;
            }
	    else { // Didn't generate any new constraints.
		keepOn = false;
	    }
        }
        else { // Don't allow to generate constraints.
            keepOn = false;
        }
    } // EOF bounding/cutting/heuristics loop
    
    //------------------------------------------------------
    // Select branching object
    //------------------------------------------------------

 //FIXME: Turn this into a function
 TERM_BRANCH:

#ifdef BLIS_DEBUG_MORE
    printf("needBranch = %d\n", needBranch);
#endif

    //FIXME: Isn't needBranch always true here?
    if (needBranch) { 
	
	bStatus = -1;
        
        while (bStatus == -1) { 
            foundSolution = false;
            if(getKnowledgeBroker()->getProcRank() == -1) {
                std::cout << "*** I AM RANK ONE: before choose:bStatus = " 
                          << bStatus << std::endl;
            }
            bStatus = selectBranchObject(model, 
                                         foundSolution, 
                                         numPassesLeft);
            --numPassesLeft;
            
            if (bStatus == -1) { 
                lpFeasible = model->resolve();
                
                //resolved = true ;
#ifdef BLIS_DEBUG_MORE
                printf("Resolve since some col fixed, Obj value %g, numRows %d, cutoff %g\n",
                       model->solver()->getObjValue(),
                       model->solver()->getNumRows(),
                       cutoff);
#endif
                if (lpFeasible) {
		    // Update new quality.
                    quality_ = model->solver()->getObjValue();
                    if (quality_ > cutoff) {
                        bStatus = -2;
                    }
                    // Check if feasible at the other branch due to random LP 
                    ipSol = model->feasibleSolution(numIntInfs, numObjInfs);
                    if (ipSol) {         
                        // IP feasible
                        model->storeSolution(BlisSolutionTypeBounding, ipSol);
                        // Update cutoff
                        cutoff = model->getCutoff();
                        setStatus(AlpsNodeStatusFathomed);
                        goto TERM_PROCESS;
                    }
		    if (msgLevel >= 100){
			printf("Final Bound:  %.3f\n", quality_);
		    }
                }
                else {
                    // Should not happen. No, it will happen when other
                    // branch is ip feasible, and cause this branch to fathom
                    // when resolving. Test enigma.
                    // assert(0);
                    bStatus = -2;
                    setStatus(AlpsNodeStatusFathomed);
                }
            }
            
            if(getKnowledgeBroker()->getProcRank() == -1) {
                std::cout << "*** I AM RANK ONE: bStatus = " << bStatus
                          << std::endl;
            }
        }
        
        assert(bStatus != -1);
        
        //----------------------------------------------------
        // If found a branching object:
        // 1. Record basis
        // 2. Record soft var bound difference. 
        // 3. Record add/del constraints.
        // NOTE: Hard var bound differences have been recorded when branch().
        //       startXXXXX have branching bounds for this node
        //----------------------------------------------------
        
        if (bStatus >= 0) {
            
#ifdef BLIS_DEBUG_MORE
            BlisBranchObjectInt *branchObject =
                dynamic_cast<BlisIntegerBranchObject *>(branchObject_);
            std::cout << "SetPregnant: branchedOn = " 
                      << model->getIntVars()[branchObject->variable()]
                      << std::endl;
#endif
            //--------------------------------------------------
            // Mark as pregnant.
            //--------------------------------------------------

            setStatus(AlpsNodeStatusPregnant);
	    
            //--------------------------------------------------
            // Save basis.
            //--------------------------------------------------

            CoinWarmStartBasis *ws = dynamic_cast<CoinWarmStartBasis*>
                (model->solver()->getWarmStart());
            BlisNodeDesc *desc = dynamic_cast<BlisNodeDesc *>(desc_);
            desc->setBasis(ws);

            //----------------------------------------------
            // Save variable/constraint bound, non-core constraints
	    // and non-core variable.
	    // The sizes of hard variable bound vectors are numCols.
	    // The sizes of soft variable bound vectors are number
	    // of modified.
            //----------------------------------------------

	    int *tempVarLBPos = model->tempVarLBPos();
	    int *tempVarUBPos = model->tempVarUBPos();
	    //int *tempConLBPos = model->tempConLBPos();
	    //int *tempConUBPos = model->tempConUBPos();
	    
	    int numModSoftColLB = 0;
	    int numModSoftColUB = 0;
	    const double *currColLB = model->solver()->getColLower();
	    const double *currColUB = model->solver()->getColUpper();
	    //const double *currRowLB = model->solver()->getRowLower();
	    //const double *currRowUB = model->solver()->getRowUpper();
	    
	    double *startColLB = model->startVarLB();
	    double *startColUB = model->startVarUB();
	    //double *startRowLB = model->startConLB();
	    //double *startRowUB = model->startConUB();

	    //BlisConstraint **tempCons = NULL;


#ifdef BLIS_DEBUG_MORE
	    // Debug survived old constraints.
	    for (k = 0; k < currNumOldCons; ++k) {
		int oldPos = oldConsPos[k];
		BlisConstraint *aCon = model->oldConstraints()[oldPos];
		assert(aCon);
		std::cout << "SAVE: DBG: oldPos=" << oldPos
			  << ", k=" << k << ", len=" << aCon->getSize()
			  << ", node=" << index_ << std::endl;
	    }
#endif

	    //----------------------------------------------
	    // Decide if save explicit decription.
	    //----------------------------------------------
	    
            BlisParams * BlisPar = model->BlisPar();
            int difference = BlisPar->entry(BlisParams::difference);
    
            if (difference == -1) {
                if (depth_ % 30 == 0 || isRoot || (phase == AlpsPhaseRampup)) {
                    explicit_ = 1;
                    //std::cout << "SAVE: node "<< index_ <<" explicitly, "
                    //  << "depth=" << depth_ << std::endl;
                }
                else {
                    explicit_ = 0;
                    //std::cout << "SAVE: node "<< index_ <<" relatively, "
                    //  << "depth=" << depth_ << std::endl;
                }
            }
            else if (difference == 0) {
                explicit_ = 1;
            }
            else {
                explicit_ = 0;
            }
	    
	    //explicit_ = 1;

	    if (explicit_ || (phase == AlpsPhaseRampup) ) {
		// NOTE: full hard bound has been stored. 
		
		int index;
		int numModify = 0;
		int numSoftVarLowers = 0;
		int numSoftVarUppers = 0;
		double value;
		
		double *fVarHardLB = new double [numCols];
		double *fVarHardUB = new double [numCols];
		int *fVarHardLBInd = new int [numCols];
		int *fVarHardUBInd = new int [numCols];
		
		double *fVarSoftLB = new double [numCols];
		double *fVarSoftUB = new double [numCols];
		int *fVarSoftLBInd = new int [numCols];
		int *fVarSoftUBInd = new int [numCols];
		
		for (k = 0; k < numCols; ++k) {
		    fVarSoftLB[k] = ALPS_DBL_MAX;
		    fVarSoftUB[k] = -ALPS_DBL_MAX;
		    fVarHardLB[k] = ALPS_DBL_MAX;
		    fVarHardUB[k] = -ALPS_DBL_MAX;
		    fVarHardLBInd[k] = k;
		    fVarHardUBInd[k] = k;
		}

		//------------------------------------------
		// Build the path to an explicit node.
		//------------------------------------------
		
		AlpsTreeNode *parent = parent_;

		// First push this node since it has branching hard bounds.
		model->leafToRootPath.push_back(this);
		BlisNodeDesc* pathDesc = NULL;
		
		if (phase != AlpsPhaseRampup) {
		    while(parent) {
			model->leafToRootPath.push_back(parent);
			if (parent->getExplicit()) {
			    // Reach an explicit node, then stop.
			    break;
			}
			else {
			    parent = parent->getParent();
			}
		    }
		}

#ifdef BLIS_DEBUG_MORE
		std::cout << "SAVE: EXP: path len = "<<model->leafToRootPath.size()
			  << std::endl;
#endif
		//------------------------------------------
		// Summarize bounds.
		//------------------------------------------
		
		for(j = static_cast<int> (model->leafToRootPath.size() - 1); j > -1; --j) {
		    
		    pathDesc = dynamic_cast<BlisNodeDesc*>
			((model->leafToRootPath.at(j))->getDesc());
		    
		    //--------------------------------------
		    // Full variable hard bounds.
		    //--------------------------------------
		    
		    numModify = pathDesc->getVars()->lbHard.numModify;
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->lbHard.posModify[k];
			value = pathDesc->getVars()->lbHard.entries[k];
			fVarHardLB[index] = value;
		    }
		    
		    numModify = pathDesc->getVars()->ubHard.numModify;
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->ubHard.posModify[k];
			value = pathDesc->getVars()->ubHard.entries[k];
			fVarHardUB[index] = value;
		    }
		    
		    //--------------------------------------
		    // Full variable soft bounds.
		    //--------------------------------------
		    
		    numModify = pathDesc->getVars()->lbSoft.numModify;
#ifdef BLIS_DEBUG_MORE
		    std::cout << "SAVE: EXP: j=" << j << ", numModify soft lb="
			      << numModify << std::endl;
#endif
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->lbSoft.posModify[k];
			value = pathDesc->getVars()->lbSoft.entries[k];
			fVarSoftLB[index] = value;
		    }
		    
		    numModify = pathDesc->getVars()->ubSoft.numModify;
#ifdef BLIS_DEBUG_MORE
		    std::cout << "SAVE: EXP: j=" << j << ", numModify soft ub="
			      << numModify << std::endl;
#endif
		    for (k = 0; k < numModify; ++k) {
			index = pathDesc->getVars()->ubSoft.posModify[k];
			value = pathDesc->getVars()->ubSoft.entries[k];
			fVarSoftUB[index] = value;
		    }

		} // EOF of for(path)

		//------------------------------------------
		// Collect modified soft bounds at this node.
		// NOTE: Do this after collecting previous soft bounds.
		//------------------------------------------
		
		numModSoftColLB = 0;
		numModSoftColUB = 0;
		for (k = 0; k < numCoreCols; ++k) {
		    if (currColLB[k] != startColLB[k]) {
			fVarSoftLB[k] = currColLB[k];
			++numModSoftColLB;
#ifdef BLIS_DEBUG_MORE
			printf("Col %d, soft lb change, start %g, curr %g\n",
			       k, startColLB[k], currColLB[k]);
#endif
			
		    }
		    if (currColUB[k] != startColUB[k]) {
			fVarSoftUB[k] = currColUB[k];
			++numModSoftColUB;
		    }
		}
		
#ifdef BLIS_DEBUG_MORE
		std::cout << "SAVE: EXP: THIS: numModSoftColLB = "<<numModSoftColLB 
			  << ", numModSoftColUB = " << numModSoftColUB << std::endl;
#endif
		
		//--------------------------------------
		// Debug if bounds are consistant.
		//--------------------------------------

#ifdef BLIS_DEBUG
		for (k = 0; k < numCols; ++k) {
		    
		    //std::cout << "EXP: COL[" << k <<"]: " 
		    //      <<"hardLB=" << fVarHardLB[k] 
		    //      <<", hardUB=" << fVarHardUB[k] << std::endl;
		    
		    // Hard lower bound should not greater than
		    // hard upper bound.
		    if (fVarHardLB[k] > fVarHardUB[k] + ALPS_GEN_TOL) {
			printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
			       k, fVarHardLB[k], fVarHardUB[k]);
			assert(0);
		    }
		    
		    if (fVarSoftLB[k] < ALPS_BND_MAX) {
			// Soft lower changed, and should not greater 
			// than its hard upper bound.
			if (fVarSoftLB[k] > fVarHardUB[k] + ALPS_GEN_TOL) {
			    printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
				   k, fVarSoftLB[k], fVarHardUB[k]);
			    assert(0);	
			}
		    }
		    
		    if (fVarSoftUB[k] > -ALPS_BND_MAX) {
			// Soft upper changed, and should not less 
			// than its hard lower bound.
			if (fVarSoftUB[k] < fVarHardLB[k] - ALPS_GEN_TOL) {
			    printf("ERROR: Col %d, \tlb %g,  \tub %g\n",
				   k, fVarHardLB[k], fVarSoftUB[k]);
			    assert(0);
			}
		    }
		}
#endif
		
		//------------------------------------------
		// Record hard variable bounds. FULL set.
		//------------------------------------------
		
		desc->assignVarHardBound(numCols,
					 fVarHardLBInd,
					 fVarHardLB,
					 numCols,
					 fVarHardUBInd,
					 fVarHardUB);
		
		//------------------------------------------
		// Recode soft variable bound. Modified.
		//------------------------------------------
		
		for (k = 0; k < numCols; ++k) {
		    if (fVarSoftLB[k] < ALPS_BND_MAX) {
			fVarSoftLBInd[numSoftVarLowers] = k;
			fVarSoftLB[numSoftVarLowers++] = fVarSoftLB[k];
		    }
		    if (fVarSoftUB[k] > -ALPS_BND_MAX) {
			fVarSoftUBInd[numSoftVarUppers] = k;
			fVarSoftUB[numSoftVarUppers++] = fVarSoftUB[k];
		    }
		}

		
#ifdef BLIS_DEBUG_MORE
		// Print soft bounds.
		std::cout << "SAVE: EXP: numSoftVarLowers=" << numSoftVarLowers
			  << ", numSoftVarUppers=" << numSoftVarUppers
			  << std::endl;
		for (k = 0; k < numSoftVarLowers; ++k) {
		    std::cout << "Col[" << fVarSoftLBInd[k] << "]: soft lb="
			      << fVarSoftLB[k] << std::endl;		    
		}
		std::cout << "------------------" << std::endl;
		for (k = 0; k < numSoftVarUppers; ++k) {
		    std::cout << "Col[" << fVarSoftUBInd[k] << "]: soft ub="
			      << fVarSoftUB[k] << std::endl;
		}
		std::cout << "------------------" << std::endl << std::endl;
#endif

		//if ( (numSoftVarUppers > 0) || (numSoftVarLowers > 0) ) {
                
                // Assign it anyway so to delete memory(fVarSoftLBInd,etc.)
                desc->assignVarSoftBound(numSoftVarLowers, 
                                         fVarSoftLBInd,
                                         fVarSoftLB,
                                         numSoftVarUppers, 
                                         fVarSoftUBInd,
                                         fVarSoftUB);

                //------------------------------------------
                // Full set of active non-core constraints.
                //------------------------------------------
                // old constraint: model->oldConstraints_, currNumOldCons.
                // new constraint: newConstraints, newNumCons.

                BcpsObject **toAddCons = new BcpsObject * [currNumOldCons + 
                                                           newNumCons];
                if (currNumOldCons > 0) {
                    // Hard copy of the survived old constraints.
                    for (k = 0; k < currNumOldCons; ++k) {
			int oldPos = oldConsPos[k];
                        BlisConstraint *aCon = model->oldConstraints()[oldPos];
                        assert(aCon);
#ifdef BLIS_DEBUG
			std::cout << "SAVE: EXP: currNumOldCons=" << currNumOldCons 
				  << ", k=" << k << ", len=" << aCon->getSize()
                                  << ", node=" << index_ << std::endl;
#endif
			
                        BlisConstraint *newCon = new BlisConstraint(*aCon);
                        toAddCons[k] = newCon;
                    }
                }
                
		if (newNumCons > 0) {
		    // Include new constraints.
		    memcpy(toAddCons + currNumOldCons, 
			   newConstraints,
			   newNumCons * sizeof(BcpsObject *));
		}

		//------------------------------------------
                // Save in description. Add first delete exiting, then add.
		// Pointers in model->oldConstraints_ can be dangling.
		// It is not safe to use model->oldConstraints_ after adding.
		
		// If this node is the root of a subtree, and before processing
		// it has a list of cuts, then model->oldConstraints_ 
		// stores pointer to the cuts when installing.
		
		// Need update model->constraints_ here OR do not be smart!
		// 1/6/06: Choose to be dumn.
		//------------------------------------------

		//------------------------------------------
		// Generating constraints,
		// also means that slack ones might been removed.
		//------------------------------------------

		int numTotal = currNumOldCons + newNumCons;
                desc->setAddedConstraints(numTotal, toAddCons);
                
#ifdef BLIS_DEBUG
		std::cout << "SAVE: EXP: currNumOldCons=" << currNumOldCons
			  << ", newNumCons=" << newNumCons
			  << std::endl;
#endif	
		
                //------------------------------------------
                // Full set of active non-core variables.
                //------------------------------------------
                // Todo.
                
		//------------------------------------------
		// Clear path vector.
		//------------------------------------------
		
		model->leafToRootPath.clear();
		assert(model->leafToRootPath.size() == 0);
	    }
	    else { // Relative.
		//------------------------------------------
		// Record soft bound changes for core vars.
		//------------------------------------------

		// Variable bound change 
		numModSoftColLB = 0;
		numModSoftColUB = 0;
		for (k = 0; k < numCoreCols; ++k) {
		    if (currColLB[k] != startColLB[k]) {
			tempVarLBPos[numModSoftColLB] = k;
			/* startColLB as a temporary storage vector */
			startColLB[numModSoftColLB] = currColLB[k];
			
#ifdef BLIS_DEBUG_MORE
			printf("Col %d, soft lb change, start %g, curr %g\n",
			       k, startColLB[k], currColLB[k]);
#endif
			
			++numModSoftColLB;
		    }
		    if (currColUB[k] != startColUB[k]) {
			tempVarUBPos[numModSoftColUB] = k;
			startColUB[numModSoftColUB] = currColUB[k];
			++numModSoftColUB;
		    }
		}

#ifdef BLIS_DEBUG_MORE
		std::cout << "SAVE: REL: numModSoftColLB = " 
                          << numModSoftColLB 
			  << ", numModSoftColUB = " 
                          << numModSoftColUB 
			  << std::endl;
#endif
            
		if (numModSoftColLB > 0 || numModSoftColUB > 0) {
#ifdef BLIS_DEBUG
		    //assert(0);
#endif
		    desc->setVarSoftBound(numModSoftColLB, 
					  tempVarLBPos,
					  startColLB,
					  numModSoftColUB, 
					  tempVarUBPos,
					  startColUB);
		}
            
		//------------------------------------------
		// TODO: Constraint bounds change.
		//------------------------------------------

#if 0
		for (k = 0; k < numCoreRows; ++k) {
		    if (currRowLB[k] != startRowLB[k]) {
			tempConLBPos[numModSoftRowLB] = k;
			startRowLB[numModSoftRowLB] = currRowLB[k];
			++numModSoftRowLB;
		    }
		    if (currRowUB[k] != startRowUB[k]) {
			tempConUBPos[numModSoftRowUB] = k;
			startRowUB[numModSoftRowUB] = currRowUB[k];
			++numModSoftRowUB;
		    }
		}
		if (numModSoftRowLB > 0 || numModSoftRowUB > 0) {
		    desc->setConSoftBound(numModSoftRowLB, 
					  tempConLBPos,
					  startRowLB,
					  numModSoftRowUB, 
					  tempConUBPos,
					  startRowUB);
		}
#endif	
            
		if (genConsHere) {
		    // NOTE: explicit_ can NOT do this if, since genConsHere maybe
		    //       false here, but there are maybe cons from parents.

		    //--------------------------------------
		    // Record add constraints.
		    //--------------------------------------

		    // Apend will copy old, then add new. 
		    // If this node has a list of cuts before pointers in 
		    // model->oldConstraints() will be kept. Safe!
		    if (newNumCons > 0) {
			desc->appendAddedConstraints(newNumCons,
                                                     newConstraints);
                    }
		    
		    //--------------------------------------
		    // Record deleted constraint positions.
		    //--------------------------------------
		    
		    int *oldLeft = new int [origNumOldCons];
		    int leftCon;
		    CoinZeroN(oldLeft, origNumOldCons);
		    
		    for (k = 0; k < currNumOldCons; ++k) {
			leftCon = oldConsPos[k];
			assert(leftCon >= 0 && leftCon < origNumOldCons);
			oldLeft[leftCon] = 1;
		    }
		    // 
		    leftCon = 0;
		    for (k = 0; k < origNumOldCons; ++k) {
			if (oldLeft[k] == 0) {
			    // Deleted. Now oldLeft stores delete position.
			    oldLeft[leftCon++] = k;
			}
                        //FIXME: clean
                        //assert(k < 15196924);
		    }
		    desc->delConstraints(leftCon, oldLeft);
		    
#ifdef BLIS_DEBUG
		    std::cout << "PROCESS: ADD: new cuts=" << newNumCons 
			      << ", numRows=" << model->solver()->getNumRows() 
			      << ", numStartRows="<< numStartRows 
                              << ", origNumStartRows="<< origNumStartRows 
                              << ", num removed=" << leftCon << std::endl;                    
#endif
		    
		}// EOF of if(genConsHere)
	    } // EOF of relative
        }
        else if (bStatus == -2) {
            
#if 0
            std::cout << "bStatus = -2, fathom this node!" << std::endl;
#endif
            //branchObject->getDown()[0], branchObject->getDown()[1]);
	    
            setStatus(AlpsNodeStatusFathomed);
        }
        else {
            throw CoinError("No branch object found", "process", 
                            "BlisTreeNode");
        }
    }
    
    //------------------------------------------------------
    // End of process()
    //------------------------------------------------------
    
 TERM_PROCESS:

    bool printCutStat = false;
    if (genConsHere) {
        if ( ((getKnowledgeBroker()->getProcType()==AlpsProcessTypeMaster)||
              (getKnowledgeBroker()->getProcType()==AlpsProcessTypeSerial)) &&
             (msgLevel >= 10) ) {
            printCutStat = true;
            //std::cout << "+++++master print cut stats"<< std::endl;
        }
        else if ( (getKnowledgeBroker()->getProcType()==AlpsProcessTypeHub) &&
                  (hubMsgLevel > 0) ) {
            printCutStat = true;
        }
        else if ((getKnowledgeBroker()->getProcType()==AlpsProcessTypeWorker)&&
		 (workerMsgLevel > 0)) {
            printCutStat = true;
        }
    }

    if (printCutStat == true) {
        printCutStat = false;
        if ( (msgLevel >= 200) || (index_ == 0 && msgLevel >= 50) ) {
            printCutStat = true;
        }
    }

    if (printCutStat) {
        int numT = model->numCutGenerators();
        for (k = 0; k < numT; ++k) {
            if ( model->cutGenerators(k)->calls() > -1) {
                model->blisMessageHandler()->message(BLIS_CUT_STAT_NODE,
                                                     model->blisMessages())
                    << index_
                    << model->cutGenerators(k)->name()
                    << model->cutGenerators(k)->calls()
                    << model->cutGenerators(k)->numConsGenerated()
                    << model->cutGenerators(k)->time()
                    << model->cutGenerators(k)->strategy()
                    << CoinMessageEol;
            }   
        }

        numT = model->numHeuristics();
        for (k = 0; k < numT; ++k) {
            if ( model->heuristics(k)->calls() > -1) {
                model->blisMessageHandler()->message(BLIS_HEUR_STAT_NODE,
                                                     model->blisMessages())
                    << index_
                    << model->heuristics(k)->name()
                    << model->heuristics(k)->calls()
                    << model->heuristics(k)->numSolutions()
                    << model->heuristics(k)->time()
                    << model->heuristics(k)->strategy()
                    << CoinMessageEol;
            }   
        }
    }

    delete [] currLpSolution;

    if (msgLevel >= 100){
	printf("Final Bound:  %.3f\n", model->solver()->getObjValue());
    }

    if (status_ == AlpsNodeStatusFathomed) {
	// Delete new cuts since no use anymore.
        for (k = 0; k < newNumCons; ++k) {
            delete newConstraints[k];
        }
    }
    delete [] newConstraints;
    delete [] oldConsPos;
    
    model->isRoot_ = false;

#ifdef BLIS_DEBUG_MORE
    // Debug survived old constraints.
    //int currNumOldCons = model->getNumOldConstraints();
    for (k = 0; k < currNumOldCons; ++k) {
	BlisConstraint *aCon = model->oldConstraints()[k];
	assert(aCon);
	std::cout << "SAVE: DBG: TERM: "
		  << "k=" << k << ", len=" << aCon->getSize()
		  << ", node=" << index_ << std::endl;
    }
#endif
    
    return returnStatus;
}

//#############################################################################

std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> >
BlisTreeNode::branch()
{

    //------------------------------------------------------
    // Change one var hard bound and record the change in nodedesc:
    // THINK: how about constraint bounds? When to update?
    // TODO: how about other SOS object, etc.?
    //------------------------------------------------------

    AlpsPhase phase = knowledgeBroker_->getPhase();

    double objVal = quality_;

    BlisNodeDesc* childDesc = NULL;  
      
    std::vector< CoinTriple<AlpsNodeDesc*, AlpsNodeStatus, double> > 
	childNodeDescs;
    
    BlisModel* model = dynamic_cast<BlisModel*>(desc_->getModel());    

    int numCols = model->getNumCols();

#ifdef BLIS_DEBUG_MORE
    // Debug survived old constraints.
    int currNumOldCons = model->getNumOldConstraints();
    for (int k = 0; k < currNumOldCons; ++k) {
	BlisConstraint *aCon = model->oldConstraints()[k];
	assert(aCon);
	std::cout << "BRANCH: DBG: "
		  << "k=" << k << ", len=" << aCon->getSize()
		  << ", node=" << index_ << std::endl;
    }
#endif

    //------------------------------------------------------
    // Get branching object. TODO: Assume integer branching object. 
    //------------------------------------------------------
    
    BlisNodeDesc* thisDesc = dynamic_cast<BlisNodeDesc*>(desc_);

    BlisBranchingObjectType type = 
	static_cast<BlisBranchingObjectType>(branchObject_->getType());

    switch(type){

     case BlisBranchingObjectTypeInt: 
	 {
	     BlisBranchObjectInt *branchObject =
		 dynamic_cast<BlisBranchObjectInt *>(branchObject_);

	     int objInd = branchObject->getObjectIndex();
	     
	     double bValue = branchObject->getValue();
	     
	     BlisObjectInt *obj = 
		 dynamic_cast<BlisObjectInt *>(model->objects(objInd));
	     int branchVar = obj->columnIndex();
	     
#ifdef BLIS_DEBUG
	     if ( (branchVar < 0) || (branchVar >= numCols) ) {
		 std::cout << "ERROR: BRANCH(): branchVar = " << branchVar 
			   << "; numCols = " << numCols  << std::endl;
		 throw CoinError("branch index is out of range", 
				 "branch", "BlisTreeNode");
	     }
#endif
	     
#ifdef BLIS_DEBUG
	     printf("BRANCH(): on %d, phase %d\n", branchVar, phase);
	     printf("DOWN: lb %g, up %g\n",
		    branchObject->getDown()[0], branchObject->getDown()[1]);
	     printf("UP  : lb %g, up %g\n",
		    branchObject->getUp()[0], branchObject->getUp()[1]);
#endif
	     
	     //======================================================
	     //------------------------------------------------------
	     // Create down-branch node description.
	     //------------------------------------------------------
	     //======================================================

	     childDesc = new BlisNodeDesc(model);
	     
	     if (phase == AlpsPhaseRampup) {
		 
		 //--------------------------------------------------
		 // Store a full description since each node will be the 
		 // root of a subtree.
		 // NOTE: this desc must be explicit during rampup.
		 //--------------------------------------------------	
		 
		 int index, k;
		 int numModify = -1;
		 double value;
		 
		 double *fVarHardLB = new double [numCols];
		 double *fVarHardUB = new double [numCols];
		 int *fVarHardLBInd = new int [numCols];
		 int *fVarHardUBInd = new int [numCols];
		 
		 double *fVarSoftLB = NULL;
		 double *fVarSoftUB = NULL;
		 int *fVarSoftLBInd = NULL;
		 int *fVarSoftUBInd = NULL;
	    
		 //--------------------------------------------------
		 // Full hard variable bounds.
		 //--------------------------------------------------
		 
		 numModify = thisDesc->getVars()->lbHard.numModify;
		 assert(numModify == numCols);
		 for (k = 0; k < numModify; ++k) {
		     index = thisDesc->getVars()->lbHard.posModify[k];
		     assert(index == k);
		     value = thisDesc->getVars()->lbHard.entries[k];
		     fVarHardLB[k] = value;
		     fVarHardLBInd[k] = index;
		 }
		 
		 numModify = thisDesc->getVars()->ubHard.numModify;
		 assert(numModify == numCols);
		 for (k = 0; k < numModify; ++k) {
		     index = thisDesc->getVars()->ubHard.posModify[k];
		     assert(index == k);
		     value = thisDesc->getVars()->ubHard.entries[k];
		     fVarHardUB[k] = value;
		     fVarHardUBInd[k] = index;
		 }
		 
		 // Branching bounds.
		 fVarHardLB[branchVar] = branchObject->getDown()[0];
		 fVarHardUB[branchVar] = branchObject->getDown()[1];
		 
		 
		 childDesc->assignVarHardBound(numCols,
					       fVarHardLBInd,
					       fVarHardLB,
					       numCols,
					       fVarHardUBInd,
					       fVarHardUB);
		 
		 //--------------------------------------------------
		 // Soft variable bounds.
		 //--------------------------------------------------
		 
		 int numSoftVarLowers = thisDesc->getVars()->lbSoft.numModify;
		 assert(numSoftVarLowers >= 0 && numSoftVarLowers <= numCols);
		 if (numSoftVarLowers > 0) {
		     fVarSoftLB = new double [numSoftVarLowers];
		     fVarSoftLBInd = new int [numSoftVarLowers];
		     for (k = 0; k < numSoftVarLowers; ++k) {
			 index = thisDesc->getVars()->lbSoft.posModify[k];
			 value = thisDesc->getVars()->lbSoft.entries[k];
			 fVarSoftLB[k] = value;
			 fVarSoftLBInd[k] = index;
		     }
		 }
		 
		 int numSoftVarUppers = thisDesc->getVars()->ubSoft.numModify;
		 assert(numSoftVarUppers >= 0 && numSoftVarUppers <= numCols);
		 if (numSoftVarUppers > 0) {
		     fVarSoftUB = new double [numSoftVarUppers];
		     fVarSoftUBInd = new int [numSoftVarUppers];
		     for (k = 0; k < numSoftVarUppers; ++k) {
			 index = thisDesc->getVars()->ubSoft.posModify[k];
			 value = thisDesc->getVars()->ubSoft.entries[k];
			 fVarSoftUB[k] = value;
			 fVarSoftUBInd[k] = index;
		     }
		 }
		 
#ifdef BLIS_DEBUG_MORE
		 // Print soft bounds.
		 std::cout << "\nBRANCH: numSoftVarLowers=" << numSoftVarLowers
			   << ", numSoftVarUppers=" << numSoftVarUppers
			   << std::endl;
		 for (k = 0; k < numSoftVarLowers; ++k) {
		     std::cout << "Col[" << fVarSoftLBInd[k] << "]: soft lb="
			       << fVarSoftLB[k] << std::endl;		    
		 }
		 std::cout << "------------------" << std::endl;
		 for (k = 0; k < numSoftVarUppers; ++k) {
		     std::cout << "Col[" << fVarSoftUBInd[k] << "]: soft ub="
			       << fVarSoftUB[k] << std::endl;
		 }
		 std::cout << "------------------" << std::endl << std::endl;
#endif
		 
		 // Assign it anyway so to transfer ownership of 
		 // memory(fVarSoftLBInd,etc.)
		 childDesc->assignVarSoftBound(numSoftVarLowers, 
					       fVarSoftLBInd,
					       fVarSoftLB,
					       numSoftVarUppers, 
					       fVarSoftUBInd,
					       fVarSoftUB);
		 
		 //--------------------------------------------------
		 // Full set of non-core constraints.
		 // NOTE: non-core constraints have been saved in description
		 //       when process() during ramp-up.
		 //--------------------------------------------------
		 
		 BcpsObject **tempCons = NULL;
		 int tempInt = 0;
		 
		 tempInt = thisDesc->getCons()->numAdd;
		 if (tempInt > 0) {
		     tempCons = new BcpsObject* [tempInt];
		     for (k = 0; k < tempInt; ++k) {
			 BlisConstraint *aCon = dynamic_cast<BlisConstraint *>
			     (thisDesc->getCons()->objects[k]);
			 
			 assert(aCon);
			 assert(aCon->getSize() > 0);
			 assert(aCon->getSize() < numCols);
			 BlisConstraint *newCon = new BlisConstraint(*aCon);
			 tempCons[k] = newCon;
		     }
		 }
		 
		 
#if 0
		 else {
		     // No cons or only root cons.
		     tempInt = model->getNumOldConstraints();
		     
		     if (tempInt > 0) {
			 tempCons = new BcpsObject* [tempInt];
		     }
		     for (k = 0; k < tempInt; ++k) {
			 BlisConstraint *aCon = model->oldConstraints()[k];                
			 assert(aCon);
			 assert(aCon->getSize() > 0);
			 assert(aCon->getSize() < numCols);
			 BlisConstraint *newCon = new BlisConstraint(*aCon);
			 tempCons[k] = newCon;
		     }
		 }
#endif
		 
#ifdef BLIS_DEBUG_MORE
		 std::cout << "BRANCH: down: tempInt=" << tempInt <<std::endl;
#endif
		 // Fresh desc, safely add.
		 childDesc->setAddedConstraints(tempInt, tempCons);
	     }else{
		 
		 //--------------------------------------------------
		 // Relative: Only need to record hard var bound change. 
		 // NOTE: soft var bound changes are Record after 
		 // selectBranchObject.
		 //--------------------------------------------------
		 
		 childDesc->setVarHardBound(1,
					    &branchVar,
					    &(branchObject->getDown()[0]),
					    1,
					    &branchVar,
					    &(branchObject->getDown()[1]));
	     }
	     
	     childDesc->setBranchedDir(-1);
	     childDesc->setBranchedInd(objInd);
	     childDesc->setBranchedVal(bValue);
	     
	     // Copy warm start.
	     CoinWarmStartBasis *ws = thisDesc->getBasis();
	     CoinWarmStartBasis *newWs = new CoinWarmStartBasis(*ws);
	     childDesc->setBasis(newWs);
	     
	     childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
						     (childDesc),
						     AlpsNodeStatusCandidate,
						     objVal));
	     
	     //======================================================
	     //------------------------------------------------------
	     // Create up-branch node description.
	     //------------------------------------------------------
	     //======================================================
	     
	     childDesc = new BlisNodeDesc(model);
	     
	     if (phase == AlpsPhaseRampup) {
		 
		 //--------------------------------------------------
		 // Store a full description since each node will be the 
		 // root of a subtree.
		 // NOTE: parent must be explicit during rampup.
		 //--------------------------------------------------	
		 
		 int index, k;
		 int numModify = -1;
		 double value;
		 
		 double *fVarHardLB = new double [numCols];
		 double *fVarHardUB = new double [numCols];
		 int *fVarHardLBInd = new int [numCols];
		 int *fVarHardUBInd = new int [numCols];
		 
		 double *fVarSoftLB = NULL;
		 double *fVarSoftUB = NULL;
		 int *fVarSoftLBInd = NULL;
		 int *fVarSoftUBInd = NULL;
		 
		 //--------------------------------------------------
		 // Full hard variable bounds.
		 //--------------------------------------------------
		 
		 numModify = thisDesc->getVars()->lbHard.numModify;
		 assert(numModify == numCols);
		 for (k = 0; k < numModify; ++k) {
		     index = thisDesc->getVars()->lbHard.posModify[k];
		     assert(index == k);
		     value = thisDesc->getVars()->lbHard.entries[k];
		     fVarHardLB[k] = value;
		     fVarHardLBInd[k] = index;
		 }
		 
		 numModify = thisDesc->getVars()->ubHard.numModify;
		 assert(numModify == numCols);
		 for (k = 0; k < numModify; ++k) {
		     index = thisDesc->getVars()->ubHard.posModify[k];
		     assert(index == k);
		     value = thisDesc->getVars()->ubHard.entries[k];
		     fVarHardUB[k] = value;
		     fVarHardUBInd[k] = index;
		 }
		 
		 // Branching bounds.
		 fVarHardLB[branchVar] = branchObject->getUp()[0];
		 fVarHardUB[branchVar] = branchObject->getUp()[1];
		 
		 childDesc->assignVarHardBound(numCols,
					       fVarHardLBInd,
					       fVarHardLB,
					       numCols,
					       fVarHardUBInd,
					       fVarHardUB);
		 
		 //--------------------------------------------------
		 // Soft variable bounds.
		 //--------------------------------------------------
		 
		 int numSoftVarLowers = thisDesc->getVars()->lbSoft.numModify;
		 assert(numSoftVarLowers >= 0 && numSoftVarLowers <= numCols);
		 if (numSoftVarLowers > 0) {
		     fVarSoftLB = new double [numSoftVarLowers];
		     fVarSoftLBInd = new int [numSoftVarLowers];
		     for (k = 0; k < numSoftVarLowers; ++k) {
			 index = thisDesc->getVars()->lbSoft.posModify[k];
			 value = thisDesc->getVars()->lbSoft.entries[k];
			 fVarSoftLB[k] = value;
			 fVarSoftLBInd[k] = index;
		     }
		 }
		 
		 int numSoftVarUppers = thisDesc->getVars()->ubSoft.numModify;
		 assert(numSoftVarUppers >= 0 && numSoftVarUppers <= numCols);
		 if (numSoftVarUppers > 0) {
		     fVarSoftUB = new double [numSoftVarUppers];
		     fVarSoftUBInd = new int [numSoftVarUppers];
		     for (k = 0; k < numSoftVarUppers; ++k) {
			 index = thisDesc->getVars()->ubSoft.posModify[k];
			 value = thisDesc->getVars()->ubSoft.entries[k];
			 fVarSoftUB[k] = value;
			 fVarSoftUBInd[k] = index;
		     }
		 }
		 
		 // Assign it anyway so to transfer ownership of 
		 // memory(fVarSoftLBInd,etc.)
		 childDesc->assignVarSoftBound(numSoftVarLowers, 
					       fVarSoftLBInd,
					       fVarSoftLB,
					       numSoftVarUppers, 
					       fVarSoftUBInd,
					       fVarSoftUB);
		 
		 //--------------------------------------------------
		 // Full set of non-core constraints.
		 // NOTE: non-core constraints have been saved in description
		 //       when process() during ramp-up.
		 //--------------------------------------------------
		 
		 BcpsObject **tempCons = NULL;
		 int tempInt = 0;
		 
		 tempInt = thisDesc->getCons()->numAdd;
		 if (tempInt > 0) {
		     tempCons = new BcpsObject* [tempInt];
		     
		     for (k = 0; k < tempInt; ++k) {
			 BlisConstraint *aCon = dynamic_cast<BlisConstraint *>
			     (thisDesc->getCons()->objects[k]);
			 
			 assert(aCon);
			 assert(aCon->getSize() > 0);
			 assert(aCon->getSize() <= numCols);
			 BlisConstraint *newCon = new BlisConstraint(*aCon);
			 tempCons[k] = newCon;
		     }
		 }
		 
#if 0
		 else {
		     // No cons or only root cons.
		     tempInt = model->getNumOldConstraints();
		     if (tempInt > 0) {
			 tempCons = new BcpsObject* [tempInt];
		     }
		     for (k = 0; k < tempInt; ++k) {
			 BlisConstraint *aCon = model->oldConstraints()[k];                
			 assert(aCon);
			 assert(aCon->getSize() > 0);
			 assert(aCon->getSize() <= numCols);
			 BlisConstraint *newCon = new BlisConstraint(*aCon);
			 tempCons[k] = newCon;
		     }
		 }
#endif
		 
#ifdef BLIS_DEBUG_MORE
		 std::cout << "BRANCH: up: tempInt=" << tempInt <<std::endl;
#endif
		 // Fresh desc, safely add.
		 childDesc->setAddedConstraints(tempInt, tempCons);
	     }else{
		 
		 //--------------------------------------------------
		 // Relative: Only need to record hard var bound change. 
		 // NOTE: soft var bound changes are Record after selectBranchObject.
		 //--------------------------------------------------
		 
		 childDesc->setVarHardBound(1,
					    &branchVar,
					    &(branchObject->getUp()[0]),
					    1,
					    &branchVar,
					    &(branchObject->getUp()[1]));
	     }
	     
	     childDesc->setBranchedDir(1);
	     childDesc->setBranchedInd(objInd);
	     childDesc->setBranchedVal(bValue);
	     
	     // Copy warm start.
	     CoinWarmStartBasis *newWs2 = new CoinWarmStartBasis(*ws);
	     childDesc->setBasis(newWs2);
	     childNodeDescs.push_back(CoinMakeTriple(static_cast<AlpsNodeDesc *>
						     (childDesc),
						     AlpsNodeStatusCandidate,
						     objVal));  
	     
	     // Change node status to branched.
	     status_ = AlpsNodeStatusBranched; 
	 }
	 
	 break;

     case BlisBranchingObjectTypeBilevel :
	 {
	     BlisBranchObjectBilevel *branchObject =
		 dynamic_cast<BlisBranchObjectBilevel *>(branchObject_);
	     std::deque<int> *branchingSet = branchObject->getBranchingSet();

	     CoinWarmStartBasis *ws = thisDesc->getBasis();
	     std::deque<int>::iterator ptr1, ptr2;
	     int size = static_cast<int>(branchingSet->size());
	     int *indices = new int[size];
	     double *values = new double[size];
	     values[0] = 1;
	     int i;
	     for (i = 0, ptr1 = branchingSet->begin(); 
		  ptr1 != branchingSet->end(); i++, ptr1++){
		 indices[i] = *ptr1;
		 values[i] = 1;
		 childDesc = new BlisNodeDesc(model);
		 childDesc->setVarHardBound(i+1, indices, values,
					    i+1, indices, values);
		 CoinWarmStartBasis *newWs = new CoinWarmStartBasis(*ws);
		 childDesc->setBasis(newWs);
		 childNodeDescs.push_back(CoinMakeTriple(
                                       static_cast<AlpsNodeDesc *>(childDesc), 
				       AlpsNodeStatusCandidate,
				       quality_));
		 values[i] = 0;
	     }
	     delete[] indices;
	     delete[] values;
	 }

	 break;
	 
     default :
	
	std::cout << "Unknown branching object type" << std::endl;
    }
    return childNodeDescs;
}

//#############################################################################

/* FIXME: need rewrite from scratch */
/* 0: find a branch var, -1 no branch var (should not happen) */

int BlisTreeNode::selectBranchObject(BlisModel *model, 
                                     bool& foundSol, 
                                     int numPassesLeft) 
{
    int bStatus = 0;
    BcpsBranchStrategy *strategy = 0;

    AlpsPhase phase = knowledgeBroker_->getPhase();
        

    if(branchObject_) {
        delete branchObject_;
        branchObject_ = NULL;
    }
    
    //------------------------------------------------------
    // Get branching strategy.
    //------------------------------------------------------

    if (phase == AlpsPhaseRampup) {
      strategy = model->rampUpBranchStrategy();
    }
    else {
      strategy = model->branchStrategy();
    }

    if (!strategy) {
        throw CoinError("No branch strategy.", "process()","BlisTreeNode");
    }

    //------------------------------------------------------
    // Create branching object candidates.
    //-----------------------------------------------------

    bStatus = strategy->createCandBranchObjects(numPassesLeft,
						model->getCutoff());   
    
    //------------------------------------------------------
    // Select the best branching objects.
    //-----------------------------------------------------
    
    if (bStatus >= 0) {
        
	branchObject_ = strategy->bestBranchObject();
        
        if (branchObject_) {
            // Move best branching object to node.

#ifdef BLIS_DEBUG_MORE
            std::cout << "SELECTBEST: Set branching obj" << std::endl;
#endif
        }
        else {
#ifdef BLIS_DEBUG
            std::cout << "ERROR: Can't find branching object" << std::endl;
#endif      
            assert(0);
        }
        
        // Set guessed solution value
        // solEstimate_ = quality_ + sumDeg;
    }
    
    if (!model->branchStrategy()) {
        delete strategy;
    }

    return bStatus;
}

//#############################################################################

int 
BlisTreeNode::bound(BcpsModel *model) 
{
    BlisLpStatus status = BlisLpStatusUnknown;

    BlisModel *m = dynamic_cast<BlisModel *>(model);   
 
    // Bounding
    m->solver()->resolve();

#if 0
    char name[50] = "";
    sprintf(name, "matrix.%i.%i", index_, m->boundingPass_);
    m->solver()->writeMps(name, "mps", 1.0);
#endif

#if 0
    int j;
    int numCols = m->solver()->getNumCols();
    const double * clb = m->solver()->getColLower();
    const double * cub = m->solver()->getColUpper();
    const double * solution = m->solver()->getColSolution();

    for (j = 0; j < numCols; ++j) {
	std::cout << "col"<< j <<": bounds ["<<clb[j]<<", "<< cub[j] << "]" 
		  << ", x = " << solution[j]
		  << std::endl;
    }
#endif

    if (m->solver()->isAbandoned()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is abandoned" << std::endl;
#endif
	status = BlisLpStatusAbandoned;
    }
    else if (m->solver()->isProvenOptimal()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is lp optimal" << std::endl;
#endif
	status = BlisLpStatusOptimal;
        BlisNodeDesc *desc = dynamic_cast<BlisNodeDesc*>(desc_);

        double objValue = m->solver()->getObjValue() *
            m->solver()->getObjSense();
        
        int dir = desc->getBranchedDir();
        if (dir != 0) {
            double objDeg = objValue - quality_;
            int objInd = desc->getBranchedInd();
            double lpX = desc->getBranchedVal();
            BlisObjectInt *intObject = 
                dynamic_cast<BlisObjectInt *>(m->objects(objInd));            
#if 0
            std::cout << "BOUND: col[" << intObject->columnIndex() 
                      << "], dir=" << dir << ", objDeg=" << objDeg
                      << ", x=" << lpX
                      << ", up=" << intObject->pseudocost().getUpCost()
                      << ", down=" << intObject->pseudocost().getDownCost()
                      << ", pre quality=" << quality_
                      << ", objValue=" << objValue
                      << std::endl;
#endif

            intObject->pseudocost().update(dir, objDeg, lpX);
            m->setSharedObjectMark(intObject->getObjectIndex());
        }

        // Update quality of this nodes.
        quality_ = objValue;
    }
    else if (m->solver()->isProvenPrimalInfeasible()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is primal inf" << std::endl;
#endif
	status = BlisLpStatusPrimalInfeasible;
    }
    else if (m->solver()->isProvenDualInfeasible()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is dual inf" << std::endl;
#endif
	status = BlisLpStatusDualInfeasible;
    }
    else if (m->solver()->isPrimalObjectiveLimitReached()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is primal limit" << std::endl;
#endif
	status = BlisLpStatusPrimalObjLim;
    }
    else if (m->solver()->isDualObjectiveLimitReached()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is dual limit" << std::endl;
#endif
	status = BlisLpStatusDualObjLim;
    }
    else if (m->solver()->isIterationLimitReached()) {
#ifdef BLIS_DEBUG
	std::cout << "BOUND: is iter limit" << std::endl;
#endif
	status = BlisLpStatusIterLim;
    }
    else {
	std::cout << "UNKNOWN LP STATUS" << std::endl;
	assert(0);
    }
    
    return status;
}

//#############################################################################

int BlisTreeNode::installSubProblem(BcpsModel *m)
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    int i, k;
    int index;
    double value;

    BlisModel *model = dynamic_cast<BlisModel *>(m);
    assert(model);
    
    BlisNodeDesc *desc = dynamic_cast<BlisNodeDesc*>(desc_);

    int numModify = 0;

    int numCoreVars = model->getNumCoreVariables();
    int numCoreCons = model->getNumCoreConstraints();

    int numCols = model->solver()->getNumCols();
    int numRows = model->solver()->getNumRows();

    //double *varSoftLB = NULL;
    //double *varSoftUB = NULL;
    double *varHardLB = NULL;
    double *varHardUB = NULL;

    //double *conSoftLB = NULL;
    //double *conSoftUB = NULL;
    //double *conHardLB = NULL;
    //double *conHardUB = NULL;
    
    double *startColLB = model->startVarLB();
    double *startColUB = model->startVarUB();
    double *startRowLB = model->startConLB();
    double *startRowUB = model->startConUB();
    
    CoinFillN(startColLB, numCoreVars, -ALPS_DBL_MAX);
    CoinFillN(startColUB, numCoreVars, ALPS_DBL_MAX);
    CoinFillN(startRowLB, numCoreCons, -ALPS_DBL_MAX);
    CoinFillN(startRowUB, numCoreCons, ALPS_DBL_MAX);

    int numOldCons = 0;
    int tempInt = 0;
    BlisConstraint *aCon = NULL;

    int nodeID = -1;
    nodeID = getIndex();
    //std::cout << "nodeID=" << nodeID << std::endl;

    AlpsPhase phase = knowledgeBroker_->getPhase();

    //======================================================
    // Restore subproblem: 
    //  1. Remove noncore var/con
    //  2. Travel back to root and correct differencing to
    //     full var/con bounds into model->startXXX
    //  3. Set col bounds
    //  4. Set row bounds (is this necessary?)
    //  5. Add contraints except cores
    //  6. Add variables except cores
    //  7. Set basis (should not need modify)
    //======================================================


    //------------------------------------------------------
    // Remove old constraints from lp solver.
    //------------------------------------------------------

    int numDelCons = numRows - numCoreCons;
    
#ifdef BLIS_DEBUG
    std::cout << "INSTALL: numDelCons = " << numDelCons << std::endl;
#endif
	
    if (numDelCons > 0) {
	int *indices = new int [numDelCons];
	if (indices == NULL) {
	    throw CoinError("Out of memory", "installSubProblem", "BlisTreeNode");
	}
	
	for (i = 0; i < numDelCons; ++i) {
	    indices[i] = numCoreCons + i;
	}
	
	model->solver()->deleteRows(numDelCons, indices);
	delete [] indices;
	indices = NULL;
    }
    
    //--------------------------------------------------------
    // Travel back to a full node, then collect diff (add/rem col/row,
    // hard/soft col/row bounds) from the node full to this node.
    //----------------------------
    // Note: if we store full set of logic/agorithm col/row, then
    //       no diff are needed for col/row
    //--------------------------------------------------------
    
    //--------------------------------------------------------
    // Collect differencing bounds. Branching bounds of this node
    // are ALSO collected.
    //--------------------------------------------------------

    BlisNodeDesc* pathDesc = NULL;
    AlpsTreeNode *parent = parent_;    
    
    /* First push this node since it has branching hard bounds. 
       NOTE: during rampup, this desc has full description when branch(). */
    model->leafToRootPath.push_back(this);
    
    if (phase != AlpsPhaseRampup) {
	while(parent) {
#ifdef BLIS_DEBUG_MORE
	    std::cout << "Parent id = " << parent->getIndex() << std::endl;
#endif     
	    model->leafToRootPath.push_back(parent);
	    if (parent->getExplicit()) {
		// Reach an explicit node, then stop.
		break;
	    }
	    else {
		parent = parent->getParent();
	    }
	}
    }
    
#ifdef BLIS_DEBUG_MORE
    std::cout << "INSTALL: path len = " << model->leafToRootPath.size()
	      << std::endl;
#endif
    
    //------------------------------------------------------
    // Travel back from this node to the explicit node to
    // collect full description.
    //------------------------------------------------------
    
    for(i = static_cast<int> (model->leafToRootPath.size() - 1); i > -1; --i) {

#ifdef BLIS_DEBUG_MORE
        if (index_ == 3487) {
            std::cout << "\n----------- NODE ------------" 
                      << model->leafToRootPath.at(i)->getIndex() << std::endl;
        }
#endif
	
	//--------------------------------------------------
	// NOTE: As away from explicit node, bounds become 
	//       tighter and tighter.
	//--------------------------------------------------
        
        pathDesc = dynamic_cast<BlisNodeDesc*>((model->leafToRootPath.at(i))->
                                               getDesc());
        
        varHardLB = pathDesc->getVars()->lbHard.entries;
        varHardUB = pathDesc->getVars()->ubHard.entries;
        
	//--------------------------------------------------      
        // Adjust bounds according to hard var lb/ub.
	// If rampup or explicit, collect hard bounds so far.
	//--------------------------------------------------
	
        numModify = pathDesc->getVars()->lbHard.numModify;
	
#ifdef BLIS_DEBUG_MORE
	std::cout << "INSTALL: numModify lb hard = " << numModify << std::endl;
#endif
	
        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->lbHard.posModify[k];
            value = pathDesc->getVars()->lbHard.entries[k];
      
#ifdef BLIS_DEBUG_MORE
	    printf("INSTALL: 1, col %d, value %g, startColLB %x\n", 
		   index, value, startColLB);
#endif      
	    // Hard bounds do NOT change according to soft bounds, so
	    // here need CoinMax.
            startColLB[index] = CoinMax(startColLB[index], value);
            
#ifdef BLIS_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 1, col %d, hard lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
#endif
        }

#ifdef BLIS_DEBUG_MORE
	std::cout << "INSTALL: numModify ub hard = " << numModify<<std::endl;
#endif

        numModify = pathDesc->getVars()->ubHard.numModify;
        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->ubHard.posModify[k];
            value = pathDesc->getVars()->ubHard.entries[k];
            startColUB[index] = CoinMin(startColUB[index], value);
	    
#ifdef BLIS_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 2, col %d, hard lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
            if (startColLB[index] > startColUB[index]) {
                    //assert(0);
            }
#endif
        }
        
        //--------------------------------------------------
        // Adjust bounds according to soft var lb/ub.
	// If rampup or explicit, collect soft bounds so far.
        //--------------------------------------------------

        numModify = pathDesc->getVars()->lbSoft.numModify;
#ifdef BLIS_DEBUG_MORE
	std::cout << "INSTALL: i=" << i << ", numModify soft lb="
		  << numModify << std::endl;
#endif
        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->lbSoft.posModify[k];
            value = pathDesc->getVars()->lbSoft.entries[k];
            startColLB[index] = CoinMax(startColLB[index], value);
            
#ifdef BLIS_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 3, col %d, soft lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
            
            if (startColLB[index] > startColUB[index]) {
		//assert(0);
            }
#endif
        }
        numModify = pathDesc->getVars()->ubSoft.numModify;
        
#ifdef BLIS_DEBUG_MORE
	std::cout << "INSTALL: i=" << i << ", numModify soft ub="
		  << numModify << std::endl;
#endif

        for (k = 0; k < numModify; ++k) {
            index = pathDesc->getVars()->ubSoft.posModify[k];
            value = pathDesc->getVars()->ubSoft.entries[k];
            startColUB[index] = CoinMin(startColUB[index], value);
            
#ifdef BLIS_DEBUG_MORE
            if (index_ == 3487) {
                printf("INSTALL: 4, col %d, soft lb %g, ub %g\n", 
                       index, startColLB[index], startColUB[index]);
            }
            
            if (startColLB[index] > startColUB[index]) {
		//assert(0);
            }
#endif
        }
        
        //--------------------------------------------------
        // TODO: Modify hard/soft row lb/ub.
        //--------------------------------------------------


        //--------------------------------------------------
        // Collect active non-core constraints at parent.
        //--------------------------------------------------

	//----------------------------------------------
	// First collect all generated cuts, then remove 
	// deleted.
	//----------------------------------------------
	
	tempInt = pathDesc->getCons()->numAdd;
	    
#ifdef BLIS_DEBUG_MORE
	std::cout << "\nINSTALL: numAdd = " << tempInt << std::endl;
#endif

	int maxOld = model->getOldConstraintsSize();
	
	for (k = 0; k < tempInt; ++k) {
	    aCon = dynamic_cast<BlisConstraint *>
		(pathDesc->getCons()->objects[k]);
	    
	    assert(aCon);
	    assert(aCon->getSize() > 0);
	    assert(aCon->getSize() < 100000);
	    
#ifdef BLIS_DEBUG_MORE
	    std::cout << "INSTALL: cut  k=" << k 
		      << ", len=" <<aCon->getSize() 
		      << ", node="<< index_ << std::endl;
#endif
	    (model->oldConstraints())[numOldCons++] = aCon;
	    
	    if (numOldCons >= maxOld) {
		// Need resize
#ifdef BLIS_DEBUG_MORE
		std::cout << "INSTALL: resize, maxOld = " 
			  << maxOld << std::endl;
#endif
		maxOld *= 2;
		BlisConstraint **tempCons = new BlisConstraint* [maxOld];
		
		memcpy(tempCons, 
		       model->oldConstraints(), 
		       numOldCons * sizeof(BlisConstraint *));
		
		model->delOldConstraints();
		model->setOldConstraints(tempCons);
		model->setOldConstraintsSize(maxOld);
	    }
	}
	    
	//----------------------------------------------
	// Remove those deleted. 
	// NOTE: model->oldConstraints_ stores all previously 
	// generated active constraints at parent.
	//----------------------------------------------
	
	tempInt = pathDesc->getCons()->numRemove;
            
	if (tempInt > 0) {
	    int tempPos;
	    int *tempMark = new int [numOldCons];
	    CoinZeroN(tempMark, numOldCons);
	    for (k = 0; k < tempInt; ++k) {
		tempPos = pathDesc->getCons()->posRemove[k];
#ifdef BLIS_DEBUG_MORE
		std::cout << "tempPos=" << tempPos 
			  << ", tempInt=" << tempInt 
			  << ", numOldCons=" << numOldCons << std::endl;
#endif
		tempMark[tempPos] = 1;
		
	    }
	    
	    tempInt = 0;
	    for (k = 0; k < numOldCons; ++k) {
		if (tempMark[k] != 1) {
		    // Survived.
		    (model->oldConstraints())[tempInt++]=
			(model->oldConstraints())[k];
		}
	    }
	    if (tempInt + pathDesc->getCons()->numRemove != numOldCons) {
		std::cout << "INSTALL: tempInt=" << tempInt
			  <<", numRemove="<<pathDesc->getCons()->numRemove
			  << ", numOldCons=" << numOldCons << std::endl;
		
		assert(0);
	    }
	    
	    // Update number of old non-core constraints.
	    numOldCons = tempInt;
	    delete [] tempMark;
	}
    } // EOF leafToRootPath.


    //--------------------------------------------------------
    // Debug variable bounds to be installed.
    //--------------------------------------------------------

#ifdef BLIS_DEBUG_MORE
    for (k = 0; k < numCols; ++k) {
        //if (index_ == -1) {
            printf("INSTALL: Col %d, \tlb %g,  \tub %g\n",
                   k, startColLB[k], startColUB[k]);
	    //}
        
        if (startColLB[k] > startColUB[k] + ALPS_GEN_TOL) {
            printf("INSTALL: Col %d, \tlb %g,  \tub %g\n",
                   k, startColLB[k], startColUB[k]);
            assert(0);
        }
    }
#endif   

    //--------------------------------------------------------
    // Clear path vector.
    //--------------------------------------------------------
    
    model->leafToRootPath.clear();
    assert(model->leafToRootPath.size() == 0);
    
    //--------------------------------------------------------
    // Adjust column bounds in lp solver  
    //--------------------------------------------------------
    
    for(i = 0; i < numCols; ++i) {
	model->solver()->setColBounds(i, startColLB[i], startColUB[i]); 
    }
    
    //--------------------------------------------------------
    // TODO: Set row bounds 
    //--------------------------------------------------------
    

    //--------------------------------------------------------
    // Add old constraints, which are collect from differencing.
    //--------------------------------------------------------
    
    // If removed cuts due to local cuts.

    model->setNumOldConstraints(numOldCons);
    
#if 0
    std::cout << "INSTALL: after collecting, numOldCons = " << numOldCons 
	      << std::endl;
#endif

    if (numOldCons > 0) {
	const OsiRowCut ** oldOsiCuts = new const OsiRowCut * [numOldCons];
	for (k = 0; k < numOldCons; ++k) {
	    OsiRowCut * acut = (model->oldConstraints()[k])->createOsiRowCut();
	    oldOsiCuts[k] = acut;
	}
	model->solver()->applyRowCuts(numOldCons, oldOsiCuts);
	for (k = 0; k < numOldCons; ++k) {
	    delete oldOsiCuts[k];
	}
	delete [] oldOsiCuts;
	oldOsiCuts = NULL;
    }
    
    //--------------------------------------------------------
    // Add parent variables, which are collect from differencing.
    //--------------------------------------------------------
    

    //--------------------------------------------------------
    // Set basis
    //--------------------------------------------------------

    CoinWarmStartBasis *pws = desc->getBasis();

    if (pws != NULL) {
	model->solver()->setWarmStart(pws);

#ifdef BLIS_DEBUG
	printf("NODE %d: set warm start\n", getIndex());
	
	numCols = model->solver()->getNumCols();
	numRows = model->solver()->getNumRows();
	int nStr = pws->getNumStructural();
	int nArt = pws->getNumArtificial();
	
	if (numCols != nStr) {
	    std::cout << "nStr=" << nStr << ", numCols=" << numCols 
		      << std::endl;
	    assert(0);
	}
	std::cout << "nArt=" << nArt << ", numRows=" << numRows 
		  << std::endl;
	if (numRows != nArt) {
	    std::cout << "nArt=" << nArt << ", numRows=" << numRows 
		      << std::endl;
	    assert(0);
	}
#endif

    }  

    return status;
}

//#############################################################################

void 
BlisTreeNode::getViolatedConstraints(BlisModel *model, 
				     const double *LpSolution, 
				     BcpsConstraintPool & conPool)
{
    int k;
    int numCons = model->constraintPoolReceive()->getNumConstraints();
    BlisConstraint *blisCon = NULL;
    std::vector<BlisConstraint *> conVector;

    // Check violation and move voilated constraints to conPool
    for (k = 0; k < numCons; ++k) {
	blisCon = dynamic_cast<BlisConstraint *>(conPool.getConstraint(k));
	
	if (blisCon->violation(LpSolution) > ALPS_SMALL_4) {    
	    conPool.addConstraint(blisCon);
	}
	else {
	    conVector.push_back(blisCon);
	}
    }

    if (numCons > 0) {
	std::cout << "Has constraints " << numCons 
		  << "; violated " << numCons-conVector.size()
		  << std::endl;
    }
    
    if ((int)conVector.size() != numCons) {
	// There violated constraints. Remove them from conPool.
	conPool.clear();
	numCons = static_cast<int>(conVector.size());
	for (k = 0; k < numCons; ++k) {
	    conPool.addConstraint(conVector[k]);
	}
	conVector.clear();
    }
}

//#############################################################################

int
BlisTreeNode::generateConstraints(BlisModel *model,BcpsConstraintPool &conPool) 
{
    int i, j, numCGs;
    BlisLpStatus status = BlisLpStatusOptimal;
    int preNumCons = 0;
    int newNumCons = 0;
    BlisCutStrategy strategy = BlisCutStrategyRoot;

    // Only autmatic stategy has depth limit.
    int maxConstraintDepth = 20;
    
    bool mustResolve = false;
    
    double genCutTime;
    
    numCGs = model->numCutGenerators();

    for (i = 0 ; i < numCGs; ++i) {
	
	//----------------------------------------------------
	// Check if call this generator.
	//----------------------------------------------------
	
	strategy =  model->cutGenerators(i)->strategy();
      
	bool useThisCutGenerator = false;
	if (strategy == BlisCutStrategyNone) {
	    useThisCutGenerator = false;
	}
	else if (strategy == BlisCutStrategyRoot) {
	    if (model->isRoot_ && (index_ == 0)) useThisCutGenerator = true;
	}
	else if (strategy == BlisCutStrategyAuto) {
	    if (depth_ < maxConstraintDepth) {
		useThisCutGenerator = true;
	    }
	}
	else if (strategy == BlisCutStrategyPeriodic) {
	    // Num of nodes is set at the beginning of process().
	    if ((model->getNumNodes()-1) %  
		model->cutGenerators(i)->cutGenerationFreq() == 0) {
	        useThisCutGenerator = true;
	    }
	}
	else {
	   throw CoinError("Unknown cut generation strategy", 
			   "generateConstraints", "BlisTreeNode");
	}
	   
#if 0
	std::cout<<"CUTGEN: " << model->cutGenerators(i)->name() 
		 <<": useThisCutGenerator ="<<useThisCutGenerator
		 << ", strategy =" << strategy 
		 << ", num of nodes =" << model->getNumNodes()
		 <<std::endl;
#endif

	//----------------------------------------------------
	// Generator constraints.
	//----------------------------------------------------

	if (useThisCutGenerator) {

	    preNumCons = conPool.getNumConstraints();
          
	    genCutTime = CoinCpuTime();

	    // Call constraint generator
	    mustResolve = model->cutGenerators(i)->generateConstraints(conPool);

	    genCutTime = CoinCpuTime() - genCutTime;
            
            // Statistics
            model->cutGenerators(i)->addTime(genCutTime);
            model->cutGenerators(i)->addCalls(1);

	    newNumCons = conPool.getNumConstraints() - preNumCons;
	    
            if (newNumCons == 0) {
                model->cutGenerators(i)->addNoConsCalls(1);
            }
            else {
                // Reset to 0
                int noCuts = model->cutGenerators(i)->noConsCalls();
                model->cutGenerators(i)->addNoConsCalls(-noCuts);
                model->cutGenerators(i)->addNumConsGenerated(newNumCons);
            }

	    if (mustResolve) {
		// TODO: Only probing will return ture.
		status = static_cast<BlisLpStatus> (bound(model));
		if (status == BlisLpStatusOptimal) {
#ifdef BLIS_DEBUG
		    std::cout << "CUTGEN: after probing, this node survived."
			      << std::endl;
#endif
		}
		else {
#ifdef BLIS_DEBUG
		    std::cout<<"CUTGEN: after probing, this node can fathomed."
			     << std::endl;
#endif
		    break;
		}
	    }
	    
	    //------------------------------------------------
	    // Modify control. 
	    // NOTE: only modify if user choose automatic.
	    //------------------------------------------------
	    
	    if (model->getCutStrategy() == BlisCutStrategyNone) {
                for (j = 0; j < numCGs; ++j) {
                    strategy =  model->cutGenerators(j)->strategy();
                    if (strategy != BlisCutStrategyNone) {
                        break;
                    }
                }
                
                if (j == numCGs) {
                    model->setCutStrategy(BlisCutStrategyNone);
                }
	    }
	}
    }

    return status;
}

//#############################################################################

BlisReturnStatus 
BlisTreeNode::applyConstraints(BlisModel *model, 
			       const double *solution,
                               BcpsConstraintPool & conPool)
{
    BlisReturnStatus status = BlisReturnStatusOk;
    int i, k;
    
    int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);

    int numRowCuts = conPool.getNumConstraints();

    int numToAdd = numRowCuts;
    int numAdded = 0;
    
    if (numRowCuts > 0) {
	BlisParams * BlisPar = model->BlisPar();
	double scaleConFactor = BlisPar->entry(BlisParams::scaleConFactor);
        
        if (numToAdd > 0) { 
	    
            BlisConstraint *blisCon = NULL;
	    
            if (msgLevel > 100) {
                printf("\nAPPLYCUT: Select cuts to be added in LP from %d candidates\n",
                       numRowCuts);
            }
            
            int numRowsNow = model->solver()->getNumRows();
            int numCols = model->solver()->getNumCols();
            CoinWarmStartBasis *ws = dynamic_cast<CoinWarmStartBasis*>
                (model->solver()->getWarmStart());
            
	    // Tranform constraints to Osi cut so that easily add them to LP.
            const OsiRowCut ** addedCuts = new const OsiRowCut * [numToAdd];
	    
            for (i = 0 ; i < numToAdd ; ++i) {
		bool keep = true;
		blisCon = dynamic_cast<BlisConstraint *>(conPool.getConstraint(i));
                
		//------------------------------------------
		// Remove:
		//  - empty cuts
		//  - dense cuts
		//  - bad scaled cuts
                //  - weak cuts
		//  - parallel cuts
		//------------------------------------------

		int length = blisCon->getSize();
		const double *elements = blisCon->getValues();
		const int *indices = blisCon->getIndices();
		
		bool check = true;
                while (check) { // while is used to turn off checking.
                    //--------------------------------------                   
                    // Empty.
                    //--------------------------------------

                    if (length <= 0) {
                        keep = false;

#if 0
                        std::cout << "APPLYCUT: A empty cut." << std::endl;
#endif
                        break;
                    }

                    //--------------------------------------
                    // Dense.
                    //--------------------------------------

                    if(length > model->getDenseConCutoff()){
                        keep = false;
                        if (msgLevel > 100) {
                            std::cout << "APPLYCUT: Discard a dense cut. length = " 
                                      << length << ", cutoff = " 
                                      << model->getDenseConCutoff() << std::endl;
                        }
                        break;
                    }

                    //--------------------------------------
                    // Compuate scale factor.
                    //--------------------------------------

                    int index;
                    double activity = 0.0;
                    
                    double maxElem = 0.0;
                    double minElem = ALPS_DBL_MAX;
                    double scaleFactor;
                    
                    for (k = 0; k < length; ++k) {
                        if (fabs(elements[k]) > maxElem) {
                            maxElem = fabs(elements[k]);
                        }
                        if (fabs(elements[k]) < minElem) {
                            minElem = fabs(elements[k]);
                        }
                        index = indices[k];
                        activity += elements[k] * solution[index];
                    }
                    if(minElem != 0.0) {
                        scaleFactor = maxElem/minElem;
                    }
                    else {
                        assert(0);
                        scaleFactor = ALPS_DBL_MAX;
                    }
                    
#ifdef BLIS_DEBUG
                    std::cout << "APPLYCUT: scaleFactor=" << scaleFactor
                              << ", maxElem=" << maxElem 
                              << ", minElem=" << minElem << std::endl;
#endif
                    if (scaleFactor > scaleConFactor) {
                        if (msgLevel > 100) {
                            std::cout<< "APPLYCUT: Discard a bad scaled cut"
                                     << std::endl;
                        }
                        keep = false;
                        break;
                    }
                    
                    //--------------------------------------
                    // Weak.
                    //--------------------------------------

                    double rowLower = CoinMax(blisCon->getLbHard(),
					      blisCon->getLbSoft());
                    double rowUpper = CoinMin(blisCon->getUbHard(),
					      blisCon->getUbSoft());
                    double violation = -9.87; // Any negative number is OK

		    if (rowLower > -ALPS_INFINITY) {
			violation = rowLower - activity;
		    }
		    if (rowUpper < ALPS_INFINITY) {
			violation = CoinMax(violation, activity-rowUpper);
		    }
                    
                    if (violation < 1.0e-6) {
                        // Found a weak cuts.
                        if (msgLevel > 100) {
                            std::cout<< "WARNING: APPLYCUT: applied a weak cut, violation="
                                     << violation << std::endl;
                        }
                        //keep = false;
                        break;
                    }

                    //--------------------------------------
                    // Parallel cuts.
                    //--------------------------------------
                    
		    bool paral = parallel(model, 
					  conPool,
					  i,
					  blisCon);
		    if (paral) {
                        if (msgLevel > 100) {
                            std::cout<< "APPLYCUT: Discard a parallel cut"
                                     << std::endl;
                        }
			keep = false;
			break;
		    }
                    
                    //--------------------------------------
                    // Check once and stop.
                    //--------------------------------------

                    check = false;
                }//while
                
		if (keep) {
                    addedCuts[numAdded++] = blisCon->createOsiRowCut();
                }
                else {
                    conPool.deleteConstraint(i);
                    --i;
                    --numToAdd;
                }
            }
	    
	    assert(numToAdd == numAdded);
	    
            if (msgLevel > 100) {
                printf("APPLYCUT: After selecting, added %d cuts to LP and discared %d cuts\n",
                       numAdded, numRowCuts - numAdded);
            }
            
            //----------------------------------------------
            // Add cuts to lp and adjust basis.
            //----------------------------------------------

	    if (numAdded > 0) {
		model->solver()->applyRowCuts(numAdded, addedCuts);
		ws->resize(numRowsNow + numToAdd, numCols);
		for (i = 0 ; i < numToAdd; ++i) { 
		    ws->setArtifStatus(numRowsNow + i,
				       CoinWarmStartBasis::basic); 
		}
		if (model->solver()->setWarmStart(ws) == false) { 
		    throw CoinError("Fail setWarmStart() after cut installation.",
				    "applyConstraints","BlisTreeNode"); 
		}

		for (k = 0; k < numAdded; ++k) {
		    delete addedCuts[k];
		}
	    }

	    delete [] addedCuts;
            delete ws;
        }   
    }
    
    return status;
}

//#############################################################################

BlisReturnStatus 
BlisTreeNode::reducedCostFix(BlisModel *model)
{ 
    int i, var;
    BlisReturnStatus status = BlisReturnStatusOk;

    int numFixedUp = 0;
    int numFixedDown = 0;
    int numTighten = 0;

    double movement;    
    double newBound;
    double boundDistance;
    double dj;

    int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);

    const double *lb = model->solver()->getColLower();
    const double *ub = model->solver()->getColUpper();
    const double *solution = model->solver()->getColSolution();
    const double *reducedCost = model->solver()->getReducedCost();
    
    double cutoff = model->getCutoff();
    
    if (cutoff >= ALPS_OBJ_MAX) return status;
    
    double lpObjValue = model->solver()->getObjValue() * 
        model->solver()->getObjSense();
    double epInt = 1.0e-5;
    
    int numIntegers = model->getNumIntObjects();
    const int *intIndices = model->getIntColIndices();
    
    for (i = 0; i < numIntegers; ++i) { 
	var = intIndices[i];
        
        dj = reducedCost[var];
        
        if (fabs(dj) < epInt) continue;
        
        boundDistance = ub[var] - lb[var];
	if (boundDistance < epInt) continue;
	
        movement = floor((cutoff - lpObjValue) / fabs(dj));
        
        if (solution[var] > ub[var] - epInt) {
            /* At upper bound */
            if (movement < boundDistance) {
                /* new lower bound. If movement is 0, then fix. */
                newBound = ub[var] - movement;
                newBound = CoinMin(newBound, ub[var]);

                if (msgLevel > 300) {
                    printf("RED-FIX: dj %g, lb %.10g, ub %.10g, newBound %.10g, movement %g\n", dj, lb[var], ub[var], newBound, movement);
                }
                
                if (movement <= ALPS_ZERO) {
                    ++numFixedUp;
                }
                else if (newBound < ub[var]){
                    ++numTighten;
                }
                model->solver()->setColLower(var, newBound);
            }
        }
        else if (solution[var] < lb[var] + epInt) {
            /* At lower bound */
            if (movement < boundDistance) {
                newBound = lb[var] + movement;
                newBound = CoinMax(newBound, lb[var]);
                
                if (msgLevel > 300) {
                    printf("RED-FIX: dj %g, lb %g, ub %g, newBound %g, movement %g\n", dj, lb[var], ub[var], newBound, movement);
                }
		
                if (movement <= ALPS_ZERO) {
                    ++numFixedDown;
                }
                else if(newBound > lb[var] ){
                    ++numTighten;
                }
                /* new upper bound. If movement is 0, then fix. */
                model->solver()->setColUpper(var, newBound);
            }
        }
    }
    
    //int change = numFixedUp + numFixedDown + numTighten;
    //model->reducedCostFixed_ += change;

    if (msgLevel > 200) {
        if (numFixedUp > 0 || numFixedDown > 0 || numTighten > 0) {
            printf("reducedCostFix: numFixedUp = %d, numFixedDown = %d, numTighten %d\n", numFixedUp, numFixedDown, numTighten);
        }
    }
    
    return status; 
}

//#############################################################################

AlpsEncoded*
BlisTreeNode::encode() const 
{
#ifdef BLIS_DEBUG
    std::cout << "BlisTreeNode::encode()--start to encode node "
	      << index_ << std::endl;
#endif

    AlpsReturnStatus status = AlpsReturnStatusOk;
    
    // NOTE: "AlpsKnowledgeTypeNode" is used as type name.
    AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeNode);
    
    // Encode decription.
    status = desc_->encode(encoded);
    
    // Encode Alps portion.
    status = encodeAlps(encoded);
    
    // Encode Bcps portion.
    status = encodeBcps(encoded);
    
    // Nothing to encode for Blis portion.
    
    return encoded;
}

//#############################################################################

AlpsKnowledge* 
BlisTreeNode::decode(AlpsEncoded& encoded) const 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    BlisTreeNode* treeNode = NULL;

    BlisModel *model = dynamic_cast<BlisModel*>(desc_->getModel());
    
    //------------------------------------------------------
    // Unpack decription.
    //------------------------------------------------------

    AlpsNodeDesc* nodeDesc = new BlisNodeDesc(model);
    status = nodeDesc->decode(encoded);
    
    //------------------------------------------------------
    // Unpack node.
    //------------------------------------------------------
    
    // Unpack Alps portion.
    treeNode = new BlisTreeNode(nodeDesc);
    nodeDesc = NULL;

    treeNode->decodeAlps(encoded);  
    
    // Unpack Bcps portion.
    int type = 0;
    encoded.readRep(type);	
    if (type == BlisBranchingObjectTypeInt) {
	// branchObject_ is simple integer.
	BlisBranchObjectInt *bo = new BlisBranchObjectInt();
	status = bo->decode(encoded);

	// Set bo in treeNode.
	treeNode->setBranchObject(bo);
        bo = NULL;
    }
    
    // Nothing to unpack for Blis portion.
    
    return treeNode;
}

//#############################################################################

void 
BlisTreeNode::convertToExplicit() 
{
#ifdef BLIS_DEBUG
    std::cout << "BLIS: convertToExplicit(); explicit_="<<explicit_ << std::endl;
#endif

    if(!explicit_) {
	
	// Convert to explicit
	explicit_ = 1;
	
	BlisModel* model = dynamic_cast<BlisModel*>(desc_->getModel());
	BlisNodeDesc *desc = dynamic_cast<BlisNodeDesc *>(desc_);
	BlisConstraint *aCon = NULL;
	
	int numCols = model->solver()->getNumCols();

	int i, k, index;
	int tempInt;

	int numModify = 0;
	int numSoftVarLowers = 0;
	int numSoftVarUppers = 0;

	double value;
		
	double *fVarHardLB = new double [numCols];
	double *fVarHardUB = new double [numCols];
	int *fVarHardLBInd = new int [numCols];
	int *fVarHardUBInd = new int [numCols];
	
	double *fVarSoftLB = new double [numCols];
	double *fVarSoftUB = new double [numCols];
	int *fVarSoftLBInd = new int [numCols];
	int *fVarSoftUBInd = new int [numCols];

	for (k = 0; k < numCols; ++k) {
	    fVarSoftLB[k] = ALPS_DBL_MAX;
	    fVarSoftUB[k] = -ALPS_DBL_MAX;
	    fVarHardLB[k] = ALPS_DBL_MAX;
	    fVarHardUB[k] = -ALPS_DBL_MAX;
	    fVarHardLBInd[k] = k;
	    fVarHardUBInd[k] = k;
	}

	int numOldCons = 0;
	int maxOld = model->getOldConstraintsSize();
	BcpsObject ** oldConstraints = new BcpsObject* [maxOld];
	
	//--------------------------------------------------
	// Travel back to a full node, then collect diff (add/rem col/row,
	// hard/soft col/row bounds) from the node full to this node.
	//--------------------------------------------------------
    
	BlisNodeDesc* pathDesc = NULL;
	AlpsTreeNode *parent = parent_;    
	
	model->leafToRootPath.push_back(this);
	
	while(parent) {
#ifdef BLIS_DEBUG_MORE
	    std::cout << "Parent id = " << parent->getIndex() << std::endl;
#endif     
	    model->leafToRootPath.push_back(parent);
	    if (parent->getExplicit()) {
		// Reach an explicit node, then stop.
		break;
	    }
	    else {
		parent = parent->getParent();
	    }
	}
    
#ifdef BLIS_DEBUG
	std::cout << "CONVERT TO EXP: path len = " << model->leafToRootPath.size()
		  << std::endl;
#endif
	
	//------------------------------------------------------
	// Travel back from this node to the explicit node to
	// collect full description.
	//------------------------------------------------------	


	for(i = static_cast<int> (model->leafToRootPath.size() - 1); i > -1; --i) {

	    pathDesc = dynamic_cast<BlisNodeDesc*>((model->leafToRootPath.at(i))->
						   getDesc());

	    //--------------------------------------
	    // Full variable hard bounds.
	    //--------------------------------------
	    
	    numModify = pathDesc->getVars()->lbHard.numModify;
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->lbHard.posModify[k];
		value = pathDesc->getVars()->lbHard.entries[k];
		fVarHardLB[index] = value;
	    }
		    
	    numModify = pathDesc->getVars()->ubHard.numModify;
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->ubHard.posModify[k];
		value = pathDesc->getVars()->ubHard.entries[k];
		fVarHardUB[index] = value;
	    }
		    
	    //--------------------------------------
	    // Full variable soft bounds.
	    //--------------------------------------
	    
	    numModify = pathDesc->getVars()->lbSoft.numModify;
#ifdef BLIS_DEBUG
	    std::cout << "CONVERT: EXP: i=" << i << ", numModify soft lb="
		      << numModify << std::endl;
#endif
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->lbSoft.posModify[k];
		value = pathDesc->getVars()->lbSoft.entries[k];
		fVarSoftLB[index] = value;
	    }
		    
	    numModify = pathDesc->getVars()->ubSoft.numModify;
#ifdef BLIS_DEBUG
	    std::cout << "CONVERT: EXP: i=" << i << ", numModify soft ub="
		      << numModify << std::endl;
#endif
	    for (k = 0; k < numModify; ++k) {
		index = pathDesc->getVars()->ubSoft.posModify[k];
		value = pathDesc->getVars()->ubSoft.entries[k];
		fVarSoftUB[index] = value;
	    }


            //----------------------------------------------
            // Collect all generated constraints, then remove deleted.
            //----------------------------------------------
            
            tempInt = pathDesc->getCons()->numAdd;
	    
#ifdef BLIS_DEBUG
            std::cout << "\nCONVERT: EXP: numAdd = " << tempInt << std::endl;
#endif
            
            for (k = 0; k < tempInt; ++k) {
                aCon = dynamic_cast<BlisConstraint *>
                    (pathDesc->getCons()->objects[k]);
                
                assert(aCon);
                assert(aCon->getSize() > 0);
                assert(aCon->getSize() < 100000);
                
#ifdef BLIS_DEBUG
		std::cout << "CONVERT: EXP: k=" << k 
			  << ", len=" <<aCon->getSize() << std::endl;
#endif
                oldConstraints[numOldCons++] = aCon;
                
                if (numOldCons >= maxOld) {
                    // Need resize
#ifdef BLIS_DEBUG
                    std::cout << "CONVERT: EXP: resize, maxOld = " 
                              << maxOld << std::endl;
#endif
                    maxOld *= 2;
                    BcpsObject **tempCons = new BcpsObject* [maxOld];
		    
                    memcpy(tempCons, 
                           oldConstraints, 
                           numOldCons * sizeof(BcpsObject *));
		    
		    delete [] oldConstraints;
		    oldConstraints = tempCons;
		    tempCons = NULL;
                }
            }
	    
            //----------------------------------------------
            // Remove those deleted. 
	    // NOTE: oldConstraints stores all previously 
	    // generated active constraints at parent.
            //----------------------------------------------
	    
            tempInt = pathDesc->getCons()->numRemove;
            
            if (tempInt > 0) {
                int tempPos;
                int *tempMark = new int [numOldCons];
                CoinZeroN(tempMark, numOldCons);
                for (k = 0; k < tempInt; ++k) {
                    tempPos = pathDesc->getCons()->posRemove[k];
#ifdef BLIS_DEBUG_MORE
		    std::cout << "tempPos=" << tempPos 
			      << ", tempInt=" << tempInt 
			      << ", numOldCons=" << numOldCons << std::endl;
#endif
                    tempMark[tempPos] = 1;
                }
		
                tempInt = 0;
                for (k = 0; k < numOldCons; ++k) {
                    if (tempMark[k] != 1) {
                        // Survived.
                        oldConstraints[tempInt++] = oldConstraints[k];
                    }
                }
		if (tempInt + pathDesc->getCons()->numRemove != numOldCons) {
		    std::cout << "INSTALL: tempInt=" << tempInt
			      <<", numRemove="<<pathDesc->getCons()->numRemove
			      << ", numOldCons=" << numOldCons << std::endl;
		    
		    assert(0);
		}

                // Update number of old non-core constraints.
                numOldCons = tempInt;
                delete [] tempMark;
            }
	    
	} // EOF for (path)
        
	//------------------------------------------
	// Record hard variable bounds. FULL set.
	//------------------------------------------
	
	desc->assignVarHardBound(numCols,
				 fVarHardLBInd,
				 fVarHardLB,
				 numCols,
				 fVarHardUBInd,
				 fVarHardUB);

	//------------------------------------------
	// Recode soft variable bound. Modified.
	//------------------------------------------
	
	for (k = 0; k < numCols; ++k) {
	    if (fVarSoftLB[k] < ALPS_BND_MAX) {
		fVarSoftLBInd[numSoftVarLowers] = k;
		fVarSoftLB[numSoftVarLowers++] = fVarSoftLB[k];
	    }
	    if (fVarSoftUB[k] > -ALPS_BND_MAX) {
		fVarSoftUBInd[numSoftVarUppers] = k;
		fVarSoftUB[numSoftVarUppers++] = fVarSoftUB[k];
	    }
	}
	// Assign it anyway so to delete memory(fVarSoftLBInd,etc.)
	desc->assignVarSoftBound(numSoftVarLowers, 
				 fVarSoftLBInd,
				 fVarSoftLB,
				 numSoftVarUppers, 
				 fVarSoftUBInd,
				 fVarSoftUB);
	
	//------------------------------------------
	// Recode added constraints.
	//------------------------------------------

	// First make a hard copy.
	for (k = 0; k < numOldCons; ++k) {
	    aCon = dynamic_cast<BlisConstraint *>(oldConstraints[k]);
	    assert(aCon);
	    
	    BlisConstraint *newCon = new BlisConstraint(*aCon);
	    oldConstraints[k] = newCon;
	}
	// Add will first delete, then add. It is safe to use in parallel.
	desc->setAddedConstraints(numOldCons, oldConstraints);

	//------------------------------------------
	// Recode deleted constraints.
	//------------------------------------------
	
	desc->delConstraints(0, NULL);
	
	//--------------------------------------------------
	// Clear path vector.
	//--------------------------------------------------
	
	model->leafToRootPath.clear();
	assert(model->leafToRootPath.size() == 0);


    } // EOF of if.
}

//#############################################################################

// Not defined yet.
void
BlisTreeNode::convertToRelative()
{
    if(explicit_) {
	

    }
}

//#############################################################################

bool 
BlisTreeNode::parallel(BlisModel *model, 
		       BcpsConstraintPool &conPool,
		       int lastNew,
		       BlisConstraint *aCon)
{
    bool parallel = false;
    int k;
    double threshold = 0.999;
    
    //------------------------------------------------------
    // Compare with old cuts
    //------------------------------------------------------
#if 0  // couse error for VRP
    int numOldCons = model->getNumOldConstraints();
    for (k = 0; k < numOldCons; ++k) {
	BlisConstraint *aCon = model->oldConstraints()[k];
	assert(aCon);
	parallel = BlisParallelCutCon(rowCut,
				      aCon,
				      threshold);
	if (parallel) return parallel;
    }
#endif    

    //------------------------------------------------------
    // Compare with new cuts
    //------------------------------------------------------

    for (k = 0; k < lastNew; ++k) {
	BlisConstraint *thisCon = 
	    dynamic_cast<BlisConstraint *>(conPool.getConstraint(k));
	parallel = BlisParallelConCon(aCon,
				      thisCon,
				      threshold);
	if (parallel) return parallel;
    }
    
    return parallel;
}

//#############################################################################

double 
BlisTreeNode::estimateSolution(BlisModel *model,
                               const double *lpSolution,
                               double lpObjValue) const
{
    // lpObjective + sum_i{downCost_i*f_i + upCost_i*(1-f_i)} 
    int k, col;
    int numInts= model->getNumIntObjects();

    double x, f, downC, upC, estimate = lpObjValue;

    BlisObjectInt *obj = NULL;
    
    for (k = 0; k < numInts; ++k) {
        obj = dynamic_cast<BlisObjectInt *>(model->objects(k));
        col = obj->columnIndex();
        x = lpSolution[col];
        f = CoinMax(0.0, x - floor(x));
        if (f > model->integerTol_) {
            downC = obj->pseudocost().getDownCost();
            upC = obj->pseudocost().getUpCost();
            estimate += (downC * f + upC * (1-f));
        }
    }
    return estimate;
}

//#############################################################################

int
BlisTreeNode::callHeuristics(BlisModel *model, bool onlyBeforeRoot)
{
    int status = 0;

    if (model->heurStrategy_ == BlisHeurStrategyNone) {
        return status;
    }

    int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);

    int foundSolution = false;
    int numCols = model->solver()->getNumCols();

    double heurObjValue = getKnowledgeBroker()->getIncumbentValue();
    double *heurSolution = new double [numCols];

    BlisSolution *bSol = NULL;

    for (int k = 0; k < model->numHeuristics(); ++k) {
	int heurStrategy = model->heuristics(k)->strategy();
        //std::cout << " call heur " << k << "; strategy = " 
        //        << heurStrategy << std::endl;

	if (heurStrategy != BlisHeurStrategyNone) {
	    if (onlyBeforeRoot) {
		// heuristics that can only be used before root.
		if (heurStrategy != BlisHeurStrategyBeforeRoot) {
		    continue;
		}
	    }
	    else {
		// regular heuristics
		if (heurStrategy == BlisHeurStrategyBeforeRoot) {
		    continue;
		}
	    }
	    
	    getKnowledgeBroker()->tempTimer().start();
	    foundSolution = false;
	    foundSolution = 
		model->heuristics(k)->searchSolution(heurObjValue,
						     heurSolution);
	    getKnowledgeBroker()->tempTimer().stop();
	    
	    model->heuristics(k)->
		addTime(getKnowledgeBroker()->tempTimer().getCpuTime());
	    
	    if (foundSolution) {
		// Check if solution from heuristic is feasible.
		bSol = model->feasibleSolutionHeur(heurSolution);
	    }
            
	    if (bSol) {
		model->heuristics(k)->addNumSolutions(1);
		int noSols = model->heuristics(k)->noSolCalls();
		model->heuristics(k)->addNoSolCalls(-noSols);
		// Store the newly found blis solution.
		model->storeSolution(BlisSolutionTypeHeuristic, bSol);
		if (onlyBeforeRoot) {
		    status = 1;
		}
		else if (quality_ >  model->getCutoff()) {
		    setStatus(AlpsNodeStatusFathomed);
		    status = 2;
		    goto TERM_HEUR;
		}
		else {
		    status = 1;
		}
		if (heurStrategy == BlisHeurStrategyBeforeRoot && 
		    msgLevel > 200) {
		    model->blisMessageHandler()->message(BLIS_HEUR_BEFORE_ROOT, 
							 model->blisMessages())
			<< (model->heuristics(k)->name())
			<< bSol->getQuality()
                        << CoinMessageEol;
		}
	    }
	    else {
		model->heuristics(k)->addNoSolCalls(1);
	    }  
	} 
    }

TERM_HEUR:

    if (heurSolution) delete [] heurSolution;

    return status;
}

//#############################################################################

