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
#include "CoinPackedVector.hpp"
#include "CoinWarmStartBasis.hpp"

#include "OsiRowCut.hpp"

#include "AlpsKnowledgeBroker.h"

#include "BlisObjectInt.h"
#include "BlisConstraint.h"
#include "BlisHelp.h"
#include "BlisModel.h"
#include "BlisSolution.h"

//#############################################################################

/** Convert a OsiRowCut to a Blis Contraint. */
BlisConstraint * BlisOsiCutToConstraint(const OsiRowCut *rowCut)
{
    int size = rowCut->row().getNumElements();
    assert(size > 0);
    
    const int *ind = rowCut->row().getIndices();
    const double *val = rowCut->row().getElements();

    double lower = rowCut->lb();
    double upper = rowCut->ub();
    
    BlisConstraint *con = new BlisConstraint(lower, upper, 
                                             lower, upper,
                                             size, ind, val);
    
    if (!con) {
        // No memory
        throw CoinError("Out of Memory", "Blis_OsiCutToConstraint", "NONE");
    }
    
    return con;
}

//#############################################################################

BlisReturnStatus
BlisStrongBranch(BlisModel *model, double objValue, int colInd, double x,
		 const double *saveLower, const double *saveUpper,
		 bool &downKeep, bool &downFinished, double &downDeg,
		 bool &upKeep, bool &upFinished, double &upDeg)
{
    BlisReturnStatus status = BlisReturnStatusOk;
    int lpStatus = 0;

    int j, numIntInfDown, numObjInfDown;

    double newObjValue;
    
    OsiSolverInterface * solver = model->solver();
    
    int numCols = solver->getNumCols();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    bool checkFeasibleInBranch =
       model->BlisPar()->entry(BlisParams::checkFeasibleInBranch);
    // Restore bounds
    int numDiff = 0;

    BlisSolution* ksol = NULL;

    int ind = model->getIntObjIndices()[colInd];
    BlisObjectInt *intObj = dynamic_cast<BlisObjectInt *>(model->objects(ind));
    
#ifdef BLIS_DEBUG_MORE
    for (j = 0; j < numCols; ++j) {
	if (saveLower[j] != lower[j]) {
	    //solver->setColLower(j, saveLower[j]);
            ++numDiff;
	}
	if (saveUpper[j] != upper[j]) {
	    //solver->setColUpper(j, saveUpper[j]);
            ++numDiff;
	}
    }
    std::cout << "BEFORE: numDiff = " << numDiff << std::endl;
#endif	 
   
    //------------------------------------------------------
    // Branching down.
    //------------------------------------------------------

    solver->setColUpper(colInd, floor(x));
    solver->solveFromHotStart();
    
    newObjValue = solver->getObjSense() * solver->getObjValue();
    downDeg = newObjValue - objValue;
    
    if (solver->isProvenOptimal()) {
	lpStatus = 0; // optimal
#ifdef BLIS_DEBUG_MORE
        printf("STRONG: COL[%d]: downDeg=%g, x=%g\n", colInd, downDeg, x);
#endif
        // Update pseudocost
        intObj->pseudocost().update(-1, downDeg, x);
        model->setSharedObjectMark(ind);        

        // Check if ip feasible
        if (checkFeasibleInBranch){
              ksol = model->feasibleSolution(numIntInfDown, numObjInfDown);
              if (ksol) {
#ifdef BLIS_DEBUG_MORE
                 printf("STRONG:Down:found a feasible solution\n");
#endif
                 
                 model->storeSolution(BlisSolutionTypeStrong, ksol);
                 downKeep = false;
              }
              else {
                 downKeep = true;
              }
              downFinished = true;
        }
    }
    else if (solver->isIterationLimitReached() && 
	     !solver->isDualObjectiveLimitReached()) {
	lpStatus = 2;      // unknown 
	downKeep = true;
	downFinished = false;
    }
    else {
        downDeg = 1.0e20;
	lpStatus = 1; // infeasible
	downKeep = false;
	downFinished = false;
    }       
            
#ifdef BLIS_DEBUG_MORE
    std::cout << "Down: lpStatus = " << lpStatus << std::endl;
#endif
    
    // restore bounds
    numDiff = 0;
    for (j = 0; j < numCols; ++j) {
	if (saveLower[j] != lower[j]) {
	    solver->setColLower(j, saveLower[j]);
            ++numDiff;
	}
	if (saveUpper[j] != upper[j]) {
	    solver->setColUpper(j, saveUpper[j]);
            ++numDiff;
	}
    }
#ifdef BLIS_DEBUG
    assert(numDiff > 0);
    //std::cout << "numDiff = " << numDiff << std::endl;
#endif	    
          
    //----------------------------------------------
    // Branching up.
    //----------------------------------------------
    
    solver->setColLower(colInd, ceil(x));
    solver->solveFromHotStart();

    newObjValue = solver->getObjSense() * solver->getObjValue();
    upDeg = newObjValue - objValue;
    
    if (solver->isProvenOptimal()) {
	lpStatus = 0; // optimal

#ifdef BLIS_DEBUG_MORE
        printf("STRONG: COL[%d]: upDeg=%g, x=%g\n", colInd, upDeg, x);
#endif

        // Update pseudocost
        intObj->pseudocost().update(1, upDeg, x);
        model->setSharedObjectMark(ind);        

        // Check if IP feasible
        if (checkFeasibleInBranch){
           ksol = model->feasibleSolution(numIntInfDown, numObjInfDown);
           if (ksol) {
#ifdef BLIS_DEBUG_MORE
              printf("STRONG:Up:found a feasible solution\n");
#endif
              
              model->storeSolution(BlisSolutionTypeStrong, ksol);
              upKeep = false;
           }
           else {
              upKeep = true;
           }
           upFinished = true;
        }
    }
    else if (solver->isIterationLimitReached()
	     &&!solver->isDualObjectiveLimitReached()) {
	lpStatus = 2; // unknown 
	upKeep = true;
	upFinished = false;
    }
    else {
	lpStatus = 1; // infeasible
	upKeep = false;
	upFinished = false;
        upDeg = 1.0e20;
    }
    
#ifdef BLIS_DEBUG_MORE
    std::cout << "STRONG: Up: lpStatus = " << lpStatus << std::endl;
#endif      
    
    // restore bounds
    for (j = 0; j < numCols; ++j) {
	if (saveLower[j] != lower[j]) {
	    solver->setColLower(j,saveLower[j]);
	}
	if (saveUpper[j] != upper[j]) {
	    solver->setColUpper(j,saveUpper[j]);
	}
    }

    return status;
}

//#############################################################################

int BlisEncodeWarmStart(AlpsEncoded *encoded, const CoinWarmStartBasis *ws)
{

    BlisReturnStatus status = BlisReturnStatusOk;
    int numCols = ws->getNumStructural();
    int numRows = ws->getNumArtificial();

    encoded->writeRep(numCols);
    encoded->writeRep(numRows);

    // Pack structural.
    int nint = (ws->getNumStructural() + 15) >> 4;
    encoded->writeRep(ws->getStructuralStatus(), nint * 4);
    
    // Pack artificial.
    nint = (ws->getNumArtificial() + 15) >> 4;
    encoded->writeRep(ws->getArtificialStatus(), nint * 4);

    return status;
}

//#############################################################################

CoinWarmStartBasis *BlisDecodeWarmStart(AlpsEncoded &encoded,
					AlpsReturnStatus *rc) 
{
    int numCols;
    int numRows;
    
    encoded.readRep(numCols);
    encoded.readRep(numRows);
    
    int tempInt;
    
    // Structural
    int nint = (numCols + 15) >> 4;
    char *structuralStatus = new char[4 * nint];
    encoded.readRep(structuralStatus, tempInt);
    assert(tempInt == nint*4);
    
    // Artificial
    nint = (numRows + 15) >> 4;
    char *artificialStatus = new char[4 * nint];
    encoded.readRep(artificialStatus, tempInt);
    assert(tempInt == nint*4);

    CoinWarmStartBasis *ws = new CoinWarmStartBasis();
    if (!ws) {
	throw CoinError("Out of memory", "BlisDecodeWarmStart", "HELP");
    }
    
    ws->assignBasisStatus(numCols, numRows, 
			  structuralStatus, artificialStatus);
    
    assert(!structuralStatus);
    assert(!artificialStatus);

    return ws;   
}

//#############################################################################

/** Compute and return a hash value of an Osi row cut. */
double BlisHashingOsiRowCut(const OsiRowCut *rowCut, 
			    const BlisModel *model)
{
    int size = rowCut->row().getNumElements();
    assert(size > 0);

    int ind, k;

    const int *indices = rowCut->row().getIndices();
    const double * randoms = model->getConRandoms();
    
    double hashValue_ = 0.0;    
    for (k = 0; k < size; ++k) {
	ind = indices[k];
	hashValue_ += randoms[ind] * ind;
    }
#ifdef BLIS_DEBUG_MORE
    std::cout << "hashValue_=" << hashValue_ << std::endl;
#endif
    return hashValue_;
}

//#############################################################################

/** Check if a row cut parallel with another row cut. */
bool BlisParallelCutCut(OsiRowCut * rowCut1,
			OsiRowCut * rowCut2,
			double threshold) 
{
    int size1 = rowCut1->row().getNumElements();
    int size2 = rowCut2->row().getNumElements();
    assert(size1 > 0 && size2 > 0);
    int i, j, k;
    bool para = false;   

    if (size1 != size2) {
        return para;
    }

    //------------------------------------------------------
    // Assume no duplicate indices.
    // rowCut->setRow() tests duplicate indices.
    //------------------------------------------------------

    rowCut1->sortIncrIndex();
    rowCut2->sortIncrIndex();

    const int *indices1 = rowCut1->row().getIndices();
    const double *elems1 = rowCut1->row().getElements();

    const int *indices2 = rowCut2->row().getIndices();
    const double *elems2 = rowCut2->row().getElements();
    
    //------------------------------------------------------
    // Compute norms.
    //------------------------------------------------------

    double norm1 = 0.0;
    double norm2 = 0.0;

    for (k = 0; k < size1; ++k) {
	norm1 += elems1[k]*elems1[k];
    }
    for (k = 0; k < size2; ++k) {
	norm2 += elems2[k]*elems2[k];
    }
    norm1 = sqrt(norm1);
    norm2 = sqrt(norm2);

    //------------------------------------------------------
    // Compute angel cut1 * cut2 
    //------------------------------------------------------

    double denorm = 0.0;
    double angle;
    int index1, index2;
    i = 0; j = 0;
    while(true) {
	index1 = indices1[i];
	index2 = indices2[j];
	if (index1 == index2) {
	    denorm += (elems1[i] * elems2[j]);
	    ++i;
	    ++j;
	    if ((i >= size1) || (j >= size2)) {
		break;
	    }
	}
	else if (index1 > index2) {
	    ++j;
	    if (j >= size2) break;
	}
	else {
	    ++i;
	    if (i >= size1) break;
	}
    }
    denorm = fabs(denorm);
    angle = denorm/(norm1 * norm2);
    assert(angle >= 0.0 && angle <= 1.000001);
    
    if (angle >= threshold) {
	para = true;
#if 0
        for (j = 0; j < size1; ++j) {
            std::cout << indices1[j] << ", " << indices2[j] << "; ";
        }
        std::cout << std::endl;
#endif
    }

    return para;
}

//#############################################################################

/** Check if a row cut parallel with a constraint. */
bool BlisParallelCutCon(OsiRowCut * rowCut,
			BlisConstraint * con,
			double threshold)
{
    bool parallel = false;

    // Convert con to row cut
    OsiRowCut * rowCut2 = con->createOsiRowCut();
    
    parallel = BlisParallelCutCut(rowCut,
				  rowCut2,
				  threshold);
    delete rowCut2;
    
    return parallel;
}

//#############################################################################

/** Check if a row cut parallel with a constraint. */
bool BlisParallelConCon(BlisConstraint * con1,
			BlisConstraint * con2,
			double threshold)
{
    bool parallel = false;

    // Convert con to row cut
    OsiRowCut * rowCut1 = con1->createOsiRowCut();
    OsiRowCut * rowCut2 = con2->createOsiRowCut();
    
    parallel = BlisParallelCutCut(rowCut1,
				  rowCut2,
				  threshold);
    delete rowCut1;
    delete rowCut2;
    
    return parallel;
}

//#############################################################################
