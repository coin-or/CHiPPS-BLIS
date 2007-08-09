/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * ALPS is distributed under the Common Public License as part of the        *
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
 * Copyright (C) 2001-2007, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include <cstdio>
#include <vector>

#include "CoinFinite.hpp"
#include "CoinTime.hpp"
#include "OsiClpSolverInterface.hpp"

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"

#include "BcpsObject.h"

#include "BlisBranchObjectInt.h"
#include "BlisBranchStrategyMaxInf.h"
#include "BlisBranchStrategyPseudo.h"
#include "BlisBranchStrategyRel.h"
#include "BlisBranchStrategyStrong.h"

#include "BlisConstraint.h"
#include "BlisHeurRound.h"
#include "BlisModel.h"
#include "BlisObjectInt.h"
#include "BlisSolution.h"
#include "BlisTreeNode.h"
#include "BlisVariable.h"

//#############################################################################
//-------------------------------------------------------------------
// Returns the greatest common denominator of two 
// positive integers, a and b, found using Euclid's algorithm 
//-------------------------------------------------------------------
static int gcd(int a, int b) 
{
    int remainder = -1;
    // make sure a<=b (will always remain so)
    if(a > b) {
	// Swap a and b
	int temp = a;
	a = b;
	b = temp;
    }
    // if zero then gcd is nonzero (zero may occur in rhs of packed)
    if (!a) {
	if (b) {
	    return b;
	} else {
	    printf("**** gcd given two zeros!!\n");
	    abort();
	}
    }
    while (remainder) {
	remainder = b % a;
	b = a;
	a = remainder;
    }
    return b;
}

//#############################################################################

void 
BlisModel::init() 
{
    processType_ = AlpsProcessTypeMaster;
    origLpSolver_ = NULL;
    presolvedLpSolver_ = NULL;
    lpSolver_ = NULL;

    numCols_ = 0;
    numRows_ = 0;
    numElems_ = 0;

    colMatrix_ = 0;

    varLB_ = NULL;
    varUB_ = NULL;
    conLB_ = NULL;
    conUB_ = NULL;

    startVarLB_ = NULL;
    startVarUB_ = NULL;
    startConLB_ = NULL;
    startConUB_ = NULL;

    tempVarLBPos_ = NULL;
    tempVarUBPos_ = NULL;
    tempConLBPos_ = NULL;
    tempConUBPos_ = NULL;
    
    objSense_ = 1.0;
    objCoef_ = NULL;

    colType_ = 0;
    numIntObjects_ = 0;
    intColIndices_ = NULL;
    intObjIndices_ = NULL;
    
    presolve_ = new BlisPresolve();
    
    numSolutions_ = 0;
    numHeurSolutions_ = 0;
    
    incumbent_ = NULL;
    
    activeNode_ = NULL;
    numStrong_ = 0;
    numBranchResolve_ = 0;
    
    cutoff_ = COIN_DBL_MAX;
    incObjValue_ = COIN_DBL_MAX;
    blisMessageHandler_ = new CoinMessageHandler();
    blisMessages_ = BlisMessage();

    objects_ = NULL;
    numObjects_ = 0;
    
    numNodes_ = 0;
    numIterations_ = 0;
    aveIterations_ = 0;
    
    branchStrategy_ = NULL;
    rampUpBranchStrategy_ = NULL;
    priority_ = NULL;
    nodeWeight_ = 1.0;

    isRoot_ = true;
    boundingPass_ = 0;
    BlisPar_ = new BlisParams;

    /// Timing
    startTime_ = 0.0;

    /// Tolerance 
    integerTol_ = 1.0e-5;
    optimalRelGap_ = 1.0e-4;
    optimalAbsGap_ = 1.0e-6;

    /// Heuristic
    heurStrategy_ = BlisHeurStrategyAuto;
    heurCallFrequency_ = 1;
    numHeuristics_ = 0;
    heuristics_ = NULL;

    /// Cons related
    cutStrategy_ = BlisCutStrategyRoot;
    cutGenerationFrequency_ = 1;
    numCutGenerators_ = 0;
    generators_ = NULL;
    constraintPool_ = NULL;
    oldConstraints_ = NULL;
    oldConstraintsSize_ = 0;
    numOldConstraints_ = 0;
    conRandoms_ = NULL;

    sharedObjectMark_ = NULL;
}

//#############################################################################

// Read from file (currently MPS format. TODO: LP format).

void
BlisModel::readInstance(const char* dataFile)
{
    int j;
    
    int msgLevel =  AlpsPar_->entry(AlpsParams::msgLevel);

    //------------------------------------------------------
    // Read in data from MPS file.
    //------------------------------------------------------
    
    CoinMpsIO *mps = new CoinMpsIO;
    
    int rc = mps->readMps(dataFile, "");
    if(rc) {
        delete mps;
        throw CoinError("Unable to read in instance",
                        "readInstance",
                        "BlisModel");
    }
    mps->messageHandler()->setLogLevel(msgLevel);

    //------------------------------------------------------
    // Get problem data.
    //------------------------------------------------------

    numCols_ = mps->getNumCols();
    numRows_ = mps->getNumRows();
    numElems_ = mps->getNumElements();

    colMatrix_ = new CoinPackedMatrix();    
    *colMatrix_ = *(mps->getMatrixByCol());
    
    varLB_ = new double [numCols_];
    varUB_ = new double [numCols_];

    conLB_ = new double [numRows_];
    conUB_ = new double [numRows_];
    
    memcpy(varLB_, mps->getColLower(), sizeof(double) * numCols_);
    memcpy(varUB_, mps->getColUpper(), sizeof(double) * numCols_);
    
    memcpy(conLB_, mps->getRowLower(), sizeof(double) * numRows_);
    memcpy(conUB_, mps->getRowUpper(), sizeof(double) * numRows_);
    
    objSense_ = BlisPar_->entry(BlisParams::objSense);
    
    objCoef_ = new double [numCols_];
    if (objSense_ > 0.0) {
        memcpy(objCoef_, mps->getObjCoefficients(), sizeof(double) * numCols_);
    }
    else {
        const double *mpsObj =  mps->getObjCoefficients();
        for (j = 0; j < numCols_; ++j) {
            objCoef_[j] = - mpsObj[j];
        }
    }    

    //------------------------------------------------------
    // Set colType_
    //------------------------------------------------------
    
    colType_ = new char [numCols_];   

    for(j = 0; j < numCols_; ++j) {
	if (mps->isContinuous(j)) {
            colType_[j] = 'C';
	}
        else {
            if (varLB_[j] == 0 && varUB_[j] == 1.0) {
		colType_[j] = 'B';
	    }
            else {
                colType_[j] = 'I';
            }
        }
    }

    //-------------------------------------------------------------
    // Create variables and constraints.
    //-------------------------------------------------------------

    createObjects();
    
    delete mps;
}

//############################################################################ 

void 
BlisModel::createObjects()
{
    int j;
    
    //------------------------------------------------------
    // Create variables and constraints.
    //------------------------------------------------------
    
#ifdef BLIS_DEBUG
    std::cout << "createObjects: numCols_ " << numCols_ 
	      << ", numRows_" << numRows_ 
	      << std::endl;
#endif

    const double *elements = colMatrix_->getElements();
    const int *indices = colMatrix_->getIndices();
    const int *lengths = colMatrix_->getVectorLengths();
    const CoinBigIndex *starts = colMatrix_->getVectorStarts();

    int beg = 0;
    
    for (j = 0; j < numCols_; ++j) {
        
        beg = starts[j];

	BlisVariable * var = new BlisVariable(varLB_[j],
					      varUB_[j], 
					      varLB_[j], 
					      varUB_[j],
                                              objCoef_[j], 
                                              lengths[j],
                                              indices + beg,
                                              elements + beg);
        
	var->setObjectIndex(j);
	var->setRepType(BCPS_CORE);
	var->setStatus(BCPS_NONREMOVALBE);
	var->setIntType(colType_[j]);
	variables_.push_back(var);
        var = NULL;
    }
    
    for (j = 0; j < numRows_; ++j) {
        BlisConstraint *con = new BlisConstraint(conLB_[j], 
                                                 conUB_[j], 
                                                 conLB_[j], 
                                                 conUB_[j]);
        con->setObjectIndex(j);
        con->setRepType(BCPS_CORE);
        con->setStatus(BCPS_NONREMOVALBE);
        constraints_.push_back(con);
        con = NULL;        
    }
    
    // Set all objects as core by default.
    numCoreVariables_ = numCols_;
    numCoreConstraints_ = numRows_;
}

//#############################################################################

/** For parallel code, only the master calls this function.
 *  1) Set colMatrix_, varLB_, varUB_, conLB_, conUB
 *     numCols_, numRows_
 *  2) Set objCoef_ and objSense_
 *  3) Set colType_ ('C', 'I', or 'B')
 *  4) Set variables_ and constraints_
 *  NOTE: Blis takes over the memory ownship of vars and cons, which 
 *        means users must NOT free vars or cons.
 */
void 
BlisModel::importModel(std::vector<BlisVariable *> vars,
                       std::vector<BlisConstraint *> cons)
{

    CoinBigIndex i, j;
    int k, size;  

    int* varIndices = NULL;
    double *varValues = NULL;

    numCols_ = static_cast<int> (vars.size());
    numRows_ = static_cast<int> (cons.size());

    varLB_ = new double [numCols_];
    varUB_ = new double [numCols_];

    conLB_ = new double [numRows_];
    conUB_ = new double [numRows_];
    
    objCoef_ = new double [numCols_];
    colType_ = new char [numCols_];
    
    // Get numElems_ and colType_
    for (i = 0; i < numCols_; ++i) {
        numElems_ += vars[i]->getSize();
        colType_[i] = vars[i]->getIntType();
    }

#if 1
    std::cout << "numCols_ = " << numCols_ 
              << "; numRows_ = " << numRows_ << std::endl;
#endif

    //------------------------------------------------------
    // Set matrix, bounds
    //------------------------------------------------------
    
    // For colMatrix_, need free memory
    CoinBigIndex * start = new CoinBigIndex [numCols_+1];
    int* indices = new int[numElems_];
    double* values = new double[numElems_];
    int* length = new int [numCols_];

    // Get varLB_, varUB_, objCoef_, and matrix from variables
    for (numElems_ = 0, i = 0; i < numCols_; ++i) {
        varLB_[i] = vars[i]->getLbHard();
        varUB_[i] = vars[i]->getUbHard();
        objCoef_[i] = vars[i]->getObjCoef();
        
        start[i] = numElems_;

        varValues = vars[i]->getValues();
        varIndices = vars[i]->getIndices();
        size = vars[i]->getSize();
        for (j = 0; j < size; ++j, ++numElems_){
            indices[numElems_] = varIndices[j];
            values[numElems_] = varValues[j];
        }
    }

    // Last position.
    start[numCols_] = numElems_;
    
    for (k = 0; k < numCols_; ++k) {
        length[k] = start[(k+1)] - start[k];
    }
    
    colMatrix_ =  new CoinPackedMatrix(true, // column-majored
                                       numRows_,
                                       numCols_,
                                       numElems_,
                                       values, 
                                       indices,
                                       start,
                                       length);


    // Get conLB_ and conUB_
    for (k = 0; k < numRows_; ++k) {
        conLB_[k] = cons[k]->getLbHard();
        conUB_[k] = cons[k]->getUbHard();
    }
        
    //------------------------------------------------------
    // Set variables_ and constraints_
    //------------------------------------------------------
    
    for (k = 0; k < numCols_; ++k) {
        variables_.push_back(vars[k]);
    }
    
    for (k = 0; k < numRows_; ++k) {
        constraints_.push_back(cons[k]);
    }
    
    //------------------------------------------------------
    // Free memory
    //------------------------------------------------------

    delete [] start;
    delete [] length;
    delete [] indices;
    delete [] values;    
}

//############################################################################ 

/** Read in Alps parameters. */
void 
BlisModel::readParameters(const int argnum, const char * const * arglist)
{    //std::cout << "Reading in ALPS parameters ..." << std::endl;
    AlpsPar_->readFromArglist(argnum, arglist);
    //std::cout << "Reading in BLIS parameters ..." << std::endl;
    BlisPar_->readFromArglist(argnum, arglist);
} 

//##############################################################################

/** Write out parameters. */
void 
BlisModel::writeParameters(std::ostream& outstream) const
{
    outstream << "\n================================================"
              <<std::endl;
    outstream << "ALPS Parameters: " << std::endl;
    AlpsPar_->writeToStream(outstream);
    outstream << "\n================================================"
              <<std::endl;
    outstream << "BLIS Parameters: " << std::endl;
    BlisPar_->writeToStream(outstream);
}

//############################################################################ 

/** Do necessary work to make model usable and presolve.
    Return success or not. 
    This function is called when constructing knowledge broker. */
bool 
BlisModel::setupSelf()
{
    int j;

    if (broker_->getMsgLevel() > 0) {

        //std::cout << "**** getProcType = " << broker_->getProcType() << std::endl;
        //if (broker_->getProcType() == AlpsProcessTypeMaster) {
        if (broker_->getProcRank() == broker_->getMasterRank()) {
            bcpsMessageHandler_->message(BCPS_S_VERSION, bcpsMessages())
                << CoinMessageEol;
            blisMessageHandler()->message(BLIS_S_VERSION, blisMessages())
                << CoinMessageEol;
        }
    }
    
    //------------------------------------------------------
    // Starting time.
    //------------------------------------------------------
    
    startTime_ = CoinCpuTime();

    processType_ = broker_->getProcType();

    //------------------------------------------------------
    // Set numIntObjects_, intColIndices_, intObjectIndices_ 
    //------------------------------------------------------

    // TODO: now integer, later other objects
    intObjIndices_ = new int [numCols_];
    memset(intObjIndices_, 0, sizeof(int) * numCols_);
    
    numIntObjects_ = 0;  
    intColIndices_ = new int [numCols_];
 
    for(j = 0; j < numCols_; ++j) {
	if (colType_[j] == 'I' || colType_[j] == 'B') {
            intColIndices_[numIntObjects_++] = j;
	}
    }
    if (numIntObjects_ == 0) {
        if (broker_->getMsgLevel() > 0) {
            bcpsMessageHandler_->message(BLIS_W_LP, blisMessages())
                << CoinMessageEol;
        }
    }   

    //------------------------------------------------------
    // Load data to LP solver.
    //------------------------------------------------------

    if (!lpSolver_) {
        // preprocessing causes this check.
        lpSolver_ = origLpSolver_;
    }
    
    lpSolver_->loadProblem(*colMatrix_,
			   varLB_, varUB_, 
			   objCoef_,
			   conLB_, conUB_);
    
    lpSolver_->setObjSense(1.0);
    lpSolver_->setInteger(intColIndices_, numIntObjects_);

    //------------------------------------------------------
    // Create integer objects and analyze objective coef.
    //------------------------------------------------------
    
    createIntgerObjects(true);   

    // Do this after loading LP.
    analyzeObjective();
    
    //------------------------------------------------------
    // Allocate memory.
    //------------------------------------------------------

    startVarLB_ = new double [numCols_];
    startVarUB_ = new double [numCols_];

    startConLB_ = new double [numRows_];
    startConUB_ = new double [numRows_];

    tempVarLBPos_ = new int [numCols_];
    tempVarUBPos_ = new int [numCols_];

    tempConLBPos_ = new int [numRows_];
    tempConUBPos_ = new int [numRows_];

    //------------------------------------------------------
    // Get parameters.
    //------------------------------------------------------
    
    integerTol_ = BlisPar_->entry(BlisParams::integerTol);
    optimalRelGap_ = BlisPar_->entry(BlisParams::optimalRelGap);
    optimalAbsGap_ = BlisPar_->entry(BlisParams::optimalAbsGap);
    
    int relibility = BlisPar_->entry(BlisParams::pseudoRelibility);
    cutoff_ =  BlisPar_->entry(BlisParams::cutoff);

    //------------------------------------------------------
    // Modify parameters.
    //------------------------------------------------------
    
    // Disable Alps message
    // AlpsPar()->setEntry(AlpsParams::msgLevel, 1);
    
#ifdef BLIS_DEBUG_MORE
    std::string problemName;
    lpSolver_->getStrParam(OsiProbName, problemName);
    printf("BLIS: setupSelf: Problem name - %s\n", problemName.c_str());
    lpSolver_->setHintParam(OsiDoReducePrint, false, OsiHintDo, 0);
#endif

    //------------------------------------------------------
    // Set branch strategy.
    //------------------------------------------------------

    int brStrategy = BlisPar_->entry(BlisParams::branchStrategy);

    if (brStrategy == BlisBranchingStrategyMaxInfeasibility) {
        // Max inf
        branchStrategy_ = new BlisBranchStrategyMaxInf(this);
    }
    else if (brStrategy == BlisBranchingStrategyPseudoCost) {
        // Pseudocost
        branchStrategy_ = new BlisBranchStrategyPseudo(this, 1);
    }
    else if (brStrategy == BlisBranchingStrategyReliability) {
        // Relibility
        branchStrategy_ = new BlisBranchStrategyRel(this, relibility);
    }
    else if (brStrategy == BlisBranchingStrategyStrong) {
        // Strong
        branchStrategy_ = new BlisBranchStrategyStrong(this);
    }
    else {
        throw CoinError("Unknown branch strategy.", "setupSelf","BlisModel");
    }

    brStrategy = BlisPar_->entry(BlisParams::branchStrategyRampUp);

    if (brStrategy == BlisBranchingStrategyMaxInfeasibility) {
        // Max inf
      rampUpBranchStrategy_ = new BlisBranchStrategyMaxInf(this);
    }
    else if (brStrategy == BlisBranchingStrategyPseudoCost) {
        // Pseudocost
        rampUpBranchStrategy_ = new BlisBranchStrategyPseudo(this, 1);
    }
    else if (brStrategy == BlisBranchingStrategyReliability) {
        // Relibility
        rampUpBranchStrategy_ = new BlisBranchStrategyRel(this, relibility);
    }
    else if (brStrategy == BlisBranchingStrategyStrong) {
        // Strong
        rampUpBranchStrategy_ = new BlisBranchStrategyStrong(this);
    }
    else {
        throw CoinError("Unknown branch strategy.", "setupSelf","BlisModel");
    }

    //------------------------------------------------------
    // Add heuristics.
    //------------------------------------------------------

    heurStrategy_ = 
       static_cast<BlisHeurStrategy> (BlisPar_->entry(BlisParams::heurStrategy));
    heurCallFrequency_ = 
       static_cast<BlisHeurStrategy> (BlisPar_->entry(BlisParams::heurCallFrequency));
    BlisHeurStrategy heurRoundStrategy = 
       static_cast<BlisHeurStrategy> (BlisPar_->entry(BlisParams::heurRoundStrategy)); 
    int callFreq;

    if (heurRoundStrategy == BlisHeurStrategyNotSet) {
        heurRoundStrategy = heurStrategy_;
	callFreq = heurCallFrequency_;
    }
    if (heurRoundStrategy != BlisHeurStrategyNone && 
	heurRoundStrategy != BlisHeurStrategyNotSet) {
        // Add rounding heuristic
        BlisHeurRound *heurRound = new BlisHeurRound(this, 
                                                     "Rounding",
                                                     heurRoundStrategy,
						     callFreq);
        addHeuristic(heurRound);
    }

    // Adjust heurStrategy
    for (j = 0; j < numHeuristics_; ++j) {
        if (heuristics_[j]->strategy() != BlisHeurStrategyNone) {
            // Doesn't matter what's the strategy, we just want to 
            // call heuristics.
            heurStrategy_ = BlisHeurStrategyAuto;
            break;
        }
    }

    //------------------------------------------------------
    // Cut generators settings.
    //------------------------------------------------------

    // Compute dense cutoff.
    
    const CoinPackedMatrix * rowMatrix = lpSolver_->getMatrixByRow();
    const int * rowLen = rowMatrix->getVectorLengths();
    double maxLen = 0.0, minLen = ALPS_DBL_MAX, sumLen = 0.0;
    double aveLen, diffLen, stdLen;
    double denseConFactor = BlisPar_->entry(BlisParams::denseConFactor);

    for (j = 0; j < numRows_; ++j) {
	if (rowLen[j] > maxLen) maxLen = rowLen[j];
	if (rowLen[j] < minLen) minLen = rowLen[j];
	sumLen += rowLen[j];
    }
    assert(numRows_ > 0);
    aveLen = sumLen / numRows_;
    sumLen = 0.0;

    for (j = 0; j < numRows_; ++j) {
	diffLen = rowLen[j] - aveLen;
	diffLen *= diffLen;
	sumLen += diffLen;
    }
    stdLen = sqrt(sumLen/numRows_);
    if (denseConFactor > 10e5) {
        denseConCutoff_ = ALPS_INT_MAX;
    }
    else {
        denseConCutoff_ = static_cast<int>(aveLen + denseConFactor*stdLen);    
        denseConCutoff_ = ALPS_MIN(numCols_/2, denseConCutoff_);
        denseConCutoff_ = ALPS_MAX(10, denseConCutoff_);
    }
    
    
#ifdef BLIS_DEBUG
    std::cout << "aveLen=" << aveLen << ", minLen=" << minLen
	      << ", maxLen=" << maxLen << ", stdLen=" << stdLen
	      << ", denseConCutoff_=" << denseConCutoff_ << std::endl;
#endif
    
    // NOTE: maxNumCons is valid only for automatic strategy.
    double cutFactor = BlisPar_->entry(BlisParams::cutFactor);

    if (cutFactor > 1.0e19) {
        // FIXME: use vector
        maxNumCons_ = 100000;
    }
    else {
        maxNumCons_ = (int)((cutFactor - 1.0) * numRows_);
    }
    
    constraintPool_ = new BcpsConstraintPool();
    oldConstraints_ = new BlisConstraint* [maxNumCons_];
    oldConstraintsSize_ = maxNumCons_;
    
    cutStrategy_ = static_cast<BlisCutStrategy> 
       (BlisPar_->entry(BlisParams::cutStrategy)); 

#ifdef BLIS_DEBUG
    std::cout << "cutStrategy_ = " << cutStrategy_ << std::endl;
#endif

    BlisCutStrategy cliqueStrategy = static_cast<BlisCutStrategy> 
       (BlisPar_->entry(BlisParams::cutCliqueStrategy));
    BlisCutStrategy fCoverStrategy = static_cast<BlisCutStrategy>
       (BlisPar_->entry(BlisParams::cutFlowCoverStrategy));
    BlisCutStrategy gomoryStrategy = static_cast<BlisCutStrategy> 
       (BlisPar_->entry(BlisParams::cutGomoryStrategy)); 
    BlisCutStrategy knapStrategy = static_cast<BlisCutStrategy>
       (BlisPar_->entry(BlisParams::cutKnapsackStrategy)); 
    BlisCutStrategy mirStrategy = static_cast<BlisCutStrategy> 
       (BlisPar_->entry(BlisParams::cutMirStrategy)); 
    BlisCutStrategy oddHoleStrategy = static_cast<BlisCutStrategy> 
       (BlisPar_->entry(BlisParams::cutOddHoleStrategy));
    BlisCutStrategy probeStrategy = static_cast<BlisCutStrategy> 
       (BlisPar_->entry(BlisParams::cutProbingStrategy));
    BlisCutStrategy twoMirStrategy = static_cast<BlisCutStrategy> 
       (BlisPar_->entry(BlisParams::cutTwoMirStrategy)); 

    int cliqueFreq = BlisPar_->entry(BlisParams::cutCliqueFreq);
    int fCoverFreq = BlisPar_->entry(BlisParams::cutFlowCoverFreq);
    int gomoryFreq = BlisPar_->entry(BlisParams::cutGomoryFreq); 
    int knapFreq = BlisPar_->entry(BlisParams::cutKnapsackFreq); 
    int mirFreq = BlisPar_->entry(BlisParams::cutMirFreq); 
    int oddHoleFreq = BlisPar_->entry(BlisParams::cutOddHoleFreq);
    int probeFreq = BlisPar_->entry(BlisParams::cutProbingFreq);
    int twoMirFreq = BlisPar_->entry(BlisParams::cutTwoMirFreq); 

    //------------------------------------------------------
    // Add cut generators.
    //------------------------------------------------------

    if (probeStrategy == BlisCutStrategyNotSet) {
        // Disable by default
        probeStrategy = BlisCutStrategyNone;
    }
    if (probeStrategy != BlisCutStrategyNone) {
        CglProbing *probing = new CglProbing;
        probing->setUsingObjective(true);
        probing->setMaxPass(1);
        probing->setMaxPassRoot(5);
        // Number of unsatisfied variables to look at
        probing->setMaxProbe(10);
        probing->setMaxProbeRoot(1000);
        // How far to follow the consequences
        probing->setMaxLook(50);
        probing->setMaxLookRoot(500);
        // Only look at rows with fewer than this number of elements
        probing->setMaxElements(200);
        probing->setRowCuts(3);
        if (probeStrategy == BlisCutStrategyPeriodic) {
	    // User freq
            addCutGenerator(probing, "Probing", probeStrategy, probeFreq); 
        }
        else {
            addCutGenerator(probing, "Probing", probeStrategy);
        }
    }

    if (cliqueStrategy == BlisCutStrategyNotSet) {
        // Only at root by default
        cliqueStrategy = BlisCutStrategyRoot;
    }
    if (cliqueStrategy != BlisCutStrategyNone) {
        CglClique *cliqueCut = new CglClique ;
        cliqueCut->setStarCliqueReport(false);
        cliqueCut->setRowCliqueReport(false);
        if (cliqueStrategy == BlisCutStrategyPeriodic) {
            addCutGenerator(cliqueCut, "Clique", cliqueStrategy, cliqueFreq);
        }
        else {
            addCutGenerator(cliqueCut, "Clique", cliqueStrategy);
        }
    }

    if (oddHoleStrategy == BlisCutStrategyNotSet) {
        // Disable by default
        oddHoleStrategy = BlisCutStrategyNone;
    }
    if (oddHoleStrategy != BlisCutStrategyNone) {
        CglOddHole *oldHoleCut = new CglOddHole;
        oldHoleCut->setMinimumViolation(0.005);
        oldHoleCut->setMinimumViolationPer(0.00002);
        // try larger limit
        oldHoleCut->setMaximumEntries(200);
        if (oddHoleStrategy == BlisCutStrategyPeriodic) {
            addCutGenerator(oldHoleCut, "OddHole",oddHoleStrategy,oddHoleFreq);
        }
        else {
             addCutGenerator(oldHoleCut, "OddHole", oddHoleStrategy);
        }
    }

    if (fCoverStrategy == BlisCutStrategyNotSet) {
         fCoverStrategy = cutStrategy_;
    }
    if (fCoverStrategy != BlisCutStrategyNone) {
        CglFlowCover *flowGen = new CglFlowCover;
        if (fCoverStrategy == BlisCutStrategyPeriodic) {
            addCutGenerator(flowGen, "Flow Cover", fCoverStrategy, fCoverFreq);
        }
        else {
            addCutGenerator(flowGen, "Flow Cover", fCoverStrategy);
        }
    }

    if (knapStrategy == BlisCutStrategyNotSet) {
        // Only at root by default
        knapStrategy = BlisCutStrategyRoot;
    }
    if (knapStrategy != BlisCutStrategyNone) {
        CglKnapsackCover *knapCut = new CglKnapsackCover;
        if (knapStrategy == BlisCutStrategyPeriodic) {
            addCutGenerator(knapCut, "Knapsack", knapStrategy, knapFreq);
        }
        else {
             addCutGenerator(knapCut, "Knapsack", knapStrategy);
        }
    }

    if (mirStrategy == BlisCutStrategyNotSet) {
        // Disable by default
        mirStrategy = BlisCutStrategyNone;
    }
    if (mirStrategy != BlisCutStrategyNone) {
        CglMixedIntegerRounding2 *mixedGen = new CglMixedIntegerRounding2;
        if (mirStrategy == BlisCutStrategyPeriodic) {
            addCutGenerator(mixedGen, "MIR", mirStrategy, mirFreq);
        }
        else {
            addCutGenerator(mixedGen, "MIR", mirStrategy);
        }
    }

    if (gomoryStrategy == BlisCutStrategyNotSet) {
        // Only at root by default
        gomoryStrategy = BlisCutStrategyRoot;
    }
    if (gomoryStrategy != BlisCutStrategyNone) {
        CglGomory *gomoryCut = new CglGomory;
        // try larger limit
        gomoryCut->setLimit(300);
        if (gomoryStrategy == BlisCutStrategyPeriodic) {
            addCutGenerator(gomoryCut, "Gomory", gomoryStrategy, gomoryFreq);
        }
        else {
            addCutGenerator(gomoryCut, "Gomory", gomoryStrategy);
        }
    }

    if (twoMirStrategy == BlisCutStrategyNotSet) {
        // Disable by default
        twoMirStrategy = BlisCutStrategyNone;
    }
    if (twoMirStrategy != BlisCutStrategyNone) {
        CglTwomir *twoMirCut =  new CglTwomir;
        if (twoMirStrategy == BlisCutStrategyPeriodic) {
            addCutGenerator(twoMirCut, "Two MIR", twoMirStrategy, twoMirFreq);
        }
        else  {
            addCutGenerator(twoMirCut, "Two MIR", twoMirStrategy);
        }
    }
    
    // Random vector
    // TODO, if generating variable, then need more space.
    conRandoms_ = new double [numCols_];
    double deseed = 12345678.0;
    double DSEED2 = 2147483647.0;
        
    for (j = 0; j < numCols_; ++j) {
        deseed *= 16807.0;
        int jseed = (int) (deseed / DSEED2);
        deseed -= (double) jseed * DSEED2;
        double random = deseed / DSEED2;
        conRandoms_[j] = random;
        //std::cout << "conRandoms_[" << j << "]="
        //      <<conRandoms_[j]<< std::endl;
    }

    // Adjust cutstrategy
    if (numCutGenerators_ > 0) {
	BlisCutStrategy strategy0 = cutGenerators(0)->strategy();
	BlisCutStrategy strategy1;
	for (j = 1; j < numCutGenerators_; ++j) {
	    strategy1 = cutGenerators(j)->strategy();
	    if (strategy1 != strategy0) {
	       // Different for each cut
	       cutStrategy_ =  BlisCutStrategyMultiple;
	       break;
	    }
	}
	if (j == numCutGenerators_) {
	    // Cut has same strategy.
	    cutStrategy_ = strategy0;
	}
        else {
            // Assume to generate cons at each node.
            cutStrategy_ = BlisCutStrategyPeriodic;
	    cutGenerationFrequency_ = 1;
        }
    }
     
#if 0
    std::cout << "AFTER: cutStrategy_ = " << cutStrategy_ << std::endl;
#endif

    return true;
}

//############################################################################ 

void 
BlisModel::preprocess()
{
    int numPasses = 10;
    double feaTol = 1.0e-6;
    bool keepIntegers = true;
    char *prohibited = 0;

    bool doPresolve = BlisPar_->entry(BlisParams::presolve);
    
    //std::cout << "Presolve = "<< doPresolve << std::endl;
    
    doPresolve = false;
    
    if (doPresolve) {
	presolvedLpSolver_ = presolve_->preprocess(*lpSolver_,
						   feaTol,
						   keepIntegers,
						   numPasses,
						   prohibited);
	lpSolver_ = presolvedLpSolver_->clone();
    }
    else {
	lpSolver_ = origLpSolver_;
    }
}

//############################################################################ 

void 
BlisModel::postprocess()
{
    
}

//#############################################################################

int 
BlisModel::storeSolution(BlisSolutionType how, BlisSolution* sol)
{
    double quality = sol->getQuality();

    // Update cutoff and lp cutoff.
    setCutoff(quality);

    ++numSolutions_;

    // Store in pool, assume minimization.
    getKnowledgeBroker()->addKnowledge(AlpsKnowledgeTypeSolution, 
                                       sol,
                                       objSense_ * sol->getQuality()); 

    // Record how the solution was found
    switch (how) {
    case BlisSolutionTypeBounding:
#ifdef BLIS_DEBUG
        std::cout << "Rounding heuristics found a better solution" 
                  <<", old cutoff = " << cutoff 
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    case BlisSolutionTypeBranching:
#ifdef BLIS_DEBUG
        std::cout << "Branching found a better solution" 
                  <<", old cutoff = " << cutoff 
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    case BlisSolutionTypeDiving:
        ++numHeurSolutions_;
#ifdef BLIS_DEBUG
        std::cout << "Branching found a better solution" 
                  <<", old cutoff = " << cutoff 
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    case BlisSolutionTypeRounding:
        ++numHeurSolutions_;
#ifdef BLIS_DEBUG
        std::cout << "Rounding heuristics found a better solution" 
                  <<", old cutoff = " << cutoff 
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    case BlisSolutionTypeStrong:
#ifdef BLIS_DEBUG
        std::cout << "Strong branching found a better solution" 
                  <<", old cutoff = " << cutoff 
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    default:
#ifdef BLIS_DEBUG
        std::cout << "Nowhere found a better solution" 
                  <<", old cutoff = " << cutoff 
                  << ", new cutoff = " << getCutoff()  << std::endl;
#endif
        break;
    }
    
    return AlpsReturnStatusOk;
}

//############################################################################ 

void 
BlisModel::createIntgerObjects(bool startAgain)
{
    assert(lpSolver_);

    if (numIntObjects_ && !startAgain && objects_) return;

    int iCol;    
    int numCols = getNumCols();
    
    const double *colLB = lpSolver_->getColLower();
    const double *colUB = lpSolver_->getColUpper();
    BlisObjectInt *intObject = NULL;

    if (intColIndices_) {
        delete [] intColIndices_;
    }
    numIntObjects_ = 0;

    for (iCol = 0; iCol < numCols; ++iCol) {
	if (lpSolver_->isInteger(iCol)) ++numIntObjects_;
    }

#if 0
    std::cout << "==== findInteger: numIntObjects_ = " 
              << numIntObjects_ << std::endl;
#endif

    double weight = BlisPar_->entry(BlisParams::pseudoWeight);
    
    int numObjects = 0;
    int iObject;
    BcpsObject ** oldObject = objects_;    

    for (iObject = 0; iObject < numObjects_; ++iObject) {
	BlisObjectInt * obj =
	    dynamic_cast <BlisObjectInt *>(oldObject[iObject]) ;
        
	if (obj) {
	    delete oldObject[iObject];
        }
	else {
	    oldObject[numObjects++] = oldObject[iObject];
        }
    }

    if (!intObjIndices_) {
        intObjIndices_ = new int [numCols];
        memset(intObjIndices_, 0, sizeof(int) * numCols);
    }
    
    objects_ = new BcpsObject * [(numIntObjects_ + numObjects)];
    intColIndices_ = new int [numIntObjects_];
    numObjects_ = numIntObjects_ + numObjects;

    // Walk the variables again, filling in the indices and creating objects 
    // for the integer variables. Initially, the objects hold the indices,
    // variable bounds and pseudocost.
    numIntObjects_ = 0;
    for (iCol = 0; iCol < numCols; ++iCol) {
	if(lpSolver_->isInteger(iCol)) {
            
	    intObject = new BlisObjectInt(numIntObjects_,
                                          iCol,
                                          colLB[iCol],
                                          colUB[iCol]);
            intObject->pseudocost().setWeight(weight);

	    intObjIndices_[iCol] = numIntObjects_;
            objects_[numIntObjects_] = intObject;
	    intColIndices_[numIntObjects_++] = iCol;
	}
    }

    if (numIntObjects_) {
        sharedObjectMark_ = new char [numIntObjects_];
        memset(sharedObjectMark_, 0, sizeof(char) * numIntObjects_);
    }
    
    // Now append other objects
    memcpy(objects_ + numIntObjects_, oldObject, numObjects*sizeof(BcpsObject *));

    // Delete old array (just array)
    delete [] oldObject;
}

//#############################################################################

bool
BlisModel::resolve()
{
    lpSolver_->resolve();
    numIterations_ += lpSolver_->getIterationCount();
    bool feasible = (lpSolver_->isProvenOptimal() &&
                     !lpSolver_->isDualObjectiveLimitReached());
    
    return feasible; 
}

//#############################################################################

// Delete all object information
void 
BlisModel::deleteObjects()
{
    delete [] priority_;
    priority_ = NULL;
    int i;
    for (i = 0; i < numObjects_; ++i) delete objects_[i];
    delete [] objects_;
    objects_ = NULL;
    numObjects_ = 0;
    createIntgerObjects(true);
}

//#############################################################################

BlisModel::~BlisModel()
{
    gutsOfDestructor();
}

//#############################################################################

void 
BlisModel::gutsOfDestructor()
{
    int i;
    bool doPresolve = BlisPar_->entry(BlisParams::presolve);
    

//    delete [] savedLpSolution_;
//    savedLpSolution_ = NULL;

    delete [] intObjIndices_;
    intObjIndices_ = NULL;

    delete [] intColIndices_;
    intColIndices_ = NULL;

    for (i = 0; i < numObjects_; ++i) delete objects_[i];
    delete [] objects_;
    objects_ = NULL;

    delete [] priority_;
    priority_ = NULL;

    delete [] colType_;
    delete colMatrix_;

    delete [] varLB_;
    delete [] varUB_;

    delete [] conLB_;
    delete [] conUB_;

    delete [] startVarLB_;
    delete [] startVarUB_;

    delete [] startConLB_;
    delete [] startConUB_;

    delete [] tempVarLBPos_;
    delete [] tempVarUBPos_;

    delete [] tempConLBPos_;
    delete [] tempConUBPos_;

    delete [] objCoef_;
    delete [] incumbent_;

    delete presolve_;
    
    if (numHeuristics_ > 0) {
#ifdef BLIS_DEBUG
        std::cout << "MODEL: distructor: numHeuristics =  " 
                  << numHeuristics_ << std::endl;
#endif
        BlisHeuristic *tempH = NULL;
        
        for (i = 0; i < numHeuristics_; ++i) {
            tempH = heuristics_[i];
            delete tempH;
        }
        
        delete [] heuristics_;
        heuristics_ = NULL;
    }

    if (generators_ != NULL) {
#ifdef BLIS_DEBUG
        std::cout << "MODEL: distructor: numCutGenerators = " 
                  << numCutGenerators_ << std::endl;
#endif
        BlisConGenerator *temp = NULL;
        for (i = 0; i < numCutGenerators_; ++i) {
            temp = generators_[i];
            delete temp;
        }
        delete [] generators_;
        generators_ = NULL;
    }

    delete constraintPool_;
    delete [] oldConstraints_;
    delete branchStrategy_;
    delete rampUpBranchStrategy_;

    delete [] conRandoms_;
    
    delete BlisPar_;
    delete blisMessageHandler_;
    if (doPresolve) {
	delete presolvedLpSolver_;
	delete lpSolver_;
    }

    delete [] sharedObjectMark_; 
    sharedObjectMark_ = NULL;
}

//#############################################################################

BlisSolution *
BlisModel::feasibleSolution(int & numIntegerInfs, int & numObjectInfs)
{
    int preferredWay, j;
    int numUnsatisfied = 0;

    bool userFeasible = true;
    double sumUnsatisfied = 0.0;
    
    BlisSolution* sol = NULL;

    for (j = 0; j < numIntObjects_; ++j) {
	const BcpsObject * object = objects_[j];
	
	double infeasibility = object->infeasibility(this, preferredWay);
	if (infeasibility) {
	    assert (infeasibility>0);
	    numUnsatisfied++;
	    sumUnsatisfied += infeasibility;
	}
    }
    numIntegerInfs = numUnsatisfied;
    for (; j < numObjects_; ++j) {
	const BcpsObject * object = objects_[j];
	double infeasibility = object->infeasibility(this, preferredWay);
	if (infeasibility) {
	    assert (infeasibility > 0);
	    numUnsatisfied++;
	    sumUnsatisfied += infeasibility;
	}
    }

    numObjectInfs = numUnsatisfied - numIntegerInfs;

    //printf("numUnsatisfied = %d\n",numUnsatisfied);

    if (!numUnsatisfied) {
        sol = userFeasibleSolution(userFeasible);
	if (!sol && userFeasible) {
	    // User doesn't provide feasibility check.
	    numBranchResolve_ = 10;
	}
    }
    
    if (!numUnsatisfied && userFeasible && !sol) {
        // User doesn't overload feasible solution function.
        sol = new BlisSolution(getNumCols(), 
                               getLpSolution(), 
                               objSense_ * getLpObjValue());
    }
    
    return sol;
}

//#############################################################################
/*
  Set branching priorities.

  Setting integer priorities looks pretty robust; the call to createIntgerObjects
  makes sure that integer objects are in place. Setting priorities for
  other objects is entirely dependent on their existence, and the routine may
  quietly fail in several directions.
*/
void 
BlisModel::passInPriorities (const int * priorities,
			     bool ifObject, 
			     int defaultValue)
{

    // FIXME: not completed.
    int i;

    createIntgerObjects(false);

    if (!priority_) {
	priority_ = new int[numObjects_];
	for (i = 0; i < numObjects_; ++i) {
            priority_[i] = defaultValue;
        }
    }

    if (priorities) {
	if (ifObject) {
	    memcpy(priority_ + numIntObjects_, 
                   priorities,
		   (numObjects_ - numIntObjects_) * sizeof(int));
        }
	else {
	    memcpy(priority_, priorities, numIntObjects_ * sizeof(int));
        }
    }
}

//#############################################################################

AlpsTreeNode * 
BlisModel::createRoot() {
    
    //-------------------------------------------------------------
    // NOTE: Root will be deleted by ALPS. Root is an explicit node.
    //-------------------------------------------------------------
    
    BlisTreeNode* root = new BlisTreeNode;
    BlisNodeDesc* desc = new BlisNodeDesc(this);
    root->setDesc(desc);

    //-------------------------------------------------------------
    // NOTE: Although original data are stored in model when reading. 
    //   Root desc still store a full copy of col and row bounds when creating.
    //   It will store soft differences after finding a branching object. 
    //   The soft difference are due to reduced cost fixing and probing.
    //   Also the added cols and rows will be stored.
    //-------------------------------------------------------------
    int k;
    
    std::vector<BcpsVariable *> vars = getVariables();
    std::vector<BcpsConstraint *> cons = getConstraints();

    int numVars = static_cast<int> (vars.size());
    int numCons = static_cast<int> (cons.size());

#ifdef BLIS_DEBUG
    std::cout << "BLIS: createRoot(): numVars=" << numVars
              << ", numCons=" << numCons 
              << "; numCoreVariables_=" << numCoreVariables_
              << ", numCoreConstraints_=" << numCoreConstraints_ << std::endl;
#endif
    
    int *varIndices1 = new int [numVars];
    int *varIndices2 = new int [numVars];
    int *varIndices3 = NULL;
    int *varIndices4 = NULL;
    double *vlhe = new double [numVars];
    double *vuhe = new double [numVars];
    double *vlse = NULL;
    double *vuse = NULL;
    
    int *conIndices1 = new int [numCons];
    int *conIndices2 = new int [numCons];
    int *conIndices3 = NULL;
    int *conIndices4 = NULL;
    double *clhe = new double [numCons];
    double *cuhe = new double [numCons];
    double *clse = NULL;
    double *cuse = NULL;

    //-------------------------------------------------------------
    // Get var bounds and indices.
    //-------------------------------------------------------------
    
    for (k = 0; k < numVars; ++k) {
        vlhe[k] = vars[k]->getLbHard();
        vuhe[k] = vars[k]->getUbHard();
        varIndices1[k] = k;
        varIndices2[k] = k;
        
#ifdef BLIS_DEBUG_MORE
        std::cout << "BLIS: createRoot(): var "<< k << ": hard: lb=" << vlhe[k]
                  << ", ub=" << vuhe[k] << std::endl;
#endif  
    }

    //-------------------------------------------------------------  
    // Get con bounds and indices.
    //-------------------------------------------------------------
    
    for (k = 0; k < numCons; ++k) {
        clhe[k] = cons[k]->getLbHard();
        cuhe[k] = cons[k]->getUbHard();
        conIndices1[k] = k;
        conIndices2[k] = k;
    }
    
    int *tempInd = NULL;
    BcpsObject **tempObj = NULL;
  
    desc->assignVars(0 /*numRem*/, tempInd,
                     0 /*numAdd*/, tempObj,
                     false, numVars, varIndices1, vlhe, /*Var hard lb*/
                     false, numVars, varIndices2, vuhe, /*Var hard ub*/
                     false, 0, varIndices3, vlse,       /*Var soft lb*/
                     false, 0, varIndices4, vuse);      /*Var soft ub*/
    desc->assignCons(0 /*numRem*/, tempInd,
                     0 /*numAdd*/, tempObj,
                     false, numCons, conIndices1, clhe, /*Con hard lb*/
                     false, numCons, conIndices2, cuhe, /*Con hard ub*/
                     false, 0,conIndices3,clse,         /*Con soft lb*/
                     false, 0,conIndices4,cuse);        /*Con soft ub*/
    
    //-------------------------------------------------------------  
    // Mark it as an explicit node.
    //-------------------------------------------------------------

    root->setExplicit(1);
    
    return root;
}

//#############################################################################

void 
BlisModel::addHeuristic(BlisHeuristic * heuristic)
{
    BlisHeuristic ** temp = heuristics_;
    heuristics_ = new BlisHeuristic * [numHeuristics_ + 1];

    memcpy(heuristics_, temp, numHeuristics_ * sizeof(BlisHeuristic *));
    delete [] temp;
    
    heuristics_[numHeuristics_++] = heuristic;
}

//#############################################################################

void 
BlisModel::addCutGenerator(CglCutGenerator * generator,
			   const char * name,
			   BlisCutStrategy strategy,
			   int freq,
			   bool normal, 
			   bool atSolution,
			   bool whenInfeasible)
{
    BlisConGenerator ** temp = generators_;

    generators_ = new BlisConGenerator * [(numCutGenerators_ + 1)];
    memcpy(generators_, temp, numCutGenerators_ * sizeof(BlisConGenerator *));
    
    generators_[numCutGenerators_++] = 
        new BlisConGenerator(this, generator, name, strategy, freq,
                             normal, atSolution, whenInfeasible);
    delete [] temp;
    temp = NULL;
}

//#############################################################################

void 
BlisModel::addCutGenerator(BlisConGenerator * generator)
{
    BlisConGenerator ** temp = generators_;

    generators_ = new BlisConGenerator * [(numCutGenerators_ + 1)];
    memcpy(generators_, temp, numCutGenerators_ * sizeof(BlisConGenerator *));
    
    generators_[numCutGenerators_++] = generator;

    delete [] temp;
    temp = NULL;
}

//#############################################################################

AlpsReturnStatus 
BlisModel::encodeBlis(AlpsEncoded *encoded) const
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    
    BlisPar_->pack(*encoded);

    encoded->writeRep(objSense_);
    
    return status;
}

//#############################################################################

AlpsReturnStatus 
BlisModel::decodeBlis(AlpsEncoded &encoded)
{
    AlpsReturnStatus status = AlpsReturnStatusOk;
    
    BlisPar_->unpack(encoded);

    encoded.readRep(objSense_);

    //------------------------------------------------------
    // Load problem to LP solver.
    //------------------------------------------------------

    std::vector<BlisVariable *> vars;
    std::vector<BlisConstraint *> cons;

    int k;
    int size = static_cast<int> (variables_.size());
    for (k = 0; k < size; ++k) {
        BlisVariable * aVar = dynamic_cast<BlisVariable *>(variables_[k]);
        vars.push_back(aVar);
    }
    
    size = static_cast<int> (constraints_.size()); 
    for (k = 0; k < size; ++k) {
        BlisConstraint * aCon = dynamic_cast<BlisConstraint *>(constraints_[k]);
        cons.push_back(aCon);
    }
    
    variables_.clear();
    constraints_.clear();

    // Import model from vars and cons.
    importModel(vars, cons);
    
    return status;
}

//#############################################################################

AlpsEncoded* 
BlisModel::encode() const 
{ 
    AlpsReturnStatus status = AlpsReturnStatusOk;

    // NOTE: "AlpsKnowledgeTypeModel" is the type name.
    AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeModel);

    status = encodeAlps(encoded);
    status = encodeBcps(encoded);
    status = encodeBlis(encoded);

    return encoded;
}

//#############################################################################

void
BlisModel::decodeToSelf(AlpsEncoded& encoded) 
{
    AlpsReturnStatus status = AlpsReturnStatusOk;

    status = decodeAlps(encoded);
    status = decodeBcps(encoded);
    status = decodeBlis(encoded);
}

//#############################################################################

AlpsEncoded* 
BlisModel::encodeKnowlegeShared()
{
    AlpsEncoded* encoded = 0;

    int k = 0;
    int numShared = 0;
    int frequency = -1, depth = -1;
    
    bool sharePseudo = true;
    bool share = false;

#if 0
    bool shareCon = false;
    bool shareVar = false;
#endif

    int phase = broker_->getPhase();
    
    if (phase == AlpsPhaseRampup) {
        sharePseudo = BlisPar_->entry(BlisParams::sharePseudocostRampUp);
    }
    else if (phase == AlpsPhaseSearch) {
        sharePseudo = BlisPar_->entry(BlisParams::sharePseudocostSearch);
        if (sharePseudo) {
            // Depth and frequency
            depth = BlisPar_->entry(BlisParams::sharePcostDepth);
            frequency =  BlisPar_->entry(BlisParams::sharePcostFrequency);
            if ( (numNodes_ % frequency != 0) || 
                 (broker_->getTreeDepth() >  depth) ) {
                sharePseudo = false;
            }
        }
    }

    if (sharePseudo) {
        for (k = 0; k < numIntObjects_; ++k) {
            if (sharedObjectMark_[k]){
                ++numShared;
            }
        }
        if (numShared) share = true;
    }
    // TODO: cuts, etc.

#if 0
    std::cout << "++++ sharePseudo =  " << sharePseudo
              << ", shareCon = " << shareCon
              << ", shareVar = " << shareVar
              << ", numShared = " << numShared 
              << ", numNodes_ = " << numNodes_
              << ", frequency = " << frequency
              << ", depth = " << depth
              << ", treeDepth = " << broker_->getTreeDepth()
              << ", share = " << share << std::endl;
#endif        

    if (share) {
        // NOTE: "AlpsKnowledgeTypeModelGen" is the type name. We don't need to
        //       register it since AlpsKnowledgeTypeModel is registered.
        
        BlisObjectInt *intObj = NULL;
        
        if (numShared > 0) {
            encoded = new AlpsEncoded(AlpsKnowledgeTypeModelGen);
            // Record how many can be shared.
            encoded->writeRep(numShared);
            for (k = 0; k < numIntObjects_; ++k) {
                if (sharedObjectMark_[k]) {
                    // Recored which object.
                    encoded->writeRep(k);
                    intObj = dynamic_cast<BlisObjectInt*>(objects_[k]);
                    (intObj->pseudocost()).encodeTo(encoded);
                }    
            }
            
            // Clear the mark for next round of sharing.
            clearSharedObjectMark();
        }
        else {
            encoded->writeRep(numShared);
        }
        
	// Make sure don't exceed large size.
	assert(encoded->size() < broker_->getLargeSize());
    }

    return encoded;
}

//#############################################################################

void 
BlisModel::decodeKnowledgeShared(AlpsEncoded& encoded)
{
    int k, objIndex, size = 0;
    
    // Encode and store pseudocost
    BlisObjectInt *intObj = NULL;
    encoded.readRep(size);
    for (k = 0; k < size; ++k) {
        encoded.readRep(objIndex);
        intObj = dynamic_cast<BlisObjectInt *>(objects_[objIndex]);
        (intObj->pseudocost()).decodeFrom(encoded);
    }

#if 0        
    // Encode and store generated constraints that should be shared.
    encoded.readRep(size);
    for (k = 0; k < size; ++k) {
        // TODO
        assert(0);
    }
    
    // Encode and store variables that should be shared.
    encoded.readRep(size);
    for (k = 0; k < size; ++k) {
        // TODO
        assert(0);
    }
#endif
}

//#############################################################################

/** Register knowledge. */
void 
BlisModel::registerKnowledge() {
    // Register model, solution, and tree node
    assert(broker_);
    broker_->registerClass(AlpsKnowledgeTypeModel, new BlisModel);
    if (broker_->getMsgLevel() > 5) {
	std::cout << "BLIS: Register Alps model." << std::endl;
    }
    
    broker_->registerClass(AlpsKnowledgeTypeNode, new BlisTreeNode(this));
    if (broker_->getMsgLevel() > 5) {
	std::cout << "BLIS: Register Alps node." << std::endl;
    }
    
    broker_->registerClass(AlpsKnowledgeTypeSolution, new BlisSolution);
    if (broker_->getMsgLevel() > 5) {
	std::cout << "BLIS: Register Alps solution." << std::endl;
    }
    
    broker_->registerClass(BcpsKnowledgeTypeConstraint, new BlisConstraint);
    if (broker_->getMsgLevel() > 5) {
	std::cout << "BLIS: Register Bcps constraint." << std::endl;
    }
    
    broker_->registerClass(BcpsKnowledgeTypeVariable, new BlisVariable);
    if (broker_->getMsgLevel() > 5) {
	std::cout << "BLIS: Register Bcps variable." << std::endl;
    }
}

//#############################################################################

/** Log of specific models. */
void 
BlisModel::modelLog() 
{
    if (processType_ != AlpsProcessTypeMaster) return;
    
    int logFileLevel = AlpsPar_->entry(AlpsParams::logFileLevel);
    int msgLevel = AlpsPar_->entry(AlpsParams::msgLevel);
    if (logFileLevel > 0) {
	std::string logfile = AlpsPar_->entry(AlpsParams::logFile);
	std::ofstream logFout(logfile.c_str(), std::ofstream::app);
	writeParameters(logFout);
    }

    if (msgLevel > 0) {
        int k;
        for (k = 0; k < numCutGenerators_; ++k) {
            if (cutGenerators(k)->calls() > 0) {
                blisMessageHandler()->message(BLIS_CUT_STAT_FINAL,
                                              blisMessages())
                    << cutGenerators(k)->name()
                    << cutGenerators(k)->calls()
                    << cutGenerators(k)->numConsGenerated()
                    << cutGenerators(k)->time()
                    << cutGenerators(k)->strategy()
                    << CoinMessageEol;
            }
        }
        for (k = 0; k < numHeuristics_; ++k) {
            if (heuristics(k)->calls() > 0) {
                blisMessageHandler()->message(BLIS_HEUR_STAT_FINAL,
                                              blisMessages())
                    << heuristics(k)->name()
                    << heuristics(k)->calls()
                    << heuristics(k)->numSolutions()
                    << heuristics(k)->time()
                    << heuristics(k)->strategy()
                    << CoinMessageEol;
            }   
        }
    }    
}

//#############################################################################

void 
BlisModel::nodeLog(AlpsTreeNode *node, bool force) 
{
    int nodeInterval = 
	broker_->getModel()->AlpsPar()->entry(AlpsParams::nodeLogInterval);
    
    int numNodesProcessed = broker_->getNumNodesProcessed();
    int numNodesLeft = broker_->updateNumNodesLeft();
    int msgLevel = broker_->getMsgLevel();
    bool printLog = false;
    
    int numCols = getNumCols();
    int numRows = getNumRows();

    //std::cout << "nodeInterval = " << nodeInterval << std::endl;
    
    AlpsTreeNode *bestNode = NULL;

    if ((msgLevel > 1) && (force||(numNodesProcessed % nodeInterval == 0))) {
        printLog = true;
    }

    if (broker_->getProcType() != AlpsProcessTypeMaster) {
        printLog = false;
    }
    
    if (msgLevel > 200) {
        printLog = true;
    }

#if 0
    std::cout << "==== Process " << broker_->getProcRank()
              << ": printLog = " << printLog 
              << ", msgLevel = " << msgLevel 
              << ", proc type = " << broker_->getProcType()
              << std::endl;
#endif
    
    if (printLog) {
        double feasBound = ALPS_OBJ_MAX;
	double relBound = ALPS_OBJ_MAX;
	double gap = ALPS_OBJ_MAX;
	double gapVal = ALPS_OBJ_MAX;
	
        if (broker_->getNumKnowledges(AlpsKnowledgeTypeSolution) > 0) {
            feasBound = (broker_->getBestKnowledge(AlpsKnowledgeTypeSolution)).second;
        }
	
        bestNode = broker_->getBestNode();
        
        if (bestNode) {
            relBound = bestNode->getQuality();
        }	

	if (numNodesProcessed == 0 ||
	    (numNodesProcessed % (nodeInterval * 30)) == 0) {
	    /* Print header. */
	    std::cout << "\n";
            std::cout << "    Node";         /*8 spaces*/
            std::cout << "      ObjValue";   /*14 spaces*/
            if (msgLevel > 2) {
                std::cout << "     Row";    /*8 Spaces*/
                std::cout << "  Column";    /*8 Spaces*/
                std::cout << "   Index";     /*8 spaces*/
                std::cout << "  Parent";     /*8 spaces*/
                std::cout << "   Depth";     /*8 spaces*/
            }

            std::cout << "  BestFeasible";
            std::cout << "     BestBound";
            std::cout << "      Gap";         /*9 spaces*/
            std::cout << "   Time";
            std::cout << "    Left";
            std::cout << std::endl;
	}
	
	if (numNodesProcessed < 10000000) {
	    printf("%8d", numNodesProcessed);
	}
	else {
	    printf("%7dK", numNodesProcessed/1000);
	}

        /* Quality */
	if (node->getStatus() == AlpsNodeStatusFathomed) {
	    printf("      Fathomed");
	}
	else {
	    printf(" %13g", node->getQuality());
	}

        if (msgLevel > 2) {
            if (numRows < 10000000) {
                printf("%8d", numRows);
            }
            else {
                printf("%7dK", numRows/1000);
            }
            if (numCols < 10000000) {
                printf("%8d", numCols);
            }
            else {
                printf("%7dK", numCols/1000);
            }
            /* This index */
            printf("%8d", node->getIndex());
            /* Paraent index */
            if (node->getParent()) {
                printf(" %7d", node->getParent()->getIndex());
            }
            else {
                printf("        ");
            }
            /* Depth */
            printf(" %7d", node->getDepth());
        }

	if (feasBound > ALPS_OBJ_MAX_LESS) {
	    printf("              ");
	}
	else {
	    printf(" %13g", feasBound);
	}

	if (relBound > ALPS_OBJ_MAX_LESS) {
	    printf("              ");
	}
	else {
	    printf(" %13g", relBound);    
	}

        /* Gap */
	if ( (feasBound < ALPS_OBJ_MAX_LESS) &&
	     (relBound < ALPS_OBJ_MAX_LESS) ) {
	    gapVal = ALPS_MAX(0, feasBound - relBound);
	    gap = 100 * gapVal / (ALPS_FABS(relBound) + 1.0);
	}
	if (gap > ALPS_OBJ_MAX_LESS) {
	    printf("         "); /* 9 spaces*/
 	}
	else {
	    if (gap < 1.0e4) {
		printf(" %7.2f%%", gap);
	    }
	    else {
		printf("% 8g", gapVal);
	    }
	}

	int solTime = static_cast<int>(broker_->timer().getCpuTime());
	if (solTime < 1000000) {
	    printf("%7d", solTime);
	}
	else {
	    solTime = static_cast<int>(solTime/3600.0);
	    printf("%6d", solTime);
	    printf("H");
	}

        /* Number of left nodes */
	if (numNodesLeft < 10000000) {
	    printf(" %7d", numNodesLeft);
	}
	else {
	    printf(" %6dK", numNodesLeft/1000);
	}

	printf("\n");
    }
}

//#############################################################################

// Borrow from CBC
void 
BlisModel::analyzeObjective()
{ 
    const double *objective = getObjCoef();
    const double *lower = getColLower();
    const double *upper = getColUpper();
    
    /*------------------------------------------------------
     * Take a first scan to see if there are unfixed continuous
     * variables in the objective.  If so, the minimum objective change 
     * could be arbitrarily small. Also pick off the maximum 
     * coefficient of an unfixed integer variable.
     * If the objective is found to contain only integer variables, 
     * set the fathoming discipline to strict.
     *------------------------------------------------------*/
    
    double maximumCost = 0.0;
    bool possibleMultiple = true;
    int iColumn;
    int numberColumns = getNumCols();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) { 
	if (upper[iColumn] > lower[iColumn]+1.0e-8) {
	    if (colType_[iColumn] == 'I' || colType_[iColumn] == 'B') {
		maximumCost = CoinMax(maximumCost,
				      fabs(objective[iColumn]));
	    }
	    else if (objective[iColumn]) {
		possibleMultiple = false;
	    } 
	}
    }
    
    //setIntParam(CbcModel::CbcFathomDiscipline,possibleMultiple);
    
    /*------------------------------------------------------
     * If a nontrivial increment is possible, try and figure 
     * it out. We're looking for gcd(c<j>) for all c<j> that 
     * are coefficients of unfixed integer variables. Since 
     * the c<j> might not be integers, try and inflate them
     * sufficiently that they look like integers (and we'll 
     * deflate the gcd later).
     * 2520.0 is used as it is a nice multiple of 2,3,5,7
     *-----------------------------------------------------*/

    if (possibleMultiple && maximumCost) { 
	int increment = 0;
	double multiplier = 2520.0;
	while (10.0*multiplier * maximumCost < 1.0e8) {
	    multiplier *= 10.0;
	}
	for (iColumn = 0; iColumn < numberColumns; iColumn++) { 
	    if (upper[iColumn] > lower[iColumn]+1.0e-8) { 
		if ( (colType_[iColumn] == 'I' || colType_[iColumn] == 'B') 
		     && objective[iColumn] ) {
		    double value = fabs(objective[iColumn])*multiplier;
		    int nearest = (int) floor(value+0.5);
		    if (fabs(value-floor(value+0.5)) > 1.0e-8) {
			increment = 0;
			break; 
		    }
		    else if (!increment) { 
			increment = nearest;
		    }
		    else { 
			increment = gcd(increment,nearest);
		    }
		} 
	    } 
	}

	/*--------------------------------------------------
	 * If the increment beats the current value for objective 
	 * change, install it.
	 *--------------------------------------------------*/

	if (increment) { 
	    double value = increment;
	    double cutoffInc = BlisPar_->entry(BlisParams::cutoffInc);

	    value /= multiplier;
	    if (value * 0.999 > cutoffInc) { 
		blisMessageHandler()->message(BLIS_CUTOFF_INC,
                                              blisMessages())
		    << value << CoinMessageEol;
		BlisPar_->setEntry(BlisParams::cutoffInc, -value*0.999);
	    } 
	} 
    }
}

//#############################################################################


