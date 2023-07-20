/*===========================================================================*
 * This file is part of the Abstract Library for Parallel Search (ALPS).     *
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

#include "BlisParams.h"

using std::make_pair;

void 
BlisParams::createKeywordList() {

  //--------------------------------------------------------
  // Create the list of keywords for parameter file reading
  //--------------------------------------------------------

  //--------------------------------------------------------
  // CharPar
  //--------------------------------------------------------
 
    keys_.push_back(make_pair(std::string("Blis_cutRampUp"),
                              AlpsParameter(AlpsBoolPar, cutRampUp)));
    
   keys_.push_back(make_pair(std::string("Blis_presolve"),
			     AlpsParameter(AlpsBoolPar, presolve)));
    
   keys_.push_back(make_pair(std::string("Blis_shareConstraints"),
			     AlpsParameter(AlpsBoolPar,
					   shareConstraints)));

   keys_.push_back(make_pair(std::string("Blis_shareVariables"),
			     AlpsParameter(AlpsBoolPar,
					   shareVariables)));

   keys_.push_back(make_pair(std::string("Blis_sharePseudocostRampUp"),
			     AlpsParameter(AlpsBoolPar,
					   sharePseudocostRampUp)));
    
   keys_.push_back(make_pair(std::string("Blis_sharePseudocostSearch"),
			     AlpsParameter(AlpsBoolPar,
					   sharePseudocostSearch)));
   keys_.push_back(make_pair(std::string("Blis_checkFeasibleInBranch"),
			     AlpsParameter(AlpsBoolPar,
					   checkFeasibleInBranch)));
    
  //--------------------------------------------------------
  // BoolArrayPar
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  // Int Parameters
  //--------------------------------------------------------

  keys_.push_back(make_pair(std::string("Blis_branchStrategy"),
			    AlpsParameter(AlpsIntPar, branchStrategy)));

  keys_.push_back(make_pair(std::string("Blis_branchStrategyRampUp"),
			    AlpsParameter(AlpsIntPar, branchStrategyRampUp)));

  keys_.push_back(make_pair(std::string("Blis_cutPass"),
			    AlpsParameter(AlpsIntPar, cutPass)));

  keys_.push_back(make_pair(std::string("Blis_quickCutPass"),
			    AlpsParameter(AlpsIntPar, quickCutPass)));
  
  keys_.push_back(make_pair(std::string("Blis_cutStrategy"),
			    AlpsParameter(AlpsIntPar, cutStrategy)));

  keys_.push_back(make_pair(std::string("Blis_cutGenerationFrequency"),
			    AlpsParameter(AlpsIntPar, 
					  cutGenerationFrequency)));

  keys_.push_back(make_pair(std::string("Blis_cutCliqueStrategy"),
			    AlpsParameter(AlpsIntPar, cutCliqueStrategy)));

  keys_.push_back(make_pair(std::string("Blis_cutGomoryStrategy"),
			    AlpsParameter(AlpsIntPar, cutGomoryStrategy)));

  keys_.push_back(make_pair(std::string("Blis_cutFlowCoverStrategy"),
			    AlpsParameter(AlpsIntPar, cutFlowCoverStrategy)));

  keys_.push_back(make_pair(std::string("Blis_cutKnapsackStrategy"),
			    AlpsParameter(AlpsIntPar, cutKnapsackStrategy)));
  
  keys_.push_back(make_pair(std::string("Blis_cutMirStrategy"),
			    AlpsParameter(AlpsIntPar, cutMirStrategy)));
  
  keys_.push_back(make_pair(std::string("Blis_cutOddHoleStrategy"),
			    AlpsParameter(AlpsIntPar, cutOddHoleStrategy)));

  keys_.push_back(make_pair(std::string("Blis_cutProbingStrategy"),
			    AlpsParameter(AlpsIntPar, cutProbingStrategy)));

  keys_.push_back(make_pair(std::string("Blis_cutTwoMirStrategy"),
			    AlpsParameter(AlpsIntPar, cutTwoMirStrategy)));

  keys_.push_back(make_pair(std::string("Blis_cutCliqueFreq"),
			    AlpsParameter(AlpsIntPar, cutCliqueFreq)));

  keys_.push_back(make_pair(std::string("Blis_cutGomoryFreq"),
			    AlpsParameter(AlpsIntPar, cutGomoryFreq)));

  keys_.push_back(make_pair(std::string("Blis_cutFlowCoverFreq"),
			    AlpsParameter(AlpsIntPar, cutFlowCoverFreq)));

  keys_.push_back(make_pair(std::string("Blis_cutKnapsackFreq"),
			    AlpsParameter(AlpsIntPar, cutKnapsackFreq)));
  
  keys_.push_back(make_pair(std::string("Blis_cutMirFreq"),
			    AlpsParameter(AlpsIntPar, cutMirFreq)));
  
  keys_.push_back(make_pair(std::string("Blis_cutOddHoleFreq"),
			    AlpsParameter(AlpsIntPar, cutOddHoleFreq)));

  keys_.push_back(make_pair(std::string("Blis_cutProbingFreq"),
			    AlpsParameter(AlpsIntPar, cutProbingFreq)));

  keys_.push_back(make_pair(std::string("Blis_cutTwoMirFreq"),
			    AlpsParameter(AlpsIntPar, cutTwoMirFreq)));

  keys_.push_back(make_pair(std::string("Blis_difference"),
			    AlpsParameter(AlpsIntPar, difference)));
  
  keys_.push_back(make_pair(std::string("Blis_heurStrategy"),
                            AlpsParameter(AlpsIntPar, heurStrategy)));
  
  keys_.push_back(make_pair(std::string("Blis_heurCallFrequencyy"),
                            AlpsParameter(AlpsIntPar, heurCallFrequency)));
  
  keys_.push_back(make_pair(std::string("Blis_heurRoundStrategy"),
			    AlpsParameter(AlpsIntPar, heurRoundStrategy)));
  
  keys_.push_back(make_pair(std::string("Blis_heurRoundFreq"),
			    AlpsParameter(AlpsIntPar, heurRoundFreq)));
  
  keys_.push_back(make_pair(std::string("Blis_lookAhead"),
			    AlpsParameter(AlpsIntPar, lookAhead)));
  
  keys_.push_back(make_pair(std::string("Blis_pseudoRelibility"),
			    AlpsParameter(AlpsIntPar, pseudoRelibility)));

  keys_.push_back(make_pair(std::string("Blis_sharePcostDepth"),
			    AlpsParameter(AlpsIntPar, sharePcostDepth)));

  keys_.push_back(make_pair(std::string("Blis_sharePcostFrequency"),
			    AlpsParameter(AlpsIntPar, sharePcostFrequency)));

  keys_.push_back(make_pair(std::string("Blis_strongCandSize"),
			    AlpsParameter(AlpsIntPar, strongCandSize)));
  
  //--------------------------------------------------------
  // Double Parameters.
  //--------------------------------------------------------
  
  keys_.push_back(make_pair(std::string("Blis_cutFactor"),
			    AlpsParameter(AlpsDoublePar, cutFactor)));
  
  keys_.push_back(make_pair(std::string("Blis_cutoff"),
			    AlpsParameter(AlpsDoublePar, cutoff)));

  keys_.push_back(make_pair(std::string("Blis_cutoffInc"),
			    AlpsParameter(AlpsDoublePar, cutoffInc)));
  
  keys_.push_back(make_pair(std::string("Blis_denseConFactor"),
			    AlpsParameter(AlpsDoublePar, denseConFactor)));

  keys_.push_back(make_pair(std::string("Blis_integerTol"),
			    AlpsParameter(AlpsDoublePar, integerTol)));

  keys_.push_back(make_pair(std::string("Blis_objSense"),
			    AlpsParameter(AlpsDoublePar, objSense)));
  
  keys_.push_back(make_pair(std::string("Blis_optimalRelGap"),
			    AlpsParameter(AlpsDoublePar, optimalRelGap)));
  
  keys_.push_back(make_pair(std::string("Blis_optimalAbsGap"),
			    AlpsParameter(AlpsDoublePar, optimalAbsGap)));
  
  keys_.push_back(make_pair(std::string("Blis_pseudoWeight"),
			    AlpsParameter(AlpsDoublePar, pseudoWeight)));
  
  keys_.push_back(make_pair(std::string("Blis_scaleConFactor"),
			    AlpsParameter(AlpsDoublePar, scaleConFactor)));

  keys_.push_back(make_pair(std::string("Blis_tailOff"),
                            AlpsParameter(AlpsDoublePar, tailOff)));
  
  //--------------------------------------------------------
  // String Parameters.
  //--------------------------------------------------------

}

//#############################################################################

void 
BlisParams::setDefaultEntries() {

  //-------------------------------------------------------------
  // Char Parameters.
  //-------------------------------------------------------------

  setEntry(cutRampUp, true);
  setEntry(presolve, false);
  setEntry(shareConstraints, false);
  setEntry(shareVariables, false);
  setEntry(sharePseudocostRampUp, true);
  setEntry(sharePseudocostSearch, false);
  setEntry(checkFeasibleInBranch, true);
  
  //-------------------------------------------------------------
  // Int Parameters.
  //-------------------------------------------------------------

  setEntry(branchStrategy, BlisBranchingStrategyPseudoCost);
  setEntry(branchStrategyRampUp, BlisBranchingStrategyPseudoCost);
  setEntry(cutStrategy, BlisCutStrategyNotSet);
  setEntry(cutGenerationFrequency, 1);
  setEntry(cutPass, 20);
  setEntry(quickCutPass, 0);
  setEntry(cutCliqueStrategy, BlisCutStrategyNotSet);
  setEntry(cutGomoryStrategy, BlisCutStrategyNotSet);
  setEntry(cutFlowCoverStrategy, BlisCutStrategyNotSet);
  setEntry(cutKnapsackStrategy, BlisCutStrategyNotSet);
  setEntry(cutMirStrategy, BlisCutStrategyNotSet);
  setEntry(cutOddHoleStrategy, BlisCutStrategyNotSet);
  setEntry(cutProbingStrategy, BlisCutStrategyNotSet);
  setEntry(cutTwoMirStrategy, BlisCutStrategyNotSet);
  setEntry(cutCliqueFreq, 1);
  setEntry(cutGomoryFreq, 1);
  setEntry(cutFlowCoverFreq, 1);
  setEntry(cutKnapsackFreq, 1);
  setEntry(cutMirFreq, 1);
  setEntry(cutOddHoleFreq, 1);
  setEntry(cutProbingFreq, 1);
  setEntry(cutTwoMirFreq, 1);
  setEntry(difference, -1);
  setEntry(heurStrategy, BlisHeurStrategyAuto);
  setEntry(heurRoundStrategy, BlisHeurStrategyNotSet);
  setEntry(heurRoundFreq, 1);
  setEntry(lookAhead, 4);
  setEntry(pseudoRelibility, 8);
  setEntry(sharePcostDepth, 30);
  setEntry(sharePcostFrequency, 100);
  setEntry(strongCandSize, 10);
  
  //-------------------------------------------------------------
  // Double Parameters
  //-------------------------------------------------------------

  setEntry(cutFactor, 4.0);
  setEntry(cutoff, ALPS_INC_MAX);
  setEntry(cutoffInc, 1.0e-6);
  setEntry(denseConFactor, 5.0);
  setEntry(integerTol, 1.0e-5);
  setEntry(objSense, 1.0);
  setEntry(optimalRelGap, 1.0e-4);
  setEntry(optimalAbsGap, 1.0e-6);
  setEntry(pseudoWeight, 0.8);
  setEntry(scaleConFactor, 1000000.0);
  setEntry(tailOff, 1e-7);
  
  //-------------------------------------------------------------
  // String Parameters
  //-------------------------------------------------------------
  
}
