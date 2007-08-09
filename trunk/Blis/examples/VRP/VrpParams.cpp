/*===========================================================================*
 * This file is part of a solver for the Vehicle Routing Problem             *
 * developed using the BiCePS Linear Integer Solver (BLIS).                  *
 *                                                                           *
 * This solver is distributed under the Common Public License as part of     * 
 * the COIN-OR repository (http://www.coin-or.org).                          *
 *                                                                           *
 * Authors: Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2007 Yan Xu and Ted Ralphs.                                 *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#include "VrpParams.h"
#include "VrpConstants.h"

using std::make_pair;

void 
VrpParams::createKeywordList() {

  //--------------------------------------------------------
  // Create the list of keywords for parameter file reading
  //--------------------------------------------------------

  //--------------------------------------------------------
  // BoolPar
  //--------------------------------------------------------

  keys_.push_back(make_pair(std::string("Vrp_doGreedy"),
			    AlpsParameter(AlpsIntPar, doGreedy)));
    
  keys_.push_back(make_pair(std::string("Vrp_doExtraInRoot"),
			    AlpsParameter(AlpsIntPar, doExtraInRoot)));

  //--------------------------------------------------------
  // BoolArrayPar
  //--------------------------------------------------------
  
  //--------------------------------------------------------
  // Int Parameters
  //--------------------------------------------------------

  keys_.push_back(make_pair(std::string("Vrp_numRoutes"),
			    AlpsParameter(AlpsIntPar, numRoutes)));

  keys_.push_back(make_pair(std::string("Vrp_verbosity"),
			    AlpsParameter(AlpsIntPar, verbosity)));

  keys_.push_back(make_pair(std::string("Vrp_greedyNumTrials"),
			    AlpsParameter(AlpsIntPar, greedyNumTrials)));

  keys_.push_back(make_pair(std::string("Vrp_whichConnectedRoutine"),
			    AlpsParameter(AlpsIntPar, whichConnectedRoutine)));

  keys_.push_back(make_pair(std::string("Vrp_maxNumCutsInShrink"),
			    AlpsParameter(AlpsIntPar, maxNumCutsInShrink)));

  //--------------------------------------------------------
  // String Parameters.
  //--------------------------------------------------------

}

//#############################################################################

void 
VrpParams::setDefaultEntries() {

  //-------------------------------------------------------------
  // Char Parameters.
  //-------------------------------------------------------------

   setEntry(doGreedy, true);

   setEntry(doExtraInRoot, false);

  //-------------------------------------------------------------
  // Int Parameters.
  //-------------------------------------------------------------

   setEntry(numRoutes, VRP_NOT_SET);

   setEntry(verbosity, 0);

   setEntry(greedyNumTrials, 5);

   setEntry(whichConnectedRoutine, BOTH);

   setEntry(maxNumCutsInShrink, 200);

  //-------------------------------------------------------------
  // Double Parameters
  //-------------------------------------------------------------

  //-------------------------------------------------------------
  // String Parameters
  //-------------------------------------------------------------
  
}

//#############################################################################
