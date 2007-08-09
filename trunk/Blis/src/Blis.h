/*===========================================================================*
 * This file is part of the BiCePS Linear Integer Solver (BLIS).             *
 *                                                                           *
 * BLIS is distributed under the Common Public License as part of the        *
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

//#############################################################################

#ifndef Blis_h_
#define Blis_h_

//#############################################################################

enum BlisLpStatus{
   BlisLpStatusOptimal,
   BlisLpStatusAbandoned,
   BlisLpStatusPrimalInfeasible,
   BlisLpStatusDualInfeasible,
   BlisLpStatusPrimalObjLim,
   BlisLpStatusDualObjLim,
   BlisLpStatusIterLim,
   BlisLpStatusUnknown
};

//#############################################################################

enum BlisReturnStatus {
   BlisReturnStatusOk = 0,
   BlisReturnStatusErrLp,
   BlisReturnStatusInfeasible,
   BlisReturnStatusUnbounded,
   BlisReturnStatusOverObjLim,
   BlisReturnStatusFeasible,
   BlisReturnStatusBranch,
   BlisReturnStatusUnknown
};

#if 0
#define BLIS_ERR_LP         100
#define BLIS_INF            200
#define BLIS_UNBOUND        201
#define BLIS_OPTIMAL          0
#define BLIS_UNKNOWN        202
#endif

//#############################################################################

enum BlisCutStrategy{
   BlisCutStrategyNotSet = -1,
   BlisCutStrategyNone = 0,
   BlisCutStrategyRoot,
   BlisCutStrategyAuto,
   BlisCutStrategyPeriodic
};

enum BlisHeurStrategy{
   BlisHeurStrategyNotSet = -1,
   BlisHeurStrategyNone = 0,
   BlisHeurStrategyPeriodic,
   BlisHeurStrategyRoot,
   BlisHeurStrategyAuto
};

#if 0
#define BLIS_NOT_SET       -555
#define BLIS_ROOT            -2
#define BLIS_AUTO            -1
#define BLIS_NONE             0
#endif

//#############################################################################

enum BlisHotStartStrategy{
   BlisHotStartBranchIncorrect,
   BlisHotStartBranchCorrect
};

//#############################################################################

enum BlisBranchingStrategy{
   BlisBranchingStrategyMaxInfeasibility,
   BlisBranchingStrategyPseudoCost,
   BlisBranchingStrategyReliability,   
   BlisBranchingStrategyStrong
};


//#############################################################################

enum BlisSolutionType {
    BlisSolutionTypeBounding,
    BlisSolutionTypeBranching,
    BlisSolutionTypeDiving,
    BlisSolutionTypeRounding,
    BlisSolutionTypeStrong
};

//#############################################################################

/** Branching object type. */
enum BlisBranchingObjectType {
    BlisBranchingObjectTypeNone = 0,
    BlisBranchingObjectTypeInt,
    BlisBranchingObjectTypeSos
};

//#############################################################################

#define BLIS_CUT_DISABLE            20

#define BLIS_HEUR_ROUND_DISABLE     1000000

#define BLIS_PSEUDO                 21

//#############################################################################

#endif
