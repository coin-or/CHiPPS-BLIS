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

//#############################################################################

#ifndef BlisHelp_h_
#define BlisHelp_h_

#include "AlpsEncoded.h"

#include "BlisConfig.h"
#include "Blis.h"

class CoinWarmStartBasis;
class OsiRowCut;
class BlisConstraint;
class BlisModel;

//#############################################################################

/** Convert a OsiRowCut to a Blis Contraint. */
BLISLIB_EXPORT BlisConstraint * BlisOsiCutToConstraint(const OsiRowCut *rowCut);

/** Strong branching on a variable colInd. */
BLISLIB_EXPORT BlisReturnStatus BlisStrongBranch(BlisModel *model, double objValue, int colInd,
                                                 double x, const double *saveLower,
                                                 const double *saveUpper, bool &downKeep,
                                                 bool &downFinished, double &downDeg,
                                                 bool &upKeep, bool &upFinished, double &upDeg);

/** Pack coin warm start into an encoded object. */
BLISLIB_EXPORT int BlisEncodeWarmStart(AlpsEncoded *encoded, const CoinWarmStartBasis *ws);

/** Unpack coin warm start from an encoded object. */
BLISLIB_EXPORT CoinWarmStartBasis *BlisDecodeWarmStart(AlpsEncoded &encoded,
                                                       AlpsReturnStatus *rc);

/** Compute and return a hash value of an Osi row cut. */
BLISLIB_EXPORT double BlisHashingOsiRowCut(const OsiRowCut *rowCut, 
                                           const BlisModel *model);

/** Check if a row cut parallel with another row cut. */
BLISLIB_EXPORT bool BlisParallelCutCut(OsiRowCut * rowCut1,
                                       OsiRowCut * rowCut2,
                                       double threshold = 1.0);

/** Check if a row cut parallel with a constraint. */
BLISLIB_EXPORT bool BlisParallelCutCon(OsiRowCut * rowCut,
                                       BlisConstraint * con,
                                       double threshold = 1.0);

/** Check if a row cut parallel with a constraint. */
BLISLIB_EXPORT bool BlisParallelConCon(BlisConstraint * con1,
                                       BlisConstraint * con2,
                                       double threshold = 1.0);


#endif
