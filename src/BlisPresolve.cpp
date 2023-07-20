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

#include "BlisPresolve.h"

//#############################################################################

OsiSolverInterface *
BlisPresolve::preprocess(OsiSolverInterface & origModel,
                         double feasibilityTolerance,
                         bool keepIntegers,
                         int numberPasses,
                         const char * prohibited)
{
    return presolvedModel(origModel,
                          feasibilityTolerance,
                          keepIntegers,
                          numberPasses,
                          prohibited);
}

//#############################################################################

void 
BlisPresolve::postprocess(bool updateStatus)
{
    postsolve(updateStatus);
}

//#############################################################################
