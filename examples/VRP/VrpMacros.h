/*===========================================================================*
 * This file is part of a solver for the Vehicle Routing Problem             *
 * developed using the BiCePS Linear Integer Solver (BLIS).                  *
 *                                                                           *
 * This solver is distributed under the Eclipse Public License as part of    * 
 * the COIN-OR repository (http://www.coin-or.org).                          *
 *                                                                           *
 * Authors: Yan Xu, Lehigh University                                        *
 *          Ted Ralphs, Lehigh University                                    *
 *                                                                           *
 * Copyright (C) 2007 Yan Xu and Ted Ralphs.                                 *
 * All Rights Reserved.                                                      *
 *===========================================================================*/

#ifndef VrpMacros_h_
#define VrpMacros_h_

#if defined(_MSC_VER) || defined(__MNO_CYGWIN) || defined(__MINGW32__)
/* Different function call in Windows */ 
#define SRANDOM(seed) srand(seed)
#define RANDOM() rand()
#else
#define SRANDOM(seed) srandom(seed)
#define RANDOM() random()
#endif

/*approximates the number of trucks necessary to service a set of customers*/
#define BINS(weight, capacity) \
((int) ceil(((double)weight)/((double)capacity)))

/*calculates the right hand side of a subtour elimination constraint*/
#define RHS(cust_num, weight, capacity) \
(cust_num-BINS(weight, capacity))

/*===========================================================================*/

#endif
