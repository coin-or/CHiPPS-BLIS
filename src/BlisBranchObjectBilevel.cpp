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
 * Copyright (C) 2001-2023, Lehigh University, Yan Xu, and Ted Ralphs.       *
 * All Rights Reserved.                                                      *
 *===========================================================================*/


#include "CoinHelperFunctions.hpp"

#include "BlisBranchObjectBilevel.h"


//#############################################################################

// Assignment operator 
BlisBranchObjectBilevel & 
BlisBranchObjectBilevel::operator=(const BlisBranchObjectBilevel& rhs)
{
    if (this != &rhs) {
	BcpsBranchObject::operator=(rhs);
	branchingSet_ = rhs.branchingSet_;
    }
    return *this;
}

//#############################################################################

/*
  Perform a branch by adjusting the bounds of the specified variable.
  Note that each arm of the branch advances the object to the next arm by
  advancing the value of direction_.
  
  Providing new values for the variable's lower and upper bounds for each
  branching direction gives a little bit of additional flexibility and will
  be easily extensible to multi-direction branching.
  Returns change in guessed objective on next branch
*/

/* As far as I can tell, this function is never called in the current
   implementation*/

double
BlisBranchObjectBilevel::branch(bool normalBranch)
{
    return 0.0;
}

//#############################################################################

// Print what would happen  
void
BlisBranchObjectBilevel::print(bool normalBranch)
{
    std::deque<int>::iterator ptr1;
    std::cout << "Branching set consists of variables ";

    for (ptr1 = branchingSet_->begin(); ptr1 != branchingSet_->end(); ptr1++){
	std::cout << " " << *ptr1;
    }
    
    std::cout << std::endl;
}

//#############################################################################
