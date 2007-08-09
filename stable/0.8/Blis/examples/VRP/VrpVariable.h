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

#ifndef VrpVariable_h_
#define VrpVariable_h_

//#############################################################################

#include "BlisVariable.h"

//#############################################################################

/** Variable class for VRP. */
class VrpVariable : public BlisVariable 
{    

 private:

   /* The endpoints of the edge */
   int ends_[2];
   int uind_;
   
public:

    VrpVariable() {
       ends_[0] = 0;
       ends_[1] = 0;
    }

    VrpVariable(int v1, int v2, int cost) {
       ends_[0] = v1 < v2 ? v1 : v2;
       ends_[1] = v1 < v2 ? v2 : v1;
       uind_ = ends_[1]*(ends_[1] - 1)/2 + ends_[0];
       int indices [2];
       double values [2];
       indices[0] = ends_[0];
       indices[1] = ends_[1];
       values[0] = values[1] = 1.0;
       setData(2, indices, values);
       setIntType('B');
       setLbHard(0.0);
       if (ends_[0]){ 
	  setUbHard(1.0);
       }else{
	  setUbHard(2.0);
       }
       setObjCoef((double) cost);
    }

    virtual ~VrpVariable() {
        //std::cout << "delete a vrp variable " << std::endl;
    }

    inline int getIndex() { return uind_; }
    inline int getv0() { return ends_[0]; }
    inline int getv1() { return ends_[1]; }
};

//#############################################################################

#endif
