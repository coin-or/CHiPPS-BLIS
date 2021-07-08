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
   
protected:

   /** Pack Vrp part into an encoded object. */
    AlpsReturnStatus encodeVrp(AlpsEncoded *encoded) {
	AlpsReturnStatus status = AlpsReturnStatusOk;

        //std::cout << "****** encodeVrp var: size_ = " << size_ << std::endl;

	encoded->writeRep(ends_[0]);
	encoded->writeRep(ends_[1]);
	encoded->writeRep(uind_);
        
	return status;
    }    

    /** Unpack Vrp part from a encode object. */
    AlpsReturnStatus decodeVrp(AlpsEncoded &encoded) {
	AlpsReturnStatus status = AlpsReturnStatusOk;

	encoded.readRep(ends_[0]);
	encoded.readRep(ends_[1]);
	encoded.readRep(uind_);
        
        //std::cout << "****** decodeVrp var: size_ = " << size_ << std::endl;

	return status;
    }
    
public:

    /** Default constructor. */
    VrpVariable() {
       ends_[0] = 0;
       ends_[1] = 0;
    }

    /** Useful constructor. */
    VrpVariable(int v1, int v2, int cost, int ub) {
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
       setUbHard((double) ub);
       setObjCoef((double) cost);
    }

    /** Destructor. */
    virtual ~VrpVariable() {
        //std::cout << "delete a vrp variable " << std::endl;
    }
  
    /** Get data  */
    /**@{*/
    inline int getIndex() { return uind_; }
    inline int getv0() { return ends_[0]; }
    inline int getv1() { return ends_[1]; }
    /**@}*/

    virtual void printDesc() {
	std::cout << "(" << getv0() << ", " << getv1() << ")";
    }
    
    /** Pack to a encode object. */
    virtual AlpsReturnStatus encode(AlpsEncoded *encoded){
	AlpsReturnStatus status;

	status = encodeBcpsObject(encoded);
	status = encodeBlis(encoded);
	status = encodeVrp(encoded);
        
	return status;
    }

    /** Decode a variable from an encoded object. */
    virtual AlpsKnowledge* decode(AlpsEncoded &encoded) const {
	AlpsReturnStatus status = AlpsReturnStatusOk;
	VrpVariable * var = new VrpVariable();    
	
	// Unpack Bcps part.
	status = var->decodeBcpsObject(encoded);
	if (status) {
	    throw CoinError("Failed to decode Bcps part of var",
			    "decode",
			    "BlisObject");
	}
	
	// Unpack Blis part.
	status = var->decodeBlis(encoded);
	if (status) {
	    throw CoinError("Failed to decode Blis part of var", 
			    "decode", 
			    "BlisObject");
	}

	// Unpack Vrp part.
	status = var->decodeVrp(encoded);
	if (status) {
	    throw CoinError("Failed to decode Vrp part of var", 
			    "decode", 
			    "BlisObject");
	}
	return var;
    }

};

//#############################################################################

#endif
