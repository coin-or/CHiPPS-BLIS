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

#ifndef VrpSolution_h_
#define VrpSolution_h_

#include "BlisSolution.h"

#include "VrpModel.h"

//#############################################################################
/** This class contains a vrp solution */
//#############################################################################

class VrpSolution : public BlisSolution {

protected:

   _node *opt_;
    
public:
    
    /** Default constructor. */
   VrpSolution() : BlisSolution(), opt_(0)
	{}

    /** Useful constructor. */
   VrpSolution(int s, const double *values, double objValue, VrpModel *vrp=0);
   
   /** Destructor. */
   virtual ~VrpSolution(){
      delete [] opt_;
   }

   /** Set opt */
   void setOpt(_node *opt) { 
       if (opt_) delete opt_;
       opt_ = opt; 
   }
   
   /** Print the solution.*/
   virtual void print(std::ostream& os) const; 
   
   /** The method that encodes the solution into a encoded object. */
   virtual AlpsEncoded* encode() const {
      AlpsEncoded* encoded = new AlpsEncoded(AlpsKnowledgeTypeSolution);
      encodeBcps(encoded);

      //Vrp part.
      int cur_vert = opt_[0].next, count = 1;  
      while (cur_vert != 0){
          cur_vert = opt_[cur_vert].next;
          ++count;
      }
      encoded->writeRep(count);
      for (int j = 0; j < count; ++j) {
          encoded->writeRep(opt_[j].next);
          encoded->writeRep(opt_[j].route);
      }

      return encoded;
   }
   
   /** The method that decodes the solution from a encoded object. */
   virtual AlpsKnowledge* decode(AlpsEncoded& encoded) const {
      VrpSolution * sol = new VrpSolution();
      sol->decodeBcps(encoded);

      // Vrp part.
      int count;
      encoded.readRep(count);
      _node *opt = new _node [count];
      for (int j = 0; j < count; ++j) {
          encoded.readRep(opt[j].next);
          encoded.readRep(opt[j].route);
      }
      sol->setOpt(opt);
      
      return sol;
   }
   
};

//#############################################################################
//#############################################################################

#endif
