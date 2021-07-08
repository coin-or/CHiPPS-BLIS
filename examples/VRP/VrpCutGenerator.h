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

#ifndef VrpCutGenerator_h_
#define VrpCutGenerator_h_

//#############################################################################

#include "BlisConGenerator.h"
#include "VrpModel.h"
#include "VrpNetwork.h"
#include "VrpMacros.h"
#ifdef DO_TSP_CUTS
extern "C"{
   #include "concorde.h"
}
#endif

class VrpModel;

//#############################################################################

class VrpCutGenerator : public BlisConGenerator  
{
private:

   VrpModel *model_;
   char **coef_list;
   int *ref_;
   double *cutVal_;
   char *cutList_;
   char *inSet_;
   
public:
    
    /** Construtors. */
   VrpCutGenerator(VrpModel *vrp=0, int vertnum = 0);
   
   /** Destructor. */
   virtual ~VrpCutGenerator() {
      delete [] ref_; ref_ = 0;
      delete [] cutVal_; cutVal_ = 0;
      delete [] cutList_; cutList_ = 0;
      delete [] inSet_; inSet_ = 0;
   }

   /** Generate cons for the client model.
       The routine returns true if reoptimisation is needed (because the 
       state of the solver interface has been modified).
   */
   virtual bool generateConstraints(BcpsConstraintPool &conPool); 

   int connectivityCuts(BcpsConstraintPool &conPool);
   
   int addVrpCut(BcpsConstraintPool &conPool, char *coef, int rhs,int type);

   void setModel(VrpModel *vrp){ model_ = vrp; }
   
   int greedyShrinking1(VrpModel *m, int max_shrink_cuts,
			BcpsConstraintPool &conPool);

   int greedyShrinking1One(VrpModel *m, int max_shrink_cuts, 
			   BcpsConstraintPool &conPool);

   int greedyShrinking6(VrpModel *m, int max_shrink_cuts, int trial_num,
			double prob, BcpsConstraintPool &conPool);

   int greedyShrinking6One(VrpModel *m, int max_shrink_cuts, int trial_num,
			   double prob, BcpsConstraintPool &conPool);

   int greedyShrinking2One(VrpModel *m, int max_shrink_cuts, 
			   BcpsConstraintPool &conPool);

#ifdef DO_TSP_CUTS
   int tspCuts(VrpModel *m, BcpsConstraintPool &conPool);
   int addTspCuts(VrpModel *m, BcpsConstraintPool &conPool,
		  CCtsp_lpcut_in **tsp_cuts, CCtsp_lpgraph *g);
#endif

};

//#############################################################################

#endif
