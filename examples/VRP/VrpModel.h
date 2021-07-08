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

#ifndef VrpModel_h_
#define VrpModel_h_

//#############################################################################

#include <vector>

#include "BlisModel.h"
#include "VrpVariable.h"
#include "VrpCommonTypes.h"
#include "VrpConstants.h"
#include "VrpParams.h"
#include "VrpCutGenerator.h"

//#############################################################################

/** Model class for VRP. */
class VrpModel : public BlisModel 
{    

   friend class VrpCutGenerator;
   friend class VrpSolution;
   
 private:

   char name_[100];
   int vertnum_;
   int edgenum_;
   int numroutes_;
   int depot_;
   int capacity_;
   int wtype_;
   int *demand_;    /*vertnum_*/
   int *posx_;      /*vertnum_*/
   int *posy_;      /*vertnum_*/
   double *coordx_; /*vertnum_*/
   double *coordy_; /*vertnum_*/
   double *coordz_; /*vertnum_*/
   double etol_;

   VrpParams *VrpPar_;
   VrpNetwork *n_;  /* Allocate when readInstance (no data filled in). */
   
   // edges_ hold the same elements as variables_ does, do not free memory.
   // For parallel, reinsert elements in variables_ to edges_
   std::vector<VrpVariable *> edges_;

protected:

   /** 1) Set colMatrix_, varLB_, varUB_, conLB_, conUB, numCols_, numRows_
    *  2) Set objCoef_ and objSense_
    *  3) Set colType_ ('C', 'I', or 'B')
    */
   void setModelData();
   
public:

    /** Default construtor. */
   VrpModel() : vertnum_(0), edgenum_(0), numroutes_(0), depot_(0),
      capacity_(0), wtype_(0), etol_(1e-5){
      demand_ = 0;
      posx_ = 0;
      posy_ = 0;
      coordx_ = 0;
      coordy_ = 0;
      coordz_ = 0;
      n_ = 0;
      VrpPar_ = new VrpParams;

      AlpsPar()->setEntry(AlpsParams::searchStrategy,
			  AlpsSearchTypeBestFirst);
      AlpsPar()->setEntry(AlpsParams::staticBalanceScheme, 1); // Spiral
      AlpsPar()->setEntry(AlpsParams::nodeLogInterval, 20);
      BlisPar()->setEntry(BlisParams::branchStrategy,
			  BlisBranchingStrategyStrong);

      BlisPar()->setEntry(BlisParams::branchStrategyRampUp,
			  BlisBranchingStrategyStrong);

      BlisPar()->setEntry(BlisParams::cutCliqueStrategy,BlisCutStrategyNone);
      BlisPar()->setEntry(BlisParams::cutFlowCoverStrategy,
			  BlisCutStrategyNone);
      BlisPar()->setEntry(BlisParams::cutGomoryStrategy,BlisCutStrategyNone);
      BlisPar()->setEntry(BlisParams::cutKnapsackStrategy,BlisCutStrategyNone);
      BlisPar()->setEntry(BlisParams::cutMirStrategy,BlisCutStrategyNone);
      BlisPar()->setEntry(BlisParams::cutOddHoleStrategy,BlisCutStrategyNone);
      BlisPar()->setEntry(BlisParams::cutProbingStrategy,BlisCutStrategyNone);
      BlisPar()->setEntry(BlisParams::cutTwoMirStrategy,BlisCutStrategyNone);
      BlisPar()->setEntry(BlisParams::heurRoundStrategy, BlisHeurStrategyNone);

      // Cuts as formulation
      BlisPar()->setEntry(BlisParams::cutFactor, ALPS_DBL_MAX);
      BlisPar()->setEntry(BlisParams::cutPass, ALPS_INT_MAX);
      BlisPar()->setEntry(BlisParams::tailOff, -1000.0);
      BlisPar()->setEntry(BlisParams::denseConFactor, ALPS_DBL_MAX);

      // Seed
      CoinSeedRandom(1234567);
   }
   
   /** Destructor. */
   virtual ~VrpModel() {
      delete [] demand_; demand_ = 0;
      delete [] posx_; posx_ = 0;
      delete [] posy_; posy_ = 0;
      delete [] coordx_; coordx_ = 0;
      delete [] coordy_; coordy_ = 0;
      delete [] coordz_; coordz_ = 0;
      delete    VrpPar_; VrpPar_ = 0;
      delete    n_; n_ = 0;
   }
   
   /** For parallel code, only the master calls this function.
    *  1) Read in the instance data
    *  2) Set colMatrix_, varLB_, varUB_, conLB_, conUB
    *     numCols_, numRows_
    *  3) Set objCoef_ and objSense_
    *  4) Set colType_ ('C', 'I', or 'B')
    *  5) Create variables and constraints
    *  6) Set numCoreVariables_ and numCoreConstraints_
    */
   virtual void readInstance(const char* dateFile);
   
   /** Read in Alps, Blis, Vrp parameters. */
   virtual void readParameters(const int argnum, const char * const *arglist);
   
   /** User's criteria for a feasible solution. 
    *  If user think the given solution is feasible then need
    *     1) set userFeasible to true, and
    *     2) return a non-null VRP solution.
    *  If user think the solution is infeasible then need
    *     1) set userFeasible to false, and
    *     2) return a null.
    */
   virtual BlisSolution * userFeasibleSolution(const double *solution,
                                               bool &userFeasible);
   
   int index (int v0, int v1) {
      return(v0 < v1 ? v1*(v1 - 1)/2 + v0 : v0*(v0 - 1)/2 + v1);
   }
   
   int computeCost(int v0, int v1); 
   
   int getNumVertices() { return vertnum_; }

   int getNumEdges() { return edgenum_; }
   
   std::vector<VrpVariable *> getEdgeList() { return edges_; }

   // Transform dense solution to a sparse vector.
   CoinPackedVector *getSolution(const double *denseSol);

   void createNet(CoinPackedVector *vec);

   /** Register knowledge. */
   virtual void registerKnowledge();  
   
   /** Pack Vrp portion of the model into an encoded object. */
   AlpsReturnStatus encodeVrp(AlpsEncoded *encoded) const;
   
   /** Unpack Vrp portion of the model from an encoded object. */
   AlpsReturnStatus decodeVrp(AlpsEncoded &encoded);  
   
   /** The method that encodes the model into an encoded object. */
   virtual AlpsEncoded* encode() const;
   
   /** The method that decodes the model from an encoded object. */
   virtual void decodeToSelf(AlpsEncoded&);
   
};

//#############################################################################

#endif
