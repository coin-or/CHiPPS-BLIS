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

#ifndef VrpHeurTSP_h_
#define VrpHeurTSP_h_

//#############################################################################

#include <vector>

#include "CoinPackedVector.hpp"

#include "BlisHeuristic.h"
#include "VrpModel.h"

//#############################################################################
#if 0
class VrpAdjList 
{
private:
    std::vector<CoinPackedVector *> list_;
    VrpModel *model_;
public:
    VrpAdjList() {}
    VrpAdjList(VrpModel *m) { 
	model_ = m; 
	createAdjList(m);
    }
    
    void createAdjList(VrpModel *model);
};

typedef struct vrp_neighbors 
{
    int n1;
    int n2;
} VrpNeighbors;

#endif

//#############################################################################

class VrpHeurTSP : public BlisHeuristic {
private:
    /** Illegal Assignment operator. */ 
    VrpHeurTSP & operator=(const VrpHeurTSP& rhs);
    
protected:
    /* Stored the predetermined next vertex to visit for vertex k if 
       the value determined_[k] greater than zero. */
    //int *determined_;
    
    /* Adjacent list of all vertices. */
    std::vector<CoinPackedVector *> adjList_;

    /** Create adjacent list for each vertex. */
    void createAdjList(VrpModel *model);

    /** TSP Tour */
    std::vector<int> tour_;
    
    /** Mark if vertices have been visited. */
    bool *visited_;

    /** The node at which this heuristic was call*/
    int preNode_;

    /** Neighbors determined from LP solution */
    //VrpNeighbor *
    int *neighbors_;

    /** Call how many time at a node. */
    int nodeCalls_;
    
    /** Edge and column relationship. Give an edge {v0, v1}, 
	edgeColMatch_[v1-1][v0] is the column index. */
    std::vector<int> *edgeColMatch_;
    
    void freeGuts() {
        if (visited_) {
            delete [] visited_;
            visited_ = NULL;
        }
	int numVertices = adjList_.size();
	for (int k = 0; k < numVertices; ++k) {
	    delete adjList_[k];
	}
        adjList_.clear();
	if (neighbors_) {
            delete [] neighbors_;
            neighbors_ = NULL;
        }
	if (edgeColMatch_) {
	    delete [] edgeColMatch_;
	    edgeColMatch_ = NULL;
	}
    }

public:
    /** Default Constructor. */
    VrpHeurTSP()
        : 
        visited_(0), preNode_(-1), 
	neighbors_(0), nodeCalls_(0), edgeColMatch_(0) {}
    
    /** Constructor with model. */
    VrpHeurTSP(VrpModel * model, const char *name,
               BlisHeurStrategy strategy, int freq)
        :
        BlisHeuristic(model, name, strategy, freq)
    { 
	visited_ = NULL;
	preNode_ = -1;
	neighbors_ = NULL;
        nodeCalls_ = 0;
	edgeColMatch_ = NULL;
        createAdjList(model);
    }

    /** Destructor. */
    ~VrpHeurTSP() 
    { 
        freeGuts();
    }
    
    /** Returns 0 if no solution, 1 if valid solution. newSolution stores
        the solution in dense format. */
    virtual bool searchSolution(double & objectiveValue, double * newSolution);
};
#endif

//#############################################################################
