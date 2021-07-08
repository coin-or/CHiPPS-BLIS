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

#include <vector>
#include "VrpHeurTSP.h"
#include "VrpModel.h"
#include "VrpVariable.h"

//#############################################################################

//------------------------------------------------------
// Given current LP solution. Rounding number closed to 1 to 1
// Let those edge be in the route, then use nearest neighbor idea
// to construct a tour.
//------------------------------------------------------
// 1. stand on an arbitrary vertex as current vertex.
// 2. find out the lightest edge connecting current vertex and 
//    an unvisited vertex V.
// 3. set current vertex be V.
// 4. mark V as visited.
// 5. if all the verteces in domain are visited then terminate.
// 6. Go to step 2.
//------------------------------------------------------
bool
VrpHeurTSP::searchSolution(double & objectiveValue, double * newSolution)
{

    VrpModel *model = dynamic_cast<VrpModel *>(model_);
    int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);

    std::vector<VrpVariable *> edges = model->getEdgeList();
    int numEdges = edges.size();

    int numNodes = model->getNumNodes();

    //------------------------------------------------------
    // Determine if call heuristic
    //------------------------------------------------------

    if (numNodes < 30) {
        // Call at most 5 times at each node
        if (preNode_ != numNodes) { // Reset count since at different node
            nodeCalls_ = 0;
        }
        if (nodeCalls_ > 4) {
            return false;
        }
    }
    else if ( (numNodes % 1 != 0) || (preNode_ == numNodes) ){
        // If nodes > 30, do every node and only once at a node
	return false;
    }

    if (msgLevel > 200) {
	std::cout << "***** Call "<< name_ << " heuristic at node " 
                  << numNodes << std::endl;
    }

    //------------------------------------------------------
    // Declaration 
    //------------------------------------------------------
    
    bool foundSolution = true;
    
    int j, k, len, tourLen;
    int currV, nextV;
    int count = 0;

    CoinPackedVector * vList = NULL;

    int numVertices = model->getNumVertices();

    double upperBound = model->getKnowledgeBroker()->getIncumbentValue();

    objectiveValue  = 0.0;
    const double *objCoef = model->getObjCoef();
    double *lpSolution = NULL;   
    
    //------------------------------------------------------
    // Start heuristic
    //------------------------------------------------------

#if 0
    for (k = 0; k < numEdges; ++k) {
        int v0 = edges[k]->getv0();
        int v1 = edges[k]->getv1();
        std::cout << "edge " << k << " : v0 = " << v0 << " ; v1 = " << v1 
                  << std::endl;
    }
#endif

    preNode_ = numNodes;
    
    memset(newSolution, 0, numEdges * sizeof(double));

    if (strategy_ == BlisHeurStrategyPeriodic) {
	lpSolution = const_cast<double *>(model_->getLpSolution());
	if (!neighbors_) {
	    neighbors_ = new int [numVertices];
	}
	for (k = 0; k < numVertices; ++k) {
	    neighbors_[k] = -1;
	}
	for (k = 0; k < numEdges; ++k) {
	    if (lpSolution[k] > 0.9) {
		double sh1 = CoinDrand48();
		double sh2 = CoinDrand48();
		int v0 = edges[k]->getv0();
		int v1 = edges[k]->getv1();
		if (sh1 > 0.5) {
		    if (neighbors_[v0] < 0 || sh2 < 0.5) {
			neighbors_[v0] = v1;
		    }
		}
		else {
		    if (neighbors_[v1] < 0 || sh2 < 0.5) {
			neighbors_[v1] = v0;		    
		    }
		}
	    }
	}
    }

    // Allocated memory and clear buffer.
    if (!visited_) {
	visited_ = new bool [numVertices];
    }
    for (k = 0; k < numVertices; ++k) {
	visited_[k] =  false;
    }
    tour_.clear();

    //------------------------------------------------------
    // 1. Stand on an arbitrary vertex as current vertex
    //------------------------------------------------------

    currV = (int)(CoinDrand48() * numVertices);
    visited_[currV] = true;

    tour_.push_back(currV);
    tourLen = tour_.size();
    
    while (tourLen < numVertices) {
	nextV = -1;
	if (neighbors_) {
	    // 2. Find nextV for pre-determined neighbor
	    nextV = neighbors_[currV];
	    if (nextV > -1 && !visited_[nextV]) {
		tour_.push_back(nextV);
		tourLen = tour_.size();
		// 3. Set current vertex be V.
		currV = nextV;
		// 4. Mark V as visited.
		visited_[nextV] = true;
	    }
	    else {  // Didn't find nextV
		nextV = -1;
	    }
	}
	if (nextV == -1) {
	    // 2. Find out the lightest edge connecting current vertex and 
	    //    an unvisited vertex V.
	    vList = adjList_[currV];
	    len = vList->getNumElements();
	    for (k = 0; k < len; ++k) {
		nextV = (vList->getIndices())[k];
		if (!visited_[nextV]) {
		    tour_.push_back(nextV);
		    tourLen = tour_.size();
		    // 3. Set current vertex be V.
		    currV = nextV;
		    // 4. Mark V as visited.
		    visited_[nextV] = true;
		    break;
		}
	    }
	}
    }
    
    if (msgLevel > 200 && foundSolution) {
    	std::cout << "Found a tour " << std::endl;
    }
    
    //------------------------------------------------------
    // Transfer the tour into a LP solution.
    //------------------------------------------------------
    
    int beg, end;
    int tourLenAddOne = tourLen + 1;
    count = 0;
    currV = tour_[0];
    for (k = 1; k < tourLenAddOne; ++k) {
        if (k < tourLen) {
            nextV = tour_[k];
        }
        else {  // edge {currV, tour_[0]}
            nextV = tour_[0];
        }
        
	// Edge is stored in the way that beginning vertex index is 
	// less than the end vertex index, like {3, 12}
	if (currV < nextV) {
	    beg = currV;
	    end = nextV;
	}
	else {
	    beg = nextV;
	    end = currV;
	}

	// Find edge(variable) {beg, end} in the edge(variable) list and 
	// set it to 1.0.
	j = edgeColMatch_[end-1][beg];
	++count;
	newSolution[j] = 1.0;
	objectiveValue += objCoef[j];
	if (objectiveValue >= upperBound) {
	    foundSolution = false;
            if (msgLevel > 200) {
                std::cout << "Discard the found tour because cost "
                          << objectiveValue << " is worse than UB "
                          << upperBound << std::endl;
            }
            break;
	}
	currV = nextV;
    }

#if 0
    if (msgLevel > 10 && foundSolution && numNodes < 2) {
	std::cout << "***** "<< name_ 
		  << " heuristic found a solution, quality " 
		  << objectiveValue << std::endl;
    }
#endif
    
    if (msgLevel > 200 && foundSolution) {

	std::cout << "Tour : len = " << tourLen 
		  << " ; cost = " << objectiveValue << std::endl;
	assert(count == tourLen);
	
	count = 0;
	for (k = 0; k < tourLen; ++k) {
	    std::cout << tour_[k]+1 << " ";
	    ++count;
	    if (count % 15 == 0) {
		std::cout << std::endl;
	    }
	}
	std::cout << std::endl;
    }

    // Add count by 1
    addCalls();

    if (strategy_ == BlisHeurStrategyBeforeRoot) {
        // Free adjacent list to solve memory
        freeGuts();
    }
    
    return foundSolution;
}

//#############################################################################

void
VrpHeurTSP::createAdjList(VrpModel *model) 
{
    int k, beg, end;
    double cost;

    VrpVariable *var = NULL;
    std::vector<VrpVariable *> edges = model->getEdgeList();
    int numVertices = model->getNumVertices();
    int numEdges = edges.size();

    assert(numEdges == model->getNumEdges());

    int msgLevel = model->AlpsPar()->entry(AlpsParams::msgLevel);
    AlpsTimer timer;

    if (msgLevel > 200) {
	timer.start();
    }

    // Allocate memory
    for (k = 0; k < numVertices; ++k) {
	CoinPackedVector * vec = new CoinPackedVector;
	adjList_.push_back(vec);
    }

    // Create adj list
    for (k = 0; k < numEdges; ++k) {
	var = edges[k];
	cost = var->getObjCoef();
        beg = var->getv0();
        end = var->getv1();
	
	// Symetric
	adjList_[beg]->insert(end, cost);
	adjList_[end]->insert(beg, cost);
    }    

    // Sort the adjacent list of each vertex incremently
    // based on weights.
    for (k = 0; k < numVertices; ++k) {
	adjList_[k]->sortIncrElement();
	if (msgLevel > 300) {
	    std::cout << "vertex " << k << " : size " 
		      << adjList_[k]->getNumElements() << std::endl;
	}
    }
    
    if (msgLevel > 200) {
	// Creating adjlist takes the half time of solving a LP.
	std::cout << "createAdjList took " << timer.getCpuTime() 
		  << " seconds." << std::endl;
	timer.start();
    }

    // Create edge and column relationship
    edgeColMatch_ = new std::vector<int> [numVertices - 1];
    int emptySpace = 0;
    int col = 0;
    int v1 = 0;
    int j;
    for (k = 1; k< numVertices; ++k) {
	emptySpace += (numVertices - k);
	v1 = k - 1;
	for (j = 0; j < k; ++j) {
	    col = k * numVertices - k + j - emptySpace;
	    edgeColMatch_[v1].push_back(col);
	    //std::cout << "{" << k << "," << j << "}: " << col << " ; ";
	}
	//std::cout << std::endl;
    }

    if (msgLevel > 200) {
	// Creating adjlist takes the half time of solving a LP.
	std::cout << "create edge col relation took " << timer.getCpuTime() 
		  << " seconds." << std::endl;
    }
}

//#############################################################################
