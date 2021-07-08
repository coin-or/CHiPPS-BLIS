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

#ifndef VrpNetwork_h_
#define VrpNetwork_h_

#include <vector>

#include "CoinPackedVector.hpp"
#include "VrpConstants.h"
#include "VrpVariable.h"

//#############################################################################

#define OTHER_END(cur_edge, v) \
        (cur_edge->data->v0 == v) ? cur_edge->data->v1 : cur_edge->data->v0

#ifndef MIN
#define MIN(x, y) (x < y ? x : y)
#endif

#ifndef MAX
#define MAX(x, y) (x > y ? x : y) 
#endif

/*-----------------------------------------------------------------------*\
| These are data tructures used in constructing the solution graph used   |
| by the cut generator to locate cuts. The graph is stored using          |
| adjacency lists                                                         |
\*-----------------------------------------------------------------------*/

typedef struct EDGE{
   int v0;      
   int v1;
   int cost;
   double weight;  
   bool scanned;
   bool tree_edge;
   bool deleted;
}edge;

typedef struct ELIST{
   struct ELIST  *next_edge; /* next edge in the edgelist */
   struct EDGE   *data;      /* the data of the edge */
   int            other_end; /* the other end of the edge */
   struct VERTEX *other;
}elist;

typedef struct VERTEX{
  int           enodenum;   /* the node number in the contracted graph */
  int           orignodenum;/* the node number in the original graph */
  struct ELIST *first; /* points to the first edge in the adjacency list */
  struct ELIST *last;  /* points to the last edge in the adjacency list */
  int           comp;  /* contains the component number if the graph is
			  disconnected */
  bool          scanned;
  int           demand; /* contains the demand for this node */
  int           degree; /* contains the degree of the node in the graph */
  int           orig_node_list_size;
  int          *orig_node_list; /* contains a list of the nodes that have been
				   contracted into this node to make a
				   "super node" */
  int           dfnumber;
  int           low;
  bool          is_art_point;
  bool          deleted; 
}vertex;

class VrpNetwork{

   friend class VrpModel;
   friend class VrpCutGenerator;
   friend class VrpSolution;
   
 private:

   int             edgenum_;     /* the number of edges in the graph */
   int             maxEdgenum_;  /* the number of edges allocated */
   int             vertnum_;     /* the number of vertices in the graph */
   bool            isIntegral_;  /* indicates whether the graph is integral or
				    not */
   int             numComps_;     /* number of components */
   struct EDGE    *edges_;       /* the list of edges in the graph */
   struct VERTEX  *verts_;       /* the list of vertices */
   double          mincut_;      /* the value of the current mincut */
   struct ELIST   *adjList_;     /* the array containing the adajacency lists
				    for each node */
   int            *compNodes_;   /* number of nodes in each component */
   int            *compDemands_; /* demand in each component */
   double         *compCuts_;    /* weight of cprresponding cut */
   int            *compMembers_; /* which component each vertex belongs to */
   int            *newDemand_;   /* the amounts of demand for each node to add
				    when the network is contracted */

 public:

   VrpNetwork() : edgenum_(0), vertnum_(0), isIntegral_(false), mincut_(0){
      adjList_ = 0;
      edges_ = 0;
      verts_ = 0;
      compNodes_ = 0;
      compDemands_ = 0;
      compCuts_ = 0;
      compMembers_ = 0;
      newDemand_ = 0;
   }

   VrpNetwork(int edgenum, int vertnum);

   virtual ~VrpNetwork() {
      gutsOfDestructor();
   }

   void createNet(CoinPackedVector *sol, int *demand,
		  std::vector<VrpVariable *> edgeList, double etol,
		  int vertnum);
   
   void computeCompNums(vertex *v, int parent_comp, int *num_comps,
			bool parent_is_art_point);

   void depthFirstSearch(vertex *v, int *count1, int *count2);
   
   int connected();
   
   int biconnected();

   void reduce_graph(double etol);

   void gutsOfDestructor();
   
};

#endif
