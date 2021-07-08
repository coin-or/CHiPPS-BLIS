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

#ifndef VrpCommonTypes_h_
#define VrpCommonTypes_h_

//#############################################################################

typedef struct _NODE{
   int next;
   int route;
}_node;

typedef struct ROUTE_DATA{
   int first;
   int last;
   int numcust;
   int weight;
   int cost;
}route_data;

typedef struct BEST_TOURS{
   int algorithm;
   double solve_time;
   int cost;
   int numroutes;
   route_data *route_info;
   _node *tour;
}best_tours;

typedef struct EDGE_DATA{
  int v0;      
  int v1;
  int cost;
}edge_data;

typedef struct SMALL_GRAPH{   /* this gets passed eg. to lin-kerninghan */
   int vertnum;               /* vertnum in the restricted (small) graph */
   int edgenum;               /* edgenum in the restricted (small) graph */
   int allocated_edgenum;
   int del_edgenum;
   edge_data *edges;       /* The data for these edges */
}small_graph;

#endif
