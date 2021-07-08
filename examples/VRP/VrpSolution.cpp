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

#include <iostream>
#include <cmath>

#include "VrpSolution.h"

VrpSolution::VrpSolution(int s, const double *values, double objValue, 
			 VrpModel *vrp) : BlisSolution(s, values, objValue)   
{
   if (!vrp) return;

   int msgLevel = vrp->AlpsPar()->entry(AlpsParams::msgLevel);

   opt_ = new _node[vrp->vertnum_];
   int cur_vert = 0, prev_vert = 0, cur_route;
   elist *cur_route_start = NULL;
   edge *edge_data;
   vertex *verts = vrp->n_->verts_;
   int cost = 0;
   
   /*construct the tour corresponding to this solution vector*/
   for (cur_route_start = verts[0].first, cur_route = 1, cost = 0,
	   edge_data = cur_route_start->data; cur_route <= vrp->numroutes_;
	cur_route++){
      edge_data = cur_route_start->data;
      edge_data->scanned = true;
      cur_vert = edge_data->v1;
      opt_[prev_vert].next = cur_vert;
      opt_[cur_vert].route = cur_route;
      prev_vert = 0;
      cost += vrp->computeCost(prev_vert, cur_vert);
      while (cur_vert){
	 if (verts[cur_vert].first->other_end != prev_vert){
	    prev_vert = cur_vert;
	    edge_data = verts[cur_vert].first->data;
	    cur_vert = verts[cur_vert].first->other_end;
	 }
	 else{
	    prev_vert = cur_vert;
	    edge_data = verts[cur_vert].last->data; /*This statement could
						      possibly be taken out to
						      speed things up a bit*/
	    cur_vert = verts[cur_vert].last->other_end;
	 }
	 opt_[prev_vert].next = cur_vert;
	 opt_[cur_vert].route = cur_route;
	 cost += vrp->computeCost(prev_vert, cur_vert);
      }
      edge_data->scanned = true;
      
      while (cur_route_start->data->scanned){
	 if (!(cur_route_start = cur_route_start->next_edge)) break;
      }
   }

   /* Display the solution (tours). */
   if (msgLevel > 4) {
       std::cout << "**** Found a better solution! Cost: " << cost << std::endl;
   }
   if (msgLevel > 200) {
       print(std::cout);
   }
}

//#############################################################################

void 
VrpSolution::print(std::ostream& os) const 
{   
   int cur_vert = opt_[0].next, prev_vert = 0, count = 0;
   
   if (opt_[0].route == 1) {
       //std::cout << std::endl << "0 ";
       std::cout << "Tour is" << std::endl;
       std::cout << "1 ";  // TSP solution start with 1
   }
   
   while (cur_vert != 0){
      if (opt_[prev_vert].route != opt_[cur_vert].route){
          std::cout << std::endl << "Route #" << opt_[cur_vert].route << ": ";
          count = 0;
      }
      if  (opt_[0].route == 1) {
          std::cout << cur_vert+1 << " ";  // Stardard TSP solution
      }
      else {
          std::cout << cur_vert << " ";    // VRP solution start with 0
      }
      count++;
      if (count > 15){
	 std::cout << std::endl;
	 count = 0;
      }
      prev_vert = cur_vert;
      cur_vert = opt_[cur_vert].next;
   }
   std::cout << std::endl;
}

//#############################################################################


