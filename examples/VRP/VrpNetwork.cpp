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

#include <cmath>
#include <iostream>

#include "VrpNetwork.h"

/*===========================================================================*/

VrpNetwork::VrpNetwork(int edgenum, int vertnum) : edgenum_(0), 
maxEdgenum_(edgenum), vertnum_(vertnum), isIntegral_(false), mincut_(0){ 

   int i;

   edges_ = new edge[maxEdgenum_];
   verts_ = new vertex[vertnum_];
   for (i = 0; i < vertnum_; i++){
      verts_[i].orig_node_list = new int[vertnum_];
   }
   adjList_ = new elist[2*maxEdgenum_];
   compNodes_ = new int[vertnum_ + 1];
   compDemands_ = new int[vertnum_ + 1];
   compCuts_ = new double[vertnum_ + 1];
   compMembers_ = new int[vertnum_ + 1];
   newDemand_ = new int[vertnum_];
}

/*===========================================================================*/

void 
VrpNetwork::createNet(CoinPackedVector *sol, int *demand,
		      std::vector<VrpVariable *> edgeList, double etol,
		      int vertnum)
{   
   edge *edges = edges_;
   elist *adjlist = adjList_;
   int nv0, nv1, i;
   double *val_low, *val_high, val_aux;
   int *ind_low, *ind_high, ind_aux;
   
   isIntegral_ = true;
   sol->sortIncrElement();
   int *indices = sol->getIndices();
   double *values = sol->getElements();
   edgenum_ = sol->getNumElements();
   for (i = 0; i < vertnum_; i++){
      delete[] verts_[i].orig_node_list;
   }
   memset(edges_, 0, edgenum_*sizeof(edge));
   memset(verts_, 0, vertnum_*sizeof(vertex));
   memset(adjList_, 0, 2*edgenum_*sizeof(elist));
   for (i = 0; i < vertnum_; i++){
      verts_[i].orig_node_list = new int[vertnum_];
   }

   for (i = 0, val_low = values, val_high = values + edgenum_ - 1, 
	   ind_low = indices, ind_high = indices + edgenum_ - 1; 
	i < (edgenum_/2); i++){
      val_aux = *val_low;
      *val_low++ = *val_high;
      *val_high-- = val_aux;
      ind_aux = *ind_low;
      *ind_low++ = *ind_high;
      *ind_high-- = ind_aux;
   }
   
   /*------------------------------------------------------------------------*\
    * set up the adjacency list
   \*------------------------------------------------------------------------*/
   
   for (i = 0; i < edgenum_; i++, values++, indices++){
      if (*values < etol) continue;
      if (fabs(floor(*values+.5) - *values) > etol){
	 isIntegral_ = false;
	 edges->weight = *values;
      }else{
	 edges->weight = floor(*values+.5);
      }
      nv0 = edges->v0 = edgeList[*indices]->getv0();
      nv1 = edges->v1 = edgeList[*indices]->getv1();
      if (!verts_[nv0].first){
	 verts_[nv0].first = verts_[nv0].last = adjlist;
	 verts_[nv0].degree++;
      }
      else{
	 verts_[nv0].last->next_edge = adjlist;
	 verts_[nv0].last = adjlist;
	 verts_[nv0].degree++;
      }
      adjlist->data = edges;
      adjlist->other_end = nv1;
      adjlist->other = verts_ + nv1;
      adjlist++;
      if (!verts_[nv1].first){
	 verts_[nv1].first = verts_[nv1].last = adjlist;
	 verts_[nv1].degree++;
      }
      else{
	 verts_[nv1].last->next_edge = adjlist;
	 verts_[nv1].last = adjlist;
	 verts_[nv1].degree++;
      }
      adjlist->data = edges;
      adjlist->other_end = nv0;
      adjlist->other = verts_ + nv0;
      adjlist++;
      
      edges++;
   }
   
   /*set the demand for each node*/
   for (i = 0; i < vertnum; i++){
      verts_[i].demand = demand[i];
      verts_[i].orignodenum = i;
   }
}

/*===========================================================================*/

void 
VrpNetwork::depthFirstSearch(vertex *v, int *count1, int *count2)
{
   register elist *e;
   bool has_child = false, has_non_tree_edge = false;
   bool is_art_point = false;
   int c, min;
   
   v->scanned = false;
   min = v->dfnumber = ++(*count1);
   for (e = v->first; e; e = e->next_edge){
      if (!e->other_end) continue;
      if (!e->other->scanned){
	 e->data->tree_edge = true;
	 depthFirstSearch(e->other, count1, count2);
	 has_child = true;
      }
      c = e->other->dfnumber;
      if (e->data->tree_edge && (c > v->dfnumber)){
	 if (min > e->other->low)
	    min = e->other->low;
	 if (e->other->low < v->dfnumber)
	    is_art_point = false;
      }else if (!e->data->tree_edge){
	 has_non_tree_edge = true;
	 if (c < v->dfnumber){
	    if (min > c){
	       min = c;
	    }
	 }
      }
   }
   v->low = min;
   if (!has_child && has_non_tree_edge) is_art_point = false;
   v->is_art_point = is_art_point;

   return;
}

/*===========================================================================*/

/*===========================================================================*\
 * Calculates the connected components of the solution graph after removing
 * the depot. Each node is assigned the number of the component in which it
 * resides. The number of nodes in each component is put in "compnodes", the
 * total demand of all customers in the component is put in "compdemands", and
 * the value of the cut induced by the component is put in "compcuts".
\*===========================================================================*/

int 
VrpNetwork::connected()
{
   int cur_node = 0, cur_comp = 0, cur_member = 0, num_nodes_to_scan = 0;
   elist *cur_edge;

   memset(compNodes_, 0, (vertnum_ + 1)*sizeof(int));
   memset(compCuts_, 0, (vertnum_ + 1)*sizeof(double));
   memset(compDemands_, 0, (vertnum_ + 1)*sizeof(int));
      
   int *nodes_to_scan = new int[vertnum_];
   
   while (true){
      for (cur_node = 1; cur_node < vertnum_; cur_node++)
	 if (!verts_[cur_node].comp){/*look for a node that hasn't been 
				       assigned to a component yet*/
	    break;
	 }
      
      if (cur_node == vertnum_) break;/*this indicates that all nodes have been
					assigned to components*/
      
      nodes_to_scan[num_nodes_to_scan++] = cur_node;/*add the first node to the
						      list of nodes to be 
						      scanned*/
      
      compMembers_[++cur_member] = cur_node;
      
      verts_[cur_node].comp = ++cur_comp;/*add the first node into the new
					   component*/
      compNodes_[cur_comp] = 1;
      verts_[cur_node].comp = cur_comp;
      compDemands_[cur_comp] = verts_[cur_node].demand;
      while(true){/*continue to execute this loop until there are no more
		    nodes to scan if there is a node to scan, then add all of
		    its neighbors in to the current component and take it off
		    the list*/
	 for (cur_node = nodes_to_scan[--num_nodes_to_scan],
		 verts_[cur_node].scanned = true,
		 cur_edge = verts_[cur_node].first, 
		 cur_comp = verts_[cur_node].comp;
	      cur_edge; cur_edge = cur_edge->next_edge){
	    if (cur_edge->other_end){
	       if (!verts_[cur_edge->other_end].comp){
		  verts_[cur_edge->other_end].comp = cur_comp;
		  compNodes_[cur_comp]++;
		  compMembers_[++cur_member] = cur_edge->other_end;
		  compDemands_[cur_comp] += verts_[cur_edge->other_end].demand;
		  nodes_to_scan[num_nodes_to_scan++] = cur_edge->other_end;
	       }
	    }
	    else{/*if this node is connected to the depot, then
		   update the value of the cut*/
	       compCuts_[cur_comp] += cur_edge->data->weight;
	    }
	 }
	 if (!num_nodes_to_scan) break;
      }
   }
   
   delete[] nodes_to_scan;

   numComps_ = cur_comp;

   return(cur_comp);
}

/*===========================================================================*/

int 
VrpNetwork::biconnected()
{
   int i;
   elist *e;
   int count1 = 0, count2 = 0;
   int num_comps = 0;
   bool is_art_point;
   
   memset(compNodes_, 0, (vertnum_ + 1)*sizeof(int));
   memset(compCuts_, 0, (vertnum_ + 1)*sizeof(double));
   memset(compDemands_, 0, (vertnum_ + 1)*sizeof(int));
      
   verts_[0].scanned = true;
   verts_[0].comp = 0;
   for (i=1; i<vertnum_; i++)
      verts_[i].scanned = false;

   for(i = 1; i < vertnum_; i++){
      if (!verts_[i].scanned){
	 is_art_point = false;
	 verts_[i].low = verts_[i].dfnumber = ++count1;
	 verts_[i].scanned = true;
	 e = verts_[i].first;
	 if (!e->other_end){
	    if (e->next_edge)
	       e = e->next_edge;
	    else
	       continue;
	 }
	 e->data->tree_edge = true;
	 depthFirstSearch(e->other, &count1, &count2);
	 is_art_point = e->other->is_art_point;
	 for(e = e->next_edge; e; e = e->next_edge){
	    if (!e->other_end) continue;
	    if (!e->other->scanned){
	       is_art_point = true;
	       e->data->tree_edge = true;
	       depthFirstSearch(e->other, &count1, &count2);
	    }
	 }
	 verts_[i].is_art_point = is_art_point;
      }
   }

   for (i=1; i<vertnum_; i++)
      verts_[i].scanned = false;

   for (i = 1; i < vertnum_; i++){
      if (!verts_[i].scanned){
	 verts_[i].scanned = true;
	 verts_[i].comp = ++num_comps;
	 for (e = verts_[i].first;e; e = e->next_edge){
	    if (!e->other_end) continue;
	    if (!e->other->scanned)
	       computeCompNums(e->other, verts_[i].comp, &num_comps,
			       verts_[i].is_art_point);
	 }
      }
   }
      

   for (i = 1; i < vertnum_; i++){
      compNodes_[verts_[i].comp]++;
      compDemands_[verts_[i].comp] += verts_[i].demand;
      for (e = verts_[i].first; e; e = e->next_edge){
	 if (e->other->comp != verts_[i].comp)
	    compCuts_[verts_[i].comp] += e->data->weight;
      }
   }

   return (num_comps);
}

/*===================================================================*/

void 
VrpNetwork::computeCompNums(vertex *v, int parent_comp, int *num_comps,
			    bool parent_is_art_point)
{
   elist *e;

   v->scanned = true;
   if (parent_is_art_point && v->is_art_point)
      v->comp = ++(*num_comps);
   else
      v->comp = parent_comp;
   for (e = v->first;e ; e = e->next_edge){
      if (!e->other_end) continue;
      if (!e->other->scanned)
	 computeCompNums(e->other, v->comp, num_comps,
			   v->is_art_point);
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * This file implements the greedy shrinking algorithm of Augerat, et al.
 * The implementation was done by Leonid Kopman.
\*===========================================================================*/

void 
VrpNetwork::reduce_graph(double etol)
{
   elist *e1, *e2, *e3;
   edge *cur_edge;
   int v1, v2, deg, count, i, k;
   int edges_deleted = 0;
   vertex *v2_pt, *third_node;
   int *demand = newDemand_;

   for (i = 0; i < vertnum_; i++){
      verts_[i].orig_node_list_size = 0;
   }

   while(true){
      edges_deleted = 0;
      for (i = 0; i < edgenum_; i++){
	 cur_edge = edges_ + i;
	 if (cur_edge->weight >= 1 - etol && cur_edge->v0 &&
	     cur_edge->v1 && !(cur_edge->deleted)){
	    cur_edge->deleted = true;
	    edgenum_--;
	    edges_deleted++;
	    v1 = (verts_[cur_edge->v0].degree ==
		  MIN(verts_[cur_edge->v0].degree, 
		      verts_[cur_edge->v1].degree))?
	          cur_edge->v0 : cur_edge->v1;
	    v2 = (v1 == cur_edge->v0) ? cur_edge->v1 : cur_edge->v0;
	    verts_[v1].deleted = true;
	    demand[v2] += demand[v1];
	    demand[v1] = 0;
	    v2_pt = verts_ + v2;
	    v2_pt->degree--;
	    if (v2_pt->first->other_end == v1){
	       v2_pt->first = v2_pt->first->next_edge;
	    }else{
	       for (e3 = v2_pt->first; e3 && e3->next_edge; e3 = e3->next_edge)
		  if (e3->next_edge->other_end == v1){
		     e3->next_edge = e3->next_edge->next_edge;
		     if (e3->next_edge == NULL) v2_pt->last = e3;
		     break;
		  }
	    }
	    
	    (v2_pt->orig_node_list)[(v2_pt->orig_node_list_size)++] = v1;
	    
	    for (k = 0; k < verts_[v1].orig_node_list_size; k++){
	       (v2_pt->orig_node_list)[(v2_pt->orig_node_list_size)++] =
		  (verts_[v1].orig_node_list)[k];
	    }
	    deg = verts_[v1].degree;
	    
	    for (e1=verts_[v1].first, count=0; e1 && (count < deg); count++ ){
	       third_node = e1->other;
	       if (third_node->orignodenum == v2){
		  e1 = e1->next_edge;
		  continue;
	       }
	       for (e2 = v2_pt->first; e2; e2 = e2->next_edge){
		  if (e2->other_end == e1->other_end ){
		     e2->data->weight += e1->data->weight;
		     e1->data->deleted = true;
		     edges_deleted++;
		     (third_node->degree)--;
		     if (third_node->first->other_end == v1){
			third_node->first=third_node->first->next_edge;
		     }else{
			for (e3 = third_node->first; e3 && e3->next_edge;
			     e3 = e3->next_edge)
			   if (e3->next_edge->other_end == v1){
			      e3->next_edge = e3->next_edge->next_edge;
			      if (e3->next_edge == NULL) third_node->last = e3;
			      break;
			   }
		     }
		     break;
		  }
	       }
	       if (e2){
		  e1 = e1->next_edge;
		  continue;
	       }
	       /* ok, so e1->other_node is not incident to v2 */
	       for (e3 = third_node->first; e3 ; e3 = e3->next_edge){
		  if (e3->other_end == v1){
		     e3->other = v2_pt;
		     e3->other_end = v2;
		     e3->data->v0 = MIN(v2, third_node->orignodenum);
		     e3->data->v1 = MAX(v2, third_node->orignodenum);
		     break;
		  }
	       }
	       v2_pt->last->next_edge = e1;
	       v2_pt->last = e1;
	       v2_pt->degree++;
	       e1=e1->next_edge;
	       v2_pt->last->next_edge = NULL;
	    }
	 }
      }
      if (!edges_deleted) break;
   }
}

/*===========================================================================*/

/*===========================================================================*\
 * Free the memory associated with the solution graph data structures
\*===========================================================================*/

void 
VrpNetwork::gutsOfDestructor()
{
  int i;
  
  if (adjList_) delete[] adjList_;
  if (verts_){
     for (i = 0; i < vertnum_; i++)
	if (verts_[i].orig_node_list)
	   delete[] verts_[i].orig_node_list;
     delete[] verts_;
  }
  if (edges_) delete[] edges_;
  if (compNodes_) delete[] compNodes_;
  if (compDemands_) delete [] compDemands_;
  if (compCuts_) delete[] compCuts_; 
  if (compMembers_) delete[] compMembers_; 
  if (newDemand_) delete[] newDemand_;
}
