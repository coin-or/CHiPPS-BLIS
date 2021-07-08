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

#include "BcpsObjectPool.h"
#include "BlisConstraint.h"
#include "BlisHelp.h"

#include "VrpModel.h"
#include "VrpCutGenerator.h"
#include "VrpMacros.h"
#include "VrpParams.h"

#define DELETE_POWER 3
#define DELETE_AND 0x07

/*===========================================================================*/

VrpCutGenerator::VrpCutGenerator(VrpModel *vrp, int vertnum)
{
   model_ = vrp;
   if (vertnum){
      ref_ = new int[model_->vertnum_];
      cutVal_ = new double[model_->vertnum_];
      cutList_ = new char[((model_->vertnum_ >> DELETE_POWER) + 1)*
		 (model_->VrpPar_->entry(VrpParams::maxNumCutsInShrink) + 1)];
      inSet_ = new char[model_->vertnum_];
   }else{
      ref_ = 0;
      cutVal_ = 0;
      cutList_ = 0;
      inSet_ = 0;
   }
   SRANDOM(1);
   setName("VRP");
}

/*===========================================================================*/

// Return if need resolve LP immediately.
// New constraints are stored in BcpsConstraintPool
bool
VrpCutGenerator::generateConstraints(BcpsConstraintPool &conPool) 
{
   int vertnum = model_->vertnum_;
   int rcnt, cur_bins = 0, i, k, max_node;
   double cur_slack = 0.0, node_cut, max_node_cut;
   int cut_size = (vertnum >> DELETE_POWER) + 1, num_cuts = 0;
   elist *cur_edge = NULL;
   VrpParams *par = model_->VrpPar_;
   int which_connected_routine = par->entry(VrpParams::whichConnectedRoutine);
   bool do_greedy = par->entry(VrpParams::doGreedy);
   double etol = model_->etol_;

   int *demand = model_->demand_;
   elist *cur_edge1 = NULL, *cur_edge2 = NULL;
   int node1 = 0, node2 = 0;
   int total_demand = demand[0], num_trials = 0;
   int type, rhs, capacity = model_->capacity_;
   VrpNetwork *n = model_->n_;

   // Get dense solution	
   const double *denseSol = model_->getLpSolution();
   // Transform it to a sparse vector.
   CoinPackedVector *sol = model_->getSolution(denseSol);
   model_->createNet(sol);
   
   if (n->isIntegral_){
      /* if the network is integral, check for connectivity */
      n->connected();
      delete sol;
      return connectivityCuts(conPool) ? true: false;
   }
   
#ifdef DO_TSP_CUTS
   if (par->tspProb){
      delete sol;
      return tspCuts(model_, conPool) ? true : false;
   }
#endif

   vertex *verts = n->verts_;

   if (which_connected_routine == BOTH) which_connected_routine = CONNECTED;
      
   int *compnodes_copy = new int[vertnum + 1];
   int *compnodes = n->compNodes_;
   double *compcuts = n->compCuts_;
   int *compdemands = n->compDemands_;

   do{
      /*------------------------------------------------------------------*\
       * Get the connected components of the solution graph without the
       * depot and see if the number of components is more than one
       \*------------------------------------------------------------------*/
      rcnt = which_connected_routine == BICONNECTED ?
	     n->biconnected() : n->connected();

      /* copy the arrays as they will be needed later */
      if (!which_connected_routine && do_greedy){
	 compnodes_copy = (int *) memcpy((char *)compnodes_copy, 
					 (char*)compnodes,
					 (vertnum + 1)*sizeof(int));
	 compnodes = compnodes_copy;
      }
      if (rcnt > 1){
	 /*---------------------------------------------------------------*\
	  * If the number of components is more then one, then check each
	  * component to see if it violates a capacity constraint
	  \*---------------------------------------------------------------*/
	 
	 char **coef_list = new char *[rcnt];
	 memset(coef_list, 0, rcnt*sizeof(char *));
	 coef_list[0] = new char[rcnt * cut_size];
	 memset(coef_list[0], 0, rcnt*cut_size*sizeof(char));
	 for(i = 1; i < rcnt; i++)
	    coef_list[i] = coef_list[0]+i*cut_size;
	 
	 for(i = 1; i < vertnum; i++)
	    (coef_list[(verts[i].comp)-1][i >> DELETE_POWER]) |=
	       (1 << (i & DELETE_AND));
	 
	 for (i = 0; i < rcnt; i++){
	    if (compnodes[i+1] < 2) continue;
	    /*check ith component to see if it violates a constraint*/
	    if (which_connected_routine == BOTH &&
		which_connected_routine == BICONNECTED && compcuts[i+1]==0)
	       continue;
	    if (compcuts[i+1] < 2*BINS(compdemands[i+1], capacity)-etol){
	       /*the constraint is violated so impose it*/
	       type = (compnodes[i+1] < vertnum/2 ?
		       SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	       rhs = (type == SUBTOUR_ELIM_SIDE ?
		      RHS(compnodes[i+1],compdemands[i+1],
			  capacity): 2*BINS(compdemands[i+1],capacity));
	       num_cuts += addVrpCut(conPool, coef_list[i], rhs, type);
	    }
	    else{/*if the constraint is not violated, then try generating a
		   violated constraint by deleting customers that don't
		   change the number of trucks required by the customers in
		   the component but decrease the value of the cut*/
	       cur_bins = BINS(compdemands[i+1], capacity);/*the current
							     number of trucks 
							     required*/
	       /*current slack in the constraint*/
	       cur_slack = (compcuts[i+1] - 2*cur_bins);
	       while (compnodes[i+1]){/*while there are still nodes in the
					component*/
		  for (max_node = 0, max_node_cut = 0, k = 1;
		       k < vertnum; k++){
		     if (verts[k].comp == i+1){
			if (BINS(compdemands[i+1]-verts[k].demand, capacity)
			    == cur_bins){
			   /*if the number of trucks doesn't decrease upon
			     deleting this customer*/
			   for (node_cut = 0, cur_edge = verts[k].first;
				cur_edge; cur_edge = cur_edge->next_edge){
			      node_cut += (cur_edge->other_end ?
					   -cur_edge->data->weight :
					   cur_edge->data->weight);
			   }
			   if (node_cut > max_node_cut){/*check whether the
							  value of the cut 
							  decrease is the best
							  seen so far*/
			      max_node = k;
			      max_node_cut = node_cut;
			   }
			}
		     }
		  }
		  if (!max_node){
		     break;
		  }
		  /*delete the customer that exhibited the greatest
		    decrease in cut value*/
		  compnodes[i+1]--;
		  compdemands[i+1] -= verts[max_node].demand;
		  compcuts[i+1] -= max_node_cut;
		  cur_slack -= max_node_cut;
		  verts[max_node].comp = 0;
		  coef_list[i][max_node >> DELETE_POWER] ^=
		     (1 << (max_node & DELETE_AND));
		  if (cur_slack < -etol){/*if the cut is now violated, impose
				       it*/
		     type = (compnodes[i+1] < vertnum/2 ?
				      SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
		     rhs = (type == SUBTOUR_ELIM_SIDE ?
			    RHS(compnodes[i+1], compdemands[i+1],
				capacity): 2*cur_bins);
		     num_cuts += addVrpCut(conPool, coef_list[i], rhs, type);
		     break;
		  }
	       }
	    }
	 }
	 delete [] coef_list[0];
	 delete [] coef_list;
      }
      which_connected_routine++;
   }while(!num_cuts && which_connected_routine == BOTH &&
	  which_connected_routine < 2);

   compnodes = n->compNodes_;
   
   if (!do_greedy){
      delete sol;
      return num_cuts ? true : false;
   }

   if (num_cuts < 10 && do_greedy){
      int numroutes = model_->numroutes_;
      char *coef = new char[cut_size];
      for (cur_edge=verts[0].first; cur_edge; cur_edge=cur_edge->next_edge){
	 for (cur_edge1 = cur_edge->other->first; cur_edge1;
	      cur_edge1 = cur_edge1->next_edge){
	    if (cur_edge1->data->weight + cur_edge->data->weight < 1 - etol)
	       continue; 
	    node1 = cur_edge->other_end; 
	    node2 = cur_edge1->other_end;
	    for (cur_edge2 = verts[node2].first; cur_edge2;
		 cur_edge2 = cur_edge2->next_edge){
	       if (!(cur_edge2->other_end) && node2){
		  if ((BINS(total_demand - demand[node1] - demand[node2],
			    capacity) > numroutes -1) &&
		      (cur_edge1->data->weight + cur_edge->data->weight +
		       cur_edge2->data->weight>2+etol)){
		     type = SUBTOUR_ELIM_ACROSS;
		     rhs =2*BINS(total_demand - demand[node1] -
				 demand[node2],capacity);
		     memset(coef, 0, cut_size);
		     for (i = 1; i <vertnum ; i++)
			if ((i != node1) && (i != node2))
			   (coef[i >> DELETE_POWER]) |= (1 << (i&DELETE_AND));
		     num_cuts += addVrpCut(conPool, coef, rhs, type);
		  }
		  break; 
	       }
	    }
	 }
      }
      delete [] coef;
   }

   //n->compNodes_ = compnodes_copy;

   if (num_cuts < 10 && do_greedy){
      
      memcpy((char *)n->newDemand_, (char *)demand, vertnum*sizeof(int));
   
      n->reduce_graph(model_->etol_);
      if (n->numComps_ > 1){
	 num_cuts += greedyShrinking1(
                     model_, par->entry(VrpParams::maxNumCutsInShrink), conPool);
      }else{
	 num_cuts += greedyShrinking1One(
                     model_, par->entry(VrpParams::maxNumCutsInShrink), conPool);
      }
   }

   if (num_cuts < 10 && do_greedy){
      if (par->entry(VrpParams::doExtraInRoot)){
	 num_trials = 
	    //level ? par->entry(VrpParams::greedyNumTrials):
	    2*par->entry(VrpParams::greedyNumTrials);
      }else{
	 num_trials = par->entry(VrpParams::greedyNumTrials);
      }
      if (n->numComps_){
	 num_cuts += greedyShrinking6(
                          model_,par->entry(VrpParams::maxNumCutsInShrink),
			  num_cuts ? num_trials : 2 * num_trials, 10.5, conPool);
      }else{
	 num_cuts += greedyShrinking6One(
                             model_, par->entry(VrpParams::maxNumCutsInShrink),
			     num_cuts ? num_trials : 2 * num_trials, 10.5, conPool);
      }
   }

   n->compNodes_ = compnodes;

   delete[] compnodes_copy;
   delete sol;
      
   return num_cuts ? true : false;
}

/*===========================================================================*/

int
VrpCutGenerator::connectivityCuts(BcpsConstraintPool &conPool)
{
   int vertnum = model_->vertnum_;
   elist *cur_route_start;
   edge *edge_data;
   int weight = 0, reduced_weight, *route;
   int cur_vert = 0, prev_vert, cust_num = 0, cur_route, rcnt;
   int i, reduced_cust_num, vert1, vert2, type, rhs;
   int cut_size = (vertnum >> DELETE_POWER) +1;
   int capacity = model_->capacity_;
   VrpNetwork *n = model_->n_;
   vertex *verts = n->verts_;
   double *compcuts = n->compCuts_;
   int *compnodes = n->compNodes_;
   int *compdemands = n->compDemands_;
   double etol = model_->etol_;
   int num_cuts = 0;

   if (!n->isIntegral_) return 0;

   /* This is a flag to tell the cut generator that the network has not been
      constructed until the next call to userFeasibleSolution();*/
   n->isIntegral_ = false;
   
   /*get the components of the solution graph without the depot to check if the
     graph is connected or not*/
   /* rcnt = n->connected(); */ /* This was previously executed */
   rcnt = n->numComps_;
   char **coef_list = new char *[rcnt];
   memset(coef_list, 0, rcnt*sizeof(char *));
   coef_list[0] = new char[rcnt * cut_size];
   memset(coef_list[0], 0, rcnt*cut_size*sizeof(char));
   for(i = 1; i < rcnt; i++)
      coef_list[i] = coef_list[0]+i*cut_size;
   for(i = 1; i < rcnt; i++)
     coef_list[i] = coef_list[0]+i*cut_size;

  for(i = 1; i < vertnum; i++)
     (coef_list[(verts[i].comp)-1][i >> DELETE_POWER]) |=
	(1 << (i & DELETE_AND));
  
   /*------------------------------------------------------------------------*\
   | For each component check to see if the cut it induces is nonzero -- each |
   | component's cut value must be either 0 or 2 since we have integrality    |
   \*------------------------------------------------------------------------*/
  
  for (i = 0; i < rcnt; i++){
     if (compcuts[i+1] < etol){/*if the cut value is zero, the graph is
				 disconnected and we have a violated cut*/
	type = (compnodes[i+1] < vertnum/2 ?
		SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	rhs = (type == SUBTOUR_ELIM_SIDE ? 
	       RHS(compnodes[i+1], compdemands[i+1], capacity) :
	       2*BINS(compdemands[i+1], capacity));
	num_cuts += addVrpCut(conPool, coef_list[i], rhs, type);
     }
  }

  delete [] coef_list[0];
  delete [] coef_list;

  /*-------------------------------------------------------------------------*\
  | if the graph is connected, check each route to see if it obeys the        |
  | capacity constraints                                                      |
  \*-------------------------------------------------------------------------*/

  int numroutes = model_->numroutes_;
  route = new int[vertnum];
  char *coef = new char[cut_size];
  for (cur_route_start = verts[0].first, cur_route = 0,
	  edge_data = cur_route_start->data; cur_route < numroutes;
       cur_route++){
     edge_data = cur_route_start->data;
     edge_data->scanned = true;
     cur_vert = edge_data->v1;
     prev_vert = weight = cust_num = 0;
     
     memset(coef, 0, cut_size*sizeof(char));

     route[0] = cur_vert;
     while (cur_vert){
	/*keep tracing around the route and whenever the addition
	  of the next customer causes a violation, impose the
	  constraint induced
	  by the set of customers seen so far on the route*/
	coef[cur_vert >> DELETE_POWER]|=(1 << (cur_vert & DELETE_AND));
	cust_num++;
	if ((weight += verts[cur_vert].demand) > capacity){
	   type = (cust_num < vertnum/2 ?
		   SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	   rhs = (type ==SUBTOUR_ELIM_SIDE ? RHS(cust_num, weight, capacity):
		  2*BINS(weight, capacity));
	   num_cuts += addVrpCut(conPool, coef, rhs, type);
	   vert1 = route[0];
	   reduced_weight = weight;
	   reduced_cust_num = cust_num;
	   while (true){
	      if ((reduced_weight -= verts[vert1].demand) > capacity){
		 reduced_cust_num--;
		 coef[vert1 >> DELETE_POWER] &=
		    ~(1 << (vert1 & DELETE_AND));
		 type = (reduced_cust_num < vertnum/2 ?
			 SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
		 rhs = (type ==SUBTOUR_ELIM_SIDE ?
			RHS(reduced_cust_num, reduced_weight, capacity):
			2*BINS(reduced_weight, capacity));
		 num_cuts += addVrpCut(conPool, coef, rhs, type);
		 vert1 = route[vert1];
	      }else{
		 break;
	      }
	   }
	   vert2 = route[0];
	   while (vert2 != vert1){
	      coef[vert2 >> DELETE_POWER] |= (1 << (vert2 & DELETE_AND));
	      vert2 = route[vert2];
	   }
	}
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
	route[prev_vert] = cur_vert;
     }
     edge_data->scanned = true;
     
     while (cur_route_start->data->scanned){/*find the next edge leading out of
					      the depot which has not yet been
					      traversed to start the next
					      route*/
	if (!(cur_route_start = cur_route_start->next_edge)) break;
     }
  }
  
  delete [] route;
  delete [] coef;
  
  for (cur_route_start = verts[0].first; cur_route_start;
       cur_route_start = cur_route_start->next_edge)
     cur_route_start->data->scanned = false;

  // return if need resolve LP immediately.
  return num_cuts;
}

/*===========================================================================*/

int
VrpCutGenerator::greedyShrinking1(VrpModel *m, 
				  int max_shrink_cuts, 
				  BcpsConstraintPool &conPool)
{
   VrpNetwork *n = m->n_;
   double set_cut_val, set_demand;
   vertex *verts = n->verts_;
   elist *e;
   int shrink_cuts = 0, i, j, k;
   char *pt, *cutpt;
   int *ipt; 
   double  *dpt;
   int vertnum = n->vertnum_;
   int truck_cap = m->capacity_, type;
   
   int max_vert = 0, set_size, begin = 1, cur_comp, end = 1, other_end;
   double maxval, weight;
   vertex *cur_nodept;
   int *compmembers = n->compMembers_;
   int *compnodes = n->compNodes_;
   int *demand = n->newDemand_;
   double etol = m->etol_;

   int rhs;
   int size = (vertnum >> DELETE_POWER) + 1;
   char *coef = new char[size];
   memset(coef, 0, size);
   memset(cutList_, 0, size * (max_shrink_cuts + 1));
   
   *inSet_ = 0;
   
   for (i = 1; i < vertnum;  i++){
      if (verts[compmembers[i]].deleted) compmembers[i] = 0;
      ref_[compmembers[i]] = i;
   }
   *ref_ = 0;
   /* ref_ is a reference array for compmembers: gives a place
      in which a vertex is listed in  compmembers */
   
   for (cur_comp = 1; cur_comp <= n->numComps_;
	begin += compnodes[cur_comp], cur_comp++){  /* for every component */
      for (i = begin, end = begin + compnodes[cur_comp]; i < end; i++){
	 if (compmembers[i] == 0) continue;
	 /* for every node as a starting one */
	 /*initialize the data structures */
	 memset(inSet_  + begin, 0, compnodes[cur_comp] * sizeof(char));
	 memset(cutVal_ + begin, 0, compnodes[cur_comp] * sizeof(double));
	 inSet_[i] = 1;
	 set_size = 1 + verts[compmembers[i]].orig_node_list_size; 
	 set_cut_val = 0;     
	 for (e = verts[compmembers[i]].first; e; e = e->next_edge){
	    if (e->other_end)
	       cutVal_[ref_[e->other_end]] = e->data->weight;
	    set_cut_val += e->data->weight;
	 }
	 set_demand = demand[compmembers[i]];  
	 
	 while(true){ 
	    if (set_cut_val < 2*(ceil(set_demand/truck_cap)) - etol &&
		set_size > 2){
	       memset(coef, 0, size*sizeof(char));
	       for (j = begin, ipt = compmembers + begin; j < end; j++, ipt++){
		  if (inSet_[j]){
		     cur_nodept = verts + (*ipt);
		     if (cur_nodept->orig_node_list_size)
			for (k = 0; k < cur_nodept->orig_node_list_size; k++)
			   (coef[(cur_nodept->orig_node_list)[k] >>
				 DELETE_POWER]) |=
			      (1 << ((cur_nodept->orig_node_list)[k] &
				     DELETE_AND));
		     (coef[(*ipt) >> DELETE_POWER]) |=
			(1 << ((*ipt) & DELETE_AND));
		  }  
	       }
	       type = (set_size < vertnum/2 ?
		       SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	       rhs =  (type == SUBTOUR_ELIM_SIDE ?
		       RHS((int)set_size,(int)set_demand,
			   (int)truck_cap):
		       2*BINS((int)set_demand, (int)truck_cap));
	       for (k = 0, cutpt = cutList_; k < shrink_cuts; k++,
		       cutpt += size)
		  if (!memcmp(coef, cutpt, size*sizeof(char)))
		     break;/* same cuts */ 
	       if (k >= shrink_cuts){ 
		  shrink_cuts += addVrpCut(conPool, coef, rhs, type);
		  memcpy(cutpt, coef, size);
	       }
	       if (shrink_cuts > max_shrink_cuts){
		  delete [] coef;
		  return(shrink_cuts);
	       }
	    } 
	    for (maxval = -1, pt = inSet_ + begin, dpt = cutVal_ + begin,
		    j = begin; j < end; pt++, dpt++, j++){
	       if (!(*pt) && *dpt > maxval){
		  maxval = cutVal_[j];
		  max_vert = j; 
	       }
	    }
	    if (maxval > 0){    /* add the vertex to the set */
	       inSet_[max_vert]=1;
	       set_size += 1+ verts[compmembers[max_vert]].orig_node_list_size;
	       set_demand += demand[compmembers[max_vert]];
	       cutVal_[max_vert] = 0;
	       for (e=verts[compmembers[max_vert]].first; e; e = e->next_edge){
		  other_end = ref_[e->other_end];
		  weight = e->data->weight;
		  set_cut_val += (inSet_[other_end]) ? (-weight):weight;
		  cutVal_[other_end] += (inSet_[other_end]) ? 0 : weight;
		  
	       }
	    }
	    else{ /* can't add anything to the set */
	       break;
	    }
	 }   
      }
   }
   
   delete [] coef;

   return(shrink_cuts);
}

/*===========================================================================*/

int
VrpCutGenerator::greedyShrinking6(VrpModel *m, 
				  int max_shrink_cuts, 
				  int trial_num, 
				  double prob,
				  BcpsConstraintPool &conPool)
{
   VrpNetwork *n = m->n_;
   double set_cut_val, set_demand;
   vertex  *verts = n->verts_;
   elist *e;
   int i, j, k, shrink_cuts = 0;
   char *pt, *cutpt;
   double *dpt;
   int vertnum = n->vertnum_, type;
   int truck_cap = m->capacity_; 
  
   int max_vert = 0, set_size, begin = 1, cur_comp, end = 1, num_trials;
   double maxval;
   double denominator=pow(2.0,31.0)-1.0;
   double r, q;
   
   int other_end;
   double weight;
   int *ipt; 
   vertex *cur_nodept;
   int *compmembers = n->compMembers_;
   int *compnodes = n->compNodes_;
   int *demand = n->newDemand_;
   double etol = m->etol_;

   int rhs;
   int size = (vertnum >> DELETE_POWER) + 1;
   char *coef = new char[size];
   memset(coef, 0, size);
   memset(cutList_, 0, size * (max_shrink_cuts +1));
   
   
   *inSet_=0;
   
   for(i = 1; i < vertnum; i++){
      if (verts[compmembers[i]].deleted) compmembers[i] = 0;
      ref_[compmembers[i]] = i;
   }
   *ref_ = 0;  
   
   /* ref_ is a reference array for compmembers: gives a place
      in which a vertex is listed in  compmembers */
   
   for (cur_comp = 1; cur_comp <= n->numComps_; begin += compnodes[cur_comp],
	   cur_comp++){
      /* for every component */
      if (compnodes[cur_comp] <= 7) continue;
      
      for (num_trials = 0; num_trials < trial_num * compnodes[cur_comp];
	num_trials++){
	 end = begin + compnodes[cur_comp];
	 /*initialize the data structures */
	 memset(inSet_ + begin, 0, compnodes[cur_comp] * sizeof(char));
	 memset(cutVal_+ begin, 0, compnodes[cur_comp] * sizeof(double));
	 
	 set_cut_val = 0;
	 set_size = 0;
	 set_demand = 0;
         for (i = begin; i < end; i++ ){
	    if (compmembers[i] == 0) continue;
	    r = (RANDOM()/denominator);
	    q = (prob/compnodes[cur_comp]);
	    if (r < q){
	       inSet_[i] = 1;
	       set_size += 1 + verts[compmembers[i]].orig_node_list_size;
	       set_demand += demand[compmembers[i]];
	       for (e = verts[compmembers[i]].first; e; e = e-> next_edge){
		  other_end = ref_[e->other_end];
		  weight = e->data->weight;
		  set_cut_val += (inSet_[other_end]) ? (-weight) : weight;
		  cutVal_[other_end] += (inSet_[other_end]) ? 0 : weight;
	       }
	    }
	 }
	 while(set_size){ 
	    if (set_cut_val < 2*(ceil(set_demand/truck_cap)) - etol &&
		set_size > 2){
	       memset(coef, 0, size*sizeof(char));
	       for (j = begin, ipt = compmembers + begin; j < end; j++, ipt++)
		  if (inSet_[j]){
		     cur_nodept = verts + (*ipt);
		     if (cur_nodept->orig_node_list_size)
			for (k = 0; k < cur_nodept->orig_node_list_size; k++)
			   (coef[(cur_nodept->orig_node_list)[k] >>
				DELETE_POWER]) |=
			      (1 << ((cur_nodept->orig_node_list)[k] &
				     DELETE_AND));
		     (coef[(*ipt) >> DELETE_POWER]) |= (1 << ((*ipt) &
							      DELETE_AND));
		  }  
	       type = (set_size < vertnum/2 ?
				SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	       rhs =  (type == SUBTOUR_ELIM_SIDE ?
		       RHS((int)set_size,(int)set_demand, (int)truck_cap):
		       2*BINS((int)set_demand, (int)truck_cap));
	       for (k = 0, cutpt = cutList_; k < shrink_cuts; k++,
		       cutpt += size)
		  if (!memcmp(coef, cutpt, size*sizeof(char))) break; 
	       if ( k >= shrink_cuts){
		  shrink_cuts += addVrpCut(conPool, coef, rhs, type);
		  memcpy(cutpt, coef, size);
	       }
	 
	       if ( shrink_cuts > max_shrink_cuts){
		  delete [] coef;
		  return(shrink_cuts);
	       }
	    } 
	    for (maxval = -1, pt = inSet_+begin, dpt = cutVal_+begin,
		    j = begin; j < end; pt++, dpt++, j++){
	       if (!(*pt) && *dpt > maxval){
		  maxval = cutVal_[j];
		  max_vert = j; 
	       }
	    }
	    if (maxval > 0){    /* add the vertex to the set */
	       inSet_[max_vert]=1;
	       set_size+=1+ verts[compmembers[max_vert]].orig_node_list_size;
	       set_demand+=demand[compmembers[max_vert]];
	       cutVal_[max_vert]=0;
	       for (e = verts[compmembers[max_vert]].first; e;
		    e = e->next_edge){
		  other_end = ref_[e->other_end];
		  weight = e->data->weight;
		  set_cut_val += (inSet_[other_end]) ? (-weight) : weight;
		  cutVal_[other_end]+=(inSet_[other_end]) ? 0 : weight;
	       }
	    }
	    else{ /* can't add anything to the set */
	       break;
	    }
	 }   
      }
   }
   
   delete [] coef;
   return shrink_cuts ? true : false;
}

/*===========================================================================*/

int
VrpCutGenerator::greedyShrinking1One(VrpModel *m, 
				     int max_shrink_cuts,
				     BcpsConstraintPool &conPool)
{
   VrpNetwork *n = m->n_; 
   double set_cut_val, set_demand;
   vertex  *verts = n->verts_;
   elist *e;
   int i, j, k, shrink_cuts = 0;
   char *pt, *cutpt;
   double  *dpt;
   int vertnum = n->vertnum_, type;
   int truck_cap = m->capacity_; 
   int max_vert = 0;
   int set_size;
   /* int flag=0; */

   double complement_demand, total_demand = verts[0].demand; 
   double complement_cut_val; 
   int complement_size; 
   double maxval;
   int other_end;
   double weight; 
   vertex *cur_nodept;
   int *demand = n->newDemand_;
   double etol = m->etol_;

   int rhs;
   int size = (vertnum >> DELETE_POWER) + 1;
   char *coef = new char[size];
   memset(coef, 0, size);
   memset(cutList_, 0, size * (max_shrink_cuts + 1));
   
   for (i = 1; i < vertnum; i++ ){
      if (verts[i].deleted) continue;/* for every node as a starting one */
      /*initialize the data structures */
      memset(inSet_, 0, vertnum*sizeof(char));
      memset(cutVal_, 0,vertnum* sizeof(double)); 
      inSet_[i] = 1;
      set_size = 1 + verts[i].orig_node_list_size; 
      set_cut_val = 0;     
      for (e= verts[i].first; e; e = e-> next_edge){
	 weight = e->data->weight;
	 cutVal_[e->other_end] = weight;
	 set_cut_val += weight;
      }
      set_demand = demand[i];  
      
      while(true){ 
	 if (set_cut_val < 2*(ceil(set_demand/truck_cap)) - etol &&
	     set_size > 2){
	    memset(coef, 0, size*sizeof(char));
	    /* printf("%d :", i); */
	    /*  printf("%d ", j); */
	    for (j = 1; j < vertnum; j++)
	       if (inSet_[j]){
		  cur_nodept = verts + j;
		  if (cur_nodept->orig_node_list_size)
		     for (k = 0; k < cur_nodept->orig_node_list_size; k++)
			(coef[(cur_nodept->orig_node_list)[k] >>
			     DELETE_POWER]) |=
			   (1 << ((cur_nodept->orig_node_list)[k] &
				  DELETE_AND));
		  (coef[j>> DELETE_POWER]) |= (1 << ( j & DELETE_AND));
	       }
	    /*  printf("%f ", set_demand);
	    printf("%f \n",set_cut_val);*/ 
	    type = (set_size < vertnum/2 ?
			     SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	    rhs =  (type == SUBTOUR_ELIM_SIDE ?
			     RHS((int)set_size,(int)set_demand,(int)truck_cap):
			     2*BINS((int)set_demand,(int)truck_cap));
	    for (k = 0, cutpt = cutList_; k < shrink_cuts; k++,
		    cutpt += size)
		  if (!memcmp(coef, cutpt, size*sizeof(char)))
		     break; /* same cuts */
	    if ( k >= shrink_cuts){
	       shrink_cuts += addVrpCut(conPool, coef, rhs, type);
	       memcpy(cutpt, coef, size);
	    }
	    
	    if ( shrink_cuts > max_shrink_cuts){
	       delete [] coef;
	       return(shrink_cuts);
	    }
	 }
	 /* check the complement */
	  
	 complement_demand = total_demand - set_demand;
	 complement_cut_val =set_cut_val- 2*(*cutVal_) + 2*m->numroutes_; 
	 complement_size = vertnum - 1 - set_size;   
	 if (complement_cut_val< 2*(ceil(complement_demand/truck_cap))-etol &&
	     complement_size > 2){
	    memset(coef, 0, size*sizeof(char));
	    for (j = 1; j < vertnum; j++)
	       if (!(inSet_[j]) && !(verts[j].deleted)){ 
		  cur_nodept = verts + j;
		  if (cur_nodept->orig_node_list_size)
		     for (k = 0; k < cur_nodept->orig_node_list_size; k++)
			(coef[(cur_nodept->orig_node_list)[k] >>
			     DELETE_POWER]) |=
			   (1 << ((cur_nodept->orig_node_list)[k] &
				  DELETE_AND));
		  (coef[j>> DELETE_POWER]) |= (1 << ( j & DELETE_AND));
	       }
	    type = (complement_size < vertnum/2 ?
			     SUBTOUR_ELIM_SIDE : SUBTOUR_ELIM_ACROSS);
	    rhs =  (type == SUBTOUR_ELIM_SIDE ?
			     RHS((int)complement_size,(int)complement_demand,
				 (int)truck_cap):
			     2*BINS((int)complement_demand,(int)truck_cap));
	    for (k=0, cutpt = cutList_; k < shrink_cuts; k++,
		    cutpt += size)
		  if (!memcmp(coef, cutpt, size*sizeof(char))) break; 
	    if ( k >= shrink_cuts){
	       shrink_cuts += addVrpCut(conPool, coef, rhs, type);
	       memcpy(cutpt, coef, size);
	    }
	 
	    if (shrink_cuts > max_shrink_cuts){
	       delete [] coef;
	       return(shrink_cuts);
	    }
	 }

	 for (maxval = -1, pt = inSet_, dpt = cutVal_,pt++, dpt++,
		 j = 1; j < vertnum; pt++, dpt++, j++){
	    if (!(*pt) && *dpt > maxval){
	       maxval = cutVal_[j];
	       max_vert = j; 
	    }
	 }
	 if (maxval > 0){    /* add the vertex to the set */
	    inSet_[max_vert] = 1;
	    set_size += 1 + verts[max_vert].orig_node_list_size ;
	    set_demand += demand[max_vert];
	    cutVal_[max_vert] = 0;
	    for (e = verts[max_vert].first; e; e = e-> next_edge){
	       other_end = e->other_end;
	       weight = e->data->weight;
	       set_cut_val += (inSet_[other_end]) ? (-weight): weight;
	       cutVal_[other_end] += weight;
	       
	    }
	 }
	 else{ /* can't add anything to the set */
	    break;
	 }
      }   
   }

   delete [] coef;
   return(shrink_cuts);
}

/*===========================================================================*/

int
VrpCutGenerator::greedyShrinking6One(VrpModel *m, 
				     int max_shrink_cuts, 
				     int trial_num, 
				     double prob,
				     BcpsConstraintPool &conPool)
{
   VrpNetwork *n = m->n_;  
   double set_cut_val, set_demand;
   vertex  *verts=n->verts_;
   elist *e;
   int i, j, k, shrink_cuts = 0;
   char *pt, *cutpt;
   double  *dpt;
   int vertnum = n->vertnum_, type;
   int truck_cap = m->capacity_; 
   
   int max_vert = 0, set_size, begin = 1, end = 1, num_trials;
   double maxval, r, q;
   double denominator = pow(2.0, 31.0) - 1.0;
   
   int other_end;
   double weight;

   double complement_demand, total_demand = verts[0].demand;
   double complement_cut_val; 
   int complement_size;
   vertex *cur_nodept; 
   /* int flag=0;*/
   int *demand = n->newDemand_;
   double etol = m->etol_;

   int rhs;
   int size = (vertnum >> DELETE_POWER) + 1;
   char *coef = new char[size];
   memset(coef, 0, size);
   memset(cutList_, 0, size * (max_shrink_cuts +1));
  
   *inSet_ = 0;
 
   for (num_trials = 0; num_trials < trial_num*vertnum ; num_trials++){
      
      /*initialize the data structures */
      memset(inSet_, 0, vertnum*sizeof(char));
      memset(cutVal_, 0,vertnum* sizeof(double)); 
      
      set_cut_val = 0;
      set_size = 0;
      set_demand = 0;
      for (i = 1 ; i < vertnum; i++ ){
	 if (verts[i].deleted) continue;
	 r = (RANDOM()/denominator);
	 q = (prob/vertnum);
	 if (r < q){
	    inSet_[i] = 1;
	    set_size += 1 + verts[i].orig_node_list_size;
	    set_demand += demand[i];
	    for (e = verts[i].first; e; e = e-> next_edge){
		other_end = e->other_end;
		weight  = e->data->weight;
		set_cut_val += (inSet_[other_end]) ? (-weight) : weight;
		cutVal_[other_end] += (inSet_[other_end]) ? 0 : weight;
	    }
	 }
      }
      while(set_size){ 
	 if (set_cut_val < 2*(ceil(set_demand/truck_cap)) - etol &&
	     set_size > 2){
	    memset(coef, 0, size*sizeof(char));
	    for (j = 1; j < vertnum; j++ )
	       if (inSet_[j]){
		  cur_nodept = verts + j;
		  if (cur_nodept->orig_node_list_size)
		     for (k = 0; k < cur_nodept->orig_node_list_size; k++)
			(coef[(cur_nodept->orig_node_list)[k] >>
			     DELETE_POWER]) |=
			   (1 << ((cur_nodept->orig_node_list)[k] &
				  DELETE_AND));
		  (coef[j>> DELETE_POWER]) |= (1 << ( j & DELETE_AND));
	       }
	    type = (set_size < vertnum/2 ?
			     SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	    rhs =  (type == SUBTOUR_ELIM_SIDE ?
			     RHS((int)set_size, (int)set_demand, (int)truck_cap):
			     2*BINS((int)set_demand, (int)truck_cap));
	    for (k = 0, cutpt = cutList_; k < shrink_cuts; k++,
		    cutpt += size)
	       if (!memcmp(coef, cutpt, size*sizeof(char))) break; 
	    if ( k >= shrink_cuts){
	       shrink_cuts += addVrpCut(conPool, coef, rhs, type);
	       memcpy(cutpt, coef, size);
	    }
	    if (shrink_cuts > max_shrink_cuts){
	       delete [] coef;
	       return(shrink_cuts);
	    }
	 }
	 
	 /* check the complement */
	 
	 complement_demand = total_demand - set_demand;
	 complement_cut_val = set_cut_val - 2*(*cutVal_) + 2*m->numroutes_; 
	 complement_size = vertnum - 1 - set_size;   
	 if (complement_cut_val< 2*(ceil(complement_demand/truck_cap))-etol &&
	     complement_size > 2){
	    memset(coef, 0, size*sizeof(char));
	    for (j = 1; j < vertnum; j++)
	       if (!(inSet_[j])&& !(verts[j].deleted)){
		  cur_nodept = verts + j;
		  if (cur_nodept->orig_node_list_size)
		     for (k = 0; k < cur_nodept->orig_node_list_size; k++)
			(coef[(cur_nodept->orig_node_list)[k] >>
			     DELETE_POWER]) |=
			   (1 << ((cur_nodept->orig_node_list)[k] &
				  DELETE_AND));
		  (coef[j>> DELETE_POWER]) |= (1 << ( j & DELETE_AND));
	       }
	    type = (complement_size  < vertnum/2 ?
			     SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	    rhs =  (type == SUBTOUR_ELIM_SIDE ?
			     RHS((int)complement_size,(int)complement_demand,
				 (int)truck_cap):
			     2*BINS((int)complement_demand,(int)truck_cap));
	    for (k = 0, cutpt = cutList_; k < shrink_cuts; k++,
		    cutpt += size)
	       if (!memcmp(coef, cutpt, size*sizeof(char))) break; 
	    if ( k >= shrink_cuts){
	       shrink_cuts += addVrpCut(conPool, coef, rhs, type);
	       memcpy(cutpt, coef, size);
	    }
	    
	    if (shrink_cuts > max_shrink_cuts){
	       delete [] coef;
	       return(shrink_cuts);
	    }
	 }
	 
	 for (maxval = -1, pt = inSet_ + begin, dpt = cutVal_ + begin,
		 j = begin; j < end; pt++, dpt++, j++){
	    if (!(*pt) && *dpt > maxval){
	       maxval = cutVal_[j];
		  max_vert = j; 
	    }
	 }
	 if (maxval > 0){    /* add the vertex to the set */
	    inSet_[max_vert] = 1;
	    set_size += 1 + verts[max_vert].orig_node_list_size ;
	    set_demand += demand[max_vert];
	    cutVal_[max_vert] = 0;
	    for (e = verts[max_vert].first; e; e = e-> next_edge){
	       other_end = e->other_end;
	       weight  = e->data->weight;
	       set_cut_val += (inSet_[other_end]) ? (-weight) : weight;
	       cutVal_[other_end] += weight;
	    }
	 }
	 else{ /* can't add anything to the set */
	    break;
	 }
      }   
   }

   delete [] coef;

   return(shrink_cuts);
}

/*===========================================================================*/

int
VrpCutGenerator::greedyShrinking2One(VrpModel *m, 
				     int max_shrink_cuts, 
				     BcpsConstraintPool &conPool)
{
   VrpNetwork *n = m->n_;  
   double set_cut_val, set_demand;
   vertex *verts = n->verts_;
   elist *e, *cur_edge1, *cur_edge2;
   int j, k, shrink_cuts = 0;
   char *pt;
   double  *dpt;
   int vertnum = n->vertnum_, type;
   int truck_cap = m->capacity_; 
   
   int max_vert = 0, set_size, begin = 1, end = 1;
   double maxval;
   
   int other_end;
   double weight;

   double complement_demand, total_demand = verts[0].demand; 
   double complement_cut_val; 
   int complement_size;
   vertex *cur_nodept;
   int *demand = n->newDemand_;
   double etol = m->etol_;

   int rhs;
   int size = (vertnum >> DELETE_POWER) + 1;
   char *coef = new char[size];
   memset(coef, 0, size);
  
   *inSet_=0;
   
   for (cur_edge1 = verts[0].first; cur_edge1;
	cur_edge1 = cur_edge1->next_edge){
      for (cur_edge2 = cur_edge1->next_edge; cur_edge2;
	   cur_edge2 = cur_edge2->next_edge){
	 
	 /*initialize the data structures */
	 memset(inSet_, 0, vertnum*sizeof(char));
	 memset(cutVal_, 0,vertnum* sizeof(double)); 
	 
	 set_cut_val = 2;
	 set_size = 2 + cur_edge1->other->orig_node_list_size +
	    cur_edge2->other->orig_node_list_size;
	 set_demand = demand[cur_edge1->other_end] +
	    demand[cur_edge2->other_end];
	 inSet_[cur_edge1->other_end] = 1;
	 
	 for (e = verts[cur_edge1->other_end].first; e; e = e-> next_edge){
	    cutVal_[e->other_end] += e->data->weight;
	 }
	 
	 inSet_[cur_edge2->other_end] = 1;
	 for (e = verts[cur_edge2->other_end].first; e; e = e-> next_edge){
	    other_end = e->other_end;
	    weight = e->data->weight;
	    set_cut_val += (inSet_[other_end]) ? (-weight) : weight;
	    cutVal_[other_end] += (inSet_[other_end]) ? 0 : weight;
	 }
	 while(set_size){ 
	    if (set_cut_val < 2*(ceil(set_demand/truck_cap)) - etol &&
		set_size > 2){
	       memset(coef, 0, size*sizeof(char));
	       for (j = 1; j < vertnum; j++ )
		  if (inSet_[j]){
		     cur_nodept = verts + j;
		     if (cur_nodept->orig_node_list_size)
			for (k = 0; k < cur_nodept->orig_node_list_size; k++)
			   (coef[(cur_nodept->orig_node_list)[k] >>
				 DELETE_POWER]) |=
			      (1 << ((cur_nodept->orig_node_list)[k] &
				     DELETE_AND));
		     (coef[j >> DELETE_POWER]) |= (1 << (j & DELETE_AND));
		  }
	       type = (set_size < vertnum/2 ?
		       SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	       rhs =  (type == SUBTOUR_ELIM_SIDE ?
		       RHS((int)set_size, (int)set_demand, (int)truck_cap):
		       2*BINS((int)set_demand, (int)truck_cap));
	       shrink_cuts += addVrpCut(conPool, coef, rhs, type);
	    }
	    
	    /* check the complement */
	    
	    complement_demand = total_demand - set_demand;
	    complement_cut_val = set_cut_val - 2*(*cutVal_) + 2*m->numroutes_; 
	    complement_size = vertnum - 1 - set_size;   
	    if (complement_cut_val<2*(ceil(complement_demand/truck_cap))-etol&&
		complement_size > 2){
	       memset(coef, 0, size*sizeof(char));
	       for (j = 1; j < vertnum; j++)
		  if (!inSet_[j]){
		     cur_nodept=verts + j;
		     if (cur_nodept->orig_node_list_size)
			for (k = 0; k < cur_nodept->orig_node_list_size; k++)
			   (coef[(cur_nodept->orig_node_list)[k] >>
				 DELETE_POWER]) |=
			      (1 << ((cur_nodept->orig_node_list)[k] &
				     DELETE_AND));
		     (coef[j>> DELETE_POWER]) |= (1 << ( j & DELETE_AND));
		  }
	       type = (complement_size  < vertnum/2 ?
		       SUBTOUR_ELIM_SIDE:SUBTOUR_ELIM_ACROSS);
	       rhs =  (type == SUBTOUR_ELIM_SIDE ?
		       RHS((int)complement_size,(int)complement_demand,
			   (int)truck_cap):
		       2*BINS((int)complement_demand,(int)truck_cap));
	       shrink_cuts += addVrpCut(conPool, coef, rhs, type);
	    }
	    
	    for (maxval = -1, pt = inSet_+begin, dpt = cutVal_+begin,
		    j = begin; j < end; pt++, dpt++, j++){
	       if (!(*pt) && *dpt > maxval){
		  maxval = cutVal_[j];
		  max_vert = j; 
	       }
	    }
	    if (maxval > 0){    /* add the vertex to the set */
	       inSet_[max_vert] = 1;
	       set_size += 1 + verts[max_vert].orig_node_list_size ;
	       set_demand += demand[max_vert];
	       cutVal_[max_vert] = 0;
	       for (e = verts[max_vert].first; e; e = e-> next_edge){
		  other_end = e->other_end;
		  weight  = e->data->weight;
		  set_cut_val += (inSet_[other_end]) ? (-weight) : weight;
		  cutVal_[other_end] += weight;
	       }
	    }
	    else{ /* can't add anything to the set */
	       break;
	    }
	 }   
      }
   }

   delete [] coef;
   
   return(shrink_cuts);
}

/*===========================================================================*/

int
VrpCutGenerator::addVrpCut(BcpsConstraintPool &conPool, 
			   char *coef, 
			   int rhs, 
			   int type)
{
   int i, nzcnt = 0, nzcnt_side = 0, nzcnt_across = 0;
   int v0, v1; 
   std::vector<VrpVariable *>edges = model_->getEdgeList();
   int *matind = NULL, *matind_across, *matind_side;
   double *matval = NULL;
   int edgenum = model_->edgenum_;
   double infinity = model_->solver()->getInfinity();
   char sense = 'A';

   switch (type){
      /*-------------------------------------------------------------------*\
       * The subtour elimination constraints are stored as a vector of
       * bits indicating which side of the cut each customer is on
       \*-------------------------------------------------------------------*/
      
    case SUBTOUR_ELIM:
      matind_side = new int[edgenum];
      matind_across = new int[edgenum];
      for (i = 0, nzcnt = 0; i < edgenum; i++){
	 v0 = edges[i]->getv0();
	 v1 = edges[i]->getv1();
	 if (coef[v0 >> DELETE_POWER] &
	     (1 << (v0 & DELETE_AND)) &&
	     (coef[v1 >> DELETE_POWER]) &
	     (1 << (v1 & DELETE_AND))){
	    matind_side[nzcnt_side++] = i;
	 }else if ((coef[v0 >> DELETE_POWER] >>
		    (v0 & DELETE_AND) & 1) ^
		   (coef[v1 >> DELETE_POWER] >>
		    (v1 & DELETE_AND) & 1)){
	    matind_across[nzcnt_across++] = i;
	 }
      }
      type = nzcnt_side < nzcnt_across ? SUBTOUR_ELIM_SIDE :
	 SUBTOUR_ELIM_ACROSS;
      switch (type){
       case SUBTOUR_ELIM_SIDE:
	 nzcnt = nzcnt_side;
	 matind = matind_side;
	 rhs = 0; /*RHS(compnodes[i+1],compdemands[i+1], capacity)*/
	 sense = 'L';
	 delete [] matind_across;
	 break;
	 
       case SUBTOUR_ELIM_ACROSS:
	 nzcnt = nzcnt_across;
	 matind = matind_across;
	 rhs = 0; /*2*BINS(compdemands[i+1], capacity)*/
	 sense = 'G';
	 delete [] matind_side;
	 break;
      }
      
      break;
      
      
    case SUBTOUR_ELIM_SIDE:
      matind = new int[edgenum];
      for (i = 0, nzcnt = 0; i < edgenum; i++){
	 v0 = edges[i]->getv0();
	 v1 = edges[i]->getv1();
	 if (coef[v0 >> DELETE_POWER] &
	     (1 << (v0 & DELETE_AND)) &&
	     (coef[v1 >> DELETE_POWER]) &
	     (1 << (v1 & DELETE_AND))){
	    matind[nzcnt++] = i;
	 }
      }
      sense = 'L';
      break;
      
    case SUBTOUR_ELIM_ACROSS:
      matind = new int[edgenum];
      for (i = 0, nzcnt = 0; i < edgenum; i++){
	 v0 = edges[i]->getv0();
	 v1 = edges[i]->getv1();
	 if ((coef[v0 >> DELETE_POWER] >>
	      (v0 & DELETE_AND) & 1) ^
	     (coef[v1 >> DELETE_POWER] >>
	      (v1 & DELETE_AND) & 1)){
	    matind[nzcnt++] = i;
	 }
      }
      sense = 'G';
      break;
      
    default:
      break;
      
   }
   
   matval = new double[nzcnt];
   for (i = nzcnt-1; i >= 0; i--)
      matval[i] = 1;

   BlisConstraint *blisCon = NULL;
   
   if (sense == 'L'){
       blisCon = new BlisConstraint(-infinity, rhs, 
				    -infinity, rhs,
				    nzcnt, matind, matval);
       blisCon->setValidRegion(BcpsValidGlobal);
   }
   else if (sense == 'G'){
       blisCon = new BlisConstraint(rhs, infinity,
				    rhs, infinity,
				    nzcnt, matind, matval);
       blisCon->setValidRegion(BcpsValidGlobal);
   }
   else{
       return 0;
   }

   conPool.addConstraint(blisCon);

   delete [] matind;
   delete [] matval;
   
   return 1;
}

/*===========================================================================*/

#ifdef DO_TSP_CUTS
int
VrpCutGenerator::tspCuts(VrpModel *m, BcpsConstraintPool &conPool)
{
   VrpNetwork *n = m->n_;  
   VrpParams *par = model_->VrpPar_;
   int edgenum = n->edgenum_;
   int vertnum = m->vertnum_;
   int verbosity = par->verbosity;

   edge *edges = n->edges_;
   CCtsp_lpcut_in *tsp_cuts = NULL;
   int *tsp_edgelist = new int[2*edgenum];
   double *tsp_x = new double[edgenum];
   int i, cutnum = 0, cuts_added = 0, rval, seed;
   CCrandstate rstate;
   CCtsp_cutselect *sel = new CCtsp_cutselect;
   CCtsp_tighten_info *stats = new CCtsp_tighten_info;
   memset (sel, 0, sizeof(CCtsp_cutselect));
   memset (stats, 0, sizeof(CCtsp_tighten_info));
   CCtsp_lpgraph g;
      
   sel->connect          = 1;
   if (par->whichTspCuts & SUBTOUR){
      sel->segments         = 1;
      sel->exactsubtour     = 1;
   }
   if (par->whichTspCuts & BLOSSOM){
      sel->fastblossom      = 1;
      sel->ghfastblossom    = 1;
      sel->exactblossom     = 0;
   }
   if (par->whichTspCuts & COMB){
      sel->blockcombs       = 1;
      sel->growcombs        = 0;
      sel->prclique         = 0;
   }
   
   for (i = 0; i < edgenum; i++, edges++){
      tsp_edgelist[i << 1] = edges->v0;
      tsp_edgelist[(i << 1) + 1] = edges->v1;
      tsp_x[i] = edges->weight;
   }

   CCtsp_init_lpgraph_struct (&g);
   CCtsp_build_lpgraph (&g, vertnum, edgenum, tsp_edgelist, (int *) NULL);   
   CCtsp_build_lpadj (&g, 0, edgenum);

   if (sel->connect){
      rval = CCtsp_connect_cuts(&tsp_cuts, &cutnum, vertnum, edgenum,
				tsp_edgelist, tsp_x);
      if (rval) {
	 fprintf(stderr, "CCtsp_connect_cuts failed\n");
	 printf("CCtsp_connect_cuts failed\n");
	 rval = 1;
      }
      if (verbosity > 3)
	 printf("Found %2d connect cuts\n", cutnum);
      if (!rval && cutnum > 0){
	 cuts_added += addTspCuts(m, conPool, &tsp_cuts, &g);
	 if (cuts_added){
	    if (verbosity > 3)
	       printf("%i connect cuts added\n", cuts_added);
	    goto CLEANUP;
	 }
      }
   }

   if (sel->segments){
      rval = CCtsp_segment_cuts(&tsp_cuts, &cutnum, vertnum, edgenum,
				tsp_edgelist, tsp_x);
      if (rval) {
	 fprintf(stderr, "CCtsp_segment_cuts failed\n");
	 printf("CCtsp_segment_cuts failed\n");
	 rval = 1;
      }
      if (verbosity > 3)
	 printf("Found %2d segment cuts\n", cutnum);
      if (!rval && cutnum > 0){
	 cuts_added += addTspCuts(m, conPool, &tsp_cuts, &g);
	 if (cuts_added){
	    if (verbosity > 3)
	       printf("%i segment cuts added\n", cuts_added);
	    goto CLEANUP;
	 }
      }
    }

   if (sel->fastblossom){
      rval = CCtsp_fastblossom(&tsp_cuts, &cutnum, vertnum, edgenum,
			       tsp_edgelist, tsp_x);
      if (rval) {
	 fprintf(stderr, "CCtsp_fastblossom failed\n");
	 printf("CCtsp_fastblossom failed\n");
	 rval = 1;
      }
      if (verbosity > 3)
	 printf("Found %2d fastblossom cuts\n", cutnum);
      if (!rval && cutnum > 0){
	 cuts_added += addTspCuts(m, conPool, &tsp_cuts, &g);
	 if (cuts_added){
	    if (verbosity > 3)
	       printf("%i fastblossom cuts added\n", cuts_added);
	    goto CLEANUP;
	 }
      }
   }

   if (sel->ghfastblossom){
      rval = CCtsp_ghfastblossom(&tsp_cuts, &cutnum, vertnum, edgenum,
				 tsp_edgelist, tsp_x);
      if (rval) {
	 fprintf(stderr, "CCtsp_ghfastblossom failed\n");
	 printf("CCtsp_ghfastblossom failed\n");
	 rval = 1;
      }
      if (verbosity > 3)
	 printf("Found %2d ghfastblossom cuts\n", cutnum);
      if (!rval && cutnum > 0){
	 cuts_added += addTspCuts(m, conPool, &tsp_cuts, &g);
	 if (cuts_added){
	    if (verbosity > 3)
	       printf("%i ghfastblossom cuts added\n", cuts_added);
	    goto CLEANUP;
	 }
      }
   }

   if (sel->blockcombs){
      rval = CCtsp_block_combs(&tsp_cuts, &cutnum, vertnum, edgenum,
			       tsp_edgelist, tsp_x, true);
      if (rval) {
	 fprintf(stderr, "CCtsp_block_combs failed\n");
	 printf("CCtsp_block_combs failed\n");
	 rval = 1;
      }
      if (verbosity > 3)
	 printf("Found %2d block combs\n", cutnum);
      if (!rval && cutnum > 0){
	 cuts_added += addTspCuts(m, conPool, &tsp_cuts, &g);
	 if (cuts_added){
	    if (verbosity > 3)
	       printf("%i block combs added\n", cuts_added);
	    goto CLEANUP;
	 }
      }
   }

   if (sel->growcombs){
      rval = CCtsp_edge_comb_grower(&tsp_cuts, &cutnum, vertnum,
				    edgenum, tsp_edgelist, tsp_x, stats);
      if (rval) {
	 fprintf(stderr, "CCtsp_edge_comb_grower failed\n");
	 printf("CCtsp_edge_comb_grower failed\n");
	 rval = 1;
      }
      if (verbosity > 3)
	 printf("Found %2d grown combs\n", cutnum);
      if (!rval && cutnum > 0){
	 cuts_added += addTspCuts(m, conPool, &tsp_cuts, &g);
	 if (cuts_added){
	    if (verbosity > 3)
	       printf("%i grown combs added\n", cuts_added);
	    goto CLEANUP;
	 }
      }
   }

   if (sel->prclique){
      rval = CCtsp_pr_cliquetree(&tsp_cuts, &cutnum, vertnum,
				 edgenum, tsp_edgelist, tsp_x, stats);
      if (rval) {
	 fprintf(stderr, "CCtsp_pr_cliquetree failed\n");
	 printf("CCtsp_pr_cliquetree failed\n");
	 rval = 1;
      }
      if (verbosity > 3)
	 printf("Found %2d PR cliquetrees\n", cutnum);
      if (!rval && cutnum > 0){
	 cuts_added += addTspCuts(m, conPool, &tsp_cuts, &g);
	 if (cuts_added){
	    if (verbosity > 3)
	       printf("%i PR cliquetrees added\n", cuts_added);
	    goto CLEANUP;
	 }
      }
   }

   if (sel->exactsubtour){
      rval = CCtsp_exact_subtours(&tsp_cuts, &cutnum, vertnum,
				  edgenum, tsp_edgelist, tsp_x);
      if (rval) {
	 fprintf(stderr, "CCtsp_exact_subtours failed\n");
	 printf("CCtsp_exact_subtours failed\n");
	 rval = 1;
      }
      if (verbosity > 3)
	 printf("Found %2d exact subtours\n", cutnum);
      if (!rval && cutnum > 0){
	 cuts_added += addTspCuts(m, conPool, &tsp_cuts, &g);
	 if (cuts_added){
	    if (verbosity > 3)
	       printf("%i exactsubtours added\n", cuts_added);
	    goto CLEANUP;
	 }
      }
   }

   if (sel->exactblossom){
      seed = (int) CCutil_real_zeit ();
      CCutil_sprand(seed, &rstate);
      rval = CCtsp_exactblossom(&tsp_cuts, &cutnum, vertnum, edgenum,
				tsp_edgelist, tsp_x, &rstate);
      if (rval) {
	 fprintf(stderr, "CCtsp_exactblossom failed\n");
	 printf("CCtsp_exactblossom failed\n");
	 rval = 1;
      }
      if (verbosity > 3)
	 printf("Found %2d exactblossoms\n", cutnum);
      if (!rval && cutnum > 0){
	 cuts_added += addTspCuts(m, conPool, &tsp_cuts, &g);
	 if (cuts_added){
	    if (verbosity > 3)
	       printf("%i exact blossoms added\n", cuts_added);
	    goto CLEANUP;
	 }
      }
   }

CLEANUP:

   delete [] stats;
   delete [] tsp_edgelist;
   delete [] tsp_x;
   delete [] sel;

   return(cuts_added);

}

/*===========================================================================*/

int
VrpCutGenerator::addTspCuts(VrpModel *m, BcpsConstraintPool &conPool, 
			    CCtsp_lpcut_in **tsp_cuts, CCtsp_lpgraph *g)
{
#if 1
   int clique_size = (m->vertnum_ >> DELETE_POWER) + 1;
   char *clique_array, *clique_set;
   int i, j, k, size, cliquecount, val, *matind, nzcnt;
   double *matval, rhs;
   BlisConstraint *blisCon = NULL;
   int num_cuts = 0;
   CCtsp_lpcut_in *tsp_cut, *tsp_cut_next;   
   int v0, v1, jj, edgenum = m->edgenum_;;
   std::vector<VrpVariable *>edges = model_->getEdgeList();
   double infinity = model_->solver()->getInfinity();

   for (tsp_cut = *tsp_cuts; tsp_cut; tsp_cut = tsp_cut->next){
      cliquecount = tsp_cut->cliquecount;
      size = cliquecount * clique_size;
      rhs = (cliquecount == 1 ? 0.0 : -((double)cliquecount)/2.0 + 1.0);
      clique_set = clique_array = new char[size];
      memset(clique_array, 0, size);
      
      for (i = 0; i < cliquecount; i++, clique_set += clique_size){
	 for(j = 0; j < tsp_cut->cliques[i].segcount; j++){
	    for(k = tsp_cut->cliques[i].nodes[j].lo;
		k <= tsp_cut->cliques[i].nodes[j].hi; k++){
	       rhs++;
	       clique_set[k >> DELETE_POWER] |= (1 << (k & DELETE_AND));
	    }
	 }
	 /*For each tooth, we want to add |T|-1 to the rhs so we have to
	   subtract off the one here. It subtracts one for the handle too
	   but that is compensated for above*/
	 rhs--;
      }
      matind = new int[cliquecount*edgenum];
      matval = new double[cliquecount*edgenum];
      for (nzcnt = 0, i = 0; i < edgenum; i++){
	 v0 = edges[i]->getv0();
	 v1 = edges[i]->getv1();
	 val = 0;
	 for (jj = 0; jj < cliquecount; jj++){
	    clique_set = clique_array + clique_size * jj;
	    if (clique_set[v0 >> DELETE_POWER] &
		(1 << (v0 & DELETE_AND)) &&
		(clique_set[v1 >> DELETE_POWER]) &
		(1 << (v1 & DELETE_AND))){
	       val += 1;
	    }
	 }
	 if (val){
	    matind[nzcnt] = i;
	    matval[nzcnt++] = val;
	 }
      }

      if (tsp_cut->sense == 'L'){
	 blisCon = new BlisConstraint(-infinity, (double) rhs, 
				      -infinity, (double) rhs, nzcnt, 
				      matind, matval);
      }else{
	 blisCon = new BlisConstraint((double) rhs, infinity, 
				      (double) rhs, infinity, nzcnt, 
				      matind, matval);
      }
	 
      conPool.addConstraint(blisCon);
      num_cuts++;

      delete [] matind;
      delete [] matval;
      delete [] clique_array;
   }

#else

   int nzlist, nzcnt;
   CCtsp_lpcut new_cut;
   int rval = 0, rhs, i, *matind;
   double *matval;
   BlisConstraint *blisCon = NULL;
   int num_cuts = 0;
   CCtsp_lpcut_in *tsp_cut;   
   double infinity = model_->solver()->getInfinity();

   for (tsp_cut = *tsp_cuts; tsp_cut; tsp_cut = tsp_cut->next){
      
      CCtsp_init_lpcut (&new_cut);
      
      new_cut.rhs         = tsp_cut->rhs;
      new_cut.sense       = tsp_cut->sense;
      new_cut.branch      = tsp_cut->branch;
      
      nzlist = CCtsp_lpcut_in_nzlist (g, tsp_cut);
      
      rval = CCtsp_copy_skeleton (&tsp_cut->skel, &new_cut.skel);
      
      rhs = new_cut.rhs;
      for (i=0; i<new_cut.modcount; i++) {
	 rhs += 2*(((int) new_cut.mods[i].mult) - 128);
      }
      
      nzcnt = 0;
      for (i=nzlist; i != -1; i = g->edges[i].coefnext) {
	 if (g->edges[i].coef) nzcnt++;
      }
      
      if (nzcnt != 0) {
	 matind = new int[nzcnt];
	 matval = new double[nzcnt];
	 for (nzcnt = 0; nzlist != -1; nzlist = i) {
	    i = g->edges[nzlist].coefnext;
	    g->edges[nzlist].coefnext = -2;
	    if (g->edges[nzlist].coef) {
	       matind[nzcnt] = nzlist;
	       matval[nzcnt] = g->edges[nzlist].coef;
	       g->edges[nzlist].coef = 0;
	       nzcnt++;
	    }
	 }
	 if (new_cut.sense == 'L'){
	    blisCon = new BlisConstraint(-infinity, (double) rhs, 
					 -infinity, (double) rhs, nzcnt, 
					 matind, matval);
	 }else{
	    blisCon = new BlisConstraint((double) rhs, infinity, 
					 (double) rhs, infinity, nzcnt, 
					 matind, matval);
	 }
	 
	 num_cuts++;
	 
	 blisCon->setValidRegion(BcpsValidGlobal);
	 
	 delete [] matind;
	 delete [] matval;
      } else {
	 printf ("WARNING: Adding an empty cut\n");
      }
   }

#endif

   return(num_cuts);

}
#endif
