/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Routines for molecular clock trees and molecular dating */


#include "rates.h"

#ifdef RWRAPPER
#include <R.h>
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Lk(t_tree *tree)
{  
  if(tree->eval_rlnL == NO) return UNLIKELY;

  tree->rates->c_lnL  = .0;
  RATES_Lk_Pre(tree->n_root,tree->n_root->v[2],NULL,tree);
  RATES_Lk_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);

  if(isnan(tree->rates->c_lnL)) assert(false);
  
  return tree->rates->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Lk_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i;
  phydbl log_dens,mu_a,mu_d,r_a,r_d,dt_a,dt_d;
  int n_a,n_d;

  log_dens = -1.;

  if(d->anc != a)
    {
      PhyML_Fprintf(stderr,"\n. d=%d d->anc=%d a=%d root=%d",d->num,d->anc->num,a->num,tree->n_root->num);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      assert(FALSE);
    }

  dt_a = -1.;
  if(a != tree->n_root) dt_a = tree->times->nd_t[a->num] - tree->times->nd_t[a->anc->num];

  if(a == tree->n_root)
    {
      tree->rates->br_r[a->num] = 1.0;
      tree->rates->nd_r[a->num] = 1.0;
    }
  
  mu_a = tree->rates->br_r[a->num];
  r_a  = tree->rates->nd_r[a->num];
  n_a  = tree->times->n_jps[a->num];
  
  dt_d = FABS(tree->times->nd_t[d->num] - tree->times->nd_t[a->num]);
  mu_d = tree->rates->br_r[d->num];
  r_d  = tree->rates->nd_r[d->num];
  n_d  = tree->times->n_jps[d->num];

  log_dens = RATES_Lk_Core(mu_a,mu_d,r_a,r_d,n_a,n_d,dt_a,dt_d,tree);
  tree->rates->c_lnL += log_dens;
  /* PhyML_Printf("\n. logR: %f [%f %f]",tree->rates->c_lnL,mu_a,mu_d); */

  if(isnan(tree->rates->c_lnL))
    {
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      MCMC_Print_Param(tree->mcmc,tree);
      Exit("\n");
    }

  tree->rates->triplet[a->num] += log_dens;

  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Lk_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Lk_Core(phydbl br_r_a, phydbl br_r_d, phydbl nd_r_a, phydbl nd_r_d, int n_a, int n_d, phydbl dt_a, phydbl dt_d, t_tree *tree)
{
  phydbl log_dens,mean,sd,min_r, max_r,cr;
  

  log_dens  = UNLIKELY;
  mean = sd = -1.;
  min_r     = tree->rates->min_rate;
  max_r     = tree->rates->max_rate;
  cr        = tree->rates->clock_r;
  
  if(br_r_d > tree->rates->max_rate) return UNLIKELY;
  if(br_r_d < tree->rates->min_rate) return UNLIKELY;
  
  /* dt_d = MAX(0.5,dt_d); // We give only one decimal when printing out node heights. It is therefore a fair approximation */

  switch(tree->rates->model_id)
    {
    case THORNE :
      {
        int err;
        phydbl log_br_r_d = log(br_r_d);
        phydbl log_br_r_a = log(br_r_a);

        /* log_dens = Log_Dnorm_Trunc(log_br_r_d,log_br_r_a,sqrt(tree->rates->nu*dt_d),log(tree->rates->min_rate),log(tree->rates->max_rate),&err); */
        log_dens = Log_Dnorm(log_br_r_d,log_br_r_a,sqrt(tree->rates->nu*dt_d),&err);
        log_dens -= log_br_r_d;

        break;
      }
    case GUINDON :
      {
	int err;	

	min_r = tree->rates->min_rate;
	max_r = tree->rates->max_rate;

	nd_r_d = log(nd_r_d*cr);
	nd_r_a = log(nd_r_a*cr);
        min_r  = log(min_r*cr);
	max_r  = log(max_r*cr);

	sd   = sqrt(tree->rates->nu*dt_d);
	mean = nd_r_a - .5*sd*sd;

        // log(density(log(Rate)))
 	log_dens = Log_Dnorm_Trunc(nd_r_d,mean,sd,min_r,max_r,&err);
        // log(density(Rate))
        log_dens -= log(exp(nd_r_d)/cr);
        
	if(err)
	  {
	    PhyML_Fprintf(stderr,"\n. Run: %d",tree->mcmc->run);
	    PhyML_Fprintf(stderr,"\n. br_r_d=%f mean=%f sd=%f min_r=%f max_r=%f dt_d=%f",br_r_d,mean,sd,min_r,max_r,dt_d);
	    PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	    Exit("\n");
	  }
	break;
      }
    case LOGNORMAL :
      {
        int err;
        phydbl log_br_r_d = log(br_r_d);
        log_dens = Log_Dnorm(log_br_r_d,-tree->rates->nu*tree->rates->nu/2.,tree->rates->nu,&err);
        log_dens -= log_br_r_d;
        break;
      }
    case STRICTCLOCK :
      {
	log_dens = 0.0;
	break;
      }

    default : 
      {
	PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("");
      }
    }

  if(isnan(log_dens))
    {
      PhyML_Fprintf(stderr,"\n. Run=%4d br_r_d=%f br_r_a=%f dt_d=%f dt_a=%f nu=%f log_dens=%G sd=%f mean=%f\n",
                    tree->mcmc->run,
                    br_r_d,br_r_a,dt_d,dt_a,tree->rates->nu,log_dens,
                    sd,mean);
      assert(false);
    }

  return log_dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Lk_Change_One_Rate(t_node *d, phydbl new_rate, t_tree *tree)
{
  tree->rates->br_r[d->num] = new_rate;
  RATES_Update_Triplet(d,tree);
  RATES_Update_Triplet(d->anc,tree);
  return(tree->rates->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl RATES_Lk_Change_One_Time(t_node *n, phydbl new_t, t_tree *tree)
{  
  if(n == tree->n_root)
    {
      PhyML_Fprintf(stderr,"\n. Moving the time of the root t_node is not permitted.");
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else
    {
      int i;
      
      tree->times->nd_t[n->num] = new_t;

      RATES_Update_Triplet(n,tree);
      
      for(i=0;i<3;i++)
	{
	  if(n->b[i] != tree->e_root) RATES_Update_Triplet(n->v[i],tree);
	  else RATES_Update_Triplet(tree->n_root,tree);
	}
    }
  return(tree->rates->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_Triplet(t_node *n, t_tree *tree)
{
  phydbl curr_triplet,new_triplet;
  phydbl dt0,dt1,dt2;
  phydbl mu1_mu0,mu2_mu0;
  phydbl mu0,mu1,mu2;
  phydbl r0,r1,r2;
  int n0,n1,n2;
  int i;
  t_node *v1,*v2;

  if(n->tax) return;

  curr_triplet = tree->rates->triplet[n->num];

  dt0 = dt1 = dt2 = -100.0;
  r0 = r1 = r2 = 0.0;

  if(n == tree->n_root)
    {
      phydbl log_dens;
      
      log_dens = 0.0;

      dt0 = tree->times->nd_t[tree->n_root->v[2]->num] - tree->times->nd_t[tree->n_root->num];
      dt1 = tree->times->nd_t[tree->n_root->v[1]->num] - tree->times->nd_t[tree->n_root->num];
      
      mu0 = tree->rates->br_r[tree->n_root->v[2]->num];
      mu1 = tree->rates->br_r[tree->n_root->v[1]->num];

      r0 = tree->rates->nd_r[tree->n_root->v[2]->num];
      r1 = tree->rates->nd_r[tree->n_root->v[1]->num];
      
      n0  = tree->times->n_jps[tree->n_root->v[2]->num];
      n1  = tree->times->n_jps[tree->n_root->v[1]->num];

      switch(tree->rates->model_id)
	{
	case COMPOUND_COR : case COMPOUND_NOCOR : 
	  {
	    log_dens  = RATES_Dmu(mu0,n0,dt0,tree->rates->nu,1./tree->rates->nu,tree->rates->lexp,0,1);
	    log_dens *= RATES_Dmu(mu1,n1,dt1,tree->rates->nu,1./tree->rates->nu,tree->rates->lexp,0,1);
	    log_dens  = log(log_dens);
	    break;
	  }
	case EXPONENTIAL : 
	  {
	    log_dens = Dexp(mu0,tree->rates->lexp) * Dexp(mu1,tree->rates->lexp);
	    log_dens = log(log_dens);
	    break;
	  }
	case LOGNORMAL :
	  {
	    log_dens = Dgamma(mu0,tree->rates->nu,1./tree->rates->nu) * Dgamma(mu1,tree->rates->nu,1./tree->rates->nu);
	    log_dens = log(log_dens);
	    break;
	  }
	case THORNE :
	  {
	    int err;
	    phydbl mean0,sd0;
	    phydbl mean1,sd1;

	    
	    sd0 = SQRT(tree->rates->nu*dt0);
	    mean0 = 1.0;

	    sd1 = SQRT(tree->rates->nu*dt1);
	    mean1 = 1.0;

	    log_dens = 
	      Log_Dnorm_Trunc(mu0,mean0,sd0,tree->rates->min_rate,tree->rates->max_rate,&err) + 
	      Log_Dnorm_Trunc(mu1,mean1,sd1,tree->rates->min_rate,tree->rates->max_rate,&err); 

	    break;
	  }
	case GUINDON :
	  {
	    PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	    Exit("\n. Not implemented yet.\n");
	    break;
	  }
	default :
	  {
	    Exit("\n. Model not implemented yet.\n");
	    break;
	  }
	}
      new_triplet = log_dens;

      if(isnan(log_dens) || isinf(FABS(log_dens)))
	{
	  PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	  MCMC_Print_Param(tree->mcmc,tree);
	  Exit("\n");
	}
    }
  else
    {
      mu0 = mu1 = mu2 = -1.;
      n0 = n1 = n2 = -1;

      mu0 = tree->rates->br_r[n->num];
      dt0 = FABS(tree->times->nd_t[n->num] - tree->times->nd_t[n->anc->num]);
      n0  = tree->times->n_jps[n->num];
      r0 = tree->rates->nd_r[n->num];

      v1 = v2 = NULL;
      for(i=0;i<3;i++)
	{
	  if((n->v[i] != n->anc) && (n->b[i] != tree->e_root))
	    {
	      if(!v1)
		{
		  v1  = n->v[i]; 
		  mu1 = tree->rates->br_r[v1->num];
		  dt1 = FABS(tree->times->nd_t[v1->num] - tree->times->nd_t[n->num]);
		  n1  = tree->times->n_jps[v1->num];
		  r1  = tree->rates->nd_r[v1->num];
		}
	      else
		{
		  v2  = n->v[i]; 
		  mu2 = tree->rates->br_r[v2->num];
		  dt2 = FABS(tree->times->nd_t[v2->num] - tree->times->nd_t[n->num]);
		  n2  = tree->times->n_jps[v2->num];
		  r2  = tree->rates->nd_r[v2->num];
		}
	    }
	}
 
      mu1_mu0 = RATES_Lk_Core(mu0,mu1,r0,r1,n0,n1,dt0,dt1,tree);
      mu2_mu0 = RATES_Lk_Core(mu0,mu2,r0,r1,n0,n2,dt0,dt2,tree);
      
      new_triplet = mu1_mu0 + mu2_mu0;
    }

  tree->rates->c_lnL = tree->rates->c_lnL + new_triplet - curr_triplet;
  tree->rates->triplet[n->num] = new_triplet;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Print_Triplets(t_tree *tree)
{
  int i;
  for(i=0;i<2*tree->n_otu-1;++i) PhyML_Printf("\n. Node %3d t=%f",i,tree->rates->triplet[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Print_Rates(t_tree *tree)
{
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[2],NULL,tree);
  RATES_Print_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Copy_Rate_Struct(t_rate *from, t_rate *to, int n_otu)
{
  int i;

  to->lexp = from->lexp;
  to->alpha = from->alpha;
  to->less_likely = from->less_likely;
  
  to->nu = from->nu;
  to->min_nu = from->min_nu;
  to->max_nu = from->max_nu;

  to->min_rate = from->min_rate;
  to->max_rate = from->max_rate;

  to->c_lnL1 = from->c_lnL1;
  to->c_lnL2 = from->c_lnL2;

  to->c_lnL = from->c_lnL;
  
  to->clock_r = from->clock_r;
  to->min_clock = from->min_clock;
  to->max_clock = from->max_clock;

  to->lbda_nu = from->lbda_nu;
  to->min_dt = from->min_dt;
  to->step_rate = from->step_rate;
  to->true_tree_size = from->true_tree_size;
  to->p_max = from->p_max;
  to->norm_fact = from->norm_fact;

  to->adjust_rates = from->adjust_rates;
  to->use_rates = from->use_rates;
  to->bl_from_rt = from->bl_from_rt;
  to->approx = from->approx;
  to->model_id = from->model_id;
  to->is_allocated = from->is_allocated;
  to->met_within_gibbs = from->met_within_gibbs;
  
  to->update_mean_l = from->update_mean_l;
  to->update_cov_l = from->update_cov_l;
  
  to->br_r_recorded = from->br_r_recorded;
  
  to->log_K_cur = from->log_K_cur;
  to->cur_comb_numb = from->cur_comb_numb;

  for(i=0;i<2*n_otu-1;++i) to->nd_r[i] = from->nd_r[i];
  for(i=0;i<2*n_otu-1;++i) to->br_r[i] = from->br_r[i];
  for(i=0;i<2*n_otu-1;++i) to->buff_br_r[i] = from->buff_br_r[i];
  for(i=0;i<2*n_otu-1;++i) to->buff_nd_r[i] = from->buff_nd_r[i];
  for(i=0;i<2*n_otu-1;++i) to->true_r[i] = from->true_r[i];
  for(i=0;i<2*n_otu-2;++i) to->dens[i] = from->dens[i];
  for(i=0;i<2*n_otu-1;++i) to->triplet[i] = from->triplet[i];
  for(i=0;i<(2*n_otu-2)*(2*n_otu-2);++i) to->cov_l[i] = from->cov_l[i];
  for(i=0;i<(2*n_otu-2)*(2*n_otu-2);++i) to->invcov[i] = from->invcov[i];
  for(i=0;i<2*n_otu-2;++i) to->mean_l[i] = from->mean_l[i];
  for(i=0;i<2*n_otu-2;++i) to->ml_l[i] = from->ml_l[i];
  for(i=0;i<2*n_otu-2;++i) to->cur_l[i] = from->cur_l[i];
  for(i=0;i<2*n_otu-3;++i) to->u_ml_l[i] = from->u_ml_l[i];
  for(i=0;i<2*n_otu-3;++i) to->u_cur_l[i] = from->u_cur_l[i];
  for(i=0;i<(2*n_otu-2)*(2*n_otu-2);++i) to->cov_r[i] = from->cov_r[i];
  for(i=0;i<2*n_otu-2;++i) to->cond_var[i] = from->cond_var[i];
  for(i=0;i<2*n_otu-2;++i) to->mean_r[i] = from->mean_r[i];
  for(i=0;i<(2*n_otu-1)*(2*n_otu-1);++i) to->lca[i] = from->lca[i];
  for(i=0;i<(2*n_otu-2)*(2*n_otu-2);++i) to->reg_coeff[i] = from->reg_coeff[i];
  for(i=0;i<(2*n_otu-2)*(6*n_otu-9);++i) to->trip_reg_coeff[i] = from->trip_reg_coeff[i];
  for(i=0;i<(2*n_otu-2)*9;++i) to->trip_cond_cov[i] = from->trip_cond_cov[i];
  for(i=0;i<2*n_otu;++i) to->_2n_vect1[i] = from->_2n_vect1[i];
  for(i=0;i<2*n_otu;++i) to->_2n_vect2[i] = from->_2n_vect2[i];
  for(i=0;i<2*n_otu;++i) to->_2n_vect3[i] = from->_2n_vect3[i];
  for(i=0;i<2*n_otu;++i) to->_2n_vect4[i] = from->_2n_vect4[i];
  for(i=0;i<2*n_otu;++i) to->_2n_vect5[i] = from->_2n_vect5[i];
  for(i=0;i<4*n_otu*n_otu;++i) to->_2n2n_vect1[i] = from->_2n2n_vect1[i];
  for(i=0;i<4*n_otu*n_otu;++i) to->_2n2n_vect2[i] = from->_2n2n_vect2[i];
  for(i=0;i<2*n_otu-1;++i) to->br_do_updt[i] = from->br_do_updt[i];
  for(i=0;i<2*n_otu-1;++i) to->cur_gamma_prior_mean[i] = from->cur_gamma_prior_mean[i];
  for(i=0;i<2*n_otu-1;++i) to->cur_gamma_prior_var[i] = from->cur_gamma_prior_var[i];
  for(i=0;i<2*n_otu-1;++i) to->n_tips_below[i] = from->n_tips_below[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Duplicate_Calib_Struct(t_tree *from, t_tree *to)
{
  int i,j;

  to->times->n_cal = from->times->n_cal;
  for(i=0;i<from->times->n_cal;i++)
    {
      to->times->a_cal[i] = Duplicate_Calib(from->times->a_cal[i]);
      for(j=0;j<from->times->a_cal[i]->clade_list_size;j++)
        Init_Target_Tip(to->times->a_cal[i]->clade_list[j],to);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Print_Rates_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{  
  if((d == tree->n_root->v[2] && d->tax) || (d == tree->n_root->v[1] && d->tax))
    PhyML_Printf("\n. a=%3d ++d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->br_r[d->num],
	   tree->times->nd_t[a->num],tree->times->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->times->nd_t[d->num]-tree->times->nd_t[a->num])*tree->rates->clock_r*tree->rates->br_r[d->num]);
  
  else if((d == tree->n_root->v[2]) || (d == tree->n_root->v[1]))
    PhyML_Printf("\n. a=%3d __d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->br_r[d->num],
	   tree->times->nd_t[a->num],tree->times->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->times->nd_t[d->num]-tree->times->nd_t[a->num])*tree->rates->clock_r*tree->rates->br_r[d->num]);
  else 
    PhyML_Printf("\n. a=%3d   d=%3d rate=%12f t_left=%12f t_rght=%12f ml=%12f l=%12f %12f",
	   a->num,d->num,
	   tree->rates->br_r[d->num],
	   tree->times->nd_t[a->num],tree->times->nd_t[d->num],
	   tree->rates->ml_l[d->num],
	   tree->rates->cur_l[d->num],
	   (tree->times->nd_t[d->num]-tree->times->nd_t[a->num])*tree->rates->clock_r*tree->rates->br_r[d->num]);
  
  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++) 
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Print_Rates_Pre(d,d->v[i],d->b[i],tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Average_Rate(t_tree *tree)
{
  int i;
  phydbl sum;
  sum = 0.0;
  for(i=0;i<2*tree->n_otu-2;++i) sum += tree->rates->br_r[i];
  return sum/(2*tree->n_otu-2);
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Average_Substitution_Rate(t_tree *tree)
{
  phydbl l,dt;
  int i;

  l = 0.0;
  dt = 0.0;
  for(i=0;i<2*tree->n_otu-1;++i)
    {
      if(tree->a_nodes[i] != tree->n_root)
        {
          dt += FABS(tree->times->nd_t[tree->a_nodes[i]->num] - tree->times->nd_t[tree->a_nodes[i]->anc->num]); 
          l  += tree->rates->cur_l[tree->a_nodes[i]->num];
        }
    }
      
  return(l/dt);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Check_Mean_Rates_True(t_tree *tree)
{
  phydbl sum;
  int i;
  
  sum = 0.0;
  for(i=0;i<2*tree->n_otu-2;++i) sum += tree->rates->true_r[i];
  return(sum/(phydbl)(2*tree->n_otu-2));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int RATES_Check_Node_Times(t_tree *tree)
{
  int err;
  err = NO;
  RATES_Check_Node_Times_Pre(tree->n_root,tree->n_root->v[2],&err,tree);
  RATES_Check_Node_Times_Pre(tree->n_root,tree->n_root->v[1],&err,tree);
  return err;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Check_Node_Times_Pre(t_node *a, t_node *d, int *err, t_tree *tree)
{
  if((tree->times->nd_t[d->num] < tree->times->nd_t[a->num]) || (FABS(tree->times->nd_t[d->num] - tree->times->nd_t[a->num]) < 1.E-20))
    {
      PhyML_Printf("\n. a->t=%f d->t=%f",tree->times->nd_t[a->num],tree->times->nd_t[d->num]);
      PhyML_Printf("\n. a->t_prior_min=%f a->t_prior_max=%f",tree->times->t_prior_min[a->num],tree->times->t_prior_max[a->num]);
      PhyML_Printf("\n. d->t_prior_min=%f d->t_prior_max=%f",tree->times->t_prior_min[d->num],tree->times->t_prior_max[d->num]);
      *err = YES;
    }
  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))	  
	  RATES_Check_Node_Times_Pre(d,d->v[i],err,tree);
    }
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Bracket_N_Jumps(int *up, int *down, phydbl param)
{
  phydbl cdf,eps,a,b,c;
  int step;

  step = 10;
  eps = 1.E-10;
  cdf = 0.0;
  c = 1;
  
  while(cdf < 1.-eps)
    {
      c = (int)FLOOR(c * step);
      cdf = Ppois(c,param);      
    }
  
  a = 0.0;
  b = (c-a)/2.;
  step = 0;
  do
    {
      step++;
      cdf = Ppois(b,param);
      if(cdf < eps) a = b;
      else 
	{
	  break;
	}
      b = (c-a)/2.;
    }
  while(step < 1000);
  
  if(step == 1000)
    {
      PhyML_Fprintf(stderr,"\n. a=%f b=%f c=%f param=%f",a,b,c,param);
      PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  *up = c;
  *down = a;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* 
   mu   : average rate of the time period dt
   dt   : time period to be considered
   a    : rate at a given time point is gamma distributed. a is the shape parameter
   b    : rate at a given time point is gamma distributed. b is the scale parameter
   lexp : the number of rate switches is Poisson distributed with parameter lexp * dt
*/ 
/* compute f(mu;dt,a,b,lexp), the probability density of mu. We need to integrate over the
   possible number of jumps (n) during the time interval dt */
phydbl RATES_Dmu(phydbl mu, int n_jumps, phydbl dt, phydbl a, phydbl b, phydbl lexp, int min_n, int jps_dens)
{
  if(n_jumps < 0) /* Marginal, i.e., the number of jumps is not fixed */
    {
      phydbl var,cumpoissprob,dens,mean,poissprob,ab2,gammadens,lexpdt;
      int n,up,down;
      
      var          = 0.0;
      cumpoissprob = 0.0;
      dens         = 0.0;
      n            = 0;
      mean         = a*b;
      ab2          = a*b*b;
      lexpdt       = lexp*dt;  
      
      RATES_Bracket_N_Jumps(&up,&down,lexpdt);
      For(n,MAX(down,min_n)-1) cumpoissprob += Dpois(n,lexpdt,NO);
      
      for(n=MAX(down,min_n);n<up+1;n++)
	{
	  poissprob    = Dpois(n,lexpdt,NO); /* probability of having n jumps */      
	  var          = (2./(n+2.))*ab2; /* var(mu|n) = var(mu|n=0) * 2 / (n+2) */
	  gammadens    = Dgamma_Moments(mu,mean,var);
	  dens         += poissprob * gammadens;
	  cumpoissprob += poissprob;
	  if(cumpoissprob > 1.-1.E-04) break;
	}
      
      if(dens < 1.E-70) dens = 1.E-70;

      return(dens);      
    }
  else /* Joint, i.e., return P(mu | dt, n_jumps) */
    {
      phydbl mean, var, density;


      mean = 1.0;
      var = (2./(n_jumps+2.))*a*b*b;

      if(jps_dens)
	density = Dgamma_Moments(mu,mean,var) * Dpois(n_jumps,dt*lexp,NO);
      else
	density = Dgamma_Moments(mu,mean,var);
      
      if(density < 1.E-70) density = 1.E-70;

      return density;
    }
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Given the times of nodes a (ta) and d (td), the shape of the gamma distribution of instantaneous
   rates, the parameter of the exponential distribution of waiting times between rate jumps and the 
   instantaneous rate at t_node a, this function works out an expected number of (amino-acids or 
   nucleotide) substitutions per site.
*/
void RATES_Expect_Number_Subst(phydbl t_beg, phydbl t_end, phydbl r_beg,  int *n_jumps, phydbl *mean_r, phydbl *r_end, t_rate *rates, t_tree *tree)
{
  phydbl curr_r, curr_t, next_t;

  switch(rates->model_id)
    {
    case COMPOUND_COR:case COMPOUND_NOCOR:
      {
	/* Compound Poisson */
	if(rates->model_id == COMPOUND_COR)
	  {
	    curr_r  = r_beg;
	    *mean_r = r_beg;
	  }
	else
	  {
	    curr_r  = Rgamma(rates->nu,1./rates->nu);;
	    *mean_r = curr_r;
	  }

	curr_t = t_beg + Rexp(rates->lexp); /* Exponentially distributed waiting times */
	next_t = curr_t;
	
	*n_jumps = 0;
	while(curr_t < t_end)
	  {
	    curr_r = Rgamma(rates->nu,1./rates->nu); /* Gamma distributed random instantaneous rate */
	    
	    (*n_jumps)++;
	    
	    next_t = curr_t + Rexp(rates->lexp);
	    
	    if(next_t < t_end)
	      {
		*mean_r = (1./(next_t - t_beg)) * (*mean_r * (curr_t - t_beg) + curr_r * (next_t - curr_t));
	      }
	    else
	      {
		*mean_r = (1./(t_end - t_beg)) * (*mean_r * (curr_t - t_beg) + curr_r * (t_end - curr_t));
	      }
	    curr_t = next_t;
	  }
	
	/*   PhyML_Printf("\n. [%3d %f %f]",*n_jumps,*mean_r,r_beg); */
	
	if(*mean_r < rates->min_rate) *mean_r = rates->min_rate;
	if(*mean_r > rates->max_rate) *mean_r = rates->max_rate;

	*r_end = curr_r;
	break;
      }
    case EXPONENTIAL:
      {
	*mean_r = Rexp(rates->nu);

	if(*mean_r < rates->min_rate) *mean_r = rates->min_rate;
	if(*mean_r > rates->max_rate) *mean_r = rates->max_rate;

	*r_end  = *mean_r;
	break;
      }
    case LOGNORMAL:
      {
	*mean_r = Rgamma(rates->nu,1./rates->nu);

	if(*mean_r < rates->min_rate) *mean_r = rates->min_rate;
	if(*mean_r > rates->max_rate) *mean_r = rates->max_rate;

	*r_end  = *mean_r;
	break;
      }
    case THORNE:
      {
	phydbl sd,mean;
	int err;
	
	sd = SQRT(rates->nu*FABS(t_beg-t_end));
	mean = r_beg;

	*mean_r = Rnorm_Trunc(mean,sd,rates->min_rate,rates->max_rate,&err);

	if(err) PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
	*r_end  = *mean_r;
	break;
      }
    default:
      {
	PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	Exit("\n. Model not implemented yet.\n");
	break;
      }
    }
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Get_Mean_Rates_Pre(t_node *a, t_node *d, t_edge *b, phydbl r_a, t_tree *tree)
{
  phydbl a_t,d_t;
  phydbl mean_r;
  int n_jumps;
  phydbl r_d;

  a_t = tree->times->nd_t[a->num];
  d_t = tree->times->nd_t[d->num];

  mean_r = -1.;
  n_jumps = -1;
  r_d = -1.;
  RATES_Expect_Number_Subst(a_t,d_t,r_a,&n_jumps,&mean_r,&r_d,tree->rates,tree);
  
  tree->rates->br_r[d->num]   = mean_r;
  tree->rates->true_r[d->num] = mean_r;
  tree->times->t_jps[d->num]  = n_jumps;


  /* Move to the next branches */
  if(d->tax) return;
  else
    {
      int i;
      
      for(i=0;i<3;i++)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      RATES_Get_Mean_Rates_Pre(d,d->v[i],d->b[i],r_d,tree);
	    }
	}
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Random_Branch_Lengths(t_tree *tree)
{
  phydbl r0;
   
  r0 = 1.0;

  tree->rates->br_r[tree->n_root->num] = r0;

  RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[2],NULL,r0,tree);
  RATES_Get_Mean_Rates_Pre(tree->n_root,tree->n_root->v[1],NULL,r0,tree);

  RATES_Check_Rates(tree);

  RATES_Update_Edge_Lengths(tree);
  RATES_Initialize_True_Rates(tree);

  tree->n_root_pos = 
    tree->rates->cur_l[tree->n_root->v[2]->num] /
    (tree->rates->cur_l[tree->n_root->v[2]->num] + tree->rates->cur_l[tree->n_root->v[1]->num]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Scale relative rates (on edges) so that they average to one */
void RATES_Normalise_Rates(t_tree *tree)
{
  phydbl scale,dt,T;
  int i;

  /* !!!!!!!!!!!!!!!!!!!!1 */
  return;
  
  scale = 0.0;
  dt = 0.0;
  T = 0.0;
  for(i=0;i<2*tree->n_otu-2;++i)
    {
      assert(tree->a_nodes[i] != tree->n_root);
      dt = fabs(tree->times->nd_t[i] - tree->times->nd_t[tree->a_nodes[i]->anc->num]);
      T+=dt;
      scale += tree->rates->br_r[i] * dt;
    }
  scale /= T;
  for(i=0;i<2*tree->n_otu-2;++i) tree->rates->br_r[i] /= scale;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Test whether rates are normalized */
void RATES_Check_Rates(t_tree *tree)
{
  phydbl scale,dt,T,eps;
  int i;

  eps = 1.E-6;
  
  scale = 0.0;
  dt = 0.0;
  T = 0.0;
  for(i=0;i<2*tree->n_otu-2;++i)
    {
      assert(tree->a_nodes[i] != tree->n_root);
      dt = fabs(tree->times->nd_t[i] - tree->times->nd_t[tree->a_nodes[i]->anc->num]);
      T+=dt;
      scale += tree->rates->br_r[i] * dt;
    }

  scale /= T;
  
  if(scale > 1.+eps || scale < 1.-eps)
    {
      PhyML_Fprintf(stderr,"\n. Relative rates are not normalised!...");
      assert(false);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Init_Triplets(t_tree *tree)
{
  int i;
  for(i=0;i<2*tree->n_otu-1;++i) tree->rates->triplet[i] = 0.0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Initialize_True_Rates(t_tree *tree)
{
  int i;
  for(i=0;i<2*tree->n_otu-2;++i) tree->rates->true_r[i] = tree->rates->br_r[i];
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Get_Rates_From_Bl(t_tree *tree)
{
  phydbl dt,cr;
  t_node *left, *rght;
  int i;

  dt = -1.0;
  cr = tree->rates->clock_r;

  if(tree->n_root)
    {
      dt = FABS(tree->times->nd_t[tree->n_root->num] - tree->times->nd_t[tree->n_root->v[2]->num]);
      tree->rates->br_r[tree->n_root->v[2]->num] = 0.5 * tree->e_root->l->v / (dt*cr);
      dt = FABS(tree->times->nd_t[tree->n_root->num] - tree->times->nd_t[tree->n_root->v[1]->num]);
      tree->rates->br_r[tree->n_root->v[1]->num] = 0.5 * tree->e_root->l->v / (dt*cr);
    }
  

  for(i=0;i<2*tree->n_otu-3;++i)
    {
      if(tree->a_edges[i] != tree->e_root)
	{
	  left = tree->a_edges[i]->left;
	  rght = tree->a_edges[i]->rght;
	  dt = FABS(tree->times->nd_t[left->num] - tree->times->nd_t[rght->num]);	  
	  
	  if(left->anc == rght) tree->rates->br_r[left->num] = tree->a_edges[i]->l->v / (dt*cr);
	  else                  tree->rates->br_r[rght->num] = tree->a_edges[i]->l->v / (dt*cr);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Lk_Jumps(t_tree *tree)
{
  int i,n_jps;
  phydbl dens,dt,lexp;
  t_node *n;

  n = NULL;
  lexp = tree->rates->lexp;
  n_jps = 0;
  dt = 0.0;
  dens = 0.0;

  for(i=0;i<2*tree->n_otu-2;++i)
    {
      n = tree->a_nodes[i];
      dt = FABS(tree->times->nd_t[n->num]-tree->times->nd_t[n->anc->num]);
      n_jps = tree->times->n_jps[n->num];
      dens += Dpois(n_jps,lexp*dt,YES);
    }

  tree->times->c_lnL_jps = dens;
  
  return dens;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Posterior_Rates(t_tree *tree)
{
  int node_num; 
  node_num = Rand_Int(0,2*tree->n_otu-3);      
  RATES_Posterior_One_Rate(tree->a_nodes[node_num]->anc,tree->a_nodes[node_num],NO,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Posterior_Times(t_tree *tree)
{
  int node_num;
  node_num = Rand_Int(tree->n_otu,2*tree->n_otu-3);
  RATES_Posterior_One_Time(tree->a_nodes[node_num]->anc,tree->a_nodes[node_num],NO,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Posterior_One_Rate(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  phydbl like_mean, like_var;
  phydbl prior_mean, prior_var;
  phydbl post_mean, post_var, post_sd;
  phydbl dt,rd,cr,cel,cvl;
  int dim;
  phydbl l_opp; /* length of the branch connected to the root, opposite to the one connected to d */
  t_edge *b;
  int i;
  t_node *v2,*v3;
  phydbl T0,T1,T2,T3;
  phydbl U0,U1,U2,U3;
  phydbl V1,V2,V3;
  phydbl r_min,r_max;
  int err;
  phydbl ratio,u;
  phydbl cvr, cer;
  phydbl nf;
  phydbl new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate;
  phydbl sd1,sd2,sd3;
  phydbl inflate_var;

  if(d == tree->n_root) return;

  dim = 2*tree->n_otu-3;
  err = NO;

  inflate_var = tree->rates->inflate_var;

  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
  else for(i=0;i<3;i++) if(d->v[i] == a) { b = d->b[i]; break; }
  
  v2 = v3 = NULL;
  if(!d->tax)
    {
      for(i=0;i<3;i++)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    if(!v2) { v2 = d->v[i]; }
	    else    { v3 = d->v[i]; }
	  }
    }

  T3 = T2 = 0.0;
  T0 = tree->times->nd_t[a->num];
  T1 = tree->times->nd_t[d->num];
  U0 = tree->rates->br_r[a->num];
  U1 = tree->rates->br_r[d->num];
  U3 = U2 = -1.0;

  if(!d->tax)
    {
      T2  = tree->times->nd_t[v2->num];
      T3  = tree->times->nd_t[v3->num];
      U2  = tree->rates->br_r[v2->num];
      U3  = tree->rates->br_r[v3->num];
    }
  

  V1 = tree->rates->nu * FABS(T1 - T0);
  V2 = tree->rates->nu * FABS(T2 - T1);
  V3 = tree->rates->nu * FABS(T3 - T1);

  dt     = T1 - T0;
  rd     = U1;
  cr     = tree->rates->clock_r;
  r_min  = tree->rates->min_rate;
  r_max  = tree->rates->max_rate;
  nf     = tree->rates->norm_fact;
  sd1    = SQRT(V1);
  sd2    = SQRT(V2);
  sd3    = SQRT(V3);


  /* Likelihood */
  cel=0.0;
  for(i=0;i<dim;i++) 
    if(i != b->num) 
	cel += tree->rates->reg_coeff[b->num*dim+i] * (tree->rates->u_cur_l[i] - tree->rates->mean_l[i]);
  cel += tree->rates->mean_l[b->num];
  cvl = tree->rates->cond_var[b->num];

  l_opp  = 0.0;
  if(a == tree->n_root)
    {
      if(d == a->v[0]) l_opp = tree->rates->cur_l[a->v[1]->num];
      else             l_opp = tree->rates->cur_l[a->v[0]->num]; 
      cel -= l_opp;
    }


  if(isnan(cvl) || isnan(cel) || cvl < .0) 
    {
      for(i=0;i<dim;i++) if(i != b->num) printf("\n. reg: %f %f %f nu=%f clock=%f",
					tree->rates->reg_coeff[b->num*dim+i],
					tree->rates->u_cur_l[i],
					tree->rates->mean_l[i],
					tree->rates->nu,tree->rates->clock_r);
      PhyML_Fprintf(stderr,"\n. cel=%f cvl=%f\n",cel,cvl); 
      PhyML_Fprintf(stderr,"\n. Warning: invalid expected and/or std. dev. values. Skipping this step.\n"); 
      Exit("\n");
    }

  /* Model rates */
  /* if(tree->mod->log_l == YES) cel = exp(cel); */
  /* like_mean = cel / (dt*cr*nf); */
  /* like_var  = cvl / POW(dt*cr*nf,2);  */

  like_mean = cel / (dt*cr*nf);
  like_var  = cvl / POW(dt*cr*nf,2);
  
  /* Prior */
  if(!d->tax)
    {
      cvr = 1./(1./V1 + 1./V2 + 1./V3);
      cer = cvr*(U0/V1 + U2/V2 + U3/V3);
    }
  else
    {
      cvr = V1;
      cer = U0;
    }
  
  if(cvr < 0.0)
    {
      PhyML_Fprintf(stderr,"\n. cvr=%f d->tax=%d V1=%f v2=%f V3=%f",cvr,d->tax,V1,V2,V3);
      PhyML_Fprintf(stderr,"\n. T0=%f T1=%f T2=%f T3=%f",T0,T1,T2,T3);
      Exit("\n");
    }

  prior_mean = cer;
  prior_var  = cvr;

  /* Posterior */
  post_mean = (prior_mean/prior_var + like_mean/like_var)/(1./prior_var + 1./like_var);
  post_var  = 1./(1./prior_var + 1./like_var);


  /* Sample according to priors */
  if(tree->eval_alnL == NO)
    {
      post_mean = prior_mean;
      post_var  = prior_var;
    }

  post_sd = SQRT(post_var);

  rd = Rnorm_Trunc(post_mean,inflate_var*post_sd,r_min,r_max,&err);

  if(err || isnan(rd))
    {
      PhyML_Fprintf(stderr,"\n");
      PhyML_Fprintf(stderr,"\n. run: %d err=%d d->tax=%d",tree->mcmc->run,err,d->tax);
      PhyML_Fprintf(stderr,"\n. rd=%f cvr=%G dt=%G cr=%G",rd,cvr,dt,cr);
      PhyML_Fprintf(stderr,"\n. prior_mean = %G prior_var = %G",prior_mean,prior_var);
      PhyML_Fprintf(stderr,"\n. like_mean = %G like_var = %G",like_mean,like_var);
      PhyML_Fprintf(stderr,"\n. post_mean = %G post_var = %G",post_mean,post_var);
      PhyML_Fprintf(stderr,"\n. clock_r = %f",tree->rates->clock_r);
      PhyML_Fprintf(stderr,"\n. T0=%f T1=%f T2=%f T3=%f",T0,T1,T2,T3);
      PhyML_Fprintf(stderr,"\n. U0=%f U1=%f U2=%f U3=%f",U0,U1,U2,U3);
      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }    

  /* !!!!!!!!!!!!!! */
/*   u = Uni(); */
/*   rd = U1 * exp(1.*(u-0.5)); */

  if(rd > r_min && rd < r_max)
    {
      
      cur_lnL_data = tree->c_lnL;
      cur_lnL_rate = tree->rates->c_lnL;
      new_lnL_data = tree->c_lnL;
      new_lnL_rate = tree->rates->c_lnL;

      tree->rates->br_r[d->num] = rd;
      RATES_Normalise_Rates(tree);
      RATES_Update_Edge_Lengths(tree);
      
      if(tree->eval_alnL) new_lnL_data = Lk(b,tree);
        /* new_lnL_rate = RATES_Lk(tree); */
      new_lnL_rate = 
	cur_lnL_rate - 
	(Log_Dnorm_Trunc(U1,U0,sd1,r_min,r_max,&err)) +
	(Log_Dnorm_Trunc(rd,U0,sd1,r_min,r_max,&err));

      if(!d->tax) 
	{
	  new_lnL_rate -= (Log_Dnorm_Trunc(U2,U1,sd2,r_min,r_max,&err) + Log_Dnorm_Trunc(U3,U1,sd3,r_min,r_max,&err));
	  new_lnL_rate += (Log_Dnorm_Trunc(U2,rd,sd2,r_min,r_max,&err) + Log_Dnorm_Trunc(U3,rd,sd3,r_min,r_max,&err));
	}

      tree->rates->c_lnL = new_lnL_rate;
      
      /* printf("\n. %f %f sd1=%f U1=%f rd=%f ra=%f a=%d d=%d [%f] [%f %f]", */
      /* 	     new_lnL_rate,RATES_Lk(tree),sd1,U1,rd,ra,a->num,d->num,tree->rates->br_r[tree->n_root->num], */
      /* 	     Log_Dnorm_Trunc(U1,ra,sd1,r_min,r_max,&err), */
      /* 	     Log_Dnorm_Trunc(rd,ra,sd1,r_min,r_max,&err)); */

      ratio = 0.0;
      /* Proposal ratio */
      ratio += (Log_Dnorm_Trunc(U1,post_mean,inflate_var*post_sd,r_min,r_max,&err) - Log_Dnorm_Trunc(rd,post_mean,inflate_var*post_sd,r_min,r_max,&err));
      /*   ratio += log(rd/U1); */
      /* Prior ratio */
      ratio += (new_lnL_rate - cur_lnL_rate);
      /* Likelihood ratio */
      ratio += (new_lnL_data - cur_lnL_data);
      
      
      ratio = exp(ratio);
      
      /*   printf("\n. R a=%3d T0=%6.1f T1=%6.1f T2=%6.1f T3=%6.1f ratio=%8f pm=%7f U1=%7.2f rd=%7.2f %f %f lr=%f %f ld=%f %f [%f]",a->num,T0,T1,T2,T3,ratio,post_mean,U1,rd, */
      /* 	 Log_Dnorm_Trunc(U1,post_mean,post_sd,r_min,r_max,&err), */
      /* 	 Log_Dnorm_Trunc(rd,post_mean,post_sd,r_min,r_max,&err), */
      /* 	 new_lnL_rate,cur_lnL_rate, */
      /* 	 new_lnL_data,cur_lnL_data, */
      /* 	 ((Pnorm(r_max,U1,sd2)-Pnorm(r_min,U1,sd2)) * */
      /* 	  (Pnorm(r_max,U1,sd3)-Pnorm(r_min,U1,sd3)))/ */
      /* 	 ((Pnorm(r_max,rd,sd2)-Pnorm(r_min,rd,sd2)) * */
      /* 	  (Pnorm(r_max,rd,sd3)-Pnorm(r_min,rd,sd3)))); */
      
      
      u = Uni();
      
      if(u > MIN(1.,ratio))
	{
	  tree->rates->br_r[d->num] = U1; /* reject */
	  tree->rates->c_lnL        = cur_lnL_rate;
	  tree->c_lnL               = cur_lnL_data;
	  RATES_Update_Edge_Lengths(tree);
	  Update_PMat_At_Given_Edge(b,tree);
	}
      else 
	{
	  tree->mcmc->acc_move[tree->mcmc->num_move_br_r+d->num]++;
	}

      RATES_Normalise_Rates(tree);
      tree->mcmc->acc_rate[tree->mcmc->num_move_br_r+d->num] = 
	(tree->mcmc->acc_move[tree->mcmc->num_move_br_r+d->num]+1.E-6)/
	(tree->mcmc->run_move[tree->mcmc->num_move_br_r+d->num]+1.E-6);
    }

  tree->mcmc->run_move[tree->mcmc->num_move_br_r+d->num]++;
      

  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  for(i=0;i<3;i++)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
	      {
		if(tree->io->lk_approx == EXACT) Update_Partial_Lk(tree,d->b[i],d);
		/* if(tree->io->lk_approx == EXACT) { tree->both_sides = YES; Lk(tree); } */
		RATES_Posterior_One_Rate(d,d->v[i],YES,tree);
	      }
	}
      
      if(tree->io->lk_approx == EXACT) Update_Partial_Lk(tree,b,d);
      /* if(tree->io->lk_approx == EXACT) { tree->both_sides = YES; Lk(tree); } */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Posterior_One_Time(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  /*
           T0 (a)
            |
            | l1,U1,b1
            |
           T1 (d)
           / \
 l2,u2,b2 /   \ L3,u3,b3
         /     \
       t2       t3
      (v2)     (v3)

   */

  phydbl l1,l2,l3;
  phydbl new_l1;
  phydbl El1,El2,El3;
  phydbl u0,r1,r2,r3;
  phydbl t0,t1,t2,t3;
  phydbl t_min, t_max;
/*   phydbl t_max_12,t_max_13; */
  phydbl bl_min, bl_max;
  phydbl t1_new;
  phydbl EX,EY;
  phydbl VX,VY;
  phydbl cr;
  int    i,j;
  phydbl *mu, *cov;
  phydbl *cond_mu, *cond_cov;
  short int *is_1;
  phydbl sig11, sig1X, sig1Y, sigXX, sigXY, sigYY;
  phydbl cov11,cov12,cov13,cov22,cov23,cov33;
  int dim;
  t_edge *b1, *b2, *b3;
  t_node *v2,*v3;
  phydbl *l2XY;
  phydbl l_opp;
  t_edge *buff_b;
  t_node *buff_n;
  int err;
  int num_1, num_2, num_3;
  phydbl nf;
  phydbl u, ratio;
  phydbl new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate;
  int num_move;
  phydbl inflate_var;

  dim = 2*tree->n_otu-3;
  num_move = tree->mcmc->num_move_times;
  inflate_var = tree->rates->inflate_var;

  if(d->tax) return;
  
  if(FABS(tree->times->t_prior_min[d->num] - tree->times->t_prior_max[d->num]) < 1.E-10) return;

  l2XY     = tree->rates->_2n_vect2;
  mu       = tree->rates->_2n_vect3;
  cov      = tree->rates->_2n2n_vect1;
  cond_mu  = tree->rates->_2n_vect1;
  cond_cov = tree->rates->_2n2n_vect2;
  is_1     = tree->rates->_2n_vect5;
  err      = NO;

  b1 = NULL;
  if(a == tree->n_root) b1 = tree->e_root;
  else for(i=0;i<3;i++) if(d->v[i] == a) { b1 = d->b[i]; break; }

  b2 = b3 = NULL;
  v2 = v3 = NULL;
  for(i=0;i<3;i++)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
	if(!v2) { v2 = d->v[i]; b2 = d->b[i]; }
	else    { v3 = d->v[i]; b3 = d->b[i]; }
      }

  t2 = tree->times->nd_t[v2->num];
  t3 = tree->times->nd_t[v3->num];

  buff_n = NULL;
  buff_b = NULL;
  if(t3 > t2)
    {
      buff_n = v2;
      v2     = v3;
      v3     = buff_n;

      buff_b = b2;
      b2     = b3;
      b3     = buff_b;
    }  
    
  t0 = tree->times->nd_t[a->num];
  t1 = tree->times->nd_t[d->num];
  t2 = tree->times->nd_t[v2->num];
  t3 = tree->times->nd_t[v3->num];
  u0 = tree->rates->br_r[a->num];
  r1 = tree->rates->br_r[d->num];
  r2 = tree->rates->br_r[v2->num];
  r3 = tree->rates->br_r[v3->num];
  l1 = tree->rates->cur_l[d->num];
  l2 = tree->rates->cur_l[v2->num];
  l3 = tree->rates->cur_l[v3->num];
  cr = tree->rates->clock_r;
  nf = tree->rates->norm_fact;

  for(i=0;i<dim;i++) is_1[i] = 0;
  is_1[b1->num] = 1;
  is_1[b2->num] = 1;
  is_1[b3->num] = 1;

/*   for(i=0;i<dim;i++)     cond_mu[i]  = 0.0; */
/*   For(i,dim*dim) cond_cov[i] = 0.0; */
/*   Normal_Conditional(tree->rates->mean_l,tree->rates->cov_l,tree->rates->u_cur_l,dim,is_1,3,cond_mu,cond_cov); */
  
/*   El1 = cond_mu[b1->num]; */
/*   El2 = cond_mu[b2->num]; */
/*   El3 = cond_mu[b3->num]; */

/*   cov11 = cond_cov[b1->num*dim+b1->num]; */
/*   cov12 = cond_cov[b1->num*dim+b2->num]; */
/*   cov13 = cond_cov[b1->num*dim+b3->num]; */
/*   cov23 = cond_cov[b2->num*dim+b3->num]; */
/*   cov22 = cond_cov[b2->num*dim+b2->num]; */
/*   cov33 = cond_cov[b3->num*dim+b3->num]; */


/*   El1 = tree->rates->u_ml_l[b1->num]; */
/*   El2 = tree->rates->u_ml_l[b2->num]; */
/*   El3 = tree->rates->u_ml_l[b3->num]; */

/*   cov11 = tree->rates->cov[b1->num*dim+b1->num]; */
/*   cov12 = tree->rates->cov[b1->num*dim+b2->num]; */
/*   cov13 = tree->rates->cov[b1->num*dim+b3->num]; */
/*   cov23 = tree->rates->cov[b2->num*dim+b3->num]; */
/*   cov22 = tree->rates->cov[b2->num*dim+b2->num]; */
/*   cov33 = tree->rates->cov[b3->num*dim+b3->num]; */

/*   PhyML_Fprintf(stderr,"\n- El1=%10f El2=%10f El3=%10f",El1,El2,El3); */
/*   PhyML_Fprintf(stderr,"\n- cov11=%10f cov12=%10f cov13=%10f cov23=%10f cov22=%10f cov33=%10f", cov11,cov12,cov13,cov23,cov22,cov33); */

  num_1 = num_2 = num_3 = -1;
  if(b1->num < b2->num && b2->num < b3->num) { num_1 = 0; num_2 = 1; num_3 = 2; }
  if(b1->num < b3->num && b3->num < b2->num) { num_1 = 0; num_2 = 2; num_3 = 1; }
  if(b2->num < b1->num && b1->num < b3->num) { num_1 = 1; num_2 = 0; num_3 = 2; }
  if(b2->num < b3->num && b3->num < b1->num) { num_1 = 2; num_2 = 0; num_3 = 1; }
  if(b3->num < b2->num && b2->num < b1->num) { num_1 = 2; num_2 = 1; num_3 = 0; }
  if(b3->num < b1->num && b1->num < b2->num) { num_1 = 1; num_2 = 2; num_3 = 0; }

  cov11 = tree->rates->trip_cond_cov[d->num * 9 + num_1 * 3 + num_1];
  cov12 = tree->rates->trip_cond_cov[d->num * 9 + num_1 * 3 + num_2];
  cov13 = tree->rates->trip_cond_cov[d->num * 9 + num_1 * 3 + num_3];
  cov23 = tree->rates->trip_cond_cov[d->num * 9 + num_2 * 3 + num_3];
  cov22 = tree->rates->trip_cond_cov[d->num * 9 + num_2 * 3 + num_2];
  cov33 = tree->rates->trip_cond_cov[d->num * 9 + num_3 * 3 + num_3];

  El1=0.0;
  for(i=0;i<dim;i++)
    if(i != b1->num && i != b2->num && i != b3->num)
      El1 += tree->rates->trip_reg_coeff[d->num * (6*tree->n_otu-9) + num_1 * dim +i] * (tree->rates->u_cur_l[i] - tree->rates->mean_l[i]);
  El1 += tree->rates->mean_l[b1->num];

  El2=0.0;
  for(i=0;i<dim;i++)
    if(i != b1->num && i != b2->num && i != b3->num)
      El2 += tree->rates->trip_reg_coeff[d->num * (6*tree->n_otu-9) + num_2 * dim +i] * (tree->rates->u_cur_l[i] - tree->rates->mean_l[i]);
  El2 += tree->rates->mean_l[b2->num];

  El3=0.0;
  for(i=0;i<dim;i++)
    if(i != b1->num && i != b2->num && i != b3->num)
      El3 += tree->rates->trip_reg_coeff[d->num * (6*tree->n_otu-9) + num_3 * dim +i] * (tree->rates->u_cur_l[i] - tree->rates->mean_l[i]);
  El3 += tree->rates->mean_l[b3->num];


/*   PhyML_Fprintf(stderr,"\n+ El1=%10f El2=%10f El3=%10f",El1,El2,El3); */
/*   PhyML_Fprintf(stderr,"\n+ cov11=%10f cov12=%10f cov13=%10f cov23=%10f cov22=%10f cov33=%10f", cov11,cov12,cov13,cov23,cov22,cov33); */


  t1_new = +1;

  t_min = MAX(t0,tree->times->t_prior_min[d->num]);
  t_max = MIN(MIN(t2,t3),tree->times->t_prior_max[d->num]);
    
  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;

  if(t_min > t_max) 
    {
      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  bl_min = bl_max = -1.0;

  l_opp = 0.0;
  if(a == tree->n_root)
    {
      l_opp = (d == a->v[0])?(tree->rates->cur_l[a->v[1]->num]):(tree->rates->cur_l[a->v[0]->num]);
      El1 -= l_opp;
    }
     
  EX = El1/r1 + El2/r2;
  EY = El1/r1 + El3/r3;
  
  VX = cov11/(r1*r1) + cov22/(r2*r2) + 2.*cov12/(r1*r2);
  VY = cov11/(r1*r1) + cov33/(r3*r3) + 2.*cov13/(r1*r3);

  mu[0] = El1;
  mu[1] = EX;
  mu[2] = EY;
  
  sig11 = cov11;
  sig1X = cov11/r1 + cov12/r2;
  sig1Y = cov11/r1 + cov13/r3;
  sigXX = VX;
  sigYY = VY;
  sigXY = cov11/(r1*r1) + cov13/(r1*r3) + cov12/(r1*r2) + cov23/(r2*r3);

  cov[0*3+0] = sig11; cov[0*3+1] = sig1X; cov[0*3+2] = sig1Y;
  cov[1*3+0] = sig1X; cov[1*3+1] = sigXX; cov[1*3+2] = sigXY;
  cov[2*3+0] = sig1Y; cov[2*3+1] = sigXY; cov[2*3+2] = sigYY;

  l2XY[0] = 0.0; /* does not matter */
  l2XY[1] = (t2-t0)*cr*nf; /* constraint 1  */
  l2XY[2] = (t3-t0)*cr*nf; /* constraint 2  */

  is_1[0] = is_1[1] = is_1[2] = 0;
  is_1[0] = 1;

  Normal_Conditional(mu,cov,l2XY,3,is_1,1,cond_mu,cond_cov);


  if(cond_cov[0*3+0] < 0.0)
    {
      PhyML_Printf("\n. a: %d d: %d",a->num,d->num);
      PhyML_Printf("\n. Conditional mean=%G var=%G",cond_mu[0],cond_cov[0*3+0]);
      PhyML_Printf("\n. t0=%G t1=%f t2=%f t3=%f l1=%G l2=%G l3=%G",t0,t1,t2,t3,l1,l2,l3);
      PhyML_Printf("\n. El1=%G El2=%G El3=%G Nu=%G r1=%G r2=%G r3=%G cr=%G",El1,El2,El3,tree->rates->nu,r1,r2,r3,cr);
      PhyML_Printf("\n. COV11=%f COV12=%f COV13=%f COV22=%f COV23=%f COV33=%f",cov11,cov12,cov13,cov22,cov23,cov33);
      PhyML_Printf("\n. constraint1: %f constraints2: %f",l2XY[1],l2XY[2]);
      PhyML_Printf("\n");
      for(i=0;i<3;i++)
      	{
      	  PhyML_Printf(". mu%d=%12lf\t",i,mu[i]);
      	  for(j=0;j<3;j++)
      	    {
      	      PhyML_Printf("%12lf ",cov[i*3+j]);
      	    }
      	  PhyML_Printf("\n");
      	}      
      cond_cov[0*3+0] = 1.E-10;
    }


  bl_min = (t_min - t0) * r1 * cr * nf;
  bl_max = (t_max - t0) * r1 * cr * nf;

  new_l1 = Rnorm_Trunc(cond_mu[0],inflate_var*SQRT(cond_cov[0*3+0]),bl_min,bl_max,&err);
/*   new_l1 = Rnorm(cond_mu[0],SQRT(cond_cov[0*3+0])); */


  if(new_l1 < bl_min) new_l1 = l1;
  if(new_l1 > bl_max) new_l1 = l1;
      
  t1_new = new_l1/(r1*cr*nf) + t0;


  if(err)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. Root ? %s",(tree->n_root==a)?("yes"):("no")); 
      PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run);
      PhyML_Printf("\n. t0=%f t1=%f t2=%f t3=%f",t0,t1,t2,t3);
      PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
      PhyML_Printf("\n. bl_min=%f bl_max=%f",bl_min,bl_max);
      PhyML_Printf("\n. cond_mu[0]=%f cond_cov[0]=%f",cond_mu[0],SQRT(cond_cov[0*3+0]));
      PhyML_Printf("\n. El1=%f El2=%f El3=%f",El1,El2,El3);
      PhyML_Printf("\n. l1=%f l2=%f l3=%f",l1,l2,l3);
      PhyML_Printf("\n. u0=%f r1=%f r2=%f r3=%f",u0,r1,r2,r3);
      PhyML_Printf("\n. COV11=%f COV22=%f",cov11,cov22);
      PhyML_Printf("\n. Clock rate = %f",tree->rates->clock_r);
      PhyML_Printf("\n. Setting t1_new to %f",t1);
      t1_new = t1;
      /* 	  Exit("\n"); */
    }
  /*     } */
  /*   else */
  /*     { */
  /*       bl_min = (t2 - t_max) * r2; */
  /*       bl_max = (t2 - t_min) * r2; */
  
  /*       l2 = Rnorm_Trunc(cond_mu[0],SQRT(cond_cov[0*3+0]),bl_min,bl_max,&err); */
  /*       t1_new = -l2/r2 + t2; */
  
  /*       if(err) */
  /* 	{ */
  /* 	  PhyML_Printf("\n"); */
  /* 	  PhyML_Printf("\n. %s %d %d",__FILE__,__LINE__,tree->mcmc->run); */
  /* 	  PhyML_Printf("\n. t0=%f t1=%f t2=%f t3=%f",t0,t1,t2,t3); */
  /* 	  PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max); */
  /* 	  PhyML_Printf("\n. bl_min=%f bl_max=%f",bl_min,bl_max); */
  /* 	  PhyML_Printf("\n. cond_mu[0]=%f cond_cov[0]=%f",cond_mu[0],SQRT(cond_cov[0*3+0])); */
  /* 	  PhyML_Printf("\n. El1=%f El2=%f El3=%f",El1,El2,El3); */
  /* 	  PhyML_Printf("\n. l1=%f l2=%f l3=%f",l1,l2,l3); */
  /* 	  PhyML_Printf("\n. u0=%f r1=%f r2=%f r3=%f",u0,r1,r2,r3); */
  /* 	  PhyML_Printf("\n. COV11=%f COV22=%f",cov11,cov22); */
  /* 	  PhyML_Printf("\n. Clock rate = %f",tree->rates->clock_r); */
  /* 	  PhyML_Printf("\n. Setting t1_new to %f",t1); */
  /* 	  t1_new = t1; */
  /* /\* 	  Exit("\n"); *\/ */
  /* 	} */
  /*     } */
  
  if(t1_new < t0)
    {
      t1_new = t0+1.E-4;
      PhyML_Printf("\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. t0 = %f t1_new = %f t1 = %f",t0,t1_new,t1);
      PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
      PhyML_Printf("\n. l1 = %f",l1);
      PhyML_Printf("\n. bl_min = %f bl_max = %f",bl_min,bl_max);
      PhyML_Printf("\n. (t1-t0)=%f (t2-t1)=%f",t1-t0,t2-t1);
      PhyML_Printf("\n. l1 = %f l2 = %f cov11=%f cov22=%f cov33=%f",l1,l2,cov11,cov22,cov33);
      PhyML_Printf("\n. clock=%G",tree->rates->clock_r);
      PhyML_Printf("\n. u0=%f r1=%f r2=%f r3=%f",u0,r1,r2,r3);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      /*       Exit("\n"); */
    }
  if(t1_new > MIN(t2,t3))
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
      PhyML_Printf("\n. t0 = %f t1_new = %f t1 = %f t2 = %f t3 = %f MIN(t2,t3)=%f",t0,t1_new,t1,t2,t3,MIN(t2,t3));
      PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
      PhyML_Printf("\n. l2 = %f",l2);
      PhyML_Printf("\n. bl_min = %f bl_max = %f",bl_min,bl_max);
      PhyML_Printf("\n. (t1-t0)=%f (t2-t1)=%f",t1-t0,t2-t1);
      PhyML_Printf("\n. l1 = %f l2 = %f cov11=%f cov22=%f cov33=%f",l1,l2,cov11,cov22,cov33);
      PhyML_Printf("\n. clock=%G",tree->rates->clock_r);
      PhyML_Printf("\n. u0=%f r1=%f r2=%f r3=%f",u0,r1,r2,r3);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      /*       Exit("\n"); */
    }
  
  if(isnan(t1_new))
    {
      PhyML_Printf("\n. run=%d",tree->mcmc->run);
      PhyML_Printf("\n. mean=%G var=%G",cond_mu[0],cond_cov[0*3+0]);
      PhyML_Printf("\n. t1=%f l1=%G r1=%G t0=%G",t1,l1,r1,t0);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      /*       Exit("\n"); */
    }
  
  if(tree->eval_alnL == YES)
    tree->times->nd_t[d->num] = t1_new;
  else
    {
      u = Uni();
      tree->times->nd_t[d->num] = u*(t_max-t_min) + t_min;
    }

  RATES_Normalise_Rates(tree);
  RATES_Update_Edge_Lengths(tree);
  
  cur_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL;
  new_lnL_data = tree->c_lnL;
  new_lnL_rate = tree->rates->c_lnL;
  
  /* !!!!!!!!!!!!!!!!! */
  /* tree->both_sides = NO; */
  /* new_lnL_data = Lk(tree); */
  
  if(tree->io->lk_approx == EXACT)
    {
      Update_PMat_At_Given_Edge(b1,tree);
      Update_PMat_At_Given_Edge(b2,tree);
      Update_PMat_At_Given_Edge(b3,tree);
      Update_Partial_Lk(tree,b1,d);
    }
  new_lnL_data = Lk(b1,tree);

  new_lnL_rate = RATES_Lk(tree);

  ratio = 0.0;

  /* Proposal ratio */
  if(tree->eval_alnL)
    ratio += (Log_Dnorm_Trunc(l1,    cond_mu[0],inflate_var*SQRT(cond_cov[0*3+0]),bl_min,bl_max,&err) - 
	      Log_Dnorm_Trunc(new_l1,cond_mu[0],inflate_var*SQRT(cond_cov[0*3+0]),bl_min,bl_max,&err));
    
  /* Prior ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);

  /* Likelihood ratio */
  ratio += (new_lnL_data - cur_lnL_data);

  /*   printf("\n* d:%d Ratio=%f l1=%f new_l1=%f mean=%f ml=%f sd=%f [%f %f]", */
  /* 	 d->num, */
  /* 	 exp(ratio), */
  /* 	 l1,new_l1, */
  /* 	 cond_mu[0], */
  /* 	 tree->rates->mean_l[b1->num], */
  /* 	 SQRT(cond_cov[0*3+0]), */
  /* 	 Log_Dnorm(l1,cond_mu[0],SQRT(cond_cov[0*3+0]),&err),Log_Dnorm(new_l1,cond_mu[0],SQRT(cond_cov[0*3+0]),&err)); */
  
  ratio = exp(ratio);


  u = Uni();
  if(u > MIN(1.,ratio))
    {
      tree->times->nd_t[d->num] = t1; /* reject */
      tree->rates->c_lnL        = cur_lnL_rate;
      tree->c_lnL               = cur_lnL_data;
      RATES_Update_Edge_Lengths(tree);
      if(tree->io->lk_approx == EXACT) 
	{
	  Update_PMat_At_Given_Edge(b1,tree);
	  Update_PMat_At_Given_Edge(b2,tree);
	  Update_PMat_At_Given_Edge(b3,tree);
	  Update_Partial_Lk(tree,b1,d);
	}
    }
  else
    {
      tree->mcmc->acc_move[num_move]++;
    }
  
  RATES_Normalise_Rates(tree);

  tree->mcmc->run_move[num_move]++;
  
  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  for(i=0;i<3;i++)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
	      {
		if(tree->io->lk_approx == EXACT) Update_Partial_Lk(tree,d->b[i],d);
		RATES_Posterior_One_Time(d,d->v[i],YES,tree);
	      }
	}
      if(tree->io->lk_approx == EXACT) Update_Partial_Lk(tree,b1,d);
    }
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Posterior_Time_Root(t_tree *tree)
{
  /* 
        Root, T0
        / \
  l1   /   \ l2
      /     \
     v1,t1   v2,t2

  */

  phydbl t1,t2;
  phydbl u0,u1,u2;
  phydbl cel,cvl;
  int i,dim;
  t_edge *b;
  t_node *root;
  phydbl t0,t0_min, t0_max;
  phydbl new_t0;
/*   phydbl t_max_01, t_max_02; */
  int err;
  phydbl cr;
  phydbl bl_min, bl_max;
  phydbl new_l;
  phydbl new_lnL_data, cur_lnL_data, cur_lnL_rate;
  phydbl u,ratio;

  dim    = 2*tree->n_otu-3;
  b      = tree->e_root;
  root   = tree->n_root;
  t0     = tree->times->nd_t[root->num];
  t1     = tree->times->nd_t[root->v[2]->num];
  t2     = tree->times->nd_t[root->v[1]->num];
  u1     = tree->rates->br_r[root->v[2]->num];
  u2     = tree->rates->br_r[root->v[1]->num];
  u0     = 1.0;
  cr     = tree->rates->clock_r;
  t0_min = -BIG;
  t0_max = MIN(t1,t2);
  
  /*   t0_max = MIN(t1 - (1./tree->rates->nu)*POW((u0-u1)/tree->rates->z_max,2), */
  /* 	       t2 - (1./tree->rates->nu)*POW((u0-u2)/tree->rates->z_max,2)); */
  
  /*   if(u1 > u0) t_max_01 = t1 - (1./tree->rates->nu)*POW((u1-u0)/(PointNormal(tree->rates->p_max*Pnorm(.0,.0,1.))),2); */
  /*   else        t_max_01 = t1 - (1./tree->rates->nu)*POW((u1-u0)/(PointNormal(tree->rates->p_max*(1.-Pnorm(.0,.0,1.)))),2); */
  
  /*   if(u2 > u0) t_max_02 = t2 - (1./tree->rates->nu)*POW((u2-u0)/(PointNormal(tree->rates->p_max*Pnorm(.0,.0,1.))),2); */
  /*   else        t_max_02 = t2 - (1./tree->rates->nu)*POW((u2-u0)/(PointNormal(tree->rates->p_max*(1.-Pnorm(.0,.0,1.)))),2); */
  
  
  /*   t_max_01 = RATES_Find_Max_Dt_Bisec(u1,u0,tree->times->t_prior_min[root->num],t1,tree->rates->nu,tree->rates->p_max,(u1 < u0)?YES:NO); */
  /*   t_max_02 = RATES_Find_Max_Dt_Bisec(u2,u0,tree->times->t_prior_min[root->num],t2,tree->rates->nu,tree->rates->p_max,(u2 < u0)?YES:NO); */
  /*   t0_max = MIN(t_max_01,t_max_02); */
    
  /*   RATES_Min_Max_Interval(u0,u1,u2,r3,t0,t2,t3,&t_min,&t_max,tree->rates->nu,tree->rates->p_max,tree); */

  t0_min = MAX(t0_min,tree->times->t_prior_min[root->num]);
  t0_max = MIN(t0_max,tree->times->t_prior_max[root->num]);

  t0_max -= tree->rates->min_dt;

  if(t0_min > t0_max) return;

  tree->times->t_prior[root->num] = Uni()*(t0_max - t0_min) + t0_min;

  u0 *= cr;
  u1 *= cr;
  u2 *= cr;

  if(FABS(tree->times->t_prior_min[root->num] - tree->times->t_prior_max[root->num]) < 1.E-10) return;

  cel=0.0;
  for(i=0;i<dim;i++) if(i != b->num) cel += tree->rates->reg_coeff[b->num*dim+i] * (tree->rates->u_cur_l[i] - tree->rates->mean_l[i]);
  cel += tree->rates->mean_l[b->num];

  cvl = tree->rates->cond_var[b->num];
  
  bl_min = u1 * (t1 - t0_max) + u2 * (t2 - t0_max); 
  bl_max = u1 * (t1 - t0_min) + u2 * (t2 - t0_min);

  if(bl_min > bl_max) return;

  new_l = Rnorm_Trunc(cel,SQRT(cvl),bl_min,bl_max,&err);

  new_t0 = (u1*t1 + u2*t2 - new_l)/(u1+u2);

  if(t0 > t1 || t0 > t2)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. Run = %d",tree->mcmc->run);
      PhyML_Printf("\n. t0=%f t1=%f t2=%f",t0,t1,t2);
      PhyML_Printf("\n. t0_min=%f t0_max=%f",t0_min,t0_max);
      PhyML_Printf("\n. new_l=%f cel=%f",new_l,cel);
      PhyML_Printf("\n. u0=%f u1=%f u2=%f",u0/cr,u1/cr,u2/cr);
      PhyML_Printf("\n. Nu = %f Clock = %f",tree->rates->nu,tree->rates->clock_r);
      PhyML_Printf("\n. Setting t0 to %f",tree->times->nd_t[root->num]);
      return;
    }

  if(t0 < t0_min || t0 > t0_max)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. Run = %d",tree->mcmc->run);
      PhyML_Printf("\n. t0=%f t1=%f t2=%f",t0,t1,t2);
      PhyML_Printf("\n. t0_min=%f t0_max=%f",t0_min,t0_max);
      PhyML_Printf("\n. u0=%f u1=%f u2=%f",u0/cr,u1/cr,u2/cr);
      PhyML_Printf("\n. Nu = %f Clock = %f",tree->rates->nu,tree->rates->clock_r);
      PhyML_Printf("\n. Setting t0 to %f",tree->times->nd_t[root->num]);
      t0 = tree->times->nd_t[root->num];      
    }

  /* Sample according to prior */
/*   tree->times->nd_t[root->num] = tree->times->t_prior[root->num]; */

  /* Sample according to posterior */
  tree->times->nd_t[root->num] = new_t0;

  RATES_Normalise_Rates(tree);
  RATES_Update_Edge_Lengths(tree);

  cur_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL;
  new_lnL_data = tree->c_lnL;
  
  new_lnL_data = Lk(NULL,tree);
    
  ratio = 0.0;
  /* Prior ratio */
  ratio += .0;
  /* Likelihood ratio */
  ratio += new_lnL_data - cur_lnL_data;

  ratio = exp(ratio);
  u = Uni();
  if(u > MIN(1.,ratio))
    {
      tree->times->nd_t[root->num] = t0; /* reject */
      tree->rates->c_lnL     = cur_lnL_rate;
      tree->c_lnL                  = cur_lnL_data;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_times]++;
    }

  RATES_Normalise_Rates(tree);
  RATES_Update_Edge_Lengths(tree);

  tree->mcmc->run_move[tree->mcmc->num_move_times]++;
  tree->mcmc->acc_rate[tree->mcmc->num_move_times] = 
    (tree->mcmc->acc_move[tree->mcmc->num_move_times]+1.E-6)/ 
    (tree->mcmc->run_move[tree->mcmc->num_move_times]+1.E-6);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_Edge_Lengths(t_tree *tree)
{  
  RATES_Update_Edge_Lengths_Pre(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);
  RATES_Update_Edge_Lengths_Pre(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
  
  RATES_Update_One_Edge_Length(tree->e_root,tree);
  
  if(tree->is_mixt_tree == YES) MIXT_RATES_Update_Edge_Lengths(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_Edge_Lengths_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  RATES_Update_One_Edge_Length(b,tree);
  
  if(d->tax) return;
  else
    {
      int i;
      for(i=0;i<3;++i) 
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Update_Edge_Lengths_Pre(d,d->v[i],d->b[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_One_Edge_Length(t_edge *b, t_tree *tree)
{
  phydbl dt,rr,ra,rd,ta,td,nu,cr;
  t_node *a, *d;

  if(b == tree->e_root)
    {
      RATES_Update_One_Edge_Length(tree->n_root->b[1],tree);
      RATES_Update_One_Edge_Length(tree->n_root->b[2],tree);

      if(tree->mod && tree->mod->log_l == YES)
        {
          tree->e_root->l->v = 
            exp(tree->n_root->b[1]->l->v) +
            exp(tree->n_root->b[2]->l->v) ;
          tree->e_root->l->v = log(tree->e_root->l->v);
        }
      else
        {
          tree->e_root->l->v = tree->n_root->b[1]->l->v + tree->n_root->b[2]->l->v;
        }
      
      tree->rates->u_cur_l[tree->e_root->num] = tree->e_root->l->v;
      tree->n_root_pos = tree->n_root->b[2]->l->v / tree->e_root->l->v;
  
      if(tree->rates->model_id == GUINDON)
        {
          phydbl t0,t1,t2;
          t_node *n0, *n1;
          
          n0 = tree->n_root->v[2];
          n1 = tree->n_root->v[1];
          t1 = tree->times->nd_t[tree->n_root->v[2]->num];
          t2 = tree->times->nd_t[tree->n_root->v[1]->num];
          t0 = tree->times->nd_t[tree->n_root->num];
          
          tree->e_root->l->v = 
            (t1-t0)/(t1+t2-2.*t0)*tree->rates->cur_gamma_prior_mean[n0->num] +
            (t2-t0)/(t1+t2-2.*t0)*tree->rates->cur_gamma_prior_mean[n1->num];
          
          tree->e_root->l_var->v = 
            pow((t1-t0)/(t1+t2-2.*t0),2)*tree->rates->cur_gamma_prior_var[n0->num] +
            pow((t2-t0)/(t1+t2-2.*t0),2)*tree->rates->cur_gamma_prior_var[n1->num];
        }      
    }
  else
    {
      if(b->left->anc == b->rght)
        {
          d = b->left;
          a = b->rght;
        }
      else
        {
          d = b->rght;
          a = b->left;
        }
      
      
      assert(a);
      assert(d);
      assert(d->anc == a);
      
      ra = rd = -1.;
      
      if(tree->rates->model_id == LOGNORMAL ||
         tree->rates->model_id == THORNE ||
         tree->rates->model_id == STRICTCLOCK)
        {
          rd = tree->rates->br_r[d->num];
          ra = tree->rates->br_r[a->num];
        }
      else if(tree->rates->model_id == GUINDON)
        {
          rd = tree->rates->nd_r[d->num];
          ra = tree->rates->nd_r[a->num];
        }
      else assert(FALSE);
      
      dt = fabs(tree->times->nd_t[d->num] - tree->times->nd_t[a->num]);
      cr = tree->rates->clock_r;
      td = tree->times->nd_t[d->num];
      ta = tree->times->nd_t[a->num];
      nu = tree->rates->nu;
      rr = -1.0;
      
      if(tree->rates->model_id == LOGNORMAL)
        {
          tree->rates->cur_l[d->num] = dt*rd*cr;
        }
      
      if(tree->rates->model_id == THORNE)
        {
          rr = (ra+rd)/2.;          
          tree->rates->cur_l[d->num] = dt*rr*cr;
        }
      
      if(tree->rates->model_id == GUINDON)
        {
          phydbl m,v;
          
          Integrated_Geometric_Brownian_Bridge_Moments(dt,ra,rd,nu,&m,&v);
          
          if(isnan(m) || isnan(v) || m < 0.0 || v < 0.0)
            {
              PhyML_Fprintf(stderr,"\n. dt: %G ra: %G rd: %G nu: %G m: %G v: %G a is root ? %d d is root ? %d",
                            dt,ra,rd,nu,m,v,
                            a==tree->n_root,
                            d==tree->n_root);
            }
          
          m *= dt*cr; // the actual rate average is m * dt. We multiply by dt in order to derive the value for the branch length
          v *= (dt*dt)*(cr*cr);
          
          tree->rates->cur_gamma_prior_mean[d->num] = m;
          tree->rates->cur_gamma_prior_var[d->num]  = v;
          
          tree->rates->cur_l[d->num] = tree->rates->cur_gamma_prior_mean[d->num]; // Required for having proper branch lengths in Write_Tree function
        }
      
      if(tree->rates->model_id == STRICTCLOCK)
        {
          tree->rates->cur_l[d->num] = dt*cr;          
        }
      
      /* printf("\n. td: %12f ta: %12f dt: %12f cr: %12f ra: %12f rd: %12f l: %12f", */
      /*        tree->times->nd_t[d->num], */
      /*        tree->times->nd_t[a->num], */
      /*        dt,cr,ra,rd,tree->rates->cur_l[d->num]); */
      
      
      if(tree->mod && tree->mod->log_l == YES) tree->rates->cur_l[d->num] = log(tree->rates->cur_l[d->num]);
      
      if(b)
        {
          b->l->v                      = tree->rates->cur_l[d->num];
          tree->rates->u_cur_l[b->num] = tree->rates->cur_l[d->num];
          b->l_var->v                  = tree->rates->cur_gamma_prior_var[d->num];
        }
      
      if(b && (isnan(b->l->v) || isnan(b->l_var->v)))
        {
          PhyML_Fprintf(stderr,"\n. dt=%G rr=%G cr=%G ra=%G rd=%G nu=%G %f %f ",dt,rr,tree->rates->clock_r,ra,rd,nu,b->l_var->v,b->l->v);	  
          PhyML_Fprintf(stderr,"\n. ta=%G td=%G ra*cr=%G rd*cr=%G sd=%G",ta,td,ra,rd,SQRT(dt*nu));
          /* assert(FALSE); */
        }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Bl_To_Bl(t_tree *tree)
{
  RATES_Bl_To_Bl_Pre(tree->n_root,tree->n_root->v[2],NULL,tree);
  RATES_Bl_To_Bl_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
  /* tree->rates->cur_l[tree->n_root->v[2]->num] = tree->a_edges[tree->e_root->num]->l->v * tree->n_root_pos; */
  /* tree->rates->cur_l[tree->n_root->v[1]->num] = tree->a_edges[tree->e_root->num]->l->v * (1. - tree->n_root_pos); */
  tree->rates->cur_l[tree->n_root->v[2]->num] = tree->a_edges[tree->e_root->num]->l->v * 0.5;
  tree->rates->cur_l[tree->n_root->v[1]->num] = tree->a_edges[tree->e_root->num]->l->v * (1. - 0.5);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Bl_To_Bl_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{

  if(b) 
    {
      tree->rates->cur_l[d->num] = b->l->v;
    }

  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++) 
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Bl_To_Bl_Pre(d,d->v[i],d->b[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Bl_To_Ml(t_tree *tree)
{
  RATES_Bl_To_Ml_Pre(tree->n_root,tree->n_root->v[2],NULL,tree);
  RATES_Bl_To_Ml_Pre(tree->n_root,tree->n_root->v[1],NULL,tree);
  tree->rates->u_ml_l[tree->e_root->num] = tree->a_edges[tree->e_root->num]->l->v;
  tree->rates->ml_l[tree->n_root->v[2]->num] = tree->rates->u_ml_l[tree->e_root->num] * tree->n_root_pos;
  tree->rates->ml_l[tree->n_root->v[1]->num] = tree->rates->u_ml_l[tree->e_root->num] * (1. - tree->n_root_pos);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Bl_To_Ml_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{

  if(b) 
    {
      tree->rates->u_ml_l[b->num] = b->l->v;
      tree->rates->ml_l[d->num]   = b->l->v;
    }

  if(d->tax) return;
  else
    {
      int i;

      for(i=0;i<3;i++) 
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Bl_To_Ml_Pre(d,d->v[i],d->b[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Get_Cov_Matrix_Rooted(phydbl *unroot_cov, t_tree *tree)
{
  int i,dim;

  dim = 2*tree->n_otu-3;

  RATES_Get_Cov_Matrix_Rooted_Pre(tree->n_root,tree->n_root->v[2],NULL,unroot_cov,tree);
  RATES_Get_Cov_Matrix_Rooted_Pre(tree->n_root,tree->n_root->v[1],NULL,unroot_cov,tree);

  for(i=0;i<dim+1;++i) tree->rates->cov_l[i*(dim+1)+tree->n_root->v[2]->num] /= 2.;
  for(i=0;i<dim+1;++i) tree->rates->cov_l[i*(dim+1)+tree->n_root->v[1]->num] /= 2.;
  for(i=0;i<dim+1;++i) tree->rates->cov_l[tree->n_root->v[2]->num*(dim+1)+i] /= 2.;
  for(i=0;i<dim+1;++i) tree->rates->cov_l[tree->n_root->v[1]->num*(dim+1)+i] /= 2.;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Get_Cov_Matrix_Rooted_Pre(t_node *a, t_node *d, t_edge *b, phydbl *cov, t_tree *tree)
{
  int i, dim;
  t_node *n;

  dim = 2*tree->n_otu-3;
  n   = NULL;

  for(i=0;i<dim;i++) 
    { 
      if(tree->a_edges[i] != tree->e_root)
	{
	  n = 
	    (tree->a_edges[i]->left->anc == tree->a_edges[i]->rght)?
	    (tree->a_edges[i]->left):
	    (tree->a_edges[i]->rght);

	  if(b)
	    {
	      tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[b->num*dim + i];
	    }
	  else
	    {
	      tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[tree->e_root->num*dim + i];
	    }
	}
      else
	{
	  n = tree->e_root->left;
	  if(b)
	    tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[b->num*dim + i];
	  else
	    tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[tree->e_root->num*dim + i];

	  n = tree->e_root->rght;
	  if(b)
	    tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[b->num*dim + i];
	  else
	    tree->rates->cov_l[d->num*(dim+1) + n->num] = cov[tree->e_root->num*dim + i];
	}
    }


  if(d->tax) return;
  else
    {
      for(i=0;i<3;i++)
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  RATES_Get_Cov_Matrix_Rooted_Pre(d,d->v[i],d->b[i],cov,tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Covariance_Mu(t_tree *tree)
{
  int i,j;
  phydbl dt,var;
  int dim;
  int lca_num;

  dim = 2*tree->n_otu-2;

  for(i=0;i<dim*dim;++i) tree->rates->cov_r[i] = 0.0;

  dt =  tree->times->nd_t[tree->n_root->v[2]->num] - tree->times->nd_t[tree->n_root->num];
  var = dt * tree->rates->nu;
  tree->rates->cov_r[tree->n_root->v[2]->num*dim+tree->n_root->v[2]->num] = var;


  dt = tree->times->nd_t[tree->n_root->v[1]->num] - tree->times->nd_t[tree->n_root->num];
  var = dt * tree->rates->nu;
  tree->rates->cov_r[tree->n_root->v[1]->num*dim+tree->n_root->v[1]->num] = var;

  RATES_Variance_Mu_Pre(tree->n_root,tree->n_root->v[2],tree);
  RATES_Variance_Mu_Pre(tree->n_root,tree->n_root->v[1],tree);

  for(i=0;i<dim;i++)
    {
      for(j=i+1;j<dim;j++)
	{
	  lca_num = tree->rates->lca[i*(dim+1)+j]->num;
	  if(lca_num < dim) 
	    {
	      tree->rates->cov_r[i*dim+j] = tree->rates->cov_r[lca_num*dim+lca_num];	    
	      tree->rates->cov_r[j*dim+i] = tree->rates->cov_r[i*dim+j];
	    }
	  else if(lca_num == dim)
	    {
	      tree->rates->cov_r[i*dim+j] = 0.0;	    
	      tree->rates->cov_r[j*dim+i] = 0.0;
	    }
	  else
	    {
	      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Exit("\n");
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Variance_Mu_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int dim;
  phydbl var0;
  phydbl dt1,var1;
  phydbl dt2,var2;
  int i;
  int dir1, dir2;

  dim = 2*tree->n_otu-2;

  if(d->tax) return;

  dir1 = dir2 = -1;
  for(i=0;i<3;i++)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  if(dir1 < 0) dir1 = i;
	  else         dir2 = i;
	}
    }


  var0 = tree->rates->cov_r[d->num*dim+d->num];

  dt1  = tree->times->nd_t[d->v[dir1]->num] - tree->times->nd_t[d->num];
  var1 = tree->rates->nu*dt1;

  dt2  = tree->times->nd_t[d->v[dir2]->num] - tree->times->nd_t[d->num];
  var2 = tree->rates->nu*dt2;

  tree->rates->cov_r[d->v[dir1]->num*dim+d->v[dir1]->num] = var0+var1;
  tree->rates->cov_r[d->v[dir2]->num*dim+d->v[dir2]->num] = var0+var2;


  for(i=0;i<3;i++)
    {
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	{
	  RATES_Variance_Mu_Pre(d,d->v[i],tree);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Fill_Lca_Table(t_tree *tree)
{
  int i,j;
  int dim;
  
  dim = 2*tree->n_otu-1;

  for(i=0;i<dim;i++)
    {
      for(j=i;j<dim;j++)
	{
	  tree->rates->lca[i*dim+j] = Find_Lca_Pair_Of_Nodes(tree->a_nodes[i],tree->a_nodes[j],tree);
	  tree->rates->lca[j*dim+i] = tree->rates->lca[i*dim+j];
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Get V(L_{i} | L_{-i}) for all i */
void RATES_Get_Conditional_Variances(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  int dim;
  t_edge *b;
  phydbl *cond_mu, *cond_cov;

  dim      = 2*tree->n_otu-3;
  a        = tree->rates->_2n_vect1;
  is_1     = tree->rates->_2n_vect5;
  b        = NULL;
  cond_mu  = tree->rates->_2n_vect2;
  cond_cov = tree->rates->_2n2n_vect1;

  for(i=0;i<dim;i++) a[i] = tree->rates->mean_l[i] * (Uni() * 0.2 + 0.9);

  for(i=0;i<dim;i++)
    {
      b = tree->a_edges[i];

      for(j=0;j<dim;j++) is_1[j] = 0;
      is_1[b->num]       = 1;

      For(j,dim*dim) cond_cov[j] = 0.0;
      for(j=0;j<dim;j++)     cond_mu[j]  = 0.0;
      Normal_Conditional(tree->rates->mean_l,tree->rates->cov_l,a,dim,is_1,1,cond_mu,cond_cov);
      
      tree->rates->cond_var[b->num] = cond_cov[b->num*dim+b->num];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Get_All_Reg_Coeff(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  int dim;
  t_edge *b;

  dim  = 2*tree->n_otu-3;
  a    = tree->rates->_2n_vect1;
  is_1 = tree->rates->_2n_vect5;
  b    = NULL;

  for(i=0;i<dim;i++) a[i] = tree->rates->mean_l[i] * (Uni() * 0.2 + 0.9);

  for(i=0;i<dim;i++)
    {
      b = tree->a_edges[i];

      for(j=0;j<dim;j++) is_1[j] = 0;
      is_1[b->num]       = 1;
      
      Get_Reg_Coeff(tree->rates->mean_l,tree->rates->cov_l,a,dim,is_1,1,tree->rates->reg_coeff+b->num*dim);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Get V(L_{i} | L_{-i}) for all i */
void RATES_Get_Trip_Conditional_Variances(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  phydbl *cond_mu, *cond_cov;
  t_node *n;
  int n_otu;

  a        = tree->rates->_2n_vect1;
  is_1     = tree->rates->_2n_vect5;
  cond_mu  = tree->rates->_2n_vect2;
  cond_cov = tree->rates->_2n2n_vect1;
  n        = NULL;
  n_otu    = tree->n_otu;

  For(i,2*n_otu-3) a[i] = tree->rates->mean_l[i] * (Uni() * 0.2 + 0.9);

  for(i=0;i<2*n_otu-2;++i)
    {
      n = tree->a_nodes[i];
      if(!n->tax)
	{
	  For(j,2*n_otu-3) is_1[j] = 0;
	  is_1[n->b[0]->num] = 1;
	  is_1[n->b[1]->num] = 1;
	  is_1[n->b[2]->num] = 1;
	  
	  Normal_Conditional_Unsorted(tree->rates->mean_l,tree->rates->cov_l,a,2*n_otu-3,is_1,3,cond_mu,cond_cov);
	  
	  for(j=0;j<9;j++) tree->rates->trip_cond_cov[n->num*9+j] = cond_cov[j];
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Get_All_Trip_Reg_Coeff(t_tree *tree)
{
  int i,j;
  short int *is_1;
  phydbl *a;
  t_node *n;
  int n_otu;

  a     = tree->rates->_2n_vect1;
  is_1  = tree->rates->_2n_vect5;
  n     = NULL;
  n_otu = tree->n_otu;

  for(i=0;i<2*n_otu-3;++i) a[i] = tree->rates->mean_l[i] * (Uni() * 0.2 + 0.9);

  for(i=0;i<2*n_otu-2;++i)
    {
      n = tree->a_nodes[i];
      if(!n->tax)
	{	  
	  For(j,2*n_otu-3) is_1[j] = 0;
	  is_1[n->b[0]->num] = 1;
	  is_1[n->b[1]->num] = 1;
	  is_1[n->b[2]->num] = 1;
	  
	  Get_Reg_Coeff(tree->rates->mean_l,tree->rates->cov_l,a,2*n_otu-3,is_1,3,tree->rates->trip_reg_coeff+n->num*(6*n_otu-9));
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Check_Lk_Rates(t_tree *tree, int *err)
{
  int i;
  phydbl u,u_anc;
  phydbl t,t_anc;
  
  *err = 0;

  For(i,2*tree->n_otu-2)
    {
      u     = tree->rates->br_r[i];
      u_anc = tree->rates->br_r[tree->a_nodes[i]->anc->num];
      t     = tree->times->nd_t[i];
      t_anc = tree->times->nd_t[tree->a_nodes[i]->anc->num];

      if(t_anc > t)
	{
	  PhyML_Printf("\n. %d %d u=%f u_anc=%f t=%f t_anc=%f",i,tree->a_nodes[i]->anc->num,u,u_anc,t,t_anc);
	  PhyML_Printf("\n. %d %d %d",tree->n_root->num,tree->n_root->v[2]->num,tree->n_root->v[1]->num);
	  *err = 1;
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Realized_Substitution_Rate(t_tree *tree)
{
  return(Tree_Length(tree)/TIMES_Tree_Length(tree));  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Expected_Tree_Length(t_tree *tree)
{
  int n;
  phydbl mean;

  n    = 0;
  mean = 0.0;
  RATES_Expected_Tree_Length_Pre(tree->n_root,tree->n_root->v[2],1.0,&mean,&n,tree);
  RATES_Expected_Tree_Length_Pre(tree->n_root,tree->n_root->v[1],1.0,&mean,&n,tree);

  if(n != 2*tree->n_otu-2)
    {
      PhyML_Fprintf(stderr,"\n. n=%d 2n-2=%d",n,2*tree->n_otu-2);
      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return mean;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Expected_Tree_Length_Pre(t_node *a, t_node *d, phydbl eranc, phydbl *mean, int *n, t_tree *tree)
{
  phydbl erdes;
  int i;
  phydbl loc_mean;
  int loc_n;
  
  
/*  erdes = u_anc - */
/*     sd*(Dnorm((tree->rates->min_rate-u_anc)/sd,.0,1.) - Dnorm((tree->rates->max_rate-u_anc)/sd,.0,1.))/ */
/*        (Pnorm((tree->rates->max_rate-u_anc)/sd,.0,1.) - Pnorm((tree->rates->min_rate-u_anc)/sd,.0,1.)); */

  /* erdes = Norm_Trunc_Mean(eranc,sd,tree->rates->min_rate,tree->rates->max_rate); */
  erdes = 1.0;

  loc_mean = *mean;
  loc_n = *n;

  loc_mean *= loc_n;
  loc_mean += erdes;
  loc_mean /= (loc_n + 1);

  *mean = loc_mean;
  *n    = loc_n + 1;


  if(d->tax) return;
  else
    for(i=0;i<3;i++)
      if(d->v[i] != a && d->b[i] != tree->e_root)
	RATES_Expected_Tree_Length_Pre(d,d->v[i],erdes,mean,n,tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Sample_Average_Rate(t_node *a, t_node *d, t_tree *tree)
{
  return(-1.);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Update_Mean_Br_Len(int iter, t_tree *tree)
{
  int i,dim;
  phydbl *mean;

  if(tree->rates->update_mean_l == NO) return;

  dim = 2*tree->n_otu-3;
  mean = tree->rates->mean_l;

  for(i=0;i<dim;i++)
    {      
     mean[i] *= (phydbl)iter;
     mean[i] += tree->a_edges[i]->l->v;
     mean[i] /= (phydbl)(iter+1);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void RATES_Update_Cov_Br_Len(int iter, t_tree *tree)
{
  int i,j,dim;
  phydbl *mean,*cov;

  if(tree->rates->update_cov_l == NO) return;

  dim = 2*tree->n_otu-3;
  mean = tree->rates->mean_l;
  cov  = tree->rates->cov_l;

  for(i=0;i<dim;i++)
    {     
      for(j=0;j<dim;j++)
	{
	  cov[i*dim+j] += mean[i]*mean[j];
	  cov[i*dim+j] *= (phydbl)tree->mcmc->run;
	  cov[i*dim+j] += tree->a_edges[i]->l->v * tree->a_edges[j]->l->v;
	  cov[i*dim+j] /= (phydbl)(tree->mcmc->run+1);
	  cov[i*dim+j] -= mean[i]*mean[j];

	  if(i == j && cov[i*dim+j] < MIN_VAR_BL) cov[i*dim+j] = MIN_VAR_BL;
	}      
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Set_Mean_L(t_tree *tree)
{
  int i;
  for(i=0;i<2*tree->n_otu-3;++i) 
    {
      tree->rates->mean_l[i] = tree->a_edges[i]->l->v;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Record_Times(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      if(tree->times->nd_t_recorded == YES)
        {
          PhyML_Fprintf(stderr,"\n. Overwriting recorded times is forbidden.\n");
          PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
          Exit("\n");
        }

      For(i,2*tree->n_otu-1) tree->times->buff_t[i] = tree->times->nd_t[i];
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Reset_Times(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  do
    {
      tree->times->nd_t_recorded = NO;
      for(i=0;i<2*tree->n_otu-1;++i) tree->times->nd_t[i] = tree->times->buff_t[i];
      tree = tree->next;
    }
  while(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Record_Rates(t_tree *tree)
{
  int i;

  if(tree->rates->br_r_recorded == YES)
    {
      PhyML_Fprintf(stderr,"\n. Overwriting recorded rates is forbidden.\n");
      PhyML_Fprintf(stderr,"\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  for(i=0;i<2*tree->n_otu-2;++i) tree->rates->buff_br_r[i] = tree->rates->br_r[i];
  for(i=0;i<2*tree->n_otu-1;++i) tree->rates->buff_nd_r[i] = tree->rates->nd_r[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Reset_Rates(t_tree *tree)
{
  int i;
  tree->rates->br_r_recorded = NO;
  for(i=0;i<2*tree->n_otu-2;++i) tree->rates->br_r[i] = tree->rates->buff_br_r[i];
  for(i=0;i<2*tree->n_otu-1;++i) tree->rates->nd_r[i] = tree->rates->buff_nd_r[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Set_Clock_And_Nu_Max(t_tree *tree)
{
  phydbl dt,nu;
  phydbl min_t;
  int i;
  phydbl step;
  phydbl l_max;
  phydbl max_clock;
  phydbl r_max;
  phydbl tune;
  phydbl pa,pb;

  if(tree->rates->model_id == THORNE || 
     tree->rates->model_id == GUINDON)
    {
      tune  = 1.05;
      r_max = tree->rates->max_rate;     
      l_max = tree->mod->l_max;
      
      min_t = .0;
      for(i=0;i<2*tree->n_otu-1;++i) if(tree->times->t_prior_min[i] < min_t) min_t = tree->times->t_prior_min[i];
      
      dt = FABS(min_t);
      max_clock = l_max / dt; 
      
      nu   = 1.E-10;
      step = 1.E-1;
      do
	{
	  do
	    {
	      nu += step;
	      pa = Dnorm(0.0,  0.0,SQRT(nu*dt)); 
	      pb = Dnorm(r_max,0.0,SQRT(nu*dt));
	    }while(pa/pb > tune);
	  nu -= step;
	  step /= 10.;
	}while(step > 1.E-10);
      
      tree->rates->max_nu    = nu;
      /* tree->rates->max_nu    = 1.0; */
      tree->rates->max_clock = max_clock;
      
      PhyML_Printf("\n. Clock rate parameter upper bound set to %f expected subst./site/time unit",tree->rates->max_clock);
      PhyML_Printf("\n. Autocorrelation parameter upper bound set to %f",tree->rates->max_nu);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Set_Birth_Rate_Boundaries(t_tree *tree)
{
  phydbl lbda;
  phydbl p_above_min,p_below_max;
  phydbl min,max;
  int assign = YES;

  min = -0.01*tree->times->t_prior_max[tree->n_root->num];
  max = -100.*tree->times->t_prior_min[tree->n_root->num];

  for(lbda = 0.0001; lbda < 10; lbda+=0.0001)
    {
      p_above_min = 1. - POW(1.-exp(-lbda*min),tree->n_otu);
      p_below_max = POW(1.-exp(-lbda*max),tree->n_otu);
 
      if(p_above_min < 1.E-10) 
	{ 
	  tree->times->birth_rate_max = lbda;
	  break;
	}
      if(p_below_max > 1.E-10 && assign==YES)
	{
	  assign = NO;
	  tree->times->birth_rate_min = lbda;
	}
    }
  
  /* tree->times->birth_rate_min = 1.E-6; */
  /* tree->times->birth_rate_max = 1.; */
  PhyML_Printf("\n. Birth rate lower bound set to %f.",tree->times->birth_rate_min);
  PhyML_Printf("\n. Birth rate upper bound set to %f.",tree->times->birth_rate_max);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl RATES_Get_Mean_Rate_In_Subtree(t_node *root, t_tree *tree)
{
  phydbl sum;
  int n;

  sum = 0.0;
  n   = 0;

  if(root->tax == NO)
    {
      if(root == tree->n_root)
	{
	  RATES_Get_Mean_Rate_In_Subtree_Pre(root,root->v[2],&sum,&n,tree);
	  RATES_Get_Mean_Rate_In_Subtree_Pre(root,root->v[1],&sum,&n,tree);
	}
      else
	{
	  int i;
	  for(i=0;i<3;i++)
	    {
	      if(root->v[i] != root->anc && root->b[i] != tree->e_root)
		{
		  RATES_Get_Mean_Rate_In_Subtree_Pre(root,root->v[i],&sum,&n,tree);
		}
	    }
	}
      return sum/(phydbl)n;
    }
  else
    {
      return 0.0;
    }  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void RATES_Get_Mean_Rate_In_Subtree_Pre(t_node *a, t_node *d, phydbl *sum, int *n, t_tree *tree)
{
  /* (*sum) += exp(tree->rates->nd_r[d->num]); */
  
  if(tree->rates->model_id == LOGNORMAL ||
     tree->rates->model_id == THORNE ||
     tree->rates->model_id == STRICTCLOCK)
    {
      (*sum) += tree->rates->br_r[d->num];
    }
  else if(tree->rates->model_id == GUINDON)
    {
      (*sum) += tree->rates->nd_r[d->num];
    }

  else assert(FALSE);

  (*n) += 1;

  if(d->tax == YES)  return;
  else
    {
      int i;
      for(i=0;i<3;++i)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      RATES_Get_Mean_Rate_In_Subtree_Pre(d,d->v[i],sum,n,tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *RATES_Get_Model_Name(int model)
{
  char *s;

  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  switch(model)
    {
    case GUINDON     : {strcpy(s,"integrated"); break;}
    case THORNE      : {strcpy(s,"autocorrelated"); break;}
    case LOGNORMAL   : {strcpy(s,"uncorrelated"); break;}
    case STRICTCLOCK : {strcpy(s,"strict clock"); break;}
    default : 
      {
	PhyML_Fprintf(stderr,"\n. Err. in file %s at line %d\n",__FILE__,__LINE__);
	Exit("\n");
	break;
     }
    }
  return s;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int RATES_Check_Edge_Length_Consistency(t_tree *mixt_tree)
{
  int i;
  t_tree *tree;

  tree = mixt_tree;
  
  do
    {
      if(tree->is_mixt_tree == YES) tree = tree->next;
      
      for(i=0;i<2*tree->n_otu-3;++i)
        {
          if(tree->a_edges[i]->left->anc == tree->a_edges[i]->rght)
            {
              if(Are_Equal(tree->rates->cur_l[tree->a_edges[i]->left->num],
                           tree->a_edges[i]->l->v,
                           1.E-5) == NO)
                {
                  PhyML_Fprintf(stderr,"\n. cur_l: %G l: %G is_root: %d",
                                tree->rates->cur_l[tree->a_edges[i]->left->num],
                                tree->a_edges[i]->l->v,
                                tree->a_edges[i] == tree->e_root);
                  return 0;
                }
            }
          
          if(tree->a_edges[i]->rght->anc == tree->a_edges[i]->left)
            {
              if(Are_Equal(tree->rates->cur_l[tree->a_edges[i]->rght->num],
                           tree->a_edges[i]->l->v,
                           1.E-5) == NO)
                {
                  PhyML_Fprintf(stderr,"\n. cur_l: %G l: %G is_root: %d",
                                tree->rates->cur_l[tree->a_edges[i]->rght->num],
                                tree->a_edges[i]->l->v,
                                tree->a_edges[i] == tree->e_root);
                  return 0;
                }
            }
        }

      tree = tree->next;
    }
  while(tree);

  return 1;
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



