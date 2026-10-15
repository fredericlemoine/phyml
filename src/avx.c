/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "assert.h"
#include "avx.h"

#if ((defined(__AVX__) || defined(__AVX2__)) && !defined(DISABLE_NATIVE))

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void AVX_Update_Eigen_Lr(t_edge *b, t_tree *tree)
{
  unsigned int i,j;
  
  unsigned const int npattern = tree->data->n_pattern;
  unsigned const int ncatg = tree->mod->ras->n_catg;
  unsigned const int ns = tree->mod->ns;
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned const int ncatgns = ncatg*ns;

  const phydbl *p_lk_left,*p_lk_rght,*pi;
  phydbl *dot_prod;
  phydbl *l_ev,*r_ev;

  __m256d *_l_ev,*_r_ev;

  l_ev         = tree->l_ev;
  _l_ev        = tree->_l_ev;
  _r_ev        = tree->_r_ev;
    
  assert(sz == 4);
  assert(tree->update_eigen_lr == YES);
  
  r_ev = tree->mod->eigen->r_e_vect;
  
  /* Copy transpose of matrix of left eigen vectors */
  for(i=0;i<ns;++i)
    for(j=0;j<ns;++j)
      l_ev[i*ns+j] = tree->mod->eigen->l_e_vect[j*ns+i];
  
  /* Load into AVX registers */
  for(i=0;i<ns;++i)
    {
      for(j=0;j<nblocks;++j)
        {
          _r_ev[i*nblocks+j] = _mm256_load_pd(r_ev + j*sz);
          _l_ev[i*nblocks+j] = _mm256_load_pd(l_ev + j*sz);
        }
      r_ev += ns;
      l_ev += ns;
    }
  
  p_lk_left = b->left->tax ? b->p_lk_tip_l : b->p_lk_left;
  p_lk_rght = b->rght->tax ? b->p_lk_tip_r : b->p_lk_rght;
  pi = tree->mod->e_frq->pi->v;
  dot_prod = tree->dot_prod;

  assert(nblocks <= LK_SIMD_MAX_NS);

  {
    const unsigned int stride_left = (b->left->tax == YES) ? ns : ncatgns;
    const unsigned int stride_rght = (b->rght->tax == YES) ? ns : ncatgns;
    const unsigned int step_left   = (b->left->tax == YES) ? 0 : ns;
    const unsigned int step_rght   = (b->rght->tax == YES) ? 0 : ns;
    const int nt = Lk_Site_Threads(tree,npattern);
    int s;

    /* Each thread uses its own scratch vectors (_prod_left, _prod_rght,
       p_lk_left_pi) and writes only the dot products of its sites. */
#pragma omp parallel for if(nt > 1) num_threads(nt) schedule(static)
    for(s=0;s<(int)npattern;++s)
      {
        __m256d _prod_left[LK_SIMD_MAX_NS],_prod_rght[LK_SIMD_MAX_NS];
        phydbl p_lk_left_pi[LK_SIMD_MAX_NS];
        const phydbl *s_lk_left,*s_lk_rght;
        phydbl *s_dot_prod;
        unsigned int c,k;

        if(tree->data->wght[s] <= SMALL) continue;

        s_lk_left  = p_lk_left + (size_t)s*stride_left;
        s_lk_rght  = p_lk_rght + (size_t)s*stride_rght;
        s_dot_prod = dot_prod + (size_t)s*ncatgns;

        for(c=0;c<ncatg;++c)
          {
            for(k=0;k<ns;++k) p_lk_left_pi[k] = s_lk_left[k] * pi[k];

            AVX_Matrix_Vect_Prod(_r_ev,p_lk_left_pi,ns,_prod_left);
            AVX_Matrix_Vect_Prod(_l_ev,s_lk_rght,ns,_prod_rght);

            for(k=0;k<nblocks;++k) _mm256_store_pd(s_dot_prod + k*sz,_mm256_mul_pd(_prod_left[k],_prod_rght[k]));

            s_dot_prod += ns;
            s_lk_left  += step_left;
            s_lk_rght  += step_rght;
          }
      }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl AVX_Lk_Core_One_Class_No_Eigen_Lr(const phydbl *p_lk_left, const phydbl *p_lk_rght, const phydbl *Pij, const phydbl *tPij, const phydbl *pi, const int ns, const int ambiguity_check, const int observed_state)
{
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;

  if(nblocks == 1)
    {
      __m256d _plk_r,_plk_l;

      if(ambiguity_check == NO) // tip case.
        {
          Pij += observed_state*ns;
          _plk_r = _mm256_mul_pd(_mm256_load_pd(Pij),_mm256_load_pd(p_lk_left));
          return pi[observed_state] * AVX_Vect_Norm(_plk_r);
        }
      else
        {
          unsigned int i;
          __m256d _pijplk,_pij;
          
          _plk_r  = _mm256_mul_pd(_mm256_load_pd(p_lk_rght),_mm256_load_pd(pi));
          _pijplk = _mm256_setzero_pd();
          
          for(i=0;i<ns;++i)
            {
              _pij = _mm256_load_pd(tPij);
              
              
#if (defined(__FMA__))
              _pijplk = _mm256_fmadd_pd(_pij,_mm256_set1_pd(p_lk_left[i]),_pijplk);
#else
              _pijplk = _mm256_add_pd(_pijplk,_mm256_mul_pd(_pij,_mm256_set1_pd(p_lk_left[i])));
#endif
              tPij += ns;
            }
          
          _plk_l = _mm256_mul_pd(_pijplk,_plk_r);      
          
          return(AVX_Vect_Norm(_plk_l));
        }
      return UNLIKELY;
    }
  else
    {
      __m256d _plk_l[nblocks],_plk_r[nblocks];
      __m256d _plk;

      /* [ Pi . Lkr ]' x Pij x Lkl */
      
      if(ambiguity_check == NO) // tip case.
        {
          unsigned int i;
          
          Pij += observed_state*ns;
          
          for(i=0;i<nblocks;++i)
            {
              _plk_l[i] = _mm256_load_pd(p_lk_left);
              _plk_r[i] = _mm256_load_pd(Pij);
              _plk_r[i] = _mm256_mul_pd(_plk_r[i],_plk_l[i]);
              p_lk_left += sz;
              Pij += sz;
            }
          
          _plk = _mm256_setzero_pd();
          for(i=0;i<nblocks;++i) _plk = _mm256_add_pd(_plk,_plk_r[i]);
          return pi[observed_state] * AVX_Vect_Norm(_plk);
        }
      else
        {
          unsigned int i,j;
          __m256d _pij[nblocks],_pijplk[nblocks];
          phydbl lk;
          
          for(i=0;i<nblocks;++i)
            {
              _plk_r[i] = _mm256_mul_pd(_mm256_load_pd(p_lk_rght),_mm256_load_pd(pi));
              p_lk_rght += sz;
              pi += sz;
            }
          
          for(i=0;i<nblocks;++i) _pijplk[i] = _mm256_setzero_pd();
          
          for(i=0;i<ns;++i)
            {
              for(j=0;j<nblocks;++j)
                {
                  _pij[j] = _mm256_load_pd(tPij);
                  tPij += sz;
                  
#if (defined(__FMA__))
                  _pijplk[j] = _mm256_fmadd_pd(_pij[j],_mm256_set1_pd(p_lk_left[i]),_pijplk[j]);
#else
                  _pijplk[j] = _mm256_add_pd(_pijplk[j],_mm256_mul_pd(_pij[j],_mm256_set1_pd(p_lk_left[i])));
#endif
                }
            }
          
          lk = 0.0;
          for(i=0;i<nblocks;++i) lk += AVX_Vect_Norm(_mm256_mul_pd(_pijplk[i],_plk_r[i]));
          return lk;
        }
      
      return UNLIKELY;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl AVX_Lk_Core_One_Class_Eigen_Lr(const phydbl *dot_prod, const phydbl *expl, const unsigned int ns)
{
  unsigned int l;
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;
  __m256d _prod[nblocks],_x;
  
  for(l=0;l<nblocks;++l)
    {
      _prod[l] = _mm256_load_pd(dot_prod); dot_prod += sz;
    }

  if(expl != NULL)
    {
      for(l=0;l<nblocks;++l)
        {
          _prod[l] = _mm256_mul_pd(_prod[l],_mm256_load_pd(expl));
          expl += sz;
        }
    }
  
  _x = _mm256_setzero_pd();
  for(l=0;l<nblocks;++l) _x = _mm256_add_pd(_x,_prod[l]);
  
  return AVX_Vect_Norm(_x);
}
 
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void AVX_Lk_dLk_Core_One_Class_Eigen_Lr(const phydbl *dot_prod, const phydbl *expl, const unsigned int ns, phydbl *lk, phydbl *dlk)
{
  unsigned int i;
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz*2;
  __m256d _x,_y,_z;

  _z = _mm256_setzero_pd();

  for(i=0;i<nblocks;++i)
    {
      _x = _mm256_blend_pd(_mm256_set1_pd(dot_prod[2*i]),
                           _mm256_set1_pd(dot_prod[2*i + 1]),12);
      
      _y = _mm256_load_pd(expl + 4*i);

#if (defined(__FMA__))
      _z = _mm256_fmadd_pd(_x,_y,_z);
#else
      _z = _mm256_add_pd(_z,_mm256_mul_pd(_x,_y));
#endif
    }
  
  *lk = ((double *)&_z)[0] + ((double *)&_z)[2];
  *dlk = ((double *)&_z)[1] + ((double *)&_z)[3];

}
 
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl AVX_Vect_Norm(__m256d _z)
{
  __m128d vlow  = _mm256_castpd256_pd128(_z);
  __m128d vhigh = _mm256_extractf128_pd(_z, 1); // high 128
  vlow  = _mm_add_pd(vlow, vhigh);     // reduce down to 128
  
  __m128d high64 = _mm_unpackhi_pd(vlow, vlow);
  return  _mm_cvtsd_f64(_mm_add_sd(vlow, high64));

      
  /* phydbl r; */
  /* __m256d _x = _mm256_hadd_pd(_z,_z); */
  /* __m256d _y = _mm256_permute2f128_pd(_x,_x,0x21); */
  /* _mm_store_sd(&r,_mm256_castpd256_pd128(_mm256_add_pd(_x,_y))); */
  /* return r; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void AVX_Update_Partial_Lk(t_tree *tree, t_edge *b, t_node *d)
{
/*
           |
           |<- b
           |
           d
          / \
         /   \
        /     \
    n_v1   n_v2
*/
  t_node *n_v1, *n_v2;
  phydbl *plk0,*plk1,*plk2;
  phydbl *Pij1,*Pij2;
  phydbl *tPij1,*tPij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  unsigned int i,j,k;
  int *p_lk_loc;
  
  const unsigned int npattern = tree->data->n_pattern;
  const unsigned int ns = tree->mod->ns;
  const unsigned int ncatg = tree->mod->ras->n_catg;

  const unsigned int ncatgns =  ncatg * ns;
  const unsigned int nsns =  ns * ns;
  
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;

  __m256d *_tPij1,*_tPij2;

  _tPij1     = tree->_tPij1;
  _tPij2     = tree->_tPij2;
    
  n_v1 = n_v2                 = NULL;
  plk0 = plk1 = plk2          = NULL;
  Pij1 = Pij2                 = NULL;
  tPij1 = tPij2               = NULL;
  sum_scale_v1 = sum_scale_v2 = NULL;
  p_lk_loc                    = NULL;

  
  Set_All_Partial_Lk(&n_v1,&n_v2,
                     &plk0,&sum_scale,&p_lk_loc,
                     &Pij1,&tPij1,&plk1,&sum_scale_v1,
                     &Pij2,&tPij2,&plk2,&sum_scale_v2,
                     d,b,tree);


  /* PhyML_Printf("\n. b: %d b->left:%d b->rght:%d d:%d [%p,%p]", */
  /*              b->num, */
  /*              b->left->num, */
  /*              b->rght->num, */
  /*              d->num, */
  /*              plk0, */
  /*              b->p_lk_rght); */
  
  // Copy transpose of transition matrices into AVX registers
  for(i=0;i<ncatg;++i)
    {
      for(j=0;j<ns;++j)
        {
          for(k=0;k<nblocks;++k)
            {
              _tPij1[k] = _mm256_load_pd(tPij1);
              _tPij2[k] = _mm256_load_pd(tPij2);
              tPij1 += sz;
              tPij2 += sz;
            }
          _tPij1 += nblocks;
          _tPij2 += nblocks;
        }
    }
  _tPij1 -= ncatg*ns*nblocks;
  _tPij2 -= ncatg*ns*nblocks;

  if(tree->mod->augmented == YES)
    {
      PhyML_Printf("\n== AVX version of the Update_Partial_Lk function does not");
      PhyML_Printf("\n== allow augmented data.");
      assert(FALSE);
    }
    
  assert(nblocks <= LK_SIMD_MAX_NS);

  {
    const unsigned int stride1 = (n_v1->tax) ? ns : ncatgns;
    const unsigned int stride2 = (n_v2->tax) ? ns : ncatgns;
    const int nt = Lk_Site_Threads(tree,npattern);
    int s;

    /* For every site in the alignment. Each thread uses its own scratch
       vectors and writes only the partial likelihoods of its sites. */
#pragma omp parallel for if(nt > 1) num_threads(nt) schedule(static)
    for(s=0;s<(int)npattern;++s)
      {
        __m256d _pmat1plk1[LK_SIMD_MAX_NS],_pmat2plk2[LK_SIMD_MAX_NS],_plk0[LK_SIMD_MAX_NS];
        const __m256d *s_tPij1,*s_tPij2;
        const phydbl *s_plk1,*s_plk2;
        phydbl *s_plk0,*site_plk0;
        int s_state_v1,s_state_v2;
        int s_ambiguity_check_v1,s_ambiguity_check_v2;
        int s_sum_scale_v1_val,s_sum_scale_v2_val;
        phydbl s_largest_p_lk;
        unsigned int c,k;

        if(tree->data->wght[s] <= SMALL) continue;

        s_state_v1 = s_state_v2 = -1;
        s_ambiguity_check_v1 = s_ambiguity_check_v2 = YES;

        if(n_v1 && n_v1->tax)
          {
            s_ambiguity_check_v1 = n_v1->c_seq->is_ambigu[s];
            if(s_ambiguity_check_v1 == NO) s_state_v1 = n_v1->c_seq->d_state[s];
          }

        if(n_v2 && n_v2->tax)
          {
            s_ambiguity_check_v2 = n_v2->c_seq->is_ambigu[s];
            if(s_ambiguity_check_v2 == NO) s_state_v2 = n_v2->c_seq->d_state[s];
          }

        s_tPij1   = _tPij1;
        s_tPij2   = _tPij2;
        site_plk0 = plk0 + (size_t)s*ncatgns;
        s_plk0    = site_plk0;
        s_plk1    = plk1 + (size_t)s*stride1;
        s_plk2    = plk2 + (size_t)s*stride2;

        for(c=0;c<ncatg;++c)
          {
            if(s_ambiguity_check_v1 == NO && s_ambiguity_check_v2 == NO)
              {
                AVX_Partial_Lk_Exex(s_tPij1,s_state_v1,
                                    s_tPij2,s_state_v2,
                                    ns,_plk0);
              }
            else if(s_ambiguity_check_v1 == YES && s_ambiguity_check_v2 == NO)
              {
                AVX_Partial_Lk_Exin(s_tPij2,s_state_v2,
                                    s_tPij1,s_plk1,_pmat1plk1,
                                    ns,_plk0);
              }
            else if(s_ambiguity_check_v1 == NO && s_ambiguity_check_v2 == YES)
              {
                AVX_Partial_Lk_Exin(s_tPij1,s_state_v1,
                                    s_tPij2,s_plk2,_pmat2plk2,
                                    ns,_plk0);
              }
            else
              {
                AVX_Partial_Lk_Inin(s_tPij1,s_plk1,_pmat1plk1,
                                    s_tPij2,s_plk2,_pmat2plk2,
                                    ns,_plk0);
              }

            for(k=0;k<nblocks;++k) _mm256_store_pd(s_plk0+sz*k,_plk0[k]);

            s_tPij1 += nsns / sz;
            s_tPij2 += nsns / sz;
            s_plk0  += ns;
            s_plk1  += (n_v1->tax) ? 0 : ns;
            s_plk2  += (n_v2->tax) ? 0 : ns;
          }

        if(tree->scaling_method == SCALE_FAST)
          {
            s_sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[s]):(0);
            s_sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[s]):(0);
            sum_scale[s] = s_sum_scale_v1_val + s_sum_scale_v2_val;

            s_largest_p_lk = -BIG;
            for(k=0;k<ncatgns;++k)
              if(site_plk0[k] > s_largest_p_lk)
                s_largest_p_lk = site_plk0[k];

            if(s_largest_p_lk < INV_TWO_TO_THE_LARGE &&
               tree->mod->augmented == NO &&
               tree->apply_lk_scaling == YES)
              {
                for(k=0;k<ncatgns;++k) site_plk0[k] *= TWO_TO_THE_LARGE;
                sum_scale[s] += LARGE;
              }
          }
      }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void AVX_Partial_Lk_Exex(const __m256d *_tPij1, const int state1, const __m256d *_tPij2, const int state2, const int ns, __m256d *plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;

  _tPij1 = _tPij1 + state1 * nblocks;
  _tPij2 = _tPij2 + state2 * nblocks;
  for(i=0;i<nblocks;++i) plk0[i] = _mm256_mul_pd(_tPij1[i],_tPij2[i]);

  /* double *x; */
  /* posix_memalign((void *)&x,BYTE_ALIGN,(size_t)4*sizeof(phydbl)); */
  
  /* _mm256_store_pd(x,plk0[0]); */
  /* for(int i=0;i<4;++i) PhyML_Printf("\n> plk0: %f",x[i]); */

  /*   _mm256_store_pd(x,_tPij1[0]); */
  /* for(int i=0;i<4;++i) PhyML_Printf("\n> Pij1: %f",x[i]); */

  /*   _mm256_store_pd(x,_tPij2[0]); */
  /* for(int i=0;i<4;++i) PhyML_Printf("\n> Pij2: %f",x[i]); */

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void AVX_Partial_Lk_Exin(const __m256d *_tPij1, const int state1, const __m256d *_tPij2, const phydbl *_plk2, __m256d *_pmat2plk2, const int ns, __m256d *_plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;
  
  _tPij1 = _tPij1 + state1 * nblocks;
  AVX_Matrix_Vect_Prod(_tPij2,_plk2,ns,_pmat2plk2);
  
  for(i=0;i<nblocks;++i) _plk0[i] = _mm256_mul_pd(_tPij1[i],_pmat2plk2[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void AVX_Partial_Lk_Inin(const __m256d *_tPij1, const phydbl *plk1, __m256d *_pmat1plk1, const __m256d *_tPij2, const phydbl *plk2, __m256d *_pmat2plk2, const int ns, __m256d *_plk0)
{
  unsigned int i;
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;

  for(i=0;i<ns;++i) if(plk1[i] > 1.0 || plk1[i] < 1.0 || plk2[i] > 1.0 || plk2[i] < 1.0) break; 

  if(i != ns)
    {      
      AVX_Matrix_Vect_Prod(_tPij1,plk1,ns,_pmat1plk1);
      AVX_Matrix_Vect_Prod(_tPij2,plk2,ns,_pmat2plk2);
      
      for(i=0;i<nblocks;++i) _plk0[i] = _mm256_mul_pd(_pmat1plk1[i],_pmat2plk2[i]);
    }
  else
    {
      for(i=0;i<nblocks;++i) _plk0[i] = _mm256_set1_pd(1.0);      
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void AVX_Matrix_Vect_Prod(const __m256d *_m_transpose, const phydbl *_v, const int ns, __m256d *_u)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;  
  unsigned int i,j;  
  __m256d _x;

  _x = _mm256_set1_pd(_v[0]);
  for(j=0;j<nblocks;++j) _u[j] = _mm256_mul_pd(_m_transpose[j],_x);
  _m_transpose = _m_transpose + nblocks;
  
  for(i=1;i<ns;++i)
    {
      _x = _mm256_set1_pd(_v[i]);
      for(j=0;j<nblocks;++j)
        {
#if (defined(__FMA__))
          _u[j] = _mm256_fmadd_pd(_m_transpose[j],_x,_u[j]);
#else
          _u[j] = _mm256_add_pd(_u[j],_mm256_mul_pd(_m_transpose[j],_x));
#endif
        }
      _m_transpose = _m_transpose + nblocks;
    }

/*   for(i=0;i<nblocks;++i) _u[i] = _mm256_setzero_pd(); */
  
/*   for(i=0;i<ns;++i) */
/*     { */
/*       _x = _mm256_set1_pd(_v[i]); */
/*       for(j=0;j<nblocks;++j) */
/*         { */
/* #if (defined(__FMA__)) */
/*           _u[j] = _mm256_fmadd_pd(_m_transpose[j],_x,_u[j]); */
/* #else */
/*           _u[j] = _mm256_add_pd(_u[j],_mm256_mul_pd(_m_transpose[j],_x)); */
/* #endif */
/*         } */
/*       _m_transpose = _m_transpose + nblocks; */
/*     } */

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

__m256d AVX_Horizontal_Add(const __m256d x[4])
{
  __m256d y[2],z[2];

  // y[0] = [x00+x01;x10+x11;x02+x03;x12+x13]
  y[0] = _mm256_hadd_pd(x[0], x[1]);
  // y[1] = [x20+x21;x30+x31;x22+x23;x32+x33]
  y[1] = _mm256_hadd_pd(x[2], x[3]);

  // z[0] = [x00+x01;x10+x11;x22+x23;x32+x33]
  /* z[0] = _mm256_blend_pd(y[0],y[1],0b1100); */
  z[0] = _mm256_blend_pd(y[0],y[1],12);
  // z[1] = [x02+x03;x12+x13;x20+x21;x30+x31]
  z[1] = _mm256_permute2f128_pd(y[0],y[1],0x21);

  return(_mm256_add_pd(z[0],z[1]));
}

#endif
