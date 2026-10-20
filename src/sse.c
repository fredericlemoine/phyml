/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "assert.h"
#include "sse.h"
#ifdef __ARM_NEON
#include "sse2neon.h"
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#if (((defined(__SSE__) || defined(__SSE2__) || defined(__SSE3__) || defined(__ARM_NEON)) && !((defined __AVX__ || defined __AVX2__))) && !defined(DISABLE_NATIVE))

void SSE_Update_Eigen_Lr(t_edge *b, t_tree *tree)
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

  __m128d *_l_ev,*_r_ev;
  
  l_ev         = tree->l_ev;
  _l_ev        = tree->_l_ev;
  _r_ev        = tree->_r_ev;
  
  assert(sz == 2);
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
          _r_ev[i*nblocks+j] = _mm_load_pd(r_ev + j*sz);
          _l_ev[i*nblocks+j] = _mm_load_pd(l_ev + j*sz);
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
        __m128d _prod_left[LK_SIMD_MAX_NS],_prod_rght[LK_SIMD_MAX_NS];
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

            SSE_Matrix_Vect_Prod(_r_ev,p_lk_left_pi,ns,_prod_left);
            SSE_Matrix_Vect_Prod(_l_ev,s_lk_rght,ns,_prod_rght);

            for(k=0;k<nblocks;++k) _mm_store_pd(s_dot_prod + k*sz,_mm_mul_pd(_prod_left[k],_prod_rght[k]));

            s_dot_prod += ns;
            s_lk_left  += step_left;
            s_lk_rght  += step_rght;
          }
      }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl SSE_Lk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, int ns)
{
  phydbl lk;
  unsigned int l;
  const unsigned sz = (int)BYTE_ALIGN / 8;
  const unsigned int nblocks = ns/sz;
  __m128d _prod[nblocks],_x;
    
  for(l=0;l<nblocks;++l) _prod[l] = _mm_load_pd(dot_prod + l*sz);
  if(expl != NULL) for(l=0;l<nblocks;++l) _prod[l] = _mm_mul_pd(_prod[l],_mm_load_pd(expl + l*sz));
  _x = _mm_setzero_pd();
  for(l=0;l<nblocks;++l) _x = _mm_add_pd(_x,_prod[l]);

#if(defined(__SSE3__))
  _x = _mm_hadd_pd(_x,_x);
#else
  PhyML_Printf("\n. SSE3 required. Try turning on the '-msse3' or '-march=native' or '-mcpu=native' option in configure.ac.");
  assert(false);
#endif
  
  _mm_store_sd(&lk,_x);
  
  return lk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SSE_Lk_dLk_Core_One_Class_Eigen_Lr(phydbl *dot_prod, phydbl *expl, unsigned int ns, phydbl *lk, phydbl *dlk)
{
  unsigned int i;
  __m128d _x,_y,_z;

  _z = _mm_setzero_pd();

  for(i=0;i<ns;++i)
    {
      _x = _mm_set1_pd(dot_prod[i]);      
      _y = _mm_load_pd(expl + 2*i);

      _z = _mm_add_pd(_z,_mm_mul_pd(_x,_y));
    }
  
  *lk = ((double *)&_z)[0];
  *dlk = ((double *)&_z)[1];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl SSE_Lk_Core_One_Class_No_Eigen_Lr(phydbl *p_lk_left, phydbl *p_lk_rght, phydbl *Pij, phydbl *tPij, phydbl *pi, int ns, int ambiguity_check, int state)
{
  phydbl lk,dum;
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;
  unsigned int i,j;
  __m128d _plk;
  __m128d _plk_l[nblocks],_plk_r[nblocks];
  
  /* [ Pi . Lkr ]' x Pij x Lkl */
  
  if(ambiguity_check == NO) // tip case.
    {
      for(i=0;i<nblocks;++i)
        {
          _plk_l[i] = _mm_load_pd(p_lk_left);      
          _plk_r[i] = _mm_load_pd(Pij + state*ns);
          p_lk_left += sz;
          Pij += sz;
        }
      
      for(i=0;i<nblocks;++i)
        {
          _plk_r[i] = _mm_mul_pd(_plk_r[i],_mm_set1_pd(pi[state]));
          _plk_r[i] = _mm_mul_pd(_plk_r[i],_plk_l[i]);
        }
      
      _plk = _mm_setzero_pd();
      lk = 0.0;
      for(i=0;i<nblocks;++i)
        {

#if(defined(__SSE3__) || defined(__SSE2__) || defined(__SSE__))
          _plk = _mm_hadd_pd(_plk_r[i],_plk_r[i]);
#else
          PhyML_Printf("\n. SSE3 required. Try turning on the '-msse3' or '-march=native' or '-mcpu=native' option in configure.ac.");
          assert(false);
#endif
          
          _mm_store_sd(&dum,_plk);
          lk += dum;
        }
      return lk;
    }
  else
    {
      __m128d _pij[nblocks],_pijplk[nblocks]; 

      for(i=0;i<nblocks;++i)
        {
          _plk_r[i] = _mm_mul_pd(_mm_load_pd(p_lk_rght),_mm_load_pd(pi));
          p_lk_rght += sz;
          pi += sz;
        }

      for(i=0;i<nblocks;++i) _pijplk[i] = _mm_setzero_pd();

      for(i=0;i<ns;++i)
        {
          for(j=0;j<nblocks;++j)
            {
              _pij[j] = _mm_load_pd(tPij);
              tPij += sz;
              
              _pijplk[j] = _mm_add_pd(_pijplk[j],_mm_mul_pd(_pij[j],_mm_set1_pd(p_lk_left[i])));
            }
        }

      lk = 0.0;
      for(i=0;i<nblocks;++i)
        {
          _plk = _mm_mul_pd(_pijplk[i],_plk_r[i]);

#if(defined(__SSE3__) || defined(__SSE2__) || defined(__SSE__))
          _plk = _mm_hadd_pd(_plk,_plk);
#else
          PhyML_Printf("\n. SSE3 required. Try turning on the '-msse3' or '-march=native' or '-mcpu=native' option in configure.ac.");
          assert(false);
#endif
          
          _mm_store_sd(&dum,_plk);
          lk += dum;
        }
      return lk;
    }
  return UNLIKELY;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Inlined copies of the kernels below, so that SSE_Update_Partial_Lk can be
   specialised for ns = 4 and ns = 20: same operations in the same order,
   with constant loop bounds. */
#if defined(__GNUC__)
#define LK_INLINE static inline __attribute__((always_inline))
#else
#define LK_INLINE static inline
#endif

LK_INLINE void SSE_Matrix_Vect_Prod_Inl(const __m128d *_m_transpose, const phydbl *_v, const int ns, __m128d *_u)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;  
  unsigned int i,j;  
  __m128d _x;


  _x = _mm_set1_pd(_v[0]);
  for(j=0;j<nblocks;++j) _u[j] = _mm_mul_pd(_m_transpose[j],_x);
  _m_transpose = _m_transpose + nblocks;
  
  for(i=1;i<ns;++i)
    {
      _x = _mm_set1_pd(_v[i]);
      for(j=0;j<nblocks;++j)
        {
          _u[j] = _mm_add_pd(_u[j],_mm_mul_pd(_m_transpose[j],_x));
        }
      _m_transpose = _m_transpose + nblocks;
    }

  /* for(i=0;i<nblocks;++i) _u[i] = _mm_setzero_pd(); */

  /* for(i=0;i<ns;++i) */
  /*   { */
  /*     _x = _mm_set1_pd(_v[i]); */
  /*     for(j=0;j<nblocks;++j) _u[j] = _mm_add_pd(_u[j],_mm_mul_pd(_m_transpose[j],_x)); */
  /*     _m_transpose = _m_transpose + nblocks; */
  /*   } */
}

LK_INLINE void SSE_Partial_Lk_Exex_Inl(const __m128d *_tPij1, const int state1, const __m128d *_tPij2, const int state2, const int ns, phydbl *plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;

  _tPij1 = _tPij1 + state1 * nblocks;
  _tPij2 = _tPij2 + state2 * nblocks;
  for(i=0;i<nblocks;++i) _mm_store_pd(plk0+sz*i,_mm_mul_pd(_tPij1[i],_tPij2[i]));
}

LK_INLINE void SSE_Partial_Lk_Exin_Inl(const __m128d *_tPij1, const int state1, const __m128d *_tPij2, const phydbl *_plk2, __m128d *_pmat2plk2, const int ns, phydbl *_plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;
  
  _tPij1 = _tPij1 + state1 * nblocks;
  SSE_Matrix_Vect_Prod_Inl(_tPij2,_plk2,ns,_pmat2plk2);
  
  for(i=0;i<nblocks;++i) _mm_store_pd(_plk0+sz*i,_mm_mul_pd(_tPij1[i],_pmat2plk2[i]));
}

LK_INLINE void SSE_Partial_Lk_Inin_Inl(const __m128d *_tPij1, const phydbl *plk1, __m128d *_pmat1plk1, const __m128d *_tPij2, const phydbl *plk2, __m128d *_pmat2plk2, const int ns, phydbl *_plk0)
{
  unsigned int i;
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  
  for(i=0;i<ns;++i) if(plk1[i] > 1.0 || plk1[i] < 1.0 || plk2[i] > 1.0 || plk2[i] < 1.0) break; 

  if(i != ns)
    {     
      SSE_Matrix_Vect_Prod_Inl(_tPij1,plk1,ns,_pmat1plk1);
      SSE_Matrix_Vect_Prod_Inl(_tPij2,plk2,ns,_pmat2plk2);
      
      for(i=0;i<nblocks;++i) _mm_store_pd(_plk0+sz*i,_mm_mul_pd(_pmat1plk1[i],_pmat2plk2[i]));
    }
  else
    {
      for(i=0;i<nblocks;++i) _mm_store_pd(_plk0+sz*i,_mm_set1_pd(1.0));      
    }
}

LK_INLINE void SSE_Update_Partial_Lk_Ns(t_tree *tree, t_edge *b, t_node *d, const unsigned int ns_c)
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
  const unsigned int ns = ns_c;
  const unsigned int ncatg = tree->mod->ras->n_catg;

  const unsigned int ncatgns =  ncatg * ns;
  const unsigned int nsns =  ns * ns;
  
  const unsigned int sz = (int)BYTE_ALIGN / 8;
  const unsigned nblocks = ns/sz;

  __m128d *_tPij1,*_tPij2;
  
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
 
  // Copy transpose of transition matrices into AVX registers
  for(i=0;i<ncatg;++i)
    {
      for(j=0;j<ns;++j)
        {
          for(k=0;k<nblocks;++k)
            {
              _tPij1[k] = _mm_load_pd(tPij1);
              _tPij2[k] = _mm_load_pd(tPij2);
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
        __m128d _pmat1plk1[LK_SIMD_MAX_NS],_pmat2plk2[LK_SIMD_MAX_NS];
        const __m128d *s_tPij1,*s_tPij2;
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
                SSE_Partial_Lk_Exex_Inl(s_tPij1,s_state_v1,
                                    s_tPij2,s_state_v2,
                                    ns,s_plk0);
              }
            else if(s_ambiguity_check_v1 == YES && s_ambiguity_check_v2 == NO)
              {
                SSE_Partial_Lk_Exin_Inl(s_tPij2,s_state_v2,
                                    s_tPij1,s_plk1,_pmat1plk1,
                                    ns,s_plk0);
              }
            else if(s_ambiguity_check_v1 == NO && s_ambiguity_check_v2 == YES)
              {
                SSE_Partial_Lk_Exin_Inl(s_tPij1,s_state_v1,
                                    s_tPij2,s_plk2,_pmat2plk2,
                                    ns,s_plk0);
              }
            else
              {
                SSE_Partial_Lk_Inin_Inl(s_tPij1,s_plk1,_pmat1plk1,
                                    s_tPij2,s_plk2,_pmat2plk2,
                                    ns,s_plk0);
              }

            // The kernels store straight into the partial likelihoods (a copy
            // from a local array became a memcpy call per site and class)

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

void SSE_Update_Partial_Lk(t_tree *tree, t_edge *b, t_node *d)
{
  // Constant ns for DNA and proteins: the loops over states are unrolled
  if(tree->mod->ns == 4)       SSE_Update_Partial_Lk_Ns(tree,b,d,4);
  else if(tree->mod->ns == 20) SSE_Update_Partial_Lk_Ns(tree,b,d,20);
  else                         SSE_Update_Partial_Lk_Ns(tree,b,d,tree->mod->ns);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SSE_Partial_Lk_Exex(const __m128d *_tPij1, const int state1, const __m128d *_tPij2, const int state2, const int ns, __m128d *plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;

  _tPij1 = _tPij1 + state1 * nblocks;
  _tPij2 = _tPij2 + state2 * nblocks;
  for(i=0;i<nblocks;++i) plk0[i] = _mm_mul_pd(_tPij1[i],_tPij2[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SSE_Partial_Lk_Exin(const __m128d *_tPij1, const int state1, const __m128d *_tPij2, const phydbl *_plk2, __m128d *_pmat2plk2, const int ns, __m128d *_plk0)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  unsigned int i;
  
  _tPij1 = _tPij1 + state1 * nblocks;
  SSE_Matrix_Vect_Prod(_tPij2,_plk2,ns,_pmat2plk2);
  
  for(i=0;i<nblocks;++i) _plk0[i] = _mm_mul_pd(_tPij1[i],_pmat2plk2[i]);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SSE_Partial_Lk_Inin(const __m128d *_tPij1, const phydbl *plk1, __m128d *_pmat1plk1, const __m128d *_tPij2, const phydbl *plk2, __m128d *_pmat2plk2, const int ns, __m128d *_plk0)
{
  unsigned int i;
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;
  
  for(i=0;i<ns;++i) if(plk1[i] > 1.0 || plk1[i] < 1.0 || plk2[i] > 1.0 || plk2[i] < 1.0) break; 

  if(i != ns)
    {     
      SSE_Matrix_Vect_Prod(_tPij1,plk1,ns,_pmat1plk1);
      SSE_Matrix_Vect_Prod(_tPij2,plk2,ns,_pmat2plk2);
      
      for(i=0;i<nblocks;++i) _plk0[i] = _mm_mul_pd(_pmat1plk1[i],_pmat2plk2[i]);
    }
  else
    {
      for(i=0;i<nblocks;++i) _plk0[i] = _mm_set1_pd(1.0);      
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SSE_Matrix_Vect_Prod(const __m128d *_m_transpose, const phydbl *_v, const int ns, __m128d *_u)
{
  unsigned const int sz = (int)BYTE_ALIGN / 8;
  unsigned const int nblocks = ns / sz;  
  unsigned int i,j;  
  __m128d _x;


  _x = _mm_set1_pd(_v[0]);
  for(j=0;j<nblocks;++j) _u[j] = _mm_mul_pd(_m_transpose[j],_x);
  _m_transpose = _m_transpose + nblocks;
  
  for(i=1;i<ns;++i)
    {
      _x = _mm_set1_pd(_v[i]);
      for(j=0;j<nblocks;++j)
        {
          _u[j] = _mm_add_pd(_u[j],_mm_mul_pd(_m_transpose[j],_x));
        }
      _m_transpose = _m_transpose + nblocks;
    }

  /* for(i=0;i<nblocks;++i) _u[i] = _mm_setzero_pd(); */

  /* for(i=0;i<ns;++i) */
  /*   { */
  /*     _x = _mm_set1_pd(_v[i]); */
  /*     for(j=0;j<nblocks;++j) _u[j] = _mm_add_pd(_u[j],_mm_mul_pd(_m_transpose[j],_x)); */
  /*     _m_transpose = _m_transpose + nblocks; */
  /*   } */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

#endif
