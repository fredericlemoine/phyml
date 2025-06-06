/*

PhyML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "models.h"

#ifdef BEAGLE
#include "beagle_utils.h"
#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Handle any number of states (>1) */
void PMat_JC69(phydbl l, int pos, phydbl *Pij, t_mod *mod)
{
  int ns;
  int i, j;

  ns = mod->ns;

  for (i = 0; i < ns; i++)
    Pij[pos + ns * i + i] =
        1. - ((ns - 1.) / ns) * (1. - exp(-ns * l / (ns - 1.)));
  for (i = 0; i < ns - 1; i++)
    for (j = i + 1; j < ns; j++)
    {
      Pij[pos + ns * i + j] = (1. / ns) * (1. - exp(-ns * l / (ns - 1.)));
      if (Pij[pos + ns * i + j] < SMALL_PIJ) Pij[pos + ns * i + j] = SMALL_PIJ;
      Pij[pos + ns * j + i] = Pij[pos + ns * i + j];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PMat_K80(phydbl l, phydbl kappa, int pos, phydbl *Pij)
{
  phydbl Ts, Tv, e1, e2, aux;
  int    i, j;
  /*0 => A*/
  /*1 => C*/
  /*2 => G*/
  /*3 => T*/

  /* Ts -> transition*/
  /* Tv -> transversion*/

  aux = -2 * l / (kappa + 2);
  e1  = (phydbl)exp(aux * 2);

  e2 = (phydbl)exp(aux * (kappa + 1));
  Tv = .25 * (1 - e1);
  Ts = .25 * (1 + e1 - 2 * e2);

  Pij[pos + 4 * 0 + 0] = Pij[pos + 4 * 1 + 1] = Pij[pos + 4 * 2 + 2] =
      Pij[pos + 4 * 3 + 3]                    = 1. - Ts - 2. * Tv;

  Pij[pos + 4 * 0 + 1] = Pij[pos + 4 * 1 + 0] = Tv;
  Pij[pos + 4 * 0 + 2] = Pij[pos + 4 * 2 + 0] = Ts;
  Pij[pos + 4 * 0 + 3] = Pij[pos + 4 * 3 + 0] = Tv;

  Pij[pos + 4 * 1 + 2] = Pij[pos + 4 * 2 + 1] = Tv;
  Pij[pos + 4 * 1 + 3] = Pij[pos + 4 * 3 + 1] = Ts;

  Pij[pos + 4 * 2 + 3] = Pij[pos + 4 * 3 + 2] = Tv;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      if (Pij[pos + 4 * i + j] < SMALL_PIJ) Pij[pos + 4 * i + j] = SMALL_PIJ;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PMat_TN93(phydbl l, t_mod *mod, int pos, phydbl *Pij)
{
  int    i, j;
  phydbl e1, e2, e3;
  phydbl a1t, a2t, bt;
  phydbl A, C, G, T, R, Y;
  phydbl kappa1, kappa2;

  A = mod->e_frq->pi->v[0];
  C = mod->e_frq->pi->v[1];
  G = mod->e_frq->pi->v[2];
  T = mod->e_frq->pi->v[3];
  R = A + G;
  Y = T + C;

  if (mod->kappa->v < .0) mod->kappa->v = 1.0e-5;

  if ((mod->whichmodel != F84) && (mod->whichmodel != TN93))
    mod->lambda->v = 1.;
  else if (mod->whichmodel == F84)
  {
    mod->lambda->v = Get_Lambda_F84(mod->e_frq->pi->v, &mod->kappa->v);
  }

  kappa2 = mod->kappa->v * 2. / (1. + mod->lambda->v);
  kappa1 = kappa2 * mod->lambda->v;

  bt = l / (2. * (A * G * kappa1 + C * T * kappa2 + R * Y));

  a1t = kappa1;
  a2t = kappa2;
  a1t *= bt;
  a2t *= bt;

  e1 = (phydbl)exp(-a1t * R - bt * Y);
  e2 = (phydbl)exp(-a2t * Y - bt * R);
  e3 = (phydbl)exp(-bt);

  /*A->A*/ Pij[pos + 4 * 0 + 0] = A + Y * A / R * e3 + G / R * e1;
  /*A->C*/ Pij[pos + 4 * 0 + 1] = C * (1 - e3);
  /*A->G*/ Pij[pos + 4 * 0 + 2] = G + Y * G / R * e3 - G / R * e1;
  /*A->T*/ Pij[pos + 4 * 0 + 3] = T * (1 - e3);

  /*C->A*/ Pij[pos + 4 * 1 + 0] = A * (1 - e3);
  /*C->C*/ Pij[pos + 4 * 1 + 1] = C + R * C / Y * e3 + T / Y * e2;
  /*C->G*/ Pij[pos + 4 * 1 + 2] = G * (1 - e3);
  /*C->T*/ Pij[pos + 4 * 1 + 3] = T + R * T / Y * e3 - T / Y * e2;

  /*G->A*/ Pij[pos + 4 * 2 + 0] = A + Y * A / R * e3 - A / R * e1;
  /*G->C*/ Pij[pos + 4 * 2 + 1] = C * (1 - e3);
  /*G->G*/ Pij[pos + 4 * 2 + 2] = G + Y * G / R * e3 + A / R * e1;
  /*G->T*/ Pij[pos + 4 * 2 + 3] = T * (1 - e3);

  /*T->A*/ Pij[pos + 4 * 3 + 0] = A * (1 - e3);
  /*T->C*/ Pij[pos + 4 * 3 + 1] = C + R * C / Y * e3 - C / Y * e2;
  /*T->G*/ Pij[pos + 4 * 3 + 2] = G * (1 - e3);
  /*T->T*/ Pij[pos + 4 * 3 + 3] = T + R * T / Y * e3 + C / Y * e2;

  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      if (Pij[pos + 4 * i + j] < SMALL_PIJ) Pij[pos + 4 * i + j] = SMALL_PIJ;

  /*   /\*A->A*\/(*Pij)[0][0] = A+Y*A/R*e3+G/R*e1;  */
  /*   /\*A->C*\/(*Pij)[0][1] = C*(1-e3); */
  /*   /\*A->G*\/(*Pij)[0][2] = G+Y*G/R*e3-G/R*e1; */
  /*   /\*A->T*\/(*Pij)[0][3] = T*(1-e3); */

  /*   /\*C->A*\/(*Pij)[1][0] = A*(1-e3); */
  /*   /\*C->C*\/(*Pij)[1][1] = C+R*C/Y*e3+T/Y*e2; */
  /*   /\*C->G*\/(*Pij)[1][2] = G*(1-e3); */
  /*   /\*C->T*\/(*Pij)[1][3] = T+R*T/Y*e3-T/Y*e2; */

  /*   /\*G->A*\/(*Pij)[2][0] = A+Y*A/R*e3-A/R*e1; */
  /*   /\*G->C*\/(*Pij)[2][1] = C*(1-e3); */
  /*   /\*G->G*\/(*Pij)[2][2] = G+Y*G/R*e3+A/R*e1; */
  /*   /\*G->T*\/(*Pij)[2][3] = T*(1-e3); */

  /*   /\*T->A*\/(*Pij)[3][0] = A*(1-e3); */
  /*   /\*T->C*\/(*Pij)[3][1] = C+R*C/Y*e3-C/Y*e2; */
  /*   /\*T->G*\/(*Pij)[3][2] = G*(1-e3); */
  /*   /\*T->T*\/(*Pij)[3][3] = T+R*T/Y*e3+C/Y*e2; */

  /*   for(i=0;i<4;i++) for(j=0;j<4;j++) */
  /*     if((*Pij)[i][j] < SMALL) (*Pij)[i][j] = SMALL; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Get_Lambda_F84(phydbl *pi, phydbl *kappa)
{
  phydbl A, C, G, T, R, Y, lambda;
  int    kappa_has_changed;

  A = pi[0];
  C = pi[1];
  G = pi[2];
  T = pi[3];
  R = A + G;
  Y = T + C;

  if (*kappa < .0) *kappa = 1.0e-5;

  kappa_has_changed = NO;

  do
  {
    lambda = (Y + (R - Y) / (2. * (*kappa))) / (R - (R - Y) / (2. * (*kappa)));

    if (lambda < .0)
    {
      *kappa += *kappa / 10.;
      kappa_has_changed = YES;
    }
  } while (lambda < .0);

  if (kappa_has_changed)
  {
    PhyML_Printf("\n. WARNING: This transition/transversion ratio\n");
    PhyML_Printf("  is impossible with these base frequencies!\n");
    PhyML_Printf("  The ratio is now set to %.3f\n", *kappa);
  }

  return lambda;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/********************************************************************/
/* void PMat_Empirical(phydbl l, t_mod *mod, phydbl ***Pij)         */
/*                                                                  */
/* Computes the substitution probability matrix                     */
/* from the initial substitution rate matrix and frequency vector   */
/* and one specific branch length                                   */
/*                                                                  */
/* input : l , branch length                                        */
/* input : mod , choosen model parameters, qmat and pi              */
/* ouput : Pij , substitution probability matrix                    */
/*                                                                  */
/* matrix P(l) is computed as follows :                             */
/* P(l) = exp(Q*t) , where :                                        */
/*                                                                  */
/*   Q = substitution rate matrix = Vr*D*inverse(Vr) , where :      */
/*                                                                  */
/*     Vr = right eigenvector matrix for Q                          */
/*     D  = diagonal matrix of eigenvalues for Q                    */
/*                                                                  */
/*   t = time interval = l / mr , where :                           */
/*                                                                  */
/*     mr = mean rate = branch length/time interval                 */
/*        = sum(i)(pi[i]*p(i->j)) , where :                         */
/*                                                                  */
/*       pi = state frequency vector                                */
/*       p(i->j) = subst. probability from i to a different state   */
/*               = -Q[ii] , as sum(j)(Q[ij]) +Q[ii] =0              */
/*                                                                  */
/* the Taylor development of exp(Q*t) gives :                       */
/* P(l) = Vr*exp(D*t)        *inverse(Vr)                           */
/*      = Vr*POW(exp(D/mr),l)*inverse(Vr)                           */
/*                                                                  */
/* for performance we compute only once the following matrixes :    */
/* Vr, inverse(Vr), exp(D/mr)                                       */
/* thus each time we compute P(l) we only have to :                 */
/* make 20 times the operation POW()                                */
/* make 2 20x20 matrix multiplications , that is :                  */
/*   16000 = 2x20x20x20 times the operation *                       */
/*   16000 = 2x20x20x20 times the operation +                       */
/*   which can be reduced to (the central matrix being diagonal) :  */
/*   8400 = 20x20 + 20x20x20 times the operation *                  */
/*   8000 = 20x20x20 times the operation +                          */
/********************************************************************/

void PMat_Empirical(phydbl l, const t_mod *mod, const int pos, phydbl *Pij,
                    phydbl *tPij)
{
  const unsigned int ns = mod->ns;
  unsigned int       i, j, k;
  const phydbl      *U, *V, *R;
  phydbl            *expt;
  phydbl            *uexpt;
  phydbl             sum;

  assert(Pij);

  expt  = mod->eigen->e_val_im;
  uexpt = mod->eigen->r_e_vect_im;
  U     = mod->eigen->r_e_vect;
  V     = mod->eigen->l_e_vect;
  R     = mod->eigen->e_val;

  for (k = 0; k < ns; ++k) expt[k] = exp(R[k] * l);

  /* multiply Vr*POW(exp(D/mr),l)*Vi into Pij */
  for (i = 0; i < ns; i++)
    for (k = 0; k < ns; k++) uexpt[i * ns + k] = U[i * ns + k] * expt[k];

  Pij = Pij + pos;
  if (tPij != NULL) tPij = tPij + pos;

  for (i = 0; i < ns; ++i)
  {
    for (j = 0; j < ns; ++j)
    {
      Pij[j] = 0.0;
      for (k = 0; k < ns; ++k)
      {
        Pij[j] += (uexpt[i * ns + k] * V[k * ns + j]);
      }
      if (Pij[j] < SMALL_PIJ) Pij[j] = SMALL_PIJ;
    }
    /* Below is required ? */
    sum = 0.0;
    for (j = 0; j < ns; ++j) sum += Pij[j];
    for (j = 0; j < ns; ++j) Pij[j] /= sum;

    Pij += ns;
  }

  Pij -= ns * ns;

  if (tPij != NULL)
  {
    for (i = 0; i < ns; ++i)
    {
      for (j = 0; j < ns; ++j)
      {
        tPij[ns * i + j] = Pij[ns * j + i];
      }
    }
  }

  /* PhyML_Printf("\n. Pmat len: %f",l); */
  /* for(i=0;i<ns;i++) */
  /*   { */
  /*     PhyML_Printf("\n"); */
  /*     for(j=0;j<ns;j++) */
  /*       { */
  /*         PhyML_Printf("%12f ",Pij[i*ns+j]); */
  /*       } */
  /*   } */
  /* Exit("\n"); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PMat_Zero_Br_Len(t_mod *mod, int pos, phydbl *Pij)
{
  int n = mod->ns;
  int i, j;

  For(i, n) For(j, n) Pij[pos + mod->ns * i + j] = .0;
  For(i, n) Pij[pos + mod->ns * i + i]           = 1.0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/*
 *Update a rate specific Transition Prob matrix for a given
 *branch-length(already adjusted with the rate prior to this function being
 *called)
 *
 *  Pij: the P-matrix that will be adjusted
 *  l  : branch length * rate
 *  pos: offset into a specific rate-category
 *
 */
void PMat(phydbl l, t_mod *mod, int pos, phydbl *Pij, phydbl *tPij)
{
  /* Warning: l is never the log of branch length here */
  if (l < 0.0)
  {
#ifdef BEAGLE
    Warn_And_Exit(TODO_BEAGLE);
#endif
    PMat_Zero_Br_Len(mod, pos, Pij);
  }
  else
  {
#ifdef BEAGLE
    // when there is no active instance (i.e. when we are building the initial
    // tree)
    if (UNINITIALIZED == mod->b_inst) PMat_Empirical(l, mod, pos, Pij, tPij);
#else
    PMat_Empirical(l, mod, pos, Pij, tPij);
#endif
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int GetDaa(phydbl *daa, phydbl *pi, char *file_name)
{
  /* Get the amino acid distance (or substitution rate) matrix
     (grantham, dayhoff, jones, etc).
  */
  FILE  *fdaa;
  int    i, j, naa;
  phydbl dmax, dmin;
  phydbl sum;
  double val;

  naa  = 20;
  dmax = .0;
  dmin = 1.E+40;

  fdaa = (FILE *)Openfile(file_name, 0);

  for (i = 0; i < naa; i++)
    for (j = 0; j < i; j++)
    {
      /* 	 if(fscanf(fdaa, "%lf", &daa[i*naa+j])) Exit("\n. err
       * aaRatefile"); */
      if (fscanf(fdaa, "%lf", &val)) Exit("\n. err aaRatefile");
      daa[i * naa + j] = (phydbl)val;
      daa[j * naa + i] = daa[i * naa + j];
      if (dmax < daa[i * naa + j]) dmax = daa[i * naa + j];
      if (dmin > daa[i * naa + j]) dmin = daa[i * naa + j];
    }

  for (i = 0; i < naa; i++)
  {
    /*        if(fscanf(fdaa,"%lf",&pi[i])!=1) Exit("\n. err aaRatefile"); */
    if (fscanf(fdaa, "%lf", &val) != 1) Exit("\n. err aaRatefile");
    pi[i] = (phydbl)val;
  }

  sum = 0.0;
  for (i = 0; i < naa; ++i) sum += pi[i];
  if (FABS(1. - sum) > 1e-4)
  {
    PhyML_Printf("\nSum of freq. = %.6f != 1 in aaRateFile\n", sum);
    exit(-1);
  }

  fclose(fdaa);

  return (0);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Update_Qmat_Generic(phydbl *rr, phydbl *pi, int ns, phydbl *qmat)
{
  int    i, j;
  phydbl sum, mr;

  for (i = 0; i < ns * ns; ++i) qmat[i] = .0;

  /* if(rr[(int)(ns*(ns-1)/2)-1] < 0.00001) */
  /*   { */
  /*     PhyML_Fprintf(stderr,"\n.
   * rr[%d]=%f",(int)(ns*(ns-1)/2)-1,rr[(int)(ns*(ns-1)/2)-1]); */
  /*     assert(FALSE); */
  /*   } */

  /* for(i=0;i<(int)(ns*(ns-1)/2);++i) */
  /*   { */
  /*     rr[i] /= rr[(int)(ns*(ns-1)/2)-1]; */
  /*   } */

  /* Fill the non-diagonal parts */
  for (i = 0; i < ns; i++)
  {
    for (j = i + 1; j < ns; j++)
    {
      qmat[i * ns + j] = rr[MIN(i, j) * ns + MAX(i, j) -
                            (MIN(i, j) + 1 + (int)POW(MIN(i, j) + 1, 2)) / 2];
      qmat[j * ns + i] = qmat[i * ns + j];
    }
  }

  /* Multiply by pi */
  for (i = 0; i < ns; i++)
  {
    for (j = 0; j < ns; j++)
    {
      qmat[i * ns + j] *= pi[j];
    }
  }

  /* Compute diagonal elements */
  mr = .0;
  for (i = 0; i < ns; i++)
  {
    sum = .0;
    for (j = 0; j < ns; j++) sum += qmat[i * ns + j];
    qmat[i * ns + i] = -sum;
    mr += sum * pi[i];
  }

  for (i = 0; i < ns * ns; i++) qmat[i] /= mr;

  return (mr);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_Qmat_GTR(phydbl *rr, phydbl *rr_val, int *rr_num, phydbl *pi,
                     phydbl *qmat)
{
  int    i;
  phydbl mr;

  for (i = 0; i < 6; i++) rr[i] = exp(rr_val[rr_num[i]]);
  
  /* for(i=0;i<6;i++) rr[i] = log(rr_val[rr_num[i]]); */
  for (i = 0; i < 6; i++)
    if (rr[i] < 0.0)
    {
      PhyML_Fprintf(stderr, "\n. rr%d: %f", i, rr[i]);
      PhyML_Fprintf(stderr, "\n. Err. in file %s at line %d (function '%s').\n",
                    __FILE__, __LINE__, __FUNCTION__);
      Exit("");
    }

  for (i = 0; i < 6; i++) rr[i] /= rr[5];
  for (i = 0; i < 6; i++)
    if (rr[i] < RR_MIN) rr[i] = RR_MIN;
  for (i = 0; i < 6; i++)
    if (rr[i] > RR_MAX) rr[i] = RR_MAX;

  qmat[0 * 4 + 1] = (rr[0] * pi[1]);
  qmat[0 * 4 + 2] = (rr[1] * pi[2]);
  qmat[0 * 4 + 3] = (rr[2] * pi[3]);

  qmat[1 * 4 + 0] = (rr[0] * pi[0]);
  qmat[1 * 4 + 2] = (rr[3] * pi[2]);
  qmat[1 * 4 + 3] = (rr[4] * pi[3]);

  qmat[2 * 4 + 0] = (rr[1] * pi[0]);
  qmat[2 * 4 + 1] = (rr[3] * pi[1]);
  qmat[2 * 4 + 3] = (rr[5] * pi[3]);

  qmat[3 * 4 + 0] = (rr[2] * pi[0]);
  qmat[3 * 4 + 1] = (rr[4] * pi[1]);
  qmat[3 * 4 + 2] = (rr[5] * pi[2]);

  qmat[0 * 4 + 0] = -(rr[0] * pi[1] + rr[1] * pi[2] + rr[2] * pi[3]);
  qmat[1 * 4 + 1] = -(rr[0] * pi[0] + rr[3] * pi[2] + rr[4] * pi[3]);
  qmat[2 * 4 + 2] = -(rr[1] * pi[0] + rr[3] * pi[1] + rr[5] * pi[3]);
  qmat[3 * 4 + 3] = -(rr[2] * pi[0] + rr[4] * pi[1] + rr[5] * pi[2]);

  /* compute diagonal terms of Q and mean rate mr = l/t */
  mr = .0;
  for (i = 0; i < 4; i++) mr += pi[i] * (-qmat[i * 4 + i]);
  for (i = 0; i < 16; i++) qmat[i] /= mr;
  /* PhyML_Printf("\n. rr[0]=%f rr[1]=%f rr[2]=%f rr[3]=%f rr[4]=%f rr[5]=%f",
   */
  /*              rr[0], */
  /*              rr[1], */
  /*              rr[2], */
  /*              rr[3], */
  /*              rr[4], */
  /*              rr[5]); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_Qmat_HKY(phydbl kappa, phydbl *pi, phydbl *qmat)
{
  int    i;
  phydbl mr;

  /* A -> C */ qmat[0 * 4 + 1] = (phydbl)(pi[1]);
  /* A -> G */ qmat[0 * 4 + 2] = (phydbl)(kappa * pi[2]);
  /* A -> T */ qmat[0 * 4 + 3] = (phydbl)(pi[3]);

  /* C -> A */ qmat[1 * 4 + 0] = (phydbl)(pi[0]);
  /* C -> G */ qmat[1 * 4 + 2] = (phydbl)(pi[2]);
  /* C -> T */ qmat[1 * 4 + 3] = (phydbl)(kappa * pi[3]);

  /* G -> A */ qmat[2 * 4 + 0] = (phydbl)(kappa * pi[0]);
  /* G -> C */ qmat[2 * 4 + 1] = (phydbl)(pi[1]);
  /* G -> T */ qmat[2 * 4 + 3] = (phydbl)(pi[3]);

  /* T -> A */ qmat[3 * 4 + 0] = (phydbl)(pi[0]);
  /* T -> C */ qmat[3 * 4 + 1] = (phydbl)(kappa * pi[1]);
  /* T -> G */ qmat[3 * 4 + 2] = (phydbl)(pi[2]);

  qmat[0 * 4 + 0] =
      (phydbl)(-(qmat[0 * 4 + 1] + qmat[0 * 4 + 2] + qmat[0 * 4 + 3]));
  qmat[1 * 4 + 1] =
      (phydbl)(-(qmat[1 * 4 + 0] + qmat[1 * 4 + 2] + qmat[1 * 4 + 3]));
  qmat[2 * 4 + 2] =
      (phydbl)(-(qmat[2 * 4 + 0] + qmat[2 * 4 + 1] + qmat[2 * 4 + 3]));
  qmat[3 * 4 + 3] =
      (phydbl)(-(qmat[3 * 4 + 0] + qmat[3 * 4 + 1] + qmat[3 * 4 + 2]));

  /* compute diagonal terms of Q and mean rate mr = l/t */
  mr = .0;
  for (i = 0; i < 4; ++i) mr += pi[i] * (-qmat[i * 4 + i]);
  for (i = 0; i < 16; i++) qmat[i] /= mr;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_Qmat_TN93(phydbl kappa, phydbl lambda, phydbl *pi, phydbl *qmat)
{
  int    i;
  phydbl mr;

  /* A -> C */ qmat[0 * 4 + 1] = (phydbl)(pi[1]);
  /* A -> G */ qmat[0 * 4 + 2] = (phydbl)(kappa * lambda * pi[2]);
  /* A -> T */ qmat[0 * 4 + 3] = (phydbl)(pi[3]);

  /* C -> A */ qmat[1 * 4 + 0] = (phydbl)(pi[0]);
  /* C -> G */ qmat[1 * 4 + 2] = (phydbl)(pi[2]);
  /* C -> T */ qmat[1 * 4 + 3] = (phydbl)(kappa * pi[3]);

  /* G -> A */ qmat[2 * 4 + 0] = (phydbl)(kappa * lambda * pi[0]);
  /* G -> C */ qmat[2 * 4 + 1] = (phydbl)(pi[1]);
  /* G -> T */ qmat[2 * 4 + 3] = (phydbl)(pi[3]);

  /* T -> A */ qmat[3 * 4 + 0] = (phydbl)(pi[0]);
  /* T -> C */ qmat[3 * 4 + 1] = (phydbl)(kappa * pi[1]);
  /* T -> G */ qmat[3 * 4 + 2] = (phydbl)(pi[2]);

  qmat[0 * 4 + 0] =
      (phydbl)(-(qmat[0 * 4 + 1] + qmat[0 * 4 + 2] + qmat[0 * 4 + 3]));
  qmat[1 * 4 + 1] =
      (phydbl)(-(qmat[1 * 4 + 0] + qmat[1 * 4 + 2] + qmat[1 * 4 + 3]));
  qmat[2 * 4 + 2] =
      (phydbl)(-(qmat[2 * 4 + 0] + qmat[2 * 4 + 1] + qmat[2 * 4 + 3]));
  qmat[3 * 4 + 3] =
      (phydbl)(-(qmat[3 * 4 + 0] + qmat[3 * 4 + 1] + qmat[3 * 4 + 2]));

  /* compute diagonal terms of Q and mean rate mr = l/t */
  mr = .0;
  for (i = 0; i < 4; ++i) mr += pi[i] * (-qmat[i * 4 + i]);
  for (i = 0; i < 16; i++) qmat[i] /= mr;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Translate_Custom_Mod_String(t_mod *mod)
{
  int i, j;

  for (i = 0; i < 6; i++) mod->r_mat->n_rr_per_cat->v[i] = 0;

  mod->r_mat->n_diff_rr = 0;

  for (i = 0; i < 6; i++)
  {
    for (j = 0; j < i; j++)
    {
      if (mod->custom_mod_string->s[i] == mod->custom_mod_string->s[j])
      {
        break;
      }
    }

    if (i == j)
    {
      mod->r_mat->rr_num->v[i] = mod->r_mat->n_diff_rr;
      mod->r_mat->n_diff_rr++;
    }
    else
    {
      mod->r_mat->rr_num->v[i] = mod->r_mat->rr_num->v[j];
    }

    mod->r_mat->n_rr_per_cat->v[mod->r_mat->rr_num->v[j]]++;
  }

  /* PhyML_Printf("\n"); */
  /* for(i=0;i<6;i++) PhyML_Printf("%d ",mod->r_mat->rr_num->v[i]); */
  /* for(i=0;i<mod->r_mat->n_diff_rr;i++) PhyML_Printf("\n. Class %d size
   * %d",i+1,mod->r_mat->n_rr_per_cat->v[i]); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

// Update rate across sites distribution settings.

int Update_RAS(t_mod *mod)
{
  phydbl sum;
  int    i;

  if (mod->ras->free_mixt_rates == NO)
  {
    DiscreteGamma(mod->ras->gamma_r_proba->v, mod->ras->gamma_rr->v,
                  mod->ras->alpha->v, mod->ras->alpha->v, mod->ras->n_catg,
                  mod->ras->gamma_median);
  }
  else
  {
    // Update class frequencies
#ifndef PHYML
    Qksort(mod->ras->gamma_r_proba_unscaled->v, NULL, 0,
           mod->ras->n_catg -
               1); // Unscaled class frequencies sorted in increasing order

    mod->ras->gamma_r_proba->v[0] =
        mod->ras->gamma_r_proba_unscaled->v[0] /
        mod->ras->gamma_r_proba_unscaled->v[mod->ras->n_catg - 1];
    for (i = 1; i < mod->ras->n_catg; i++)
    {
      mod->ras->gamma_r_proba->v[i] =
          (mod->ras->gamma_r_proba_unscaled->v[i] -
           mod->ras->gamma_r_proba_unscaled->v[i - 1]) /
          mod->ras->gamma_r_proba_unscaled->v[mod->ras->n_catg - 1];
    }
#endif

#ifdef PHYML
    sum = 0.0;
    for (i = 0; i < mod->ras->n_catg; i++)
      sum += exp(mod->ras->gamma_r_proba_unscaled->v[i]);
    for (i = 0; i < mod->ras->n_catg; i++)
      mod->ras->gamma_r_proba->v[i] =
          exp(mod->ras->gamma_r_proba_unscaled->v[i]) / sum;
#endif

    // Update class rates
    if (mod->ras->normalise_rr == YES)
    {
#ifdef PHYML
      sum = .0;
      for (i = 0; i < mod->ras->n_catg; i++)
        sum += mod->ras->gamma_r_proba->v[i] *
               exp(mod->ras->gamma_rr_unscaled->v[i]);
      for (i = 0; i < mod->ras->n_catg; i++)
        mod->ras->gamma_rr->v[i] = exp(mod->ras->gamma_rr_unscaled->v[i]) / sum;
#else
      sum = .0;
      for (i = 0; i < mod->ras->n_catg; i++)
        sum += mod->ras->gamma_r_proba->v[i] *
               fabs(mod->ras->gamma_rr_unscaled->v[i]);
      for (i = 0; i < mod->ras->n_catg; i++)
        mod->ras->gamma_rr->v[i] =
            fabs(mod->ras->gamma_rr_unscaled->v[i]) / sum;
#endif
    }
    else
    {
      for (i = 0; i < mod->ras->n_catg; i++)
        mod->ras->gamma_rr->v[i] = exp(mod->ras->gamma_rr_unscaled->v[i]);
    }

    /* printf("\n"); */
    /* for(i=0;i<mod->ras->n_catg;i++) */
    /*   printf("\nx %3d %12f %12f xx %12f %12f", */
    /*          mod->ras->normalise_rr, */
    /*          mod->ras->gamma_r_proba->v[i], */
    /*          mod->ras->gamma_rr->v[i], */
    /*          mod->ras->gamma_r_proba_unscaled->v[i], */
    /*          mod->ras->gamma_rr_unscaled->v[i]); */
  }
#ifdef BEAGLE
  if (UNINITIALIZED != mod->b_inst) update_beagle_ras(mod);
#endif

  /* for(i=0;i<mod->ras->n_catg;i++) PhyML_Printf("\n. REALW%d: %12f [%12f]
   * (%12f;%12f--%2d)", */
  /*                                              i, */
  /*                                              mod->ras->gamma_r_proba->v[i],
   */
  /*                                              mod->ras->gamma_r_proba_unscaled->v[i],
   */
  /*                                              mod->ras->gamma_rr->v[i], */
  /*                                              mod->ras->gamma_rr_unscaled->v[i],
   */
  /*                                              mod->ras->normalise_rr); */

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Update_Efrq(t_mod *mod)
{
  unsigned int i;
  phydbl       sum;

  if (mod->e_frq->type == ML)
  {
    for (i = 0; i < mod->ns; ++i)
      mod->e_frq->pi->v[i] = exp(mod->e_frq->pi_unscaled->v[i]);
    sum = 0.0;
    for (i = 0; i < mod->ns; ++i) sum += mod->e_frq->pi->v[i];
    for (i = 0; i < mod->ns; ++i)
      mod->e_frq->pi->v[i] = mod->e_frq->pi->v[i] / sum;

#ifdef BEAGLE
    if (UNINITIALIZED != mod->b_inst) update_beagle_efrqs(mod);
#endif
  }

  for (i = 0; i < mod->ns; ++i)
    if (mod->e_frq->pi->v[i] < E_FRQ_MIN) mod->e_frq->pi->v[i] = E_FRQ_MIN;
  for (i = 0; i < mod->ns; ++i)
    if (mod->e_frq->pi->v[i] > E_FRQ_MAX) mod->e_frq->pi->v[i] = E_FRQ_MAX;

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Set_Model_Parameters(t_mod *mod)
{
  if (!Update_RAS(mod)) return 0;
  if (!Update_Efrq(mod)) return 0;
  if (!Update_Eigen(mod)) return 0;
  if (!Update_Boundaries(mod)) return 0;
  if (mod->is_mixt_mod == YES) MIXT_Set_Model_Parameters(mod);
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Update_Boundaries(t_mod *mod)
{
  int i;

  if (mod->kappa->v > TSTV_MAX) mod->kappa->v = TSTV_MAX;
  if (mod->kappa->v < TSTV_MIN) mod->kappa->v = TSTV_MIN;

  if (mod->ras->alpha->v < ALPHA_MIN) mod->ras->alpha->v = ALPHA_MIN;
  if (mod->ras->alpha->v > ALPHA_MAX) mod->ras->alpha->v = ALPHA_MAX;

  if (mod->ras->free_mixt_rates == YES)
  {
    for (i = 0; i < mod->ras->n_catg; i++)
    {
      if (mod->ras->gamma_rr_unscaled->v[i] < GAMMA_RR_UNSCALED_MIN)
        mod->ras->gamma_rr_unscaled->v[i] = GAMMA_RR_UNSCALED_MIN;

      if (mod->ras->gamma_rr_unscaled->v[i] > GAMMA_RR_UNSCALED_MAX)
        mod->ras->gamma_rr_unscaled->v[i] = GAMMA_RR_UNSCALED_MAX;

      if (mod->ras->gamma_r_proba_unscaled->v[i] < GAMMA_R_PROBA_UNSCALED_MIN)
        mod->ras->gamma_r_proba_unscaled->v[i] = GAMMA_R_PROBA_UNSCALED_MIN;

      if (mod->ras->gamma_r_proba_unscaled->v[i] > GAMMA_R_PROBA_UNSCALED_MAX)
        mod->ras->gamma_r_proba_unscaled->v[i] = GAMMA_R_PROBA_UNSCALED_MAX;
    }
  }

  if (mod->whichmodel == CUSTOM || mod->whichmodel == GTR)
  {
    for (i = 0; i < 6; i++)
      if (mod->r_mat->rr_val->v[i] < UNSCALED_RR_MIN)
        mod->r_mat->rr_val->v[i] = UNSCALED_RR_MIN;
    for (i = 0; i < 6; i++)
      if (mod->r_mat->rr_val->v[i] > UNSCALED_RR_MAX)
        mod->r_mat->rr_val->v[i] = UNSCALED_RR_MAX;

    for (i = 0; i < 6; i++)
      if (mod->r_mat->rr->v[i] < RR_MIN) mod->r_mat->rr->v[i] = RR_MIN;
    for (i = 0; i < 6; i++)
      if (mod->r_mat->rr->v[i] > RR_MAX) mod->r_mat->rr->v[i] = RR_MAX;
  }

  for (i = 0; i < mod->ns; i++)
  {
    if (mod->e_frq->pi_unscaled->v[i] < UNSCALED_E_FRQ_MIN)
      mod->e_frq->pi_unscaled->v[i] = UNSCALED_E_FRQ_MIN;

    if (mod->e_frq->pi_unscaled->v[i] > UNSCALED_E_FRQ_MAX)
      mod->e_frq->pi_unscaled->v[i] = UNSCALED_E_FRQ_MAX;

    if (mod->e_frq->pi->v[i] < E_FRQ_MIN) mod->e_frq->pi->v[i] = E_FRQ_MIN;

    if (mod->e_frq->pi->v[i] > E_FRQ_MAX) mod->e_frq->pi->v[i] = E_FRQ_MAX;
  }

  if (mod->r_mat_weight->v < R_MAT_WEIGHT_MIN)
    mod->r_mat_weight->v = R_MAT_WEIGHT_MIN;
  if (mod->r_mat_weight->v > R_MAT_WEIGHT_MAX)
    mod->r_mat_weight->v = R_MAT_WEIGHT_MAX;

  if (mod->e_frq_weight->v < E_FRQ_WEIGHT_MIN)
    mod->e_frq_weight->v = E_FRQ_WEIGHT_MIN;
  if (mod->e_frq_weight->v > E_FRQ_WEIGHT_MAX)
    mod->e_frq_weight->v = E_FRQ_WEIGHT_MAX;

  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int Update_Eigen(t_mod *mod)
{
  int    result, n_iter;
  phydbl scalar, eig_res;
  int    i;

  if (mod->update_eigen == YES)
  {
    // Update the Q-matrix first before computing the eigen decomposition
    // (because eigen values/vectors is computed based on the Q-matrix)
    if (mod->io->datatype == NT)
    {
      if (mod->whichmodel == GTR)
        Update_Qmat_GTR(mod->r_mat->rr->v, mod->r_mat->rr_val->v,
                        mod->r_mat->rr_num->v, mod->e_frq->pi->v,
                        mod->r_mat->qmat->v);
      else if (mod->whichmodel == CUSTOM)
        Update_Qmat_GTR(mod->r_mat->rr->v, mod->r_mat->rr_val->v,
                        mod->r_mat->rr_num->v, mod->e_frq->pi->v,
                        mod->r_mat->qmat->v);
      else if (mod->whichmodel == HKY85)
        Update_Qmat_HKY(mod->kappa->v, mod->e_frq->pi->v, mod->r_mat->qmat->v);
      else if (mod->whichmodel == TN93)
        Update_Qmat_TN93(mod->kappa->v, mod->lambda->v, mod->e_frq->pi->v,
                         mod->r_mat->qmat->v);
      else /* Any other nucleotide-based t_mod */
        Update_Qmat_HKY(mod->kappa->v, mod->e_frq->pi->v, mod->r_mat->qmat->v);
    }
    else if (mod->io->datatype == AA)
    {
      Update_Qmat_Generic(mod->r_mat->rr_val->v, mod->e_frq->pi->v, mod->ns,
                          mod->r_mat->qmat->v);
    }
    else if (mod->io->datatype == GENERIC)
    {
      Update_Qmat_Generic(mod->r_mat->rr_val->v, mod->e_frq->pi->v, mod->ns,
                          mod->r_mat->qmat->v);
    }

    // Now compute the eigen vectors and values
    scalar = 1.0;
    n_iter = 0;
    result = 0;

    for (i = 0; i < mod->ns * mod->ns; ++i)
    {
      mod->r_mat->qmat_buff->v[i] = mod->r_mat->qmat->v[i];
    }

    /* For(i,mod->eigen->size*mod->eigen->size) PhyML_Printf("\n. qmat
     * [%d,%d]=%f",i/mod->eigen->size,i%mod->eigen->size,mod->r_mat->qmat->v[i]);
     */
    /* for(i=0;i<6;i++) PhyML_Printf("\n. rr: %f rr_Val:
     * %f",mod->r_mat->rr_val->v[i],mod->r_mat->rr->v[i]); */
    /* for(i=0;i<4;i++) PhyML_Printf("\n. f: %f
     * %f",mod->e_frq->pi->v[i],mod->e_frq->pi_unscaled->v[i]); */
    /* PhyML_Printf("\n"); */

    /* compute eigenvectors/values */
    /*       if(!EigenRealGeneral(mod->eigen->size,mod->r_mat->qmat,mod->eigen->e_val,
     */
    /* 			  mod->eigen->e_val_im, mod->eigen->r_e_vect, */
    /* 			  mod->eigen->space_int,mod->eigen->space)) */

    eig_res =
        Eigen(1, mod->r_mat->qmat_buff->v, mod->eigen->size, mod->eigen->e_val,
              mod->eigen->e_val_im, mod->eigen->r_e_vect,
              mod->eigen->r_e_vect_im, mod->eigen->space);
    if (eig_res == 0)
    {
      /* compute inverse(Vr) into Vi */
      for (i = 0; i < mod->ns * mod->ns; ++i)
        mod->eigen->l_e_vect[i] = mod->eigen->r_e_vect[i];
      while (!Matinv(mod->eigen->l_e_vect, mod->eigen->size, mod->eigen->size,
                     YES))
      {
        PhyML_Printf("\n. Trying Q<-Q*scalar and then Root<-Root/scalar to fix "
                     "this...\n");
        scalar += scalar / 3.;
        For(i, mod->eigen->size * mod->eigen->size)
            PhyML_Printf("\n. qmat [%d,%d]=%f", i / mod->eigen->size,
                         i % mod->eigen->size, mod->r_mat->qmat->v[i]);
        For(i, mod->eigen->size * mod->eigen->size)
            mod->r_mat->qmat_buff->v[i] = mod->r_mat->qmat->v[i];
        For(i, mod->eigen->size * mod->eigen->size)
            mod->r_mat->qmat_buff->v[i] *= scalar;
        result =
            Eigen(1, mod->r_mat->qmat_buff->v, mod->eigen->size,
                  mod->eigen->e_val, mod->eigen->e_val_im, mod->eigen->r_e_vect,
                  mod->eigen->r_e_vect_im, mod->eigen->space);
        if (result == -1)
        {
          PhyML_Fprintf(stderr, "\n. Eigenvalues/vectors computation did not "
                                "converge: computation cancelled.");
          Exit("\n");
        }
        else if (result == 1)
        {
          PhyML_Fprintf(
              stderr,
              "\n. Complex eigenvalues/vectors: computation cancelled.");
          Exit("\n");
        }

        for (i = 0; i < mod->eigen->size * mod->eigen->size; ++i)
          mod->eigen->l_e_vect[i] = mod->eigen->r_e_vect[i];
        n_iter++;
        if (n_iter > 100)
        {
          PhyML_Fprintf(stderr, "\n. Cannot work out eigen vectors.");
          Exit("\n");
        }
      }
      for (i = 0; i < mod->eigen->size; i++) mod->eigen->e_val[i] /= scalar;

#ifdef BEAGLE
      // Recall that BEAGLE is initialized *after* all the model parameters are
      // set IOW, this function may be called before BEAGLE is initialized
      // ("chicken-egg")
      if (UNINITIALIZED != mod->b_inst) update_beagle_eigen(mod);
#endif
    }
    else if (eig_res == -1)
    {
      PhyML_Fprintf(stderr, "\n. kappa: %f", mod->kappa->v);
      for (int i = 0; i < mod->ns; ++i)
      {
        for (int j = 0; j < mod->ns; ++j)
        {
          PhyML_Fprintf(stderr, "\n. Q[%d,%d]=%g", i, j,
                        mod->r_mat->qmat->v[i * mod->ns + j]);
        }
      }

      PhyML_Fprintf(stderr, "\n. Eigenvalues/vectors computation does not "
                            "converge. Computation cancelled.");
      return 0;
    }
    else if (eig_res == 1)
    {
      PhyML_Fprintf(stderr, "\n. WARNING: imaginary eigenvectors not null.");
    }
  }
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Switch_From_M4mod_To_Mod(t_mod *mod)
{
  int i;

  mod->use_m4mod = 0;
  mod->ns        = mod->m4mod->n_o;
  for (i = 0; i < mod->ns; i++) mod->e_frq->pi->v[i] = mod->m4mod->o_fq[i];
  mod->eigen->size = mod->ns;
  Set_Update_Eigen(YES, mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PMat_MGF_Gamma(phydbl mu, phydbl sigsq, const t_mod *mod, const int pos,
                    phydbl *Pij, phydbl *tPij)
{
  const unsigned int ns = mod->ns;
  unsigned int       i, j, k;
  phydbl            *uexpt, *imbd;

  Pij = Pij + pos;
  if (tPij != NULL) tPij = tPij + pos;

  uexpt = mod->eigen->r_e_vect_im;
  imbd  = mod->eigen->e_val_im;

  for (i = 0; i < ns; i++)
    imbd[i] = POW(1. - mod->eigen->e_val[i] * sigsq / mu, -mu * mu / sigsq);
  for (i = 0; i < ns; i++)
    for (k = 0; k < ns; k++)
      uexpt[i * ns + k] = mod->eigen->r_e_vect[i * ns + k] * imbd[k];

  for (i = 0; i < ns; i++)
    for (k = 0; k < ns; k++) Pij[ns * i + k] = .0;

  for (i = 0; i < ns; i++)
  {
    for (j = 0; j < ns; j++)
    {
      for (k = 0; k < ns; k++)
      {
        Pij[j] += (uexpt[i * ns + k] * mod->eigen->l_e_vect[k * ns + j]);
      }
      if (Pij[j] < SMALL_PIJ) Pij[ns * i + j] = SMALL_PIJ;
    }
    Pij += ns;
  }

  Pij -= ns * ns;

  if (tPij != NULL)
  {
    for (i = 0; i < ns; ++i)
    {
      for (j = 0; j < ns; ++j)
      {
        tPij[ns * i + j] = Pij[ns * j + i];
      }
    }
  }

  /* PhyML_Printf("\n. MGF Pmat"); */
  /* for(i=0;i<ns;i++) */
  /*   { */
  /*     PhyML_Printf("\n"); */
  /*     for(j=0;j<ns;j++) */
  /*       { */
  /*         PhyML_Printf("%12f ",Pij[i*ns+j]); */
  /*       } */
  /*   } */
  /* Exit("\n"); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Switch_From_Mod_To_M4mod(t_mod *mod)
{
  int i;
  mod->use_m4mod = 1;
  mod->ns        = mod->m4mod->n_o * mod->m4mod->n_h;
  for (i = 0; i < mod->ns; i++)
    mod->e_frq->pi->v[i] = mod->m4mod->o_fq[i % mod->m4mod->n_o] *
                           mod->m4mod->h_fq->v[i / mod->m4mod->n_o];
  mod->eigen->size = mod->ns;
  Set_Update_Eigen(YES, mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl General_Dist(phydbl *F, t_mod *mod, eigen *eigen_struct)
{
  phydbl *pi, *mod_pi;
  int     i, j, k;
  phydbl  dist;
  phydbl *F_phydbl;

  /* TO DO : call eigen decomposition function for all nt models */

  F_phydbl = (phydbl *)mCalloc(eigen_struct->size * eigen_struct->size,
                               sizeof(phydbl));
  pi       = (phydbl *)mCalloc(eigen_struct->size, sizeof(phydbl));
  mod_pi   = (phydbl *)mCalloc(eigen_struct->size, sizeof(phydbl));

  for (i = 0; i < mod->ns; i++) mod_pi[i] = mod->e_frq->pi->v[i];

  for (i = 0; i < eigen_struct->size; i++)
  {
    for (j = 0; j < eigen_struct->size; j++)
    {
      pi[i] +=
          (F[eigen_struct->size * i + j] + F[eigen_struct->size * j + i]) / 2.;
    }
  }

  Make_Symmetric(&F, eigen_struct->size);
  Divide_Mat_By_Vect(&F, mod->e_frq->pi->v, eigen_struct->size);

  /* Eigen decomposition of pi^{-1} x F */
  for (i = 0; i < eigen_struct->size; i++)
    for (j = 0; j < eigen_struct->size; j++)
      F_phydbl[eigen_struct->size * i + j] = F[eigen_struct->size * i + j];

  if (Eigen(1, F_phydbl, mod->eigen->size, mod->eigen->e_val,
            mod->eigen->e_val_im, mod->eigen->r_e_vect, mod->eigen->r_e_vect_im,
            mod->eigen->space))
  {
    for (i = 0; i < mod->ns; i++) mod->e_frq->pi->v[i] = mod_pi[i];
    Update_Qmat_GTR(mod->r_mat->rr->v, mod->r_mat->rr_val->v,
                    mod->r_mat->rr_num->v, mod->e_frq->pi->v,
                    mod->r_mat->qmat->v);
    Free(pi);
    Free(mod_pi);
    return -1.;
  }

  /* Get the left eigen vector of pi^{-1} x F */
  For(i, eigen_struct->size * eigen_struct->size) eigen_struct->l_e_vect[i] =
      eigen_struct->r_e_vect[i];
  if (!Matinv(eigen_struct->l_e_vect, eigen_struct->size, eigen_struct->size,
              YES))
  {
    for (i = 0; i < mod->ns; i++) mod->e_frq->pi->v[i] = mod_pi[i];
    Update_Qmat_GTR(mod->r_mat->rr->v, mod->r_mat->rr_val->v,
                    mod->r_mat->rr_num->v, mod->e_frq->pi->v,
                    mod->r_mat->qmat->v);
    Free(pi);
    Free(mod_pi);
    return -1.;
  }

  /* log of eigen values */
  for (i = 0; i < eigen_struct->size; i++)
  {
    /*       if(eigen_struct->e_val[i] < 0.0) eigen_struct->e_val[i] = 0.0001;
     */
    eigen_struct->e_val[i] = (phydbl)log(eigen_struct->e_val[i]);
  }

  /* Matrix multiplications log(pi^{-1} x F) */
  for (i = 0; i < eigen_struct->size; i++)
    for (j = 0; j < eigen_struct->size; j++)
      eigen_struct->r_e_vect[eigen_struct->size * i + j] =
          eigen_struct->r_e_vect[eigen_struct->size * i + j] *
          eigen_struct->e_val[j];
  for (i = 0; i < eigen_struct->size; i++)
    for (j = 0; j < eigen_struct->size; j++) F[eigen_struct->size * i + j] = .0;
  for (i = 0; i < eigen_struct->size; i++)
    for (j = 0; j < eigen_struct->size; j++)
      for (k = 0; k < eigen_struct->size; k++)
        F[eigen_struct->size * i + j] +=
            eigen_struct->r_e_vect[eigen_struct->size * i + k] *
            eigen_struct->l_e_vect[eigen_struct->size * k + j];

  /* Trace */
  dist = .0;
  for (i = 0; i < eigen_struct->size; i++)
    dist += F[eigen_struct->size * i + i];

  dist /= -4.;

  Free(pi);
  Free(mod_pi);
  Free(F_phydbl);

  if (isnan(dist)) return -1.;
  return dist;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl GTR_Dist(phydbl *F, phydbl alpha, eigen *eigen_struct)
{
  phydbl *pi;
  int     i, j, k;
  phydbl  dist;
  phydbl *F_phydbl;

  pi       = (phydbl *)mCalloc(eigen_struct->size, sizeof(phydbl));
  F_phydbl = (phydbl *)mCalloc(eigen_struct->size * eigen_struct->size,
                               sizeof(phydbl));

  /*   /\* Waddell and Steel's example *\/ */
  /*   F[4*0+0] = 1415./4898.; F[4*0+1] = 8./4898.;    F[4*0+2] = 55./4898.;
   * F[4*0+3] = 2./4898.; */
  /*   F[4*1+0] = 4./4898.;    F[4*1+1] = 1371./4898.; F[4*1+2] = 1./4898.;
   * F[4*1+3] = 144./4898.; */
  /*   F[4*2+0] = 73./4898.;   F[4*2+1] = 0./4898.;    F[4*2+2] = 578./4898.;
   * F[4*2+3] = 0./4898.; */
  /*   F[4*3+0] = 3./4898.;    F[4*3+1] = 117./4898.;  F[4*3+2] = 1./4898.;
   * F[4*3+3] = 1126./4898.; */

  for (i = 0; i < eigen_struct->size; i++)
  {
    for (j = 0; j < eigen_struct->size; j++)
    {
      pi[i] +=
          (F[eigen_struct->size * i + j] + F[eigen_struct->size * j + i]) / 2.;
    }
  }

  /* /\*   Jukes and Cantor correction *\/ */
  /*   sum = .0; */
  /*   for(i=0;i<eigen_struct->size;i++) sum += F[eigen_struct->size*i+i]; */
  /*   sum = 1.-sum; */
  /*   For(i,eigen_struct->size*eigen_struct->size) F[i] = sum/12.; */
  /*   for(i=0;i<eigen_struct->size;i++) F[eigen_struct->size*i+i] =
   * (1.-sum)/4.; */
  /*   for(i=0;i<eigen_struct->size;i++) pi[i] = 1./(phydbl)eigen_struct->size;
   */

  Make_Symmetric(&F, eigen_struct->size);
  Divide_Mat_By_Vect(&F, pi, eigen_struct->size);

  /* Eigen decomposition of pi^{-1} x F */
  for (i = 0; i < eigen_struct->size; i++)
    for (j = 0; j < eigen_struct->size; j++)
      F_phydbl[eigen_struct->size * i + j] = F[eigen_struct->size * i + j];
  if (Eigen(1, F_phydbl, eigen_struct->size, eigen_struct->e_val,
            eigen_struct->e_val_im, eigen_struct->r_e_vect,
            eigen_struct->r_e_vect_im, eigen_struct->space))
  {
    Free(pi);
    return -1.;
  }

  /* Get the left eigen vector of pi^{-1} x F */
  For(i, eigen_struct->size * eigen_struct->size) eigen_struct->l_e_vect[i] =
      eigen_struct->r_e_vect[i];
  if (!Matinv(eigen_struct->l_e_vect, eigen_struct->size, eigen_struct->size,
              YES))
  {
    Free(pi);
    return -1.;
  }

  /* Equation (3) + inverse of the moment generating function for the gamma
   * distribution (see Waddell & Steel, 1997) */
  for (i = 0; i < eigen_struct->size; i++)
  {
    if (eigen_struct->e_val[i] < 0.0)
    {
      eigen_struct->e_val[i] = 0.0001;
    }
    if (alpha < .0)
      eigen_struct->e_val[i] = (phydbl)log(eigen_struct->e_val[i]);
    else
      eigen_struct->e_val[i] =
          alpha * (1. - (phydbl)POW(eigen_struct->e_val[i], -1. / alpha));
  }

  /* Matrix multiplications pi x log(pi^{-1} x F) */
  for (i = 0; i < eigen_struct->size; i++)
    for (j = 0; j < eigen_struct->size; j++)
      eigen_struct->r_e_vect[eigen_struct->size * i + j] =
          eigen_struct->r_e_vect[eigen_struct->size * i + j] *
          eigen_struct->e_val[j];
  for (i = 0; i < eigen_struct->size; i++)
    for (j = 0; j < eigen_struct->size; j++) F[eigen_struct->size * i + j] = .0;
  for (i = 0; i < eigen_struct->size; i++)
    for (j = 0; j < eigen_struct->size; j++)
      for (k = 0; k < eigen_struct->size; k++)
        F[eigen_struct->size * i + j] +=
            eigen_struct->r_e_vect[eigen_struct->size * i + k] *
            eigen_struct->l_e_vect[eigen_struct->size * k + j];
  for (i = 0; i < eigen_struct->size; i++)
    for (j = 0; j < eigen_struct->size; j++)
      F[eigen_struct->size * i + j] *= pi[i];

  /* Trace */
  dist = .0;
  for (i = 0; i < eigen_struct->size; i++)
    dist -= F[eigen_struct->size * i + i];

  Free(pi);
  Free(F_phydbl);

  if (isnan(dist)) return -1.;
  return dist;
}
