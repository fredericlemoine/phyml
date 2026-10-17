/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Protein models with one replacement matrix per site rate category:

   S.Q. Le, C.C. Dang, O. Gascuel. Modeling protein evolution with several
   amino acid replacement matrices depending on site rates. Molecular
   Biology and Evolution 29(10):2921-2936, 2012.

   LG4M (equation 5): four equally weighted discrete Gamma categories, the
   k-th of which uses matrix M_k (M1 for the slowest rates, M4 for the
   fastest); one free parameter, alpha.

   LG4X (equation 6): four categories of free weights w_k and rates rho_k,
   with sum_k w_k = 1 and sum_k w_k rho_k = 1; category k uses matrix X_k.
   Six free parameters. Following the paper, the starting tree, rates and
   weights (0.25) are first estimated with LG4M.

   Each matrix has its own exchangeabilities and frequencies and is
   normalised to one substitution per unit time. The models have no +F
   and no +I variants. */

#ifndef LG4_H
#define LG4_H

#include "utilities.h"

#define LG4_N_CATG 4

/* Models (-m LG4M, -m LG4X) */
#define LG4_M 1
#define LG4_X 2

typedef struct __LG4_Mod
{
  int    model;                   /* LG4_M or LG4_X */
  eigen *eig[LG4_N_CATG];         /* eigen system of each normalised matrix */
  phydbl (*pi)[20]; /* equilibrium frequencies of each matrix, 32-byte
                       aligned for the SIMD kernels */
} t_lg4;

t_lg4      *LG4_Make(int model);
void        LG4_Free(t_lg4 *g);
t_lg4      *LG4_Copy(const t_lg4 *g);
void        LG4_PMat(phydbl l, const t_lg4 *g, int catg, int ns, int pos,
                     phydbl *Pij, phydbl *tPij);
void        LG4_Set_Rates(t_mod *mod, int model);
void        LG4_Switch(t_mod *mod, int model);
void        LG4_Print(FILE *fp, const t_tree *tree);
const char *LG4_Model_Name(int model);

#endif
