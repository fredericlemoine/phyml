/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* SMS (Smart Model Selection): selects the substitution model before the
   analysis (option --sms).

   V. Lefort, J.-E. Longueville, O. Gascuel. SMS: Smart Model Selection in
   PhyML. Molecular Biology and Evolution 34(9):2422-2424, 2017. */

#ifndef SMS_H
#define SMS_H

#include "utilities.h"

/* Rate and frequency decorations of a model, combined as bit flags */
#define SMS_G 1 /* discrete Gamma */
#define SMS_I 2 /* invariable sites */
#define SMS_R 4 /* FreeRate */
#define SMS_F 8 /* ML equilibrium frequencies */

typedef struct __Sms_Fit
{
  int       model;      /* whichmodel */
  int       deco;       /* SMS_G | SMS_I | SMS_R | SMS_F */
  short int fixed_pinv; /* YES: p-inv was fixed, not estimated */
  int       K;          /* number of free parameters */
  phydbl    lnL, aic, bic;

  /* Estimates, used as starting values of the final analysis */
  int     ns, n_catg;
  phydbl  alpha, pinv, kappa, lambda;
  phydbl  rr_val[6], rr[6];
  phydbl  pi[20], pi_unscaled[20];
  phydbl *gamma_rr_unscaled, *gamma_r_proba_unscaled;

  char *tree; /* tree after the fit (Step 1 only) */

  struct __Sms_Fit *next;
} t_sms_fit;

void SMS_Select_Model(calign *cdata, option *io);
void SMS_Set_Parameters(t_mod *mod, const t_sms_fit *fit);
void SMS_Free_Fits(t_sms_fit *fit);
char *SMS_Model_Name(int model, int deco, int datatype);
int SMS_Free_Params(int datatype, int model, int deco, int n_branches,
                    int n_catg);
phydbl SMS_Chi2_Distance(const phydbl *f1, const phydbl *f2);
int SMS_Sort_AA_Models(calign *cdata, int *models, int *n_residues);

#endif
