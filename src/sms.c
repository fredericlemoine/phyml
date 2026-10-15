/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* SMS (Smart Model Selection), after the SMS 2.1.1 scripts of V. Lefort and
   J.-E. Longueville (https://gite.lirmm.fr/atgc/sms). The candidate models
   are fitted on a fixed tree (branch lengths and model parameters are
   optimised) and compared with AIC or BIC, in four steps:

   1. Fixed tree: BioNJ tree (or -u tree) with GTR+G+I+F (DNA) or LG
      (proteins), branch lengths and parameters optimised.
   2. Rate decoration of the Step 1 matrix: +G, +I, +G+I, +R (all with +F
      for DNA), or +G, +G+I, +R for proteins.
   3. Matrix with the Step 2 decoration. DNA: GTR, then TN93, HKY85 and K80
      while each simpler model is better. Proteins: the matrices ordered by
      the chi-square distance between their frequencies and those of the
      alignment, with and without +F; one of the two is only tested until
      the other has been better three times in a row. For proteins, a p-inv
      from Step 2 is kept fixed at this step.
   4. Final decoration of the Step 3 matrix: all the decorations above, with
      and without +F for proteins.

   A fit is computed once and reused by later steps. */

#include "sms.h"
#include "bionj.h"
#include "free.h"
#include "init.h"
#include "io.h"
#include "lk.h"
#include "make.h"
#include "models.h"
#include "optimiz.h"
#include "pars.h"
#include "spr.h"

#ifdef MPI
#include "mpi_boot.h"
#endif

/* Proteins: alignments with fewer residues than this are 'small' */
#define SMS_THRESHOLD 5000

#define SMS_N_AA_MODELS 23

/* Protein matrices tested at Step 3, in the order used to break ties */
static const struct
{
  int         model;
  const char *name;
  int (*init)(phydbl *daa, phydbl *pi);
} sms_aa_models[SMS_N_AA_MODELS] = {
    {LG, "LG", Init_Qmat_LG},
    {WAG, "WAG", Init_Qmat_WAG},
    {JTT, "JTT", Init_Qmat_JTT},
    {MTREV, "MtREV", Init_Qmat_MtREV},
    {DAYHOFF, "Dayhoff", Init_Qmat_Dayhoff},
    {DCMUT, "DCMut", Init_Qmat_DCMut},
    {RTREV, "RtREV", Init_Qmat_RtREV},
    {CPREV, "CpREV", Init_Qmat_CpREV},
    {VT, "VT", Init_Qmat_VT},
    {BLOSUM62, "Blosum62", Init_Qmat_Blosum62},
    {MTMAM, "MtMam", Init_Qmat_MtMam},
    {MTART, "MtArt", Init_Qmat_MtArt},
    {HIVW, "HIVw", Init_Qmat_HIVw},
    {HIVB, "HIVb", Init_Qmat_HIVb},
    {AB, "AB", Init_Qmat_AB},
    {FLU, "FLU", Init_Qmat_FLU},
    {MTZOA, "MtZoa", Init_Qmat_MtZoa},
    {QBIRD, "Q.bird", Init_Qmat_Q_Bird},
    {QINSECT, "Q.insect", Init_Qmat_Q_Insect},
    {QMAMMAL, "Q.mammal", Init_Qmat_Q_Mammal},
    {QPFAM, "Q.pfam", Init_Qmat_Q_Pfam},
    {QPLANT, "Q.plant", Init_Qmat_Q_Plant},
    {QYEAST, "Q.yeast", Init_Qmat_Q_Yeast},
};

typedef struct
{
  calign    *cdata;
  option    *io;
  int        criterion;
  int        n_catg;     /* categories of +G and +R */
  int        n_branches; /* 2n-3 */
  int        topo_model; /* matrix of Step 1 */
  char      *tree;       /* fixed tree, from Step 1 */
  t_sms_fit *fits;       /* all fits, in the order they were computed */
} t_sms;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static const char *SMS_Matrix_Name(int model)
{
  int i;

  switch (model)
  {
  case GTR: return "GTR";
  case TN93: return "TN93";
  case HKY85: return "HKY85";
  case K80: return "K80";
  }
  for (i = 0; i < SMS_N_AA_MODELS; i++)
    if (sms_aa_models[i].model == model) return sms_aa_models[i].name;
  return "?";
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Decoration as written by SMS: +R, +G, +I, +F (no +F for DNA, where the
   frequencies are always estimated). */

static void SMS_Deco_Name(int deco, int datatype, char *s)
{
  s[0] = '\0';
  if (deco & SMS_R) strcat(s, "+R");
  if (deco & SMS_G) strcat(s, "+G");
  if (deco & SMS_I) strcat(s, "+I");
  if ((deco & SMS_F) && datatype == AA) strcat(s, "+F");
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

char *SMS_Model_Name(int model, int deco, int datatype)
{
  char *s;

  s = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  strcpy(s, SMS_Matrix_Name(model));
  SMS_Deco_Name(deco, datatype, s + strlen(s));
  return s;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Number of free parameters, as counted by SMS: branch lengths, exchangeability
   and frequency parameters of the DNA matrices, rate decoration, and 19
   frequencies for proteins with +F. */

int SMS_Free_Params(int datatype, int model, int deco, int n_branches,
                    int n_catg)
{
  int k;

  k = n_branches;
  if (datatype == NT)
  {
    if (model == GTR) k += 8;
    else if (model == TN93) k += 5;
    else if (model == HKY85) k += 4;
    else if (model == K80) k += 1;
  }
  else if (deco & SMS_F)
    k += 19;
  if (deco & SMS_G) k += 1;
  if (deco & SMS_R) k += 2 * n_catg - 2;
  if (deco & SMS_I) k += 1;
  return k;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Chi-square distance between two amino-acid frequency vectors. Amino acids
   absent from both are ignored. */

phydbl SMS_Chi2_Distance(const phydbl *f1, const phydbl *f2)
{
  int    i;
  phydbl d, s;

  d = 0.0;
  for (i = 0; i < 20; i++)
  {
    s = f1[i] + f2[i];
    if (s > 0.0) d += (f1[i] - f2[i]) * (f1[i] - f2[i]) / s;
  }
  return sqrt(d);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Fill 'models' with the protein matrices sorted by increasing chi-square
   distance between their frequencies and the alignment's, and set
   'n_residues' to the number of the 20 amino acids in the alignment (gaps
   and ambiguities are not counted). Returns the number of matrices. */

int SMS_Sort_AA_Models(calign *cdata, int *models, int *n_residues)
{
  const char *aa = "ARNDCQEGHILKMFPSTWYV";
  phydbl      count[20], f_aln[20], f_mod[20], dist[SMS_N_AA_MODELS], total;
  phydbl     *daa;
  int         i, j, k, idx[SMS_N_AA_MODELS];
  char        c;

  for (k = 0; k < 20; k++) count[k] = 0.0;
  for (j = 0; j < cdata->n_pattern; j++)
    for (i = 0; i < cdata->n_otu; i++)
    {
      c = (char)toupper((unsigned char)cdata->c_seq[i]->state[j]);
      for (k = 0; k < 20; k++)
        if (aa[k] == c)
        {
          count[k] += cdata->wght[j];
          break;
        }
    }

  total = 0.0;
  for (k = 0; k < 20; k++) total += count[k];
  *n_residues = (int)(total + 0.5);
  for (k = 0; k < 20; k++) f_aln[k] = (total > 0.0) ? count[k] / total : 0.0;

  daa = (phydbl *)mCalloc(400, sizeof(phydbl));
  for (i = 0; i < SMS_N_AA_MODELS; i++)
  {
    sms_aa_models[i].init(daa, f_mod);
    dist[i] = SMS_Chi2_Distance(f_mod, f_aln);
  }
  Free(daa);

  // Stable sort: ties keep the order of sms_aa_models
  for (i = 0; i < SMS_N_AA_MODELS; i++)
  {
    for (j = i; j > 0 && dist[idx[j - 1]] > dist[i]; j--) idx[j] = idx[j - 1];
    idx[j] = i;
  }
  for (i = 0; i < SMS_N_AA_MODELS; i++) models[i] = sms_aa_models[idx[i]].model;

  return SMS_N_AA_MODELS;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Set a model as the PhyML command line would with -m, -f, -c, -a, -v and
   --freerates. A p-inv >= 0 is fixed, a negative one is estimated. */

static void SMS_Configure_Model(t_mod *mod, option *io, int model, int deco,
                                int n_catg, phydbl fixed_pinv)
{
  mod->ns         = (io->datatype == AA) ? 20 : 4;
  mod->whichmodel = model;
  Set_Model_Name(mod);

  if (mod->r_mat) Free_Rmat(mod->r_mat);
  mod->r_mat = (t_rmat *)Make_Rmat(mod->ns);
  Init_Rmat(mod->r_mat);
  Make_Custom_Model(mod);
  if (model == GTR) mod->r_mat->optimize = YES;

  if (mod->e_frq) Free_Efrq(mod->e_frq);
  mod->e_frq = (t_efrq *)Make_Efrq(mod->ns);
  Init_Efrq(NULL, mod->e_frq);
  if (io->datatype == NT || (deco & SMS_F))
  {
    mod->e_frq->type                  = ML;
    mod->e_frq->pi->optimize          = YES;
    mod->e_frq->pi_unscaled->optimize = YES;
  }
  else
  {
    mod->e_frq->type                  = MODEL;
    mod->e_frq->pi->optimize          = NO;
    mod->e_frq->pi_unscaled->optimize = NO;
  }

  mod->kappa->optimize =
      (model == K80 || model == HKY85 || model == TN93) ? YES : NO;
  mod->lambda->optimize = YES;

  // Make_RAS_Complete only allocates the rate arrays when they are NULL:
  // free them, as the number of categories may change (next data set)
  if (mod->ras->gamma_r_proba->v)
  {
    Free(mod->ras->gamma_r_proba->v);
    Free(mod->ras->gamma_r_proba_unscaled->v);
    Free(mod->ras->gamma_rr->v);
    Free(mod->ras->gamma_rr_unscaled->v);
    Free(mod->ras->skip_rate_cat);
    mod->ras->gamma_r_proba->v          = NULL;
    mod->ras->gamma_r_proba_unscaled->v = NULL;
    mod->ras->gamma_rr->v               = NULL;
    mod->ras->gamma_rr_unscaled->v      = NULL;
    mod->ras->skip_rate_cat             = NULL;
  }
  mod->ras->n_catg          = (deco & (SMS_G | SMS_R)) ? n_catg : 1;
  mod->ras->alpha->optimize = (deco & SMS_G) ? YES : NO;
  mod->ras->free_mixt_rates = (deco & SMS_R) ? YES : NO;
  if (deco & SMS_I)
  {
    mod->ras->invar = YES;
    if (fixed_pinv >= 0.0)
    {
      mod->ras->pinvar->v        = fixed_pinv;
      mod->ras->pinvar->optimize = NO;
    }
    else
      mod->ras->pinvar->optimize = YES;
  }
  else
  {
    mod->ras->invar            = NO;
    mod->ras->pinvar->v        = 0.0;
    mod->ras->pinvar->optimize = NO;
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Fit a model on the fixed tree (on the starting tree at Step 1), as
   'phyml -o lr' would. */

static t_sms_fit *SMS_Fit(t_sms *sms, int model, int deco, phydbl fixed_pinv)
{
  option    *io = sms->io;
  t_mod     *mod, *io_mod;
  t_tree    *tree;
  t_sms_fit *fit;
  short int  quiet;
  int        i;

  // The optimiser draws the order of some parameters at random: each fit
  // starts from the seed, as a separate 'phyml --r_seed' run would (and
  // all MPI processes get the same fits)
  srand(io->r_seed);

  mod = Make_Model_Basic();
  Set_Defaults_Model(mod);
  mod->s_opt = Make_Optimiz();
  Set_Defaults_Optimiz(mod->s_opt);
  mod->s_opt->opt_topo = NO;
  mod->io              = io;
  SMS_Configure_Model(mod, io, model, deco, sms->n_catg, fixed_pinv);
  Make_Model_Complete(mod);

  // Some routines read io->mod rather than tree->mod
  io_mod    = io->mod;
  quiet     = io->quiet;
  io->mod   = mod;
  io->quiet = YES;

  Init_Model(sms->cdata, mod, io);
  Set_Model_Parameters(mod);

  if (sms->tree)
  {
    char *s = sms->tree;
    tree    = Read_Tree(&s);
  }
  else if (io->in_tree == 2)
  {
    long pos = ftell(io->fp_in_tree);
    tree     = Read_User_Tree(sms->cdata, mod, io);
    fseek(io->fp_in_tree, pos, SEEK_SET);
  }
  else
    tree = Dist_And_BioNJ(sms->cdata, mod, io);

  tree->mod          = mod;
  tree->io           = io;
  tree->data         = sms->cdata;
  tree->n_root       = NULL;
  tree->e_root       = NULL;
  tree->n_tot_bl_opt = 0;

  Set_Both_Sides(YES, tree);
  Connect_CSeqs_To_Nodes(tree->data, io, tree);
  Make_Tree_For_Pars(tree);
  Make_Tree_For_Lk(tree);
  Make_Spr(tree);
  Br_Len_Not_Involving_Invar(tree);
  Unscale_Br_Len_Multiplier_Tree(tree);

  Set_Update_Eigen(YES, mod);
  Lk(NULL, tree);
  Set_Update_Eigen(NO, mod);

  Round_Optimize(tree, ROUND_MAX);

  Set_Both_Sides(YES, tree);
  Lk(NULL, tree);

  fit             = (t_sms_fit *)mCalloc(1, sizeof(t_sms_fit));
  fit->model      = model;
  fit->deco       = deco;
  fit->fixed_pinv = ((deco & SMS_I) && fixed_pinv >= 0.0) ? YES : NO;
  fit->lnL        = tree->c_lnL;
  fit->K   = SMS_Free_Params(io->datatype, model, deco, sms->n_branches,
                             sms->n_catg);
  fit->aic = 2. * fit->K - 2. * fit->lnL;
  fit->bic = fit->K * log(sms->cdata->init_len) - 2. * fit->lnL;

  fit->ns     = mod->ns;
  fit->n_catg = mod->ras->n_catg;
  fit->alpha  = mod->ras->alpha->v;
  fit->pinv   = (mod->ras->invar == YES) ? mod->ras->pinvar->v : 0.0;
  fit->kappa  = mod->kappa->v;
  fit->lambda = mod->lambda->v;
  if (mod->ns == 4)
    for (i = 0; i < 6; i++)
    {
      fit->rr_val[i] = mod->r_mat->rr_val->v[i];
      fit->rr[i]     = mod->r_mat->rr->v[i];
    }
  for (i = 0; i < mod->ns; i++)
  {
    fit->pi[i]          = mod->e_frq->pi->v[i];
    fit->pi_unscaled[i] = mod->e_frq->pi_unscaled->v[i];
  }
  fit->gamma_rr_unscaled = (phydbl *)mCalloc(fit->n_catg, sizeof(phydbl));
  fit->gamma_r_proba_unscaled = (phydbl *)mCalloc(fit->n_catg, sizeof(phydbl));
  for (i = 0; i < fit->n_catg; i++)
  {
    fit->gamma_rr_unscaled[i]      = mod->ras->gamma_rr_unscaled->v[i];
    fit->gamma_r_proba_unscaled[i] = mod->ras->gamma_r_proba_unscaled->v[i];
  }

  if (sms->tree == NULL) // Step 1: the tree is kept for the next fits
  {
    Check_Br_Lens(tree);
    Br_Len_Involving_Invar(tree);
    Rescale_Br_Len_Multiplier_Tree(tree);
    // Written as main() writes the tree file (rooted, then read and written
    // again): the next fits then visit the branches in the same order as
    // in the SMS scripts, which read that file
    Get_Best_Root_Position(tree);
    {
      char   *s, *p;
      t_tree *dum;
      s         = Write_Tree(tree);
      p         = s;
      dum       = Read_Tree(&p);
      dum->data = tree->data;
      dum->mod  = mod;
      dum->io   = io;
      Connect_CSeqs_To_Nodes(tree->data, io, dum);
      Insert_Duplicates(dum);
      fit->tree = Write_Tree(dum);
      Free_Tree(dum);
      Free(s);
    }
  }

  Free_Best_Spr(tree);
  Free_Spr_List_One_Edge(tree);
  Free_Spr_List_All_Edge(tree);
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Tree(tree);
  Free_Model_Complete(mod);
  Free_Optimiz(mod->s_opt);
  Free_Model_Basic(mod);

  io->mod   = io_mod;
  io->quiet = quiet;

  return fit;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static phydbl SMS_Crit(const t_sms *sms, const t_sms_fit *fit)
{
  return (sms->criterion == SMS_AIC) ? fit->aic : fit->bic;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static const char *SMS_Crit_Name(int criterion)
{
  return (criterion == SMS_AIC) ? "AIC" : "BIC";
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static t_sms_fit *SMS_Find(t_sms *sms, int model, int deco,
                           short int fixed_pinv)
{
  t_sms_fit *fit;

  for (fit = sms->fits; fit; fit = fit->next)
    if (fit->model == model && fit->deco == deco &&
        fit->fixed_pinv == fixed_pinv)
      return fit;
  return NULL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* The fit of a model, computed if needed. A p-inv >= 0 is fixed (for +I). */

static t_sms_fit *SMS_Get(t_sms *sms, int model, int deco, phydbl fixed_pinv)
{
  t_sms_fit *fit, **last;
  short int  fixed;
  char      *name;

  fixed = ((deco & SMS_I) && fixed_pinv >= 0.0) ? YES : NO;

  // The Step 1 matrix was fitted with an estimated p-inv: reuse that fit
  if (fixed == YES && model == sms->topo_model &&
      (fit = SMS_Find(sms, model, deco, NO)))
    return fit;
  if ((fit = SMS_Find(sms, model, deco, fixed))) return fit;

  fit = SMS_Fit(sms, model, deco, fixed ? fixed_pinv : -1.0);
  for (last = &sms->fits; *last; last = &(*last)->next)
    ;
  *last = fit;

  if (sms->io->quiet == NO)
  {
    name = SMS_Model_Name(model, deco, sms->io->datatype);
    PhyML_Printf("\n\t%-18s %-7s lnL = %14.4f   K = %4d   %s = %14.4f", name,
                 fit->fixed_pinv ? "(p-inv)" : "", fit->lnL, fit->K,
                 SMS_Crit_Name(sms->criterion), SMS_Crit(sms, fit));
    Free(name);
  }

  return fit;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Best of several decorations of one matrix (the first one wins ties). */

static t_sms_fit *SMS_Best_Deco(t_sms *sms, int model, const int *decos,
                                int n_decos)
{
  t_sms_fit *fit, *best;
  int        i;

  best = NULL;
  for (i = 0; i < n_decos; i++)
  {
    fit = SMS_Get(sms, model, decos[i], -1.0);
    if (best == NULL || SMS_Crit(sms, fit) < SMS_Crit(sms, best)) best = fit;
  }
  return best;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static t_sms_fit *SMS_Best_Matrix_NT(t_sms *sms, int deco)
{
  t_sms_fit *gtr, *tn93, *hky, *k80;

  gtr  = SMS_Get(sms, GTR, deco, -1.0);
  tn93 = SMS_Get(sms, TN93, deco, -1.0);
  if (SMS_Crit(sms, tn93) >= SMS_Crit(sms, gtr)) return gtr;
  hky = SMS_Get(sms, HKY85, deco, -1.0);
  if (SMS_Crit(sms, hky) >= SMS_Crit(sms, tn93)) return tn93;
  k80 = SMS_Get(sms, K80, deco, -1.0);
  if (SMS_Crit(sms, k80) >= SMS_Crit(sms, hky)) return hky;
  return k80;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static t_sms_fit *SMS_Best_Matrix_AA(t_sms *sms, int deco, phydbl pinv)
{
  t_sms_fit *fit, *alt, *best;
  int        models[SMS_N_AA_MODELS], n_models, n_residues, i, m;
  int        deco_all, deco_alt, all_is_better;

  n_models = SMS_Sort_AA_Models(sms->cdata, models, &n_residues);

  if (n_residues < SMS_THRESHOLD)
  {
    // Small alignment: without +F first, the farthest matrices first
    deco_all = deco;
    deco_alt = deco | SMS_F;
  }
  else
  {
    deco_all = deco | SMS_F;
    deco_alt = deco;
  }

  if (sms->io->quiet == NO)
    PhyML_Printf("\n\t(%d residues: %s alignment)", n_residues,
                 n_residues < SMS_THRESHOLD ? "small" : "large");

  best          = NULL;
  all_is_better = 0;
  for (i = 0; i < n_models; i++)
  {
    m   = (n_residues < SMS_THRESHOLD) ? models[n_models - 1 - i] : models[i];
    fit = SMS_Get(sms, m, deco_all, pinv);
    if (best == NULL || SMS_Crit(sms, fit) < SMS_Crit(sms, best)) best = fit;

    // The other decoration is no longer tested once deco_all has been
    // better three times in a row
    if (all_is_better < 3)
    {
      alt = SMS_Get(sms, m, deco_alt, pinv);
      if (SMS_Crit(sms, fit) < SMS_Crit(sms, alt))
        all_is_better++;
      else
        all_is_better = 0;
      if (SMS_Crit(sms, alt) < SMS_Crit(sms, best)) best = alt;
    }
  }
  return best;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Write all fits, sorted by criterion, except the fits of the selected
   matrix with a fixed p-inv (Step 4 estimated it). */

static void SMS_Write_Fits(t_sms *sms, const t_sms_fit *best)
{
  option     *io = sms->io;
  t_sms_fit **fits, *fit, *tmp;
  FILE       *fp;
  char       *file, deco[20];
  int         n, i, j;

#ifdef MPI
  if (Global_myRank != 0) return;
#endif

  n = 0;
  for (fit = sms->fits; fit; fit = fit->next) n++;
  fits = (t_sms_fit **)mCalloc(n, sizeof(t_sms_fit *));
  n    = 0;
  for (fit = sms->fits; fit; fit = fit->next)
  {
    if (fit->fixed_pinv == YES && fit->model == best->model) continue;
    for (j = n; j > 0 && SMS_Crit(sms, fits[j - 1]) > SMS_Crit(sms, fit); j--)
      fits[j] = fits[j - 1];
    fits[j] = fit;
    n++;
  }

  file = (char *)mCalloc(strlen(io->in_align_file) + strlen(io->run_id_string) +
                             20,
                         sizeof(char));
  strcpy(file, io->in_align_file);
  strcat(file, "_phyml_sms");
  if (io->append_run_ID)
  {
    strcat(file, "_");
    strcat(file, io->run_id_string);
  }
  strcat(file, ".csv");

  fp = fopen(file, "w");
  if (fp == NULL)
  {
    PhyML_Fprintf(stderr, "\n. Cannot write the SMS results to '%s'.\n", file);
    Exit("\n");
  }
  fprintf(fp, "Model,Decoration,K,lnL,AIC,BIC,FixedPinv\n");
  for (i = 0; i < n; i++)
  {
    tmp = fits[i];
    SMS_Deco_Name(tmp->deco, io->datatype, deco);
    fprintf(fp, "%s,%s,%d,%.6f,%.6f,%.6f,%s\n", SMS_Matrix_Name(tmp->model),
            deco, tmp->K, tmp->lnL, tmp->aic, tmp->bic,
            tmp->fixed_pinv ? "yes" : "no");
  }
  fclose(fp);

  if (io->quiet == NO)
    PhyML_Printf("\n\n. The %d fits are in '%s'.", n, Basename(file));

  Free(fits);
  Free(file);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void SMS_Free_Fits(t_sms_fit *fit)
{
  t_sms_fit *next;

  while (fit)
  {
    next = fit->next;
    if (fit->tree) Free(fit->tree);
    Free(fit->gamma_rr_unscaled);
    Free(fit->gamma_r_proba_unscaled);
    Free(fit);
    fit = next;
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

static void SMS_Print_Step(t_sms *sms, int step, const char *what,
                           const t_sms_fit *fit)
{
  char *name;

  if (sms->io->quiet == YES) return;
  name = SMS_Model_Name(fit->model, fit->deco, sms->io->datatype);
  PhyML_Printf("\n  Step %d: %s %s (%s = %.4f)\n", step, what, name,
               SMS_Crit_Name(sms->criterion), SMS_Crit(sms, fit));
  Free(name);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Select the model, write the fits to *_phyml_sms.csv, set io->mod to the
   selected model and keep its estimates in io->sms_fit. */

void SMS_Select_Model(calign *cdata, option *io)
{
  static const int decos_nt[] = {SMS_F, SMS_G | SMS_F, SMS_I | SMS_F,
                                 SMS_G | SMS_I | SMS_F, SMS_R | SMS_F};
  static const int decos_aa_step2[] = {SMS_G, SMS_G | SMS_I, SMS_R};
  static const int decos_aa_step4[] = {0,
                                       SMS_F,
                                       SMS_G,
                                       SMS_G | SMS_F,
                                       SMS_G | SMS_I,
                                       SMS_G | SMS_I | SMS_F,
                                       SMS_R,
                                       SMS_R | SMS_F};
  t_sms      sms;
  t_sms_fit *fit, *best, **prev;
  phydbl     pinv;
  char      *name;

  sms.cdata      = cdata;
  sms.io         = io;
  sms.criterion  = io->sms_criterion;
  sms.n_catg     = io->sms_n_catg;
  sms.n_branches = 2 * cdata->n_otu - 3;
  sms.topo_model = (io->datatype == NT) ? GTR : LG;
  sms.tree       = NULL;
  sms.fits       = NULL;

  if (io->quiet == NO)
    PhyML_Printf("\n\n. Selecting the substitution model with SMS (%s)...\n",
                 SMS_Crit_Name(sms.criterion));

  // Step 1: fixed tree
  fit = SMS_Get(&sms, sms.topo_model,
                (io->datatype == NT) ? (SMS_G | SMS_I | SMS_F) : 0, -1.0);
  sms.tree = fit->tree;
  SMS_Print_Step(&sms, 1,
                 (io->in_tree == 2) ? "tree from the -u tree with"
                                    : "tree from BioNJ with",
                 fit);

  // Step 2: rate decoration
  if (io->datatype == NT)
    best = SMS_Best_Deco(&sms, sms.topo_model, decos_nt, 5);
  else
    best = SMS_Best_Deco(&sms, sms.topo_model, decos_aa_step2, 3);
  pinv = (best->deco & SMS_I) ? best->pinv : -1.0;
  SMS_Print_Step(&sms, 2, "decoration", best);

  // Step 3: matrix
  if (io->datatype == NT)
    best = SMS_Best_Matrix_NT(&sms, best->deco);
  else
    best = SMS_Best_Matrix_AA(&sms, best->deco, pinv);
  SMS_Print_Step(&sms, 3, "matrix", best);

  // Step 4: final decoration
  if (io->datatype == NT)
    best = SMS_Best_Deco(&sms, best->model, decos_nt, 5);
  else
    best = SMS_Best_Deco(&sms, best->model, decos_aa_step4, 8);
  SMS_Print_Step(&sms, 4, "final decoration", best);

  SMS_Write_Fits(&sms, best);

  if (io->quiet == NO)
  {
    name = SMS_Model_Name(best->model, best->deco, io->datatype);
    PhyML_Printf("\n\n. Selected model: %s (%s = %.4f, lnL = %.4f).\n", name,
                 SMS_Crit_Name(sms.criterion), SMS_Crit(&sms, best), best->lnL);
    Free(name);
  }

  // The analysis uses the selected model, starting from its estimates
  SMS_Configure_Model(io->mod, io, best->model, best->deco, sms.n_catg, -1.0);
  io->mod->ras->alpha->v  = best->alpha; // shown by Print_Settings if fixed
  io->mod->ras->pinvar->v = best->pinv;
  io->mod->kappa->v       = best->kappa;
  io->mod->lambda->v      = best->lambda;
  if (io->mod->s_opt->opt_subst_param == NO)
  {
    io->mod->ras->alpha->optimize  = NO;
    io->mod->kappa->optimize       = NO;
    io->mod->lambda->optimize      = NO;
    io->mod->ras->pinvar->optimize = NO;
    io->mod->r_mat->optimize       = NO;
  }

  // The analysis draws the same random numbers as without --sms
#ifdef MPI
  srand(io->r_seed + Global_myRank);
#else
  srand(io->r_seed);
#endif

  for (prev = &sms.fits; *prev != best; prev = &(*prev)->next)
    ;
  *prev      = best->next;
  best->next = NULL;
  if (best->tree)
  {
    Free(best->tree);
    best->tree = NULL;
  }
  SMS_Free_Fits(sms.fits);
  SMS_Free_Fits(io->sms_fit);
  io->sms_fit = best;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Start (or fix, if not optimised) the parameters of the selected model at
   the SMS estimates. Call after Init_Model, before Set_Model_Parameters. */

void SMS_Set_Parameters(t_mod *mod, const t_sms_fit *fit)
{
  int i;

  if (fit == NULL || mod->ns != fit->ns || mod->ras->n_catg != fit->n_catg)
    return;

  mod->kappa->v      = fit->kappa;
  mod->lambda->v     = fit->lambda;
  mod->ras->alpha->v = fit->alpha;
  if (mod->ras->invar == YES) mod->ras->pinvar->v = fit->pinv;
  if (mod->whichmodel == GTR)
    for (i = 0; i < 6; i++)
    {
      mod->r_mat->rr_val->v[i] = fit->rr_val[i];
      mod->r_mat->rr->v[i]     = fit->rr[i];
    }
  if (mod->e_frq->type == ML)
    for (i = 0; i < mod->ns; i++)
    {
      mod->e_frq->pi->v[i]          = fit->pi[i];
      mod->e_frq->pi_unscaled->v[i] = fit->pi_unscaled[i];
    }
  if (mod->ras->free_mixt_rates == YES)
    for (i = 0; i < mod->ras->n_catg; i++)
    {
      mod->ras->gamma_rr_unscaled->v[i]      = fit->gamma_rr_unscaled[i];
      mod->ras->gamma_r_proba_unscaled->v[i] = fit->gamma_r_proba_unscaled[i];
    }
}
