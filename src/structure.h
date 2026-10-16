/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Protein models accounting for solvent accessibility and secondary
   structure:

   S.Q. Le, O. Gascuel. Accounting for solvent accessibility and secondary
   structure in protein phylogenetics is clearly beneficial. Systematic
   Biology 59(3):277-287, 2010.

   Sites belong to the classes of a partition: EX (buried/exposed), EHO
   (extended/helix/other) or EX_EHO (their 6 combinations), each with its
   own replacement matrix and global rate. Four ways of combining them:
   MIX (mixture, classes unknown), PART (each site uses the matrix of its
   annotated class), CONF/LG and CONF/MIX (an annotation is trusted with
   probability chi; otherwise the site follows LG, or the mixture). A single
   discrete Gamma distribution of rates is used throughout.

   In the likelihood engine, each component (class matrix, or LG) times
   each Gamma category is one rate category: its transition matrix uses the
   component's eigen system, and the site likelihood weights the
   categories with site-specific weights. */

#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "utilities.h"

#define STRUCT_N_MATRICES 11
#define STRUCT_MAX_COMP 7 /* 6 classes + LG */

/* Partitions (-m) */
#define STRUCT_EX 1
#define STRUCT_EHO 2
#define STRUCT_EX_EHO 3

/* Models (--struct_model) */
#define STRUCT_MIX 1
#define STRUCT_PART 2
#define STRUCT_CONF_LG 3
#define STRUCT_CONF_MIX 4

typedef struct __Struct_Mod
{
  int partition; /* STRUCT_EX, STRUCT_EHO or STRUCT_EX_EHO */
  int model;     /* STRUCT_MIX, STRUCT_PART, STRUCT_CONF_LG or STRUCT_CONF_MIX */
  int n_class;   /* 2, 3 or 6 */
  int n_comp;    /* n_class, plus LG for CONF/LG */
  int n_gamma;   /* Gamma categories (-c) */

  phydbl  rate[STRUCT_MAX_COMP];   /* global rate of each component */
  phydbl  pi[STRUCT_MAX_COMP][20]; /* equilibrium frequencies */
  eigen  *eig[STRUCT_MAX_COMP];    /* eigen system of each normalised matrix */
  phydbl *gamma_rr;                /* the n_gamma Gamma rates (mean 1) */
  phydbl *gamma_p;

  phydbl    chi;          /* confidence in the annotations (CONF models) */
  short int opt_chi;      /* YES: chi is estimated */
  phydbl    props[6];     /* class proportions P_C (MIX, CONF/MIX) */
  phydbl    props_log[6]; /* unconstrained values, props = softmax */
  short int opt_props;    /* YES: proportions are estimated */

  int   n_pattern;
  short *site_class; /* class of each site pattern, -1 if unknown */
  int   *class_count; /* number of annotated sites per class */
} t_struct;

t_struct *STRUCT_Make(int partition, int model, int n_gamma);
void      STRUCT_Free(t_struct *s);
t_struct *STRUCT_Copy(const t_struct *s);
void      STRUCT_Set_Site_Classes(t_struct *s, const short *site_class,
                                  calign *cdata);
void      STRUCT_Set_Props(t_struct *s, const phydbl *props);
void      STRUCT_Update_Props(t_struct *s);
void      STRUCT_Update_RAS(t_mod *mod);
void      STRUCT_PMat(phydbl l, const t_struct *s, int catg, int ns, int pos,
                      phydbl *Pij, phydbl *tPij);
phydbl    STRUCT_Weight(const t_struct *s, int site, int catg);
const phydbl *STRUCT_Pi(const t_struct *s, int catg);
eigen    *STRUCT_Eigen(const t_struct *s, int catg);
short    *STRUCT_Read_File(const char *file, int partition, int n_site);
void      STRUCT_Optimize(t_tree *tree, int verbose);
phydbl    STRUCT_Mean_Rate(const t_tree *tree);
phydbl    STRUCT_Mean_Rate_Data(const t_struct *s, const calign *data);
int       STRUCT_Free_Params(const t_struct *s);
void      STRUCT_Print(FILE *fp, const t_tree *tree);
const char *STRUCT_Model_Name(int partition, int model);
const char *STRUCT_Class_Name(int partition, int class_num);

#endif
