/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Entry points of libphyml (see phyml.h): the PhyML analysis, which was
   main() until the library was split from the program. */

#include <setjmp.h>
#include <unistd.h>

#include "alrt.h"
#include "bionj.h"
#include "eigen.h"
#include "error.h"
#include "free.h"
#include "help.h"
#include "invitee.h"
#include "lk.h"
#include "mixt.h"
#include "models.h"
#include "optimiz.h"
#include "pars.h"
#include "simu.h"
#include "sms.h"
#include "structure.h"
#include "lg4.h"
#include "spr.h"
#include "utilities.h"
#include "phyml.h"

#ifdef MPI
#include "mpi_boot.h"
#endif

#ifdef BEAGLE
#include "beagle_utils.h"
#endif

#if (defined PHYML)

/* Structure model (-m EX, EHO, EX_EHO): one rate category per component
   and Gamma category. */

static void Setup_Structure_Model(calign *cdata, option *io)
{
  t_struct *s;

  s = STRUCT_Make(io->struct_partition, io->struct_model, io->struct_n_gamma);
  if (io->struct_site_class && io->struct_model != STRUCT_MIX)
    STRUCT_Set_Site_Classes(s, io->struct_site_class, cdata);

  if (io->struct_chi >= 0.0)
  {
    s->chi     = io->struct_chi;
    s->opt_chi = NO;
  }
  if (io->struct_n_props > 0)
  {
    STRUCT_Set_Props(s, io->struct_props);
    s->opt_props = NO;
  }
  if (io->struct_fix_props == YES) s->opt_props = NO;
  if (io->mod->s_opt->opt_subst_param == NO) s->opt_chi = s->opt_props = NO;

  io->mod->smod        = s;
  io->mod->ras->n_catg = s->n_comp * s->n_gamma;
}

/* Tree search, or optimisation on the fixed topology. */

static void Search_Or_Optimize(t_tree *tree)
{
  if (tree->mod->s_opt->opt_topo)
  {
    Global_Spr_Search(tree);
    if (tree->n_root) Add_Root(tree->a_edges[0], tree);
  }
  else
  {
#ifdef BEAGLE
    tree->b_inst = create_beagle_instance(tree, tree->io->quiet, tree->io);
#endif
    if (tree->mod->s_opt->opt_subst_param ||
        tree->mod->s_opt->opt_bl_one_by_one)
      Round_Optimize(tree, ROUND_MAX);
  }
}

// Truncates the tree file to pos and writes from there (also in --append
// mode, where every write goes to the end of the file)
static void Truncate_Tree_File(FILE *fp, long pos)
{
  fflush(fp);
  if (ftruncate(fileno(fp), (off_t)pos) != 0 || fseek(fp, pos, SEEK_SET) != 0)
  {
    PhyML_Fprintf(stderr, "\n. Err.: could not rewrite the tree file.");
    Exit("\n");
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int PhyML_Main(int argc, char **argv)
{
  calign *cdata;
  option *io;
  t_tree *tree;
  int     num_data_set;
  int     num_tree, num_rand_tree;
  t_mod  *mod;
  time_t  t_beg, t_end;
  phydbl  best_lnL;
  int     r_seed;
  char   *most_likely_tree = NULL;
  int     orig_random_input_tree, orig_n_rand_starts;
  long    tree_file_start = 0, tree_pos = 0;

#ifdef MPI
  int rc;
  rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS)
  {
    PhyML_Fprintf(stderr, "\n. Err. starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  if (MPI_Comm_size(MPI_COMM_WORLD, &Global_numTask) != MPI_SUCCESS)
    MPI_Abort(MPI_COMM_WORLD, rc);
  if (MPI_Comm_rank(MPI_COMM_WORLD, &Global_myRank) != MPI_SUCCESS)
    MPI_Abort(MPI_COMM_WORLD, rc);
  PhyML_Fprintf(stdout, "\n\n. Running the analysis on %d CPU%s.",
                Global_numTask, Global_numTask > 1 ? "s." : ".");
#endif

#ifdef QUIET
  setvbuf(stdout, NULL, _IOFBF, 2048);
#endif

  tree     = NULL;
  mod      = NULL;
  best_lnL = UNLIKELY;

  // Errors are given the code of what the program is doing (error.h)
  PhyML_Error_Scope(PHYML_ERR_USAGE);
  io = (option *)Get_Input(argc, argv);
  PhyML_Error_Scope(PHYML_ERR);
  if (!io)
    return (0);
  else if (io->use_xml == YES)
  {
    Free(io);
    return (0);
  }

  r_seed = (io->r_seed < 0) ? (time(NULL)) : (io->r_seed);
#ifdef MPI
  srand(r_seed + Global_myRank);
#else
  srand(r_seed);
#endif
  io->r_seed = r_seed;

  if (io->in_tree == 2)
    Test_Multiple_Data_Set_Format(io);
  else
    io->n_trees = 1;

  if (io->n_trees == 0 && io->in_tree == 2)
  {
    PhyML_Printf("\n== Err.: the input tree file does not provide a tree in "
                 "valid format.");
    Exit("\n");
  }

  if ((io->n_data_sets > 1) && (io->n_trees > 1))
  {
    io->n_data_sets = MIN(io->n_trees, io->n_data_sets);
    io->n_trees     = MIN(io->n_trees, io->n_data_sets);
  }

  // The last random start of each data set switches random_input_tree off
  // and adds one start from BioNJ: each data set starts from these values
  orig_random_input_tree = io->mod->s_opt->random_input_tree;
  orig_n_rand_starts     = io->mod->s_opt->n_rand_starts;

  // Where this run's trees start in the tree file (its end with --append)
  if (io->fp_out_tree && fseek(io->fp_out_tree, 0, SEEK_END) == 0)
    tree_file_start = ftell(io->fp_out_tree);

  for (num_data_set = 0; num_data_set < io->n_data_sets; num_data_set++)
  {
    best_lnL = UNLIKELY;
    PhyML_Error_Scope(PHYML_ERR_INPUT);
    Get_Seq(io);
    PhyML_Error_Scope(PHYML_ERR);

    cdata = NULL;
    if (io->data)
    {
      if (io->n_data_sets > 1)
        PhyML_Printf("\n. Data set [#%d]\n", num_data_set + 1);
      // Structure annotations: a site pattern must not mix classes, so the
      // sites are not compressed
      if (io->struct_partition && io->struct_file &&
          io->struct_model != STRUCT_MIX)
      {
        if (io->struct_site_class) Free(io->struct_site_class);
        io->struct_site_class = STRUCT_Read_File(
            io->struct_file, io->struct_partition, io->data[0]->len);
        io->colalias = NO;
      }

      cdata = Compact_Data(io->data, io);
      Free_Seq(io->data, cdata->n_otu);

      // Sets io->mod to the selected model
      if (io->sms == YES) SMS_Select_Model(cdata, io);

      if (io->struct_partition) Setup_Structure_Model(cdata, io);
      if (io->lg4_model) io->mod->lg4 = LG4_Make(io->lg4_model);
    }

    Make_Model_Complete(io->mod);
    Set_Model_Name(io->mod);
    Print_Settings(io);
    mod                               = io->mod;
    io->mod->s_opt->random_input_tree = orig_random_input_tree;
    io->mod->s_opt->n_rand_starts     = orig_n_rand_starts;

    if (cdata && io->sms_only == YES)
    {
      Free_Calign(cdata);
    }
    else if (cdata)
    {
      for (num_tree = (io->n_trees == 1) ? (0) : (num_data_set);
           num_tree < io->n_trees; num_tree++)
      {
        if (io->mod->s_opt->random_input_tree == NO)
          io->mod->s_opt->n_rand_starts = 1;

        if (orig_random_input_tree == YES && io->n_trees > 1)
        {
          PhyML_Printf("\n== Cannot combine random starting trees with "
                       "multiple input trees.");
          Exit("\n");
        }

        for (num_rand_tree = 0; num_rand_tree < io->mod->s_opt->n_rand_starts;
             num_rand_tree++)
        {
          if ((io->mod->s_opt->random_input_tree) &&
              (io->mod->s_opt->topo_search != NNI_MOVE))
            if (!io->quiet)
              PhyML_Printf("\n\n. [Random start %3d/%3d]", num_rand_tree + 1,
                           io->mod->s_opt->n_rand_starts);

          Init_Model(cdata, mod, io);
          SMS_Set_Parameters(mod, io->sms_fit);
          if (mod->lg4) LG4_Set_Rates(mod, mod->lg4->model);
          Set_Model_Parameters(mod);

#ifdef M4
          if (io->mod->use_m4mod) M4_Init_Model(mod->m4mod, cdata, mod);
#endif

          // Make the initial tree
          switch (io->in_tree)
          {
          case 0:
          case 1:
          {
            tree = Dist_And_BioNJ(cdata, mod, io);
            if (io->print_mat_and_exit == YES)
            {
              Print_Mat(ML_Dist(cdata, mod));
              PhyML_Stop(-1); // asked to stop here: not an error
            }
            break;
          }
          case 2:
          {
            PhyML_Error_Scope(PHYML_ERR_INPUT);
            tree = Read_User_Tree(cdata, mod, io);
            PhyML_Error_Scope(PHYML_ERR);
            break;
          }
          }

          if (io->mod->s_opt->opt_topo == YES)
            Remove_Duplicates(cdata, io, tree);

          if (io->fp_in_constraint_tree != NULL)
          {
            char *s;

            PhyML_Printf("\n. Reading constraint tree file...");

            io->cstr_tree = Read_Tree_File_Phylip(io->fp_in_constraint_tree);

            if (io->cstr_tree == NULL)
              Exit("\n. Constraint tree not found...\n");

            if (io->cstr_tree->n_root != NULL)
            {
              PhyML_Printf("\n== The constraint tree file must be unrooted");
              Exit("\n");
            }

            s = Add_Taxa_To_Constraint_Tree(io->fp_in_constraint_tree, cdata);
            fflush(NULL);
            Free_Tree(tree);
            tree        = Read_Tree(&s);
            io->in_tree = 2;
            Free(s);
            Check_Constraint_Tree_Taxa_Names(io->cstr_tree, cdata);
            Alloc_Bip(io->cstr_tree);
            Get_Bip(io->cstr_tree->a_nodes[0], io->cstr_tree->a_nodes[0]->v[0],
                    io->cstr_tree);
            if (tree->has_branch_lengths == NO)
              Add_BioNJ_Branch_Lengths(tree, cdata, mod, NULL);
          }

          if (!tree) continue;

          time(&t_beg);
          time(&(tree->t_beg));

          tree->mod          = mod;
          tree->io           = io;
          tree->data         = cdata;
          tree->n_root       = NULL;
          tree->e_root       = NULL;
          tree->n_tot_bl_opt = 0;

          Set_Both_Sides(YES, tree);

          if ((!num_data_set) && (!num_tree) && (!num_rand_tree))
            Check_Memory_Amount(tree);

          if (io->cstr_tree && !Check_Topo_Constraints(tree, io->cstr_tree))
          {
            PhyML_Printf("\n\n== The initial tree does not satisfy the "
                         "topological constraint.");
            PhyML_Printf("\n== Please use the user input tree option with an "
                         "adequate tree topology.");
            Exit("\n");
          }

          Connect_CSeqs_To_Nodes(tree->data, tree->io, tree);
          Make_Tree_For_Pars(tree);
          Make_Tree_For_Lk(tree);
          Make_Spr(tree);
          Br_Len_Not_Involving_Invar(tree);
          Unscale_Br_Len_Multiplier_Tree(tree);

#ifdef BEAGLE
          if (mod->bootstrap == YES)
          {
            PhyML_Printf("\n== PhyML-BEAGLE does not support bootstrap "
                         "analysis yet... ");
            Exit("\n");
          }
          if (mod->ras->invar == YES)
          {
            PhyML_Printf("\n== PhyML-BEAGLE does not support invariant site "
                         "models yet... ");
            Exit("\n");
          }
#endif

          if (tree->io->print_json_trace == YES)
            JSON_Tree_Io(tree, tree->io->fp_out_json_trace);

          Set_Update_Eigen(YES, tree->mod);
          Lk(NULL, tree);
          Set_Update_Eigen(NO, tree->mod);

          PhyML_Printf("\n. Init log-likelihood: %f", tree->c_lnL);

          // LG4X (Le, Dang & Gascuel 2012): "the starting tree, rates, and
          // weights (=0.25) are first estimated with LG4M", then the six
          // rate and weight parameters are optimised
          if (tree->mod->lg4 && tree->mod->lg4->model == LG4_X &&
              tree->mod->s_opt->opt_subst_param == YES &&
              (tree->mod->s_opt->opt_topo ||
               tree->mod->s_opt->opt_bl_one_by_one))
          {
            int alpha_opt = tree->mod->ras->alpha->optimize;
            tree->mod->ras->alpha->optimize = YES;
            LG4_Switch(tree->mod, LG4_M);
            Set_Update_Eigen(YES, tree->mod);
            Lk(NULL, tree);
            Set_Update_Eigen(NO, tree->mod);
            PhyML_Printf("\n. LG4X: first fitting LG4M (log-likelihood %f)",
                         tree->c_lnL);
            Search_Or_Optimize(tree);
            tree->mod->ras->alpha->optimize = alpha_opt;
            LG4_Switch(tree->mod, LG4_X);
            Set_Update_Eigen(YES, tree->mod);
            Lk(NULL, tree);
            Set_Update_Eigen(NO, tree->mod);
            PhyML_Printf("\n. LG4X: starting from the LG4M tree and rates "
                         "(log-likelihood %f)", tree->c_lnL);
          }

          Search_Or_Optimize(tree);

        
          /* if(tree->mod->gamma_mgf_bl)
           * Best_Root_Position_IL_Model(tree); */

          Set_Both_Sides(YES, tree);
          Lk(NULL, tree);
          Pars(NULL, tree);
          Get_Tree_Size(tree);
          PhyML_Printf("\n\n. Log likelihood of the current tree: %.*f.",
                       DECIMAL_DIG, tree->c_lnL);

          if (tree->io->ancestral == YES) Ancestral_Sequences(tree, YES);

          Check_Br_Lens(tree);
          Br_Len_Involving_Invar(tree);
          Rescale_Br_Len_Multiplier_Tree(tree);

          if (!tree->n_root) Get_Best_Root_Position(tree);

          /* Print the tree estimated using the current random (or
           * BioNJ) starting tree */
          /* if(io->mod->s_opt->n_rand_starts > 1) */
          if (orig_random_input_tree == YES)
          {
            Print_Tree(io->fp_out_trees, tree);
            fflush(NULL);
          }

          /* Record the most likely tree in a string of characters */
          if (tree->c_lnL > best_lnL)
          {
            best_lnL = tree->c_lnL;
            if (most_likely_tree) Free(most_likely_tree);
            most_likely_tree = Write_Tree(tree);

            time(&t_end);

            Print_Fp_Out(
                io->fp_out_stats, t_beg, t_end, tree, io, num_data_set + 1,
                (orig_random_input_tree == YES) ? (num_rand_tree) : (num_tree),
                (num_rand_tree == io->mod->s_opt->n_rand_starts - 1) ? (YES)
                                                                     : (NO),
                io->precision);

            if (tree->io->print_site_lnl) Print_Site_Lk(tree, io->fp_out_lk);
          }

          /* Start from BioNJ tree */
          if ((num_rand_tree == io->mod->s_opt->n_rand_starts - 1) &&
              (tree->mod->s_opt->random_input_tree))
          {
            /* Do one more iteration in the loop, but don't randomize the tree
             */
            tree->mod->s_opt->n_rand_starts++;
            tree->mod->s_opt->random_input_tree = NO;
          }
#ifdef BEAGLE
          finalize_beagle_instance(tree);
#endif
          Free_Best_Spr(tree);
          Free_Spr_List_One_Edge(tree);
          Free_Spr_List_All_Edge(tree);
          Free_Tree_Pars(tree);
          Free_Tree_Lk(tree);
          Free_Tree(tree);
        } // Tree done

        // The tree is written now, so that the file holds it during the
        // bootstrap, and rewritten with its supports below. With one data
        // set, the tree of each input tree replaces the previous one.
        if (io->n_data_sets == 1)
          tree_pos = tree_file_start;
        else
        {
          fseek(io->fp_out_tree, 0, SEEK_END);
          tree_pos = ftell(io->fp_out_tree);
        }
        Truncate_Tree_File(io->fp_out_tree, tree_pos);
        if (most_likely_tree)
          PhyML_Fprintf(io->fp_out_tree, "%s\n", most_likely_tree);

        /* Launch bootstrap analysis */
        if (io->do_boot || io->do_tbe || io->do_bayesboot)
        {
          if (!io->quiet)
            PhyML_Printf(
                "\n\n. Launch bootstrap analysis on the most likely tree...");

#ifdef MPI
          Bcast_String(&most_likely_tree);
          if (!io->quiet)
            PhyML_Printf("\n\n. The bootstrap analysis will use %d CPU%c.",
                         Global_numTask, Global_numTask > 1 ? 's' : '\0');
#endif

          most_likely_tree =
              Bootstrap_From_String(most_likely_tree, cdata, mod, io);

          PhyML_Printf("\n\n. Completed the bootstrap analysis succesfully.");
          fflush(NULL);
        }
        else if (io->ratio_test != NO)
        {
          /* Launch aLRT */
          most_likely_tree = aLRT_From_String(most_likely_tree, cdata, mod, io);
        }

        /* Print the most likely tree in the output file */
        if (!io->quiet)
          PhyML_Printf("\n\n. Printing the most likely tree in file '%s'.",
                       Basename(io->out_tree_file));
        Truncate_Tree_File(io->fp_out_tree, tree_pos);

        t_tree *dum;
        dum       = Read_Tree(&most_likely_tree);
        dum->data = cdata;
        dum->mod  = mod;
        dum->io   = io;
        Connect_CSeqs_To_Nodes(cdata, io, dum);
        Insert_Duplicates(dum);
        // Structure models: branch lengths in expected substitutions per
        // site, i.e. multiplied by the model's mean rate (Le & Gascuel 2010)
        if (mod->smod)
        {
          phydbl f = STRUCT_Mean_Rate_Data(mod->smod, cdata);
          for (int e = 0; e < 2 * dum->n_otu - 3; ++e) dum->a_edges[e]->l->v *= f;
        }
        Free(most_likely_tree);
        most_likely_tree = Write_Tree(dum);
        Free_Tree(dum);

        PhyML_Fprintf(io->fp_out_tree, "%s\n", most_likely_tree);

        if (io->n_trees > 1 && io->n_data_sets > 1) break;
      }
      Free_Calign(cdata);
    }
    else
    {
      PhyML_Printf("\n== No data was found.\n");
      PhyML_Printf("\n== Err. in file %s at line %d\n", __FILE__, __LINE__);
      Exit("\n");
    }
    Free_Model_Complete(mod);
    if (mod->smod)
    {
      STRUCT_Free(mod->smod);
      mod->smod = NULL;
    }
    if (mod->lg4)
    {
      LG4_Free(mod->lg4);
      mod->lg4 = NULL;
    }
    SMS_Free_Fits(io->sms_fit);
    io->sms_fit = NULL;
  }

  if (most_likely_tree) Free(most_likely_tree);

  if (mod->s_opt->n_rand_starts > 1)
    PhyML_Printf("\n. Best log likelihood: %f\n", best_lnL);

  Free_Optimiz(mod->s_opt);
  Free_Model_Basic(mod);

  if (io->fp_in_constraint_tree) fclose(io->fp_in_constraint_tree);
  if (io->fp_in_align) fclose(io->fp_in_align);
  if (io->fp_in_tree) fclose(io->fp_in_tree);
  if (io->fp_out_lk) fclose(io->fp_out_lk);
  if (io->fp_out_tree) fclose(io->fp_out_tree);
  if (io->fp_out_trees) fclose(io->fp_out_trees);
  if (io->fp_out_stats) fclose(io->fp_out_stats);
  if (io->fp_out_trace) fclose(io->fp_out_trace);
  if (io->fp_out_json_trace) fclose(io->fp_out_json_trace);

  if (io->fp_in_constraint_tree != NULL) Free_Tree(io->cstr_tree);
  Free_Input(io);

  time(&t_end);
  Print_Time_Info(t_beg, t_end);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* PhyML_Main with a recovery point: an error comes back here with its code
   instead of ending the process (see error.h). */

int PhyML_Run(int argc, char **argv)
{
  jmp_buf env;
  int     code;

  if (PhyML_Error_Has_Recovery()) return PHYML_ERR_BUSY;

  PhyML_Error_Clear();
  optind = 1; // getopt starts again from the first argument

  // What the run allocates and opens is recorded from here on, and released
  // below: nothing is left of a run, whether it fails or not
  PhyML_Track_Begin();

  code = setjmp(env);
  if (code == 0)
  {
    PhyML_Error_Set_Recovery(&env);
    PhyML_Main(argc, argv);
    code = PHYML_OK;
  }
  else if (code == PHYML_STOP_REQUESTED)
    code = PHYML_OK;

  PhyML_Error_Set_Recovery(NULL);
  PhyML_Error_Scope(PHYML_ERR);
  fflush(NULL);
  PhyML_Track_End(NULL, NULL);

  return code;
}

#endif
