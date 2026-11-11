

// PhyML:  a program that  computes maximum likelihood phylogenies from
// DNA or AA homologous sequences.

// Copyright (C) Stephane Guindon. Oct 2003 onward.

// All parts of the source except where indicated are distributed under
// the GNU public licence. See http://www.opensource.org for details.

#include "alrt.h"
#include "bionj.h"
#include "eigen.h"
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

/* The analysis is in the library (libphyml.c). Without a recovery point
   (PhyML_Run), an error ends the process with status 1, as before. */
int main(int argc, char **argv) { return PhyML_Main(argc, argv); }

#elif (M4)
#include "m4.h"
int main(int argc, char **argv)
{
  M4_main(argc, argv);
  return 1;
}

#elif (PART)
#include "mg.h"
int main(int argc, char **argv)
{
  PART_main(argc, argv);
  return 1;
}

/* #elif(PHYTIME) */
/* #include "times.h" */
/* int main(int argc, char **argv) */
/* { */
/*   TIMES_main(argc, argv); */
/*   return 1; */
/* } */

#elif (PHYCONT)
#include "continuous.h"
int main(int argc, char **argv)
{
  CONT_main(argc, argv);
  return 1;
}

#elif (RF)
int main(int argc, char **argv)
{
  // option *io;
  // int     r_seed;

  // io = (option *)Get_Input(argc, argv);
  // if (!io) return (0);

  // r_seed = (io->r_seed < 0) ? (time(NULL)) : (io->r_seed);
  // srand(r_seed);
  // io->r_seed = r_seed;

  // Get_Seq(io);
  // Shuffle_Sites(io->mod->ras->pinvar->v, io->data, io->n_otu);
  // Print_Seq(stdout, io->data, io->n_otu);

  /* t_tree *tree1, *tree2; */
  /* FILE *fp_tree1, *fp_tree2; */
  /* int i,j; */

  /* fp_tree1 = (FILE *)fopen(argv[1],"r"); */
  /* fp_tree2 = (FILE *)fopen(argv[2],"r"); */

  /* tree1 = Read_Tree_File_Phylip(fp_tree1); */
  /* tree2 = Read_Tree_File_Phylip(fp_tree2); */

  /* Prune_Tree(tree1,tree2); */
  /* Prune_Tree(tree2,tree1); */

  /* PhyML_Printf("%s\n%s", */
  /*              Write_Tree(tree1), */
  /*              Write_Tree(tree2)); */

  /* Match_Nodes_In_Small_Tree(tree1,tree2); */

  /* For(i,2*tree1->n_otu-2) */
  /*   { */
  /*     printf("\n. Node %d in tree1 matches node %d in
   * tree2",i,(tree1->noeud[i]->match_node)?(tree1->noeud[i]->match_node->num):(-1));
   */
  /*   } */

  t_tree *tree1, *tree2;
  FILE *fp_tree1, *fp_tree2;
  int i,j,rf,n_edges,n_common,bip_size;
  phydbl thresh;
  t_edge *b; 

  fp_tree1 = (FILE *)fopen(argv[1],"r");
  fp_tree2 = (FILE *)fopen(argv[2],"r");

  tree1 = Read_Tree_File_Phylip(fp_tree1); 
  tree2 = Read_Tree_File_Phylip(fp_tree2); 

  // Get_Rid_Of_Prefix('_',tree1);
  // Find_Common_Tips(tree1,tree2); 

  Alloc_Bip(tree1); 
  Alloc_Bip(tree2); 

  Get_Bip(tree1->a_nodes[0],tree1->a_nodes[0]->v[0],tree1); 
  Get_Bip(tree2->a_nodes[0],tree2->a_nodes[0]->v[0],tree2);

  PhyML_Printf("\n. RFxLENGTH=%f",
               .5 * (Compare_Bip(tree1, tree2, NO, TREE_COMP_RF_PLUS_LENGTH, -1.) + Compare_Bip(tree2, tree1, NO, TREE_COMP_RF_PLUS_LENGTH, -1.)));
  PhyML_Printf("\n. PLAIN RF=%f",.5 * (Compare_Bip(tree1, tree2, NO, TREE_COMP_RF_PLAIN, -1.) + Compare_Bip(tree2, tree1, NO, TREE_COMP_RF_PLAIN, -1.)));
  // PhyML_Printf("\n. INT LENGTHS=%f",.5 * (Compare_Bip(tree1, tree2, NO, TREE_COMP_LENGTH_SHARED_EDGES) + Compare_Bip(tree2, tree1, NO, TREE_COMP_LENGTH_SHARED_EDGES)));
  // PhyML_Printf("\n. EXT LENGTHS=%f",.5 * (Compare_Bip(tree1, tree2, NO, TREE_COMP_LENGTH_EXTERNAL_EDGES) + Compare_Bip(tree2, tree1, NO, TREE_COMP_LENGTH_EXTERNAL_EDGES)));
  PhyML_Printf("\n. INT LENGTHS=%f",(Compare_Bip(tree1, tree2, NO, TREE_COMP_LENGTH_SHARED_EDGES, -1.)));
  PhyML_Printf("\n. EXT LENGTHS=%f",(Compare_Bip(tree1, tree2, NO, TREE_COMP_LENGTH_EXTERNAL_EDGES, -1.)));

  // for (i = 0; i < 2 * tree1->n_otu - 3; ++i) tree1->a_edges[i]->bip_score = 0;
  // for (i = 0; i < 2 * tree2->n_otu - 3; ++i) tree2->a_edges[i]->bip_score = 0;

  /*   rf = 0; */
  /*   n_edges = 0; */

  /*   /\* First tree *\/ */
  /*   For(i,2*tree1->n_otu-3)  */
  /*     { */
  /*       /\* Consider the branch only if the corresponding bipartition has
   * size > 1 *\/ */
  /*       b = tree1->a_edges[i]; */
  /*       bip_size = MIN(b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]);
   */

  /*       if(bip_size > 1) */
  /* 	{ */
  /* 	  /\* with non-zero length *\/ */
  /* 	  if(tree1->a_edges[i]->l > thresh)   */
  /* 	    { */
  /* 	      n_edges++; */
  /* 	      /\* This t_edge is not found in tree2 *\/ */
  /* 	      if(!tree1->a_edges[i]->bip_score) rf++; ; */
  /* 	    } */
  /* 	} */
  /*     } */

  /*   /\* Second tree *\/ */
  /*   For(i,2*tree2->n_otu-3)  */
  /*     { */
  /*       b = tree2->a_edges[i]; */
  /*       bip_size = MIN(b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]);
   */

  /*       if(bip_size > 1) */
  /* 	{ */
  /* 	  if(tree2->a_edges[i]->l > thresh)   */
  /* 	    { */
  /* 	      n_edges++; */
  /* 	      /\* This t_edge is not found in tree1 *\/ */
  /* 	      if(!tree2->a_edges[i]->bip_score) rf++; ; */
  /* 	    } */
  /* 	} */
  /*     } */

  /*   if(!n_edges) */
  /*     { */
  /*       Exit("\n. No comparable internal edges were found.\n"); */
  /*     } */
  /*   else */
  /*     { */
  /*       PhyML_Printf("\n. Robinson and Foulds distance:
   * %f.",(double)rf/(n_edges)); */
  /* /\*       PhyML_Printf("\n. %d internal edges were processed (%d in the
   * first tree, %d in the second).\n",n_edges,n_edges_t1,n_edges-n_edges_t1);
   * *\/ */
  /*       PhyML_Printf("\n"); */
  /*     } */

  return 1;
}

#elif (TIPORDER)
#include "tiporder.h"
int main(int argc, char **argv)
{
  TIPO_main(argc, argv);
  return 1;
}

#elif (TEST)
#include "xml.h"
int main(int argc, char **argv)
{

  //   /* Prediction using linear extrapolation of velocities estimated at the
  //   tips. Takes as input */
  //   /* the XML file used for running the MCMC analysis (so as to get the
  //   coordinates + sampling dates), */
  //   /* the list of sampled trees, the time interval for which prediction will
  //   be made and the burnin */
  //   /* proportion */

  //    option *io;
  //    t_tree *tree;
  //    xml_node *root;
  //    phydbl date_old,date_recent,burnin;
  //    short int first;
  //    int stepsize,sampsize;

  //    date_recent = atof(argv[4]);
  //    date_old = atof(argv[5]);
  //    burnin = atof(argv[6]);
  //    sampsize = 1000;
  //    first = YES;

  //    PhyML_Printf("\n. Time interval considered:
  //    [%f,%f]",date_old,date_recent); assert(date_old < date_recent);

  //    io = (option *)Get_Input(argc,argv);
  //    if(!io) return(0);

  //    Read_User_Tree(NULL,NULL,io);

  //    root = XML_Load_File(io->fp_in_xml);

  //    stepsize = io->treelist->list_size*(1. - burnin) / sampsize;

  //    for(int i=0;i<io->treelist->list_size;i++)
  //      {
  //        if(i > (int)(burnin*io->treelist->list_size))
  //          {
  //            PhyML_Fprintf(stderr,"\n. Processing tree %d",i+1);
  //            tree = io->treelist->tree[i];

  //            tree->times = TIMES_Make_Time_Struct(tree->n_otu);
  //            TIMES_Init_Time_Struct(tree->times,NULL,tree->n_otu);

  //            tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  //            RATES_Init_Rate_Struct(tree->rates,NULL,tree->n_otu);

  //            tree->mmod = PHYREX_Make_Migrep_Model(tree->n_otu,2);
  //            tree->mmod->n_dim = 2;
  //            PHYREX_Set_Default_Migrep_Mod(tree->n_otu,tree->mmod);

  //            XML_Read_Calibration(root,tree);
  //            MIXT_Chain_Cal(tree);

  //            TIMES_Randomize_Tip_Times_Given_Calibrations(tree); // Topology
  //            is unchanged tree->rates->clock_r = 1.0;
  //            TIMES_Bl_To_Times(1,tree);
  //            Update_Ancestors(tree->n_root,tree->n_root->v[2],tree->n_root->b[2],tree);
  //            Update_Ancestors(tree->n_root,tree->n_root->v[1],tree->n_root->b[1],tree);

  //            PHYREX_Make_And_Connect_Tip_Disks(tree);
  //            PHYREX_Tree_To_Ldsk(tree);

  //            Node_Labels_To_Velocities(tree);
  //            Node_Labels_To_Locations(tree);

  //            if(first == YES) PhyML_Printf("\n. XXX Tree\t Time\t Tax\t
  //            Longitude\t Latitude\t NextLongitude\t NextLatitude"); first =
  //            NO;

  //            for(int j=0;j<2*tree->n_otu-1;++j)
  //              {
  //                if(tree->a_nodes[j]->tax == YES)
  //                  {
  //                    phydbl survival_duration,survival_rate,t,delta_t;
  //                    int k,K;

  //                    survival_rate = 1.0;
  //                    K = 4;
  //                    delta_t = (date_recent - date_old) / K;

  //                    survival_duration = Rexp(survival_rate);

  //                    if(tree->times->nd_t[j] + survival_duration > date_old)
  //                      {
  //                        t = date_old;
  //                        k = 0;
  //                        do
  //                          {
  //                            PhyML_Printf("\n. XXX %d\t %f\t %d\t %f\t %f\t
  //                            %f\t %f",
  //                                         i+1,
  //                                         tree->times->nd_t[j],
  //                                         tree->a_nodes[j]->tax,
  //                                         tree->a_nodes[j]->ldsk->coord->lonlat[0],
  //                                         tree->a_nodes[j]->ldsk->coord->lonlat[1],
  //                                         tree->a_nodes[j]->ldsk->coord->lonlat[0]
  //                                         + (k*delta_t) *
  //                                         tree->a_nodes[j]->ldsk->veloc->deriv[0],
  //                                         tree->a_nodes[j]->ldsk->coord->lonlat[1]
  //                                         + (k*delta_t) *
  //                                         tree->a_nodes[j]->ldsk->veloc->deriv[1]);
  //                            k++;
  //                            t += delta_t;
  //                          }
  //                        while(t < MIN(date_recent,tree->times->nd_t[j] +
  //                        survival_duration));
  //                      }
  //                  }
  //              }

  //            PHYREX_Free_Ldsk_Struct(tree);
  //            RATES_Free_Rates(tree->rates);
  //            TIMES_Free_Times(tree->times);
  //            Free_Mmod(tree->mmod);
  //            Free_Tree(tree);

  //          }
  //      }

  //    Exit("\n");

  // /* Select samples in the 801 sequence WNV data set within a user defined
  // date range ./test --xml=../WNV_RRW_tree_1.xml 2000 1990 1.0 */ option *io;
  // xml_node
  // *beast_root,*beast_taxa,*beast_taxon,*beast_date,*beast_loc,*beast_seq,**valid_taxa;
  // xml_node *phyrex_root,*nd,*ndnd;
  // int i,n_selected,n_selected_max,r_seed;
  // FILE *fp_xml,*fp_coord,*fp_seq;
  // char *filename,*dum,*xml_filename;
  // phydbl date_recent, date_f;
  // int *permut;
  // phydbl *select_proba,sum;

  // valid_taxa = (xml_node **)mCalloc(800,sizeof(xml_node *));
  // select_proba = (phydbl *)mCalloc(800,sizeof(phydbl));

  // filename = (char *)mCalloc(100,sizeof(char));
  // xml_filename = (char *)mCalloc(100,sizeof(char));
  // r_seed = time(NULL);
  // PhyML_Printf("\n seed: %d",r_seed);
  // srand(r_seed);

  // io = (option *)Get_Input(argc,argv);
  // if(!io) return(0);

  // date_recent = atof(argv[2]);
  // n_selected_max = atoi(argv[3]);

  // dum = (char *)mCalloc(100,sizeof(char));
  // sprintf(dum,"%s%d%s","wnv_config_",r_seed,".xml");
  // strcpy(xml_filename,dum);
  // Free(dum);

  // PhyML_Printf("\n. Upper limit of time: %f",date_recent);

  // strcpy(filename,"coord.txt");
  // fp_coord = Openfile(filename,WRITE);

  // PhyML_Fprintf(fp_coord,"traits lat long\n");
  // PhyML_Fprintf(fp_coord,"|NorthEast| 1000 1000\n");
  // PhyML_Fprintf(fp_coord,"|SouthWest| -1000 -1000");

  // strcpy(filename,"seq.txt");
  // fp_seq = Openfile(filename,WRITE);

  // dum = (char *)mCalloc(100,sizeof(char));
  // sprintf(dum,"%s%d","wnv_predict_",r_seed);
  // phyrex_root = Generate_PhyREX_XMLObj("ibm",dum,"seq.txt","coord.txt");
  // Free(dum);

  // beast_root = XML_Load_File(io->fp_in_xml);

  // beast_taxa = XML_Search_Node_Name("taxa",NO,beast_root);

  // i = 1;
  // n_selected = 0;
  // beast_taxon = beast_taxa->child;
  // do
  //   {
  //     beast_date = XML_Search_Node_Name("date",NO,beast_taxon);
  //     assert(beast_date);

  //     date_f = atof(XML_Get_Attribute_Value(beast_date,"value"));

  //     if(date_f < date_recent)
  //       {
  //         valid_taxa[n_selected] = beast_taxon;
  //         select_proba[n_selected] = exp(-fabs(date_f-date_recent));
  //         n_selected++;
  //       }

  //     beast_taxon = beast_taxon->next;
  //   }
  // while(beast_taxon);

  // PhyML_Fprintf(fp_seq,"%d 10302\n",MIN(n_selected,n_selected_max));
  // permut = Permutate(MIN(n_selected,n_selected_max));

  // sum = 0.0;
  // for(i=0;i<n_selected;++i) sum += select_proba[i];
  // for(i=0;i<n_selected;++i) select_proba[i] /= sum;

  // for(i=0;i<n_selected;++i) PhyML_Printf("\n. Target %s (%d) prob:
  // %f",XML_Get_Attribute_Value(valid_taxa[i],"id"),i,select_proba[i]);

  // for(i=0;i<MIN(n_selected,n_selected_max);++i)
  //   {
  //     permut[i] = Sample_i_With_Proba_pi(select_proba,n_selected);

  //     beast_date = XML_Search_Node_Name("date",NO,valid_taxa[permut[i]]);
  //     assert(beast_date);
  //     date_f = atof(XML_Get_Attribute_Value(beast_date,"value"));

  //     nd = XML_Add_Node(phyrex_root,"clade");
  //     dum = (char *)mCalloc(100,sizeof(char));
  //     sprintf(dum,"%s%d","clad",i+1);
  //     nd->attr = XML_Make_Attribute(NULL,"id",dum);
  //     ndnd = XML_Add_Node(nd,"taxon");
  //     ndnd->attr =
  //     XML_Make_Attribute(NULL,"value",XML_Get_Attribute_Value(valid_taxa[permut[i]],"id"));

  //     nd = XML_Add_Node(phyrex_root,"calibration");
  //     dum = (char *)mCalloc(100,sizeof(char));
  //     sprintf(dum,"%s%d","cal",i+1);
  //     nd->attr = XML_Make_Attribute(NULL,"id",dum);
  //     Free(dum);

  //     ndnd = XML_Add_Node(nd,"lower");
  //     XML_Set_Node_Value(ndnd,XML_Get_Attribute_Value(beast_date,"value"));

  //     ndnd = XML_Add_Node(nd,"upper");
  //     XML_Set_Node_Value(ndnd,XML_Get_Attribute_Value(beast_date,"value"));

  //     ndnd = XML_Add_Node(nd,"appliesto");
  //     dum = (char *)mCalloc(100,sizeof(char));
  //     sprintf(dum,"%s%d","clad",i+1);
  //     ndnd->attr = XML_Make_Attribute(NULL,"clade.id",dum);
  //     Free(dum);

  //     PhyML_Fprintf(fp_coord,"\n");
  //     PhyML_Fprintf(fp_coord,"
  //     %s",XML_Get_Attribute_Value(valid_taxa[permut[i]],"id"));

  //     beast_loc =
  //     XML_Search_Node_Generic("attr","name","latitude",NO,valid_taxa[permut[i]]);
  //     assert(beast_loc);
  //     PhyML_Fprintf(fp_coord,"\t %s",beast_loc->value);

  //     beast_loc =
  //     XML_Search_Node_Generic("attr","name","longitude",NO,valid_taxa[permut[i]]);
  //     assert(beast_loc);
  //     PhyML_Fprintf(fp_coord,"\t %s",beast_loc->value);

  //     beast_seq =
  //     XML_Search_Node_Attribute_Value("idref",XML_Get_Attribute_Value(valid_taxa[permut[i]],"id"),NO,beast_root);
  //     assert(beast_seq);
  //     PhyML_Printf("\n. Found match with %s (%d) prob:
  //     %f",XML_Get_Attribute_Value(valid_taxa[permut[i]],"id"),permut[i],select_proba[permut[i]]);
  //     PhyML_Fprintf(fp_seq,"%s",XML_Get_Attribute_Value(valid_taxa[permut[i]],"id"));
  //     PhyML_Fprintf(fp_seq,"\t\t%s\n",beast_seq->value);

  //     select_proba[permut[i]] = 0.0;
  //     sum = 0.0;
  //     for(int i=0;i<n_selected;++i) sum += select_proba[i];
  //     for(int i=0;i<n_selected;++i) select_proba[i] /= sum;

  //   }

  // PhyML_Printf("\n. Number of selected sequences:
  // %d",MIN(n_selected,n_selected_max));

  // fp_xml = Openfile(xml_filename,WRITE);
  // XML_Write_XML_Graph(fp_xml,phyrex_root);
  // fclose(fp_xml);

  // fclose(fp_coord);
  // fclose(fp_seq);

  // Philippe Lemey data preprocessisng steps
  option *io;
  int     i;
  char   *year, *lat, *lon;
  t_tree *tree;

  io = (option *)Get_Input(argc, argv);
  if (!io) return (0);

  tree = Read_User_Tree(NULL, NULL, io);

  PhyML_Printf("\n. n_otus: %d", tree->n_otu);
  for (i = 0; i < tree->n_otu; ++i)
  {
    PhyML_Printf("\n<clade id=\"clad%d\">", i + 1);
    PhyML_Printf("\n\t<taxon value=\"%s\"/>", tree->a_nodes[i]->name);
    PhyML_Printf("\n</clade>");
    PhyML_Printf("\n<calibration id=\"cal%d\">", i + 1);
    PhyML_Printf("\n\t<lower>-%f</lower>",
                 atof(strrchr(tree->a_nodes[i]->name, '|') + 1));
    PhyML_Printf("\n\t<upper>-%f</upper>",
                 atof(strrchr(tree->a_nodes[i]->name, '|') + 1));
    PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>", i + 1);
    PhyML_Printf("\n</calibration>");
  }

  PhyML_Printf("\n\n%d\t 4", tree->n_otu);
  for (i = 0; i < tree->n_otu; ++i)
  {
    PhyML_Printf("\n%s\tATGC", tree->a_nodes[i]->name);
  }
  Exit("\n");

  //    Get_Seq(io);

  //    for(i=0;i<io->n_otu;i++)
  //       {
  //         sscanf(io->data[i]->name,"%d",&year);
  //         PhyML_Printf("\n%s\t%d",io->data[i]->name,year);
  //       }

  for (i = 0; i < io->n_otu; i++)
  {
    //    sscanf(io->data[i]->name,"%[^_]_%[^_]_%lf",a,b,&year);
    //    PhyML_Printf("\n. a: %s b: %s",a,b);
    strtok(io->data[i]->name, "_");
    char *a = strtok(NULL, "_");
    year    = strtok(NULL, "_");
    lat     = strtok(NULL, "_");
    lon     = strtok(NULL, "_");

    PhyML_Printf("\n\t <taxon id=\"%s_%s_%s\">", io->data[i]->name, a, year);
    PhyML_Printf(
        "\n\t\t <date value=\"%s\" direction=\"forwards\" units=\"years\"/>",
        year);
    PhyML_Printf("\n\t\t <attr name=\"lat\">");
    PhyML_Printf("\n\t\t\t %s", lat);
    PhyML_Printf("\n\t\t </attr>");
    PhyML_Printf("\n\t\t <attr name=\"lon\">");
    PhyML_Printf("\n\t\t\t %s", lon);
    PhyML_Printf("\n\t\t </attr>");
    PhyML_Printf("\n\t\t <attr name=\"location\">");
    PhyML_Printf("\n\t\t\t %s %s", lat, lon);
    PhyML_Printf("\n\t\t </attr>");
    PhyML_Printf("\n\t </taxon>");

    // PhyML_Printf("\n<clade id=\"clad%d\">",i+1);
    // PhyML_Printf("\n\t<taxon value=\"%s\"/>",io->data[i]->name);
    // PhyML_Printf("\n</clade>");
    // PhyML_Printf("\n<calibration id=\"cal%d\">",i+1);
    // PhyML_Printf("\n\t<lower>%d</lower>",year);
    // PhyML_Printf("\n\t<upper>%d</upper>",year);
    // PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>",i+1);
    // PhyML_Printf("\n</calibration>");
  }

  //  FILE *fp;
  //  char *name, *date;
  //  int i;

  //  name = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  //  date = (char *)mCalloc(T_MAX_NAME, sizeof(char));

  //  i = 0;
  //  fp = Openfile(argv[1], READ);
  //  do
  //  {
  //      if (fscanf(fp, "%s", name) == EOF)
  //          break;
  //      if (fscanf(fp, "%s", date) == EOF)
  //          break;
  //      PhyML_Printf("\n<clade id=\"clad%d\">", i + 1);
  //      PhyML_Printf("\n\t<taxon value=\"%s\"/>", name);
  //      PhyML_Printf("\n</clade>");
  //      PhyML_Printf("\n<calibration id=\"cal%d\">", i + 1);
  //      PhyML_Printf("\n\t<lower>%s</lower>", date);
  //      PhyML_Printf("\n\t<upper>%s</upper>", date);
  //      PhyML_Printf("\n\t<appliesto clade.id=\"clad%d\"/>", i + 1);
  //      PhyML_Printf("\n</calibration>");
  //      ++i;
  //  } while (1);
}

#elif (INVITEE)
#include "invitee.h"
int main(int argc, char **argv)
{
  // My_Function(argc, argv);
  // PhyTime_XML(argc, argv);
  // Get_Input(argc, argv);
  // return 1;
}

#elif (GEO)
#include "geo.h"
int main(int argc, char **argv)
{
  GEO_Main(argc, argv);
  return 1;
}

#elif (defined PHYREX || PHYREXSIM)
#include "phyrex.h"
int main(int argc, char **argv)
{
  PHYREX_Main(argc, argv);
  return 1;
}

#elif (PHYTIME)
#include "date.h"
int main(int argc, char **argv)
{
  DATE_Main(argc, argv);
  return 1;
}

#elif (CHECKPOINT)
#include "checkpoint.h"
int main(int argc, char **argv)
{
  CHECK_Main(argc, argv);
  return 1;
}

#elif (EVOLVE)
#include "evolve.h"
int main(int argc, char **argv)
{
  EVOLVE_Main(argc, argv);
  return 1;
}

#endif
