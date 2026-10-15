/* Unit tests for Newick parsing and writing (io.c), bipartitions and
   Robinson-Foulds comparison (utilities.c), transfer distances (tbe.c) and
   BioNJ (bionj.c). */

#include "tinytest.h"
#include "utilities.h"
#include "bionj.h"
#include "free.h"
#include "io.h"
#include "make.h"

static t_tree *parse(const char *newick)
{
  char *s = (char *)mCalloc((int)strlen(newick) + 1, sizeof(char));
  strcpy(s, newick);
  return Read_Tree(&s);
}

static t_node *tip(t_tree *t, const char *name)
{
  for (int i = 0; i < t->n_otu; i++)
    if (!strcmp(t->a_nodes[i]->name, name)) return t->a_nodes[i];
  tt_fail(__FILE__, __LINE__, "no tip named %s", name);
  return NULL;
}

static void prepare_bip(t_tree *t)
{
  Alloc_Bip(t);
  Get_Bip(t->a_nodes[0], t->a_nodes[0]->v[0], t);
}

static double rf(const char *a, const char *b)
{
  t_tree *t1 = parse(a), *t2 = parse(b);
  Match_Tip_Numbers(t1, t2);
  prepare_bip(t1);
  prepare_bip(t2);
  return Compare_Bip(t1, t2, NO, TREE_COMP_RF_PLAIN, -1.);
}

#define T5 "((A:0.1,B:0.2):0.05,C:0.3,(D:0.4,E:0.5):0.06);"

/* ---------------------------------------------------------- parsing */

TEST(newick, reads_taxa_and_branch_lengths)
{
  t_tree *t = parse(T5);
  CHECK_EQ_INT(t->n_otu, 5);
  CHECK(t->n_root == NULL);
  CHECK_NEAR(tip(t, "A")->b[0]->l->v, 0.1, 1e-12);
  CHECK_NEAR(tip(t, "E")->b[0]->l->v, 0.5, 1e-12);
  CHECK_NEAR(Get_Tree_Size(t), 1.61, 1e-9);
}

TEST(newick, tips_are_numbered_first_and_internal_nodes_have_degree_three)
{
  t_tree *t = parse(T5);
  for (int i = 0; i < t->n_otu; i++) CHECK(t->a_nodes[i]->tax == YES);
  for (int i = t->n_otu; i < 2 * t->n_otu - 2; i++)
  {
    CHECK(t->a_nodes[i]->tax == NO);
    for (int j = 0; j < 3; j++) CHECK(t->a_nodes[i]->v[j] != NULL);
  }
}

TEST(newick, write_round_trip)
{
  t_tree *t = parse(T5);
  char   *s = Write_Tree(t);
  CHECK_STR_EQ(s, "((A:0.10000000,B:0.20000000):0.05000000,C:0.30000000,"
                  "(D:0.40000000,E:0.50000000):0.06000000);");
  t_tree *u = parse(s);
  CHECK_EQ_INT(u->n_otu, 5);
  CHECK_NEAR(Get_Tree_Size(u), Get_Tree_Size(t), 1e-9);
}

TEST(newick, rooted_tree_is_recognised)
{
  t_tree *t = parse("((A:0.1,B:0.2):0.05,(C:0.3,(D:0.4,E:0.5):0.06):0.07);");
  CHECK_EQ_INT(t->n_otu, 5);
  CHECK(t->n_root != NULL);
  CHECK_STR_EQ(Write_Tree(t), "((A:0.10000000,B:0.20000000):0.05000000,"
                              "(C:0.30000000,(D:0.40000000,E:0.50000000):"
                              "0.06000000):0.07000000);");
}

/* Clean_Multifurcation gives the edges it creates the length of a sibling
   edge instead of zero: (A:.1,B:.2,C:.3,D:.4,E:.5) is written back as
   (((A:.1,B:.2):.3,C:.3):.4,D:.4,E:.5), inflating the tree length. */
KNOWN_BUG(newick, multifurcating_root_is_resolved,
          "io.c Clean_Multifurcation copies sibling lengths onto new edges")
{
  t_tree *t = parse("(A:0.1,B:0.2,C:0.3,D:0.4,E:0.5);");
  CHECK_EQ_INT(t->n_otu, 5);
  CHECK_NEAR(Get_Tree_Size(t), 1.5, 1e-9);
}

/* R_rtree (io.c:247-256) frees the subtree string, allocates an empty
   buffer and parses it again without copying the resolved subtrees into it. */
KNOWN_BUG(newick, internal_multifurcation_is_resolved_with_zero_length_edges,
          "io.c:247 R_rtree re-parses an empty buffer after Clean_Multifurcation")
{
  t_tree *t = parse("((A:0.1,B:0.2,C:0.3):0.1,D:0.4,E:0.5);");
  CHECK_EQ_INT(t->n_otu, 5);
  CHECK_NEAR(Get_Tree_Size(t), 1.6, 1e-6);
}

TEST(newick, support_labels_are_accepted)
{
  t_tree *t = parse("((A:0.1,B:0.2)95:0.05,C:0.3,(D:0.4,E:0.5)0.87:0.06);");
  CHECK_EQ_INT(t->n_otu, 5);
  CHECK_NEAR(Get_Tree_Size(t), 1.61, 1e-9);
}

/* Read_Node_Name (io.c:327) copies the name until ':' or '[', so a tip
   without a branch length runs past the end of its buffer. */
KNOWN_BUG(newick, tree_without_branch_lengths,
          "io.c:327 Read_Node_Name overflows when a tip has no branch length")
{
  t_tree *t = parse("((A,B),C,(D,E));");
  CHECK_EQ_INT(t->n_otu, 5);
  CHECK(t->has_branch_lengths == NO);
}

/* Read_Tree skips [comments] when counting taxa, but Sub_Trees does not
   handle them and the parse fails. */
KNOWN_BUG(newick, comments_are_skipped,
          "io.c Sub_Trees does not skip [comments] in Newick strings")
{
  t_tree *t = parse("((A:0.1,B:0.2)[a, comment]:0.05,C:0.3,(D:0.4,E:0.5):0.06);");
  CHECK_EQ_INT(t->n_otu, 5);
}

TEST(newick, read_from_file)
{
  FILE *fp = tmpfile();
  fputs(T5 "\n", fp);
  rewind(fp);
  t_tree *t = Read_Tree_File_Phylip(fp);
  CHECK(t != NULL);
  CHECK_EQ_INT(t->n_otu, 5);
  fclose(fp);
}

/* The label skipper in Read_Tree (io.c:39) walks past the end of the
   string when a '[' is never closed. The parser should reject the tree. */
KNOWN_BUG(newick, unterminated_comment_is_rejected,
          "io.c:39 scans past the end of the string on an unclosed '['")
{
  t_tree *t = parse("((A:0.1,B:0.2)[oops:0.05,C:0.3,(D:0.4,E:0.5):0.06);");
  CHECK(t == NULL);
}

/* -------------------------------------------------- bipartitions, RF */

TEST(bipartitions, rf_of_identical_trees_is_zero)
{
  CHECK_NEAR(rf(T5, T5), 0.0, 0);
}

TEST(bipartitions, rf_ignores_taxon_order_and_rooting)
{
  CHECK_NEAR(rf(T5, "((E:1,D:1):1,(B:1,A:1):1,C:1);"), 0.0, 0);
  CHECK_NEAR(rf(T5, "(((A:1,B:1):1,C:1):1,(D:1,E:1):1);"), 0.0, 0);
}

TEST(bipartitions, rf_counts_missing_splits)
{
  /* {A,B}|{C,D,E} vs {A,C}|{B,D,E}: one of the two splits differs. */
  CHECK_NEAR(rf(T5, "((A:1,C:1):1,B:1,(D:1,E:1):1);"), 1.0, 0);
  /* No shared non-trivial split. */
  CHECK_NEAR(rf(T5, "((A:1,D:1):1,B:1,(C:1,E:1):1);"), 2.0, 0);
}

TEST(bipartitions, rf_is_symmetric)
{
  const char *a = "(((A:1,B:1):1,C:1):1,D:1,((E:1,F:1):1,G:1):1);";
  const char *b = "(((A:1,G:1):1,C:1):1,D:1,((E:1,F:1):1,B:1):1);";
  CHECK_NEAR(rf(a, b), rf(b, a), 0);
  CHECK(rf(a, b) > 0);
}

TEST(bipartitions, bip_sizes_sum_to_n)
{
  t_tree *t = parse(T5);
  prepare_bip(t);
  for (int i = 0; i < 2 * t->n_otu - 3; i++)
  {
    t_edge *b = t->a_edges[i];
    CHECK_EQ_INT(b->left->bip_size[b->l_r] + b->rght->bip_size[b->r_l], t->n_otu);
  }
}

TEST(bipartitions, bip_score_marks_shared_edges)
{
  t_tree *t1 = parse(T5), *t2 = parse("((A:1,C:1):1,B:1,(D:1,E:1):1);");
  int     shared_internal = 0;
  Match_Tip_Numbers(t1, t2);
  prepare_bip(t1);
  prepare_bip(t2);
  Compare_Bip(t1, t2, NO, TREE_COMP_RF_PLAIN, -1.);
  for (int i = 0; i < 2 * t1->n_otu - 3; i++)
  {
    t_edge *b = t1->a_edges[i];
    if (b->left->tax || b->rght->tax)
      CHECK_EQ_INT(b->bip_score, 1); /* trivial splits are always shared */
    else
      shared_internal += b->bip_score;
  }
  CHECK_EQ_INT(shared_internal, 1);
}

TEST(bipartitions, brlen_cutoff_ignores_short_edges_of_second_tree)
{
  t_tree *t1 = parse(T5), *t2 = parse("((A:1,B:1):0.001,C:1,(D:1,E:1):1);");
  int     shared = 0;
  Match_Tip_Numbers(t1, t2);
  prepare_bip(t1);
  prepare_bip(t2);
  Compare_Bip(t1, t2, NO, TREE_COMP_RF_PLAIN, 0.01);
  for (int i = 0; i < 2 * t1->n_otu - 3; i++)
    if (!t1->a_edges[i]->left->tax && !t1->a_edges[i]->rght->tax)
      shared += t1->a_edges[i]->bip_score;
  CHECK_EQ_INT(shared, 1);
}

/* ------------------------------------------------ transfer distance */

static void tdist(const char *ref, const char *other, t_tree **out)
{
  t_tree *t1 = parse(ref), *t2 = parse(other);
  Match_Tip_Numbers(t1, t2);
  prepare_bip(t1);
  prepare_bip(t2);
  Compare_Bip_Distance(t1, t2);
  *out = t1;
}

TEST(tbe, identical_trees_have_zero_transfer_distance)
{
  t_tree *t;
  tdist(T5, T5, &t);
  for (int i = 0; i < 2 * t->n_otu - 3; i++) CHECK_NEAR(t->a_edges[i]->tdist_score, 0.0, 0);
}

TEST(tbe, moving_one_taxon_costs_one_transfer)
{
  t_tree *t;
  /* Swap B and C: the split {A,B} is one transfer away from {A,C}. */
  tdist(T5, "((A:1,C:1):1,B:1,(D:1,E:1):1);", &t);
  int ones = 0;
  for (int i = 0; i < 2 * t->n_otu - 3; i++)
  {
    double d = t->a_edges[i]->tdist_score;
    CHECK(d == 0.0 || d == 1.0);
    ones += d == 1.0;
  }
  CHECK_EQ_INT(ones, 1);
}

TEST(tbe, transfer_distance_bounded_by_split_size)
{
  const char *a = "(((A:1,B:1):1,(C:1,D:1):1):1,E:1,((F:1,G:1):1,(H:1,I:1):1):1);";
  const char *b = "(((A:1,F:1):1,(C:1,H:1):1):1,E:1,((B:1,G:1):1,(D:1,I:1):1):1);";
  t_tree     *t;
  tdist(a, b, &t);
  for (int i = 0; i < 2 * t->n_otu - 3; i++)
  {
    t_edge *e = t->a_edges[i];
    int     p = MIN(e->left->bip_size[e->l_r], e->rght->bip_size[e->r_l]);
    CHECK(e->tdist_score >= 0.0);
    CHECK(e->tdist_score <= p - 1);
  }
}

/* ------------------------------------------------------------ BioNJ */

/* Path lengths from `from` to every node, by depth-first search. */
static void dfs_dist(t_node *a, t_node *from, double acc, double *out, int n)
{
  if (a->tax && a->num < n) out[a->num] = acc;
  for (int i = 0; i < 3; i++)
    if (a->v[i] && a->v[i] != from) dfs_dist(a->v[i], a, acc + a->b[i]->l->v, out, n);
}

/* Runs BioNJ on the additive distances of `newick` and returns the tree. */
static t_tree *bionj_from_tree(const char *newick)
{
  t_tree *truth = parse(newick);
  int     n     = truth->n_otu;

  matrix *mat   = Make_Mat(n);

  mat->n_otu    = n;
  mat->r        = n;
  mat->curr_int = n;
  mat->method   = 1;
  for (int i = 0; i < n; i++)
  {
    strcpy(mat->name[i], truth->a_nodes[i]->name);
    mat->on_off[i] = 1;
  }
  for (int i = 0; i < n; i++)
  {
    double row[64];
    dfs_dist(truth->a_nodes[i], NULL, 0.0, row, n);
    for (int j = 0; j < n; j++) mat->dist[i][j] = (i == j) ? 0.0 : row[j];
  }

  mat->tree = Make_Tree_From_Scratch(n, NULL);
  for (int i = 0; i < n; i++)
  {
    mat->tree->a_nodes[i]->name = (char *)mCalloc(T_MAX_NAME, sizeof(char));
    strcpy(mat->tree->a_nodes[i]->name, truth->a_nodes[i]->name);
  }
  Bionj(mat);
  return mat->tree;
}

TEST(bionj, recovers_tree_from_additive_distances)
{
  const char *truth = "(((A:0.1,B:0.2):0.15,(C:0.12,D:0.3):0.1):0.05,E:0.25,"
                      "((F:0.1,G:0.05):0.2,H:0.4):0.1);";
  t_tree *est = bionj_from_tree(truth);
  char   *s   = Write_Tree(est);
  CHECK_NEAR(rf(truth, s), 0.0, 0);
  CHECK_NEAR(Get_Tree_Size(est), Get_Tree_Size(parse(truth)), 1e-6);
}

TEST(bionj, recovers_branch_lengths)
{
  const char *truth = "((A:0.1,B:0.2):0.05,C:0.3,(D:0.4,E:0.5):0.06);";
  t_tree *est = bionj_from_tree(truth);
  t_tree *ref = parse(truth);
  const char *names[] = {"A", "B", "C", "D", "E"};
  for (int i = 0; i < 5; i++)
    CHECK_NEAR(tip(est, names[i])->b[0]->l->v, tip(ref, names[i])->b[0]->l->v, 1e-6);
}
