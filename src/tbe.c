/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Transfer index of the branches of a reference tree in a bootstrap tree
   (Lemoine et al. 2018, Nature 556:452), in O(n) memory.

   The transfer distance between a reference split L|L' and a bootstrap
   split C|C' is min(|L xor C|, n - |L xor C|); the transfer index of a
   reference branch is its minimum over the bootstrap splits (trivial ones
   included, so it is at most p-1, p = |L| the light side).

   As in Lutteropp, Kozlov & Stamatakis (2020, Bioinformatics 36:2280) and
   gotree (github.com/evolbioinfo/gotree), the computation avoids the
   n x n tables of Brehelin, Gascuel & Martin (2008), which needed three
   (2n)^2 matrices per bootstrap tree (about 380 MB for 4,000 taxa):

   1. Splits present in the bootstrap tree are found by hashing (128-bit
      XOR of fixed random keys per taxon): their transfer index is 0.
   2. For the others, only the light side is used: each of its p taxa walks
      up the rooted bootstrap tree until it meets a node already visited,
      which gives the union U of their paths to the root. The number of
      light taxa below each node of U is summed in post-order. The edges
      above the nodes of U give |L xor C| = p + |C| - 2|L inter C|; any
      other edge lies in a subtree without light taxa, whose best distance
      is n - p - |C| for the largest such subtree, a child of a node of U.
      The work is O(|U| log |U|) instead of O(n), and stops at distance 1.

   Memory is a few arrays of 2n entries. The keys come from a fixed hash
   of the taxon number, not from rand(), so the random stream of the
   bootstrap is unchanged. */

#include "tbe.h"
#include "utilities.h"

typedef struct
{
  uint64_t a, b;
} tbe_key;

/* splitmix64: fixed, well mixed 64-bit keys */
static uint64_t TBE_Mix(uint64_t x)
{
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

static int TBE_Key_Cmp(const void *x, const void *y)
{
  const tbe_key *u = (const tbe_key *)x, *v = (const tbe_key *)y;
  if (u->a != v->a) return (u->a < v->a) ? -1 : 1;
  if (u->b != v->b) return (u->b < v->b) ? -1 : 1;
  return 0;
}

static int TBE_Int_Cmp(const void *x, const void *y)
{
  int u = *(const int *)x, v = *(const int *)y;
  return (u > v) - (u < v);
}

/* A split and its complement get the same key */
static tbe_key TBE_Canonical(uint64_t a, uint64_t b, uint64_t tot_a,
                             uint64_t tot_b)
{
  tbe_key k1 = {a, b}, k2 = {a ^ tot_a, b ^ tot_b};
  return (TBE_Key_Cmp(&k1, &k2) <= 0) ? k1 : k2;
}

/* A tree rooted at an internal node, nodes indexed by their ->num (tips:
   0..n-1, the taxon numbers made common by Match_Tip_Numbers). The tree is
   read through v[], as unrooted: a root node n_root, if any, is not
   reached. */
typedef struct
{
  int      n, n_idx;
  int      root;
  t_node **node;  /* node of each index */
  int     *par;   /* parent index, -1 at the root */
  int     *order; /* indexes in post-order */
  int     *pos;   /* position of each index in order[] */
  int     *size;  /* number of tips below */
  int     *lo;    /* tips below: leaf[lo..hi-1] */
  int     *hi;
  int     *leaf;  /* tip indexes, in DFS order */
  uint64_t *ha, *hb; /* clade hashes */
} tbe_rooted;

static void TBE_Root(t_tree *tree, tbe_rooted *r, const uint64_t *key_a,
                     const uint64_t *key_b)
{
  int  i, k, n = tree->n_otu, top, n_order, n_leaf;
  int *stack, *next;

  r->n     = n;
  r->n_idx = 2 * n;
  r->node  = (t_node **)mCalloc(r->n_idx, sizeof(t_node *));
  r->par   = (int *)mCalloc(r->n_idx, sizeof(int));
  r->order = (int *)mCalloc(r->n_idx, sizeof(int));
  r->pos   = (int *)mCalloc(r->n_idx, sizeof(int));
  r->size  = (int *)mCalloc(r->n_idx, sizeof(int));
  r->lo    = (int *)mCalloc(r->n_idx, sizeof(int));
  r->hi    = (int *)mCalloc(r->n_idx, sizeof(int));
  r->leaf  = (int *)mCalloc(n, sizeof(int));
  r->ha    = (uint64_t *)mCalloc(r->n_idx, sizeof(uint64_t));
  r->hb    = (uint64_t *)mCalloc(r->n_idx, sizeof(uint64_t));
  stack    = (int *)mCalloc(r->n_idx, sizeof(int));
  next     = (int *)mCalloc(r->n_idx, sizeof(int));

  /* root: an internal node with three neighbours */
  i = n;
  while (!tree->a_nodes[i]->v[0] || !tree->a_nodes[i]->v[1] ||
         !tree->a_nodes[i]->v[2])
    i++;
  r->root = tree->a_nodes[i]->num;
  assert(r->root >= 0 && r->root < r->n_idx);

  /* iterative depth-first search: tips in DFS order, nodes in post-order */
  top      = 0;
  n_order  = 0;
  n_leaf   = 0;
  stack[0] = r->root;
  r->node[r->root] = tree->a_nodes[i];
  r->par[r->root]  = -1;
  next[r->root]    = 0;
  r->lo[r->root]   = 0;
  while (top >= 0)
  {
    int     u  = stack[top];
    t_node *nu = r->node[u];

    if (nu->tax)
    {
      r->lo[u]           = n_leaf;
      r->leaf[n_leaf++]  = u;
      r->hi[u]           = n_leaf;
      r->size[u]         = 1;
      r->ha[u]           = key_a[u];
      r->hb[u]           = key_b[u];
      r->pos[u]          = n_order;
      r->order[n_order++] = u;
      top--;
      continue;
    }

    /* next child of u not yet visited */
    for (k = next[u]; k < 3; k++)
      if (nu->v[k] && (r->par[u] < 0 || nu->v[k] != r->node[r->par[u]]))
        break;

    if (k < 3)
    {
      int c   = nu->v[k]->num;
      next[u] = k + 1;
      assert(c >= 0 && c < r->n_idx);
      r->node[c] = nu->v[k];
      r->par[c]  = u;
      next[c]    = 0;
      r->lo[c]   = n_leaf;
      stack[++top] = c;
    }
    else
    {
      /* all children done */
      r->size[u] = 0;
      r->ha[u] = r->hb[u] = 0;
      for (k = 0; k < 3; k++)
        if (nu->v[k] && (r->par[u] < 0 || nu->v[k] != r->node[r->par[u]]))
        {
          int c = nu->v[k]->num;
          r->size[u] += r->size[c];
          r->ha[u] ^= r->ha[c];
          r->hb[u] ^= r->hb[c];
        }
      r->hi[u]            = n_leaf;
      r->pos[u]           = n_order;
      r->order[n_order++] = u;
      top--;
    }
  }
  assert(n_leaf == n);

  Free(stack);
  Free(next);
}

static void TBE_Free_Rooted(tbe_rooted *r)
{
  Free(r->node);
  Free(r->par);
  Free(r->order);
  Free(r->pos);
  Free(r->size);
  Free(r->lo);
  Free(r->hi);
  Free(r->leaf);
  Free(r->ha);
  Free(r->hb);
}

/* Transfer index of the split whose light side is the p taxa of light[] in
   the rooted bootstrap tree b. stamp/cnt/upos are scratch arrays of b->n_idx
   entries; *mark is a fresh value for stamp[]. */
static int TBE_Transfer_Index(const tbe_rooted *b, const int *light, int p,
                              int *stamp, int mark, int *cnt, int *upos)
{
  int i, k, m = 0, n = b->n, best = p - 1;

  /* union of the paths from the light taxa to the root */
  for (i = 0; i < p; i++)
  {
    int u = light[i];
    stamp[u]  = mark;
    cnt[u]    = 1;
    upos[m++] = b->pos[u];
    for (u = b->par[u]; u >= 0 && stamp[u] != mark; u = b->par[u])
    {
      stamp[u]  = mark;
      cnt[u]    = 0;
      upos[m++] = b->pos[u];
    }
  }

  /* children before parents */
  qsort(upos, m, sizeof(int), TBE_Int_Cmp);

  for (i = 0; i < m && best > 1; i++)
  {
    int     u  = b->order[upos[i]];
    t_node *nu = b->node[u];
    int     x, d;

    /* subtrees without light taxa hanging from u */
    if (!nu->tax)
      for (k = 0; k < 3; k++)
        if (nu->v[k] && (b->par[u] < 0 || nu->v[k] != b->node[b->par[u]]))
        {
          int c = nu->v[k]->num;
          if (stamp[c] != mark)
          {
            x = p + b->size[c];
            d = MIN(x, n - x);
            if (d < best) best = d;
          }
        }

    if (b->par[u] < 0) continue; /* the root has no edge above it */

    x = p + b->size[u] - 2 * cnt[u];
    d = MIN(x, n - x);
    if (d < best) best = d;
    cnt[b->par[u]] += cnt[u];
  }

  return best;
}

/* Adds to the tdist_score of each branch of ref_tree its transfer index in
   boot_tree. Call Match_Tip_Numbers(ref_tree, boot_tree) first. */
void TBE_Add_Transfer_Index(t_tree *ref_tree, t_tree *boot_tree)
{
  int         i, n = ref_tree->n_otu, n_keys, mark;
  uint64_t   *key_a, *key_b, tot_a = 0, tot_b = 0;
  tbe_rooted  ref, boot;
  tbe_key    *keys;
  int        *stamp, *cnt, *upos, *light;

  assert(boot_tree->n_otu == n);

  key_a = (uint64_t *)mCalloc(2 * n, sizeof(uint64_t));
  key_b = (uint64_t *)mCalloc(2 * n, sizeof(uint64_t));
  for (i = 0; i < n; i++)
  {
    key_a[i] = TBE_Mix(2 * (uint64_t)i + 1);
    key_b[i] = TBE_Mix(2 * (uint64_t)i + 2 + 0x5851F42D4C957F2DULL);
    tot_a ^= key_a[i];
    tot_b ^= key_b[i];
  }

  for (i = 0; i < n; i++)
  {
    assert(ref_tree->a_nodes[i]->tax && ref_tree->a_nodes[i]->num >= 0 &&
           ref_tree->a_nodes[i]->num < n);
    assert(boot_tree->a_nodes[i]->tax && boot_tree->a_nodes[i]->num >= 0 &&
           boot_tree->a_nodes[i]->num < n);
  }

  TBE_Root(ref_tree, &ref, key_a, key_b);
  TBE_Root(boot_tree, &boot, key_a, key_b);

  /* splits of the bootstrap tree */
  keys   = (tbe_key *)mCalloc(2 * n, sizeof(tbe_key));
  n_keys = 0;
  for (i = 0; i < 2 * n; i++)
    if (boot.node[i] && i != boot.root)
      keys[n_keys++] = TBE_Canonical(boot.ha[i], boot.hb[i], tot_a, tot_b);
  qsort(keys, n_keys, sizeof(tbe_key), TBE_Key_Cmp);

  stamp = (int *)mCalloc(2 * n, sizeof(int));
  cnt   = (int *)mCalloc(2 * n, sizeof(int));
  upos  = (int *)mCalloc(2 * n, sizeof(int));
  light = (int *)mCalloc(n, sizeof(int));
  mark  = 0;

  for (i = 0; i < 2 * n - 3; i++)
  {
    t_edge *e = ref_tree->a_edges[i];
    int     c, s, p, j, k;
    tbe_key key;

    /* the side of e away from the root of ref */
    c = (ref.par[e->rght->num] == e->left->num) ? e->rght->num : e->left->num;
    assert(ref.par[c] >= 0);
    s = ref.size[c];
    p = MIN(s, n - s);
    if (p <= 1) continue; /* terminal branch: always in boot_tree */

    key = TBE_Canonical(ref.ha[c], ref.hb[c], tot_a, tot_b);
    if (bsearch(&key, keys, n_keys, sizeof(tbe_key), TBE_Key_Cmp)) continue;

    /* light side: the clade below c, or its complement */
    k = 0;
    if (s <= n - s)
      for (j = ref.lo[c]; j < ref.hi[c]; j++) light[k++] = ref.leaf[j];
    else
    {
      for (j = 0; j < ref.lo[c]; j++) light[k++] = ref.leaf[j];
      for (j = ref.hi[c]; j < n; j++) light[k++] = ref.leaf[j];
    }
    assert(k == p);

    e->tdist_score +=
        TBE_Transfer_Index(&boot, light, p, stamp, ++mark, cnt, upos);
  }

  Free(key_a);
  Free(key_b);
  Free(keys);
  Free(stamp);
  Free(cnt);
  Free(upos);
  Free(light);
  TBE_Free_Rooted(&ref);
  TBE_Free_Rooted(&boot);
}
