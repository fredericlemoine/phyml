/* Unit tests for small helpers in utilities.c and stats.c. */

#include "tinytest.h"
#include "utilities.h"
#include "stats.h"

TEST(states, nucleotides_map_to_indices)
{
  char c[1];
  const char *nt = "ACGTU";
  int expected[] = {0, 1, 2, 3, 3};
  for (int i = 0; i < 5; i++)
  {
    c[0] = nt[i];
    CHECK_EQ_INT(Assign_State(c, NT, 1), expected[i]);
    CHECK_EQ_INT(Is_Ambigu(c, NT, 1), NO);
  }
}

TEST(states, nucleotide_ambiguities)
{
  const char *amb = "RYKMSWBDHVN-?X";
  for (int i = 0; amb[i]; i++)
  {
    char c[1] = {amb[i]};
    CHECK_EQ_INT(Is_Ambigu(c, NT, 1), YES);
    CHECK_EQ_INT(Assign_State(c, NT, 1), -1);
  }
}

TEST(states, amino_acids_map_to_paml_order)
{
  const char *aa = "ARNDCQEGHILKMFPSTWYV";
  for (int i = 0; i < 20; i++)
  {
    char c[1] = {aa[i]};
    CHECK_EQ_INT(Assign_State(c, AA, 1), i);
    CHECK_EQ_INT(Is_Ambigu(c, AA, 1), NO);
  }
  char x[1] = {'X'};
  CHECK_EQ_INT(Is_Ambigu(x, AA, 1), YES);
}

TEST(strings, uppercase)
{
  char s[] = "acgtn-?";
  for (int i = 0; s[i]; i++) Uppercase(&s[i]);
  CHECK_STR_EQ(s, "ACGTN-?");
}

TEST(strings, basename)
{
  char a[] = "/tmp/dir/file.phy", b[] = "file.phy", c[] = "/tmp/dir/";
  CHECK_STR_EQ(Basename(a), "file.phy");
  CHECK_STR_EQ(Basename(b), "file.phy");
  CHECK_STR_EQ(Basename(c), "dir");
}

TEST(numbers, are_equal)
{
  CHECK(Are_Equal(1.0, 1.0 + 1e-12, 1e-9));
  CHECK(!Are_Equal(1.0, 1.1, 1e-9));
}

TEST(linear_algebra, matrix_inverse)
{
  phydbl x[4] = {4, 7, 2, 6};
  CHECK_EQ_INT(Matinv(x, 2, 2, NO), 1); /* 1 = success */
  CHECK_NEAR(x[0], 0.6, 1e-12);
  CHECK_NEAR(x[1], -0.7, 1e-12);
  CHECK_NEAR(x[2], -0.2, 1e-12);
  CHECK_NEAR(x[3], 0.4, 1e-12);
}

TEST(linear_algebra, matrix_inverse_times_matrix_is_identity)
{
  phydbl a[9] = {2, -1, 0, -1, 2, -1, 0, -1, 2}, inv[9];
  memcpy(inv, a, sizeof(a));
  Matinv(inv, 3, 3, NO);
  phydbl *p = Matrix_Mult(a, inv, 3, 3, 3, 3);
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) CHECK_NEAR(p[i * 3 + j], i == j ? 1.0 : 0.0, 1e-12);
  Free(p);
}
