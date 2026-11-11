/* libphyml: PhyML_Run returns error codes instead of ending the process
   (src/phyml.h, src/error.h). The unit tests run from test/, each in its
   own process. */

#include <stdlib.h>
#include <unistd.h>

#include "tinytest.h"
#include "utilities.h"
#include "io.h"
#include "error.h"
#include "phyml.h"

static char lib_dir[256];

/* Copies data/<name> to a new temporary directory and returns its path
   there (PhyML writes its output files next to the alignment). */
static const char *Stage(const char *name)
{
  static char path[512];
  char        src[512], line[4096];
  FILE       *in, *out;
  size_t      n;

  if (lib_dir[0] == '\0')
  {
    const char *tmp = getenv("TMPDIR");
    snprintf(lib_dir, sizeof(lib_dir), "%s/phyml_lib_XXXXXX", tmp ? tmp : "/tmp");
    if (!mkdtemp(lib_dir)) tt_fail(__FILE__, __LINE__, "mkdtemp");
  }
  snprintf(src, sizeof(src), "data/%s", name);
  snprintf(path, sizeof(path), "%s/%s", lib_dir, name);
  in  = fopen(src, "r");
  out = fopen(path, "w");
  if (!in || !out) tt_fail(__FILE__, __LINE__, "cannot stage %s", src);
  while ((n = fread(line, 1, sizeof(line), in)) > 0) fwrite(line, 1, n, out);
  fclose(in);
  fclose(out);
  return path;
}

static char *Slurp(const char *path)
{
  static char buf[1 << 16];
  FILE       *fp = fopen(path, "r");
  size_t      n;

  if (!fp) tt_fail(__FILE__, __LINE__, "cannot read %s", path);
  n      = fread(buf, 1, sizeof(buf) - 1, fp);
  buf[n] = '\0';
  fclose(fp);
  return buf;
}

/* PhyML's own output is not part of the test report */
static void Silence(void)
{
  if (!freopen("/dev/null", "w", stdout) || !freopen("/dev/null", "w", stderr))
    tt_fail(__FILE__, __LINE__, "freopen");
}

static int Run(const char *aln, const char *extra)
{
  char *argv[12];
  int   argc = 0;

  argv[argc++] = "phyml";
  argv[argc++] = "-i";
  argv[argc++] = (char *)aln;
  argv[argc++] = "--r_seed";
  argv[argc++] = "3";
  argv[argc++] = "-b";
  argv[argc++] = "0";
  if (extra) argv[argc++] = (char *)extra;
  argv[argc] = NULL;
  return PhyML_Run(argc, argv);
}

TEST(lib, missing_alignment_is_a_usage_error)
{
  Silence();
  int rc = Run("data/no_such_file.phy", NULL);
  CHECK_EQ_INT(rc, PHYML_ERR_USAGE);
  CHECK_EQ_INT(PhyML_Error_Code(), PHYML_ERR_USAGE);
  CHECK(strstr(PhyML_Error_Message(), "does not exist") != NULL);
}

TEST(lib, unknown_option_is_a_usage_error)
{
  Silence();
  CHECK_EQ_INT(Run(Stage("dna_small.phy"), "--no_such_option"), PHYML_ERR_USAGE);
}

TEST(lib, malformed_alignment_is_an_input_error)
{
  Silence();
  CHECK_EQ_INT(Run(Stage("bad_short_seq.phy"), NULL), PHYML_ERR_INPUT);
  CHECK(PhyML_Error_Message()[0] != '\0');
}

TEST(lib, malformed_user_tree_is_an_input_error)
{
  char  tree[512];
  char *argv[] = {"phyml", "-i", NULL, "-u", tree, "-b", "0", NULL};

  Silence();
  snprintf(tree, sizeof(tree), "%s", Stage("bad_tree.nwk"));
  argv[2] = strdup(Stage("dna_small.phy"));
  CHECK_EQ_INT(PhyML_Run(7, argv), PHYML_ERR_INPUT);
}

TEST(lib, run_succeeds_and_writes_the_tree)
{
  char out[600];

  Silence();
  const char *aln = Stage("dna_small.phy");
  CHECK_EQ_INT(Run(aln, NULL), PHYML_OK);
  CHECK_EQ_INT(PhyML_Error_Code(), PHYML_OK);
  snprintf(out, sizeof(out), "%s_phyml_tree.txt", aln);
  CHECK(strchr(Slurp(out), '(') != NULL);
}

TEST(lib, run_after_an_error_and_repeated_runs_agree)
{
  char out[600], first[1 << 16];

  Silence();
  CHECK(Run("data/no_such_file.phy", NULL) != PHYML_OK);

  const char *aln = strdup(Stage("dna_small.phy"));
  snprintf(out, sizeof(out), "%s_phyml_tree.txt", aln);

  CHECK_EQ_INT(Run(aln, NULL), PHYML_OK);
  CHECK_STR_EQ(PhyML_Error_Message(), "");
  strcpy(first, Slurp(out));
  CHECK_EQ_INT(Run(aln, NULL), PHYML_OK);
  CHECK_STR_EQ(first, Slurp(out));
}

TEST(lib, internal_error_comes_back_with_its_code)
{
  jmp_buf      env;
  volatile int code;

  Silence();
  PhyML_Error_Clear();
  code = setjmp(env);
  if (code == 0)
  {
    PhyML_Error_Set_Recovery(&env);
    Generic_Exit(__FILE__, __LINE__, __FUNCTION__);
    CHECK(0); // not reached
  }
  PhyML_Error_Set_Recovery(NULL);
  CHECK_EQ_INT(code, PHYML_ERR_INTERNAL);
  CHECK(strstr(PhyML_Error_Message(), "finished prematurely") != NULL);
}

TEST(lib, exit_uses_the_code_of_the_scope)
{
  jmp_buf      env;
  volatile int code;

  Silence();
  PhyML_Error_Clear();
  code = setjmp(env);
  if (code == 0)
  {
    PhyML_Error_Set_Recovery(&env);
    CHECK_EQ_INT(PhyML_Error_Scope(PHYML_ERR_INPUT), PHYML_ERR);
    Exit("\n. bad input\n");
  }
  PhyML_Error_Set_Recovery(NULL);
  CHECK_EQ_INT(code, PHYML_ERR_INPUT);
  CHECK(strstr(PhyML_Error_Message(), "bad input") != NULL);
}

TEST(lib, stop_request_is_not_an_error)
{
  jmp_buf      env;
  volatile int code;

  PhyML_Error_Clear();
  code = setjmp(env);
  if (code == 0)
  {
    PhyML_Error_Set_Recovery(&env);
    PhyML_Stop(-1);
  }
  PhyML_Error_Set_Recovery(NULL);
  CHECK_EQ_INT(code, PHYML_STOP_REQUESTED);
  CHECK_EQ_INT(PhyML_Error_Code(), PHYML_OK);
}

TEST(lib, error_text_keeps_the_last_characters)
{
  jmp_buf env;
  int     i;

  Silence();
  PhyML_Error_Clear();
  PhyML_Error_Set_Recovery(&env); // capture is on while a recovery point is set
  for (i = 0; i < 2000; i++) PhyML_Fprintf(stderr, "line %d\n", i);
  PhyML_Fprintf(stderr, "the end");
  PhyML_Error_Set_Recovery(NULL);
  CHECK(strlen(PhyML_Error_Message()) < 4096);
  CHECK(strstr(PhyML_Error_Message(), "line 1999\nthe end") != NULL);
  CHECK(strstr(PhyML_Error_Message(), "line 0\n") == NULL);
}

TEST(lib, error_names)
{
  CHECK_STR_EQ(PhyML_Error_Name(PHYML_OK), "PHYML_OK");
  CHECK_STR_EQ(PhyML_Error_Name(PHYML_ERR_MEMORY), "PHYML_ERR_MEMORY");
  CHECK_STR_EQ(PhyML_Error_Name(1234), "PHYML_ERR_UNKNOWN");
}

TEST(lib, empty_allocation_is_not_a_failure)
{
  void *p = mCalloc(0, sizeof(int));
  CHECK(p != NULL);
  Free(p);
}
