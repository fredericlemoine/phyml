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

TEST(lib, failed_assert_comes_back_as_internal_error)
{
  jmp_buf      env;
  volatile int code;
  int          zero = 0;

  Silence();
  PhyML_Error_Clear();
  code = setjmp(env);
  if (code == 0)
  {
    PhyML_Error_Set_Recovery(&env);
    assert(zero == 1);
    CHECK(0); // not reached
  }
  PhyML_Error_Set_Recovery(NULL);
  CHECK_EQ_INT(code, PHYML_ERR_INTERNAL);
  CHECK(strstr(PhyML_Error_Message(), "Assertion failed: (zero == 1)") != NULL);
}

/* What a run that fails had allocated and opened is released. */
TEST(lib, failed_run_releases_its_blocks_and_files)
{
  jmp_buf      env;
  volatile int code;
  long         n_blocks = -1, n_files = -1;
  char        *argv[]   = {"phyml", "-i", NULL, "-b", "0", NULL};

  Silence();
  argv[2] = strdup(Stage("bad_short_seq.phy"));
  optind  = 1;
  PhyML_Error_Clear();
  PhyML_Track_Begin();
  code = setjmp(env);
  if (code == 0)
  {
    PhyML_Error_Set_Recovery(&env);
    PhyML_Main(5, argv);
  }
  PhyML_Error_Set_Recovery(NULL);
  CHECK_EQ_INT(code, PHYML_ERR_INPUT);
  PhyML_Track_End(&n_blocks, &n_files);
  CHECK(n_blocks > 10); // the options, the model, the sequences read so far
  CHECK(n_files >= 3);  // the alignment, the tree and the statistics files

  // Nothing is recorded outside a run
  PhyML_Track_Begin();
  PhyML_Track_End(&n_blocks, &n_files);
  CHECK_EQ_INT(n_blocks, 0);
  CHECK_EQ_INT(n_files, 0);
}

/* 400 failed runs would exhaust the file descriptors (256 by default on
   macOS, 1,024 on Linux with 3 files per run) if their files stayed open. */
TEST(lib, many_failed_runs_do_not_use_up_file_descriptors)
{
  int   i, fd_before, fd_after;
  char *aln;

  Silence();
  aln       = strdup(Stage("bad_short_seq.phy"));
  fd_before = dup(0);
  close(fd_before);
  for (i = 0; i < 400; i++) CHECK_EQ_INT(Run(aln, NULL), PHYML_ERR_INPUT);
  fd_after = dup(0);
  close(fd_after);
  CHECK_EQ_INT(fd_after, fd_before);
  CHECK_EQ_INT(Run(Stage("dna_small.phy"), NULL), PHYML_OK);
}

/* Blocks freed or reallocated during a run are no longer recorded. */
TEST(lib, tracking_follows_free_and_realloc)
{
  long  n_blocks = -1, n_files = -1;
  int   i;
  void *p[2000];
  FILE *fp;

  PhyML_Track_Begin();
  for (i = 0; i < 2000; i++) p[i] = mCalloc(1 + i % 7, sizeof(int));
  for (i = 0; i < 2000; i += 2) Free(p[i]);
  for (i = 1; i < 2000; i += 4) p[i] = mRealloc(p[i], 5000, sizeof(int));
  fp = fopen("/dev/null", "r");
  CHECK(fp != NULL);
  fclose(fp);
  fp = fopen("/dev/null", "r");
  PhyML_Track_End(&n_blocks, &n_files);
  CHECK_EQ_INT(n_blocks, 1000);
  CHECK_EQ_INT(n_files, 1);
}

/* Error messages printed just before an exit go to stderr, so they are in
   the error text (they were printed on stdout). */
TEST(lib, option_error_text_is_returned)
{
  Silence();
  CHECK_EQ_INT(Run(Stage("dna_small.phy"), "-mXYZ"), PHYML_ERR_USAGE);
  CHECK(strstr(PhyML_Error_Message(), "model name is incorrect") != NULL);
}

TEST(lib, version_and_help_are_not_errors)
{
  char *version[] = {"phyml", "--version", NULL};
  char *help[]    = {"phyml", "--help", NULL};

  Silence();
  CHECK_EQ_INT(PhyML_Run(2, version), PHYML_OK);
  CHECK_EQ_INT(PhyML_Run(2, help), PHYML_OK);
  CHECK_STR_EQ(PhyML_Error_Message(), "");
}

/* A run with threads gives the same tree as without, also when PhyML is
   called as a library, where each parallel region sets an error guard. */
TEST(lib, threaded_run_matches_serial_run)
{
  char  out[600], tree1[1 << 16];
  char *argv[] = {"phyml", "-i", NULL, "-u", NULL, "-m", "GTR", "-o", "lr",
                  "-b", "0", "--r_seed", "5", "--threads", NULL, NULL};

  Silence();
  argv[2] = strdup(Stage("threads_dna.phy"));
  argv[4] = strdup(Stage("threads_dna.true.nwk"));
  snprintf(out, sizeof(out), "%s_phyml_tree.txt", argv[2]);

  argv[14] = "1";
  CHECK_EQ_INT(PhyML_Run(15, argv), PHYML_OK);
  strcpy(tree1, Slurp(out));
  argv[14] = "4";
  CHECK_EQ_INT(PhyML_Run(15, argv), PHYML_OK);
  CHECK_STR_EQ(Slurp(out), tree1);
}

#ifdef _OPENMP
#include <omp.h>

/* An error in one thread of a parallel region ends that thread's work; the
   region ends normally and the error is raised after it. */
/* Not a local of the test: it is changed between setjmp and longjmp */
static int n_done;

static int Parallel_Region_With_Error(int failing_thread, int use_assert)
{
  n_done = 0;
  // A shared counter, not a reduction variable: a local changed after the
  // error point may be stored before it (its value is then indeterminate
  // after the longjmp, as GCC showed)
#pragma omp parallel num_threads(4)
  {
    PHYML_THREAD_TRY
    if (omp_get_thread_num() == failing_thread)
    {
      if (use_assert)
        assert(failing_thread < 0);
      else
        Exit("\n. Error in a thread.\n");
    }
#pragma omp atomic
    n_done++;
    PHYML_THREAD_END
  }
  PhyML_Thread_Errors_Raise();
  return 0;
}

TEST(lib, error_in_a_thread_comes_back_after_the_region)
{
  jmp_buf      env;
  volatile int code, round;
  int          n_threads = 0;

  Silence();
  for (round = 0; round < 2; round++)
  {
    PhyML_Error_Clear();
    PhyML_Error_Scope(PHYML_ERR_INPUT);
    n_done = -1;
    code   = setjmp(env);
    if (code == 0)
    {
      PhyML_Error_Set_Recovery(&env);
      Parallel_Region_With_Error(2, round == 1);
      CHECK(0); // not reached
    }
    PhyML_Error_Set_Recovery(NULL);
    CHECK_EQ_INT(code, round == 0 ? PHYML_ERR_INPUT : PHYML_ERR_INTERNAL);
    CHECK_EQ_INT(n_done, 3); // the other threads finished their work
    CHECK(strstr(PhyML_Error_Message(), round == 0 ? "Error in a thread"
                                                   : "Assertion failed") != NULL);
    CHECK_EQ_INT(omp_get_level(), 0); // the region was left normally
  }

  // The OpenMP runtime still gives the next region its threads
#pragma omp parallel num_threads(4)
  {
#pragma omp atomic
    n_threads++;
  }
  CHECK_EQ_INT(n_threads, 4);
}

/* An error in an inactive region (one thread) is also raised after it. */
TEST(lib, error_in_an_inactive_region_comes_back_after_it)
{
  jmp_buf      env;
  volatile int code;
  int          nt = 1;

  Silence();
  PhyML_Error_Clear();
  code = setjmp(env);
  if (code == 0)
  {
    PhyML_Error_Set_Recovery(&env);
#pragma omp parallel if (nt > 1) num_threads(nt)
    {
      PHYML_THREAD_TRY
      Exit("\n. Error in an inactive region.\n");
      PHYML_THREAD_END
    }
    PhyML_Thread_Errors_Raise();
    CHECK(0);
  }
  PhyML_Error_Set_Recovery(NULL);
  CHECK_EQ_INT(code, PHYML_ERR);
  CHECK_EQ_INT(omp_get_level(), 0);
}
#endif
