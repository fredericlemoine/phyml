/* Runner for the tinytest framework: see tinytest.h. */

#include "tinytest.h"

#include <signal.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sys/wait.h>
#include <unistd.h>

#define TT_MAX_TESTS 1024
#define TT_FAIL_CODE 97
#define TT_TIMEOUT_S 60

typedef struct
{
  const char *group, *name, *known_bug;
  tt_fn       fn;
} tt_test;

static tt_test tt_tests[TT_MAX_TESTS];
static int     tt_n_tests;

void tt_register(const char *group, const char *name, tt_fn fn, const char *known_bug)
{
  if (tt_n_tests == TT_MAX_TESTS)
  {
    fprintf(stderr, "tinytest: too many tests\n");
    exit(2);
  }
  tt_tests[tt_n_tests++] = (tt_test){group, name, known_bug, fn};
}

void tt_fail(const char *file, int line, const char *fmt, ...)
{
  va_list ap;
  fprintf(stderr, "    %s:%d: ", file, line);
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);
  fprintf(stderr, "\n");
  fflush(NULL);
  _exit(TT_FAIL_CODE);
}

static int tt_cmp(const void *a, const void *b)
{
  const tt_test *x = a, *y = b;
  int c = strcmp(x->group, y->group);
  return c ? c : strcmp(x->name, y->name);
}

/* Returns 1 when the test passed. */
static int tt_run_one(const tt_test *t, char *why, size_t why_len)
{
  int   status;
  pid_t pid;

  fflush(NULL);
  pid = fork();
  if (pid == 0)
  {
    /* PhyML prints a lot; keep the report readable. */
    if (!freopen("/dev/null", "w", stdout)) _exit(2);
    alarm(TT_TIMEOUT_S);
    srand(1);
    t->fn();
    fflush(NULL);
    _exit(0);
  }
  if (pid < 0 || waitpid(pid, &status, 0) < 0)
  {
    snprintf(why, why_len, "could not run test");
    return 0;
  }
  if (WIFEXITED(status) && WEXITSTATUS(status) == 0) return 1;
  if (WIFEXITED(status) && WEXITSTATUS(status) == TT_FAIL_CODE)
    snprintf(why, why_len, "check failed");
  else if (WIFEXITED(status))
    snprintf(why, why_len, "code under test called exit(%d)", WEXITSTATUS(status));
  else if (WIFSIGNALED(status) && WTERMSIG(status) == SIGALRM)
    snprintf(why, why_len, "timed out after %d s", TT_TIMEOUT_S);
  else if (WIFSIGNALED(status))
    snprintf(why, why_len, "killed by signal %d (%s)", WTERMSIG(status),
             strsignal(WTERMSIG(status)));
  return 0;
}

int main(int argc, char **argv)
{
  int  i, n_pass = 0, n_fail = 0, n_xfail = 0, n_xpass = 0, n_run = 0;
  char why[256];

  qsort(tt_tests, tt_n_tests, sizeof(tt_test), tt_cmp);

  for (i = 0; i < tt_n_tests; i++)
  {
    const tt_test *t = &tt_tests[i];
    char           full[256];
    int            ok;

    snprintf(full, sizeof(full), "%s.%s", t->group, t->name);
    if (argc > 1 && !strstr(full, argv[1])) continue;
    n_run++;

    why[0] = '\0';
    ok     = tt_run_one(t, why, sizeof(why));

    if (t->known_bug == NULL)
    {
      if (ok) n_pass++, printf("PASS   %s\n", full);
      else n_fail++, printf("FAIL   %s: %s\n", full, why);
    }
    else
    {
      if (ok) n_xpass++, printf("XPASS  %s (known bug seems fixed: %s)\n", full, t->known_bug);
      else n_xfail++, printf("XFAIL  %s: %s [%s]\n", full, t->known_bug, why);
    }
    fflush(stdout);
  }

  printf("\n%d tests: %d passed, %d failed, %d expected failures (known bugs), "
         "%d unexpected passes\n",
         n_run, n_pass, n_fail, n_xfail, n_xpass);
  return (n_fail || n_xpass) ? 1 : 0;
}
