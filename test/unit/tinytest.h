/*
  Minimal unit-test framework for the PhyML test suite.

  TEST(group, name) { ... }        defines a test, registered automatically.
  KNOWN_BUG(group, name, "why")    a test of correct behaviour that currently
                                   fails: reported XFAIL, or XPASS once fixed.

  Each test runs in its own child process, so a PhyML routine that calls
  exit(), aborts on an assert or crashes only fails that test.

  Checks: CHECK(cond), CHECK_EQ_INT(a, b), CHECK_NEAR(a, b, tol),
          CHECK_STR_EQ(a, b). A failed check reports and ends the test.
*/

#ifndef TINYTEST_H
#define TINYTEST_H

#include <math.h>
#include <stdio.h>
#include <string.h>

typedef void (*tt_fn)(void);

void tt_register(const char *group, const char *name, tt_fn fn, const char *known_bug);
void tt_fail(const char *file, int line, const char *fmt, ...);

#define TT_REGISTER(group, name, bug)                                      \
  static void tt_##group##_##name(void);                                   \
  __attribute__((constructor)) static void tt_reg_##group##_##name(void)   \
  {                                                                        \
    tt_register(#group, #name, tt_##group##_##name, bug);                  \
  }                                                                        \
  static void tt_##group##_##name(void)

#define TEST(group, name) TT_REGISTER(group, name, NULL)
#define KNOWN_BUG(group, name, why) TT_REGISTER(group, name, why)

#define CHECK(cond)                                                        \
  do {                                                                     \
    if (!(cond)) tt_fail(__FILE__, __LINE__, "CHECK(%s)", #cond);          \
  } while (0)

#define CHECK_EQ_INT(a, b)                                                 \
  do {                                                                     \
    long long tt_a = (long long)(a), tt_b = (long long)(b);                \
    if (tt_a != tt_b)                                                      \
      tt_fail(__FILE__, __LINE__, "%s == %s (%lld vs %lld)", #a, #b, tt_a, \
              tt_b);                                                       \
  } while (0)

#define CHECK_NEAR(a, b, tol)                                              \
  do {                                                                     \
    double tt_a = (double)(a), tt_b = (double)(b), tt_t = (double)(tol);   \
    if (!(fabs(tt_a - tt_b) <= tt_t))                                      \
      tt_fail(__FILE__, __LINE__, "%s ~= %s (%.12g vs %.12g, tol %g)", #a, \
              #b, tt_a, tt_b, tt_t);                                       \
  } while (0)

#define CHECK_STR_EQ(a, b)                                                 \
  do {                                                                     \
    const char *tt_a = (a), *tt_b = (b);                                   \
    if (!tt_a || !tt_b || strcmp(tt_a, tt_b))                              \
      tt_fail(__FILE__, __LINE__, "%s == %s (\"%s\" vs \"%s\")", #a, #b,   \
              tt_a ? tt_a : "(null)", tt_b ? tt_b : "(null)");             \
  } while (0)

#endif
