/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Error handling of libphyml.

   PhyML reports an error by printing a message and calling Exit(),
   Warn_And_Exit() or Generic_Exit(): there is no error return path in the
   1,200 call sites. These functions now end in PhyML_Raise(), which

   - returns to the recovery point of PhyML_Run() (phyml.h), if there is
     one, with an error code, so that the caller of the library gets the
     code instead of losing its process;
   - ends the process with status 1 otherwise, as the phyml program always
     did.

   The code of an error is the one given to PhyML_Raise(), or for Exit()
   the code of the current scope (PhyML_Error_Scope), which says what the
   program is doing: reading the command line, reading input files, ... */

#ifndef ERROR_H
#define ERROR_H

#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>

#include "phyml.h"

#if defined(__GNUC__)
#define PHYML_NORETURN __attribute__((noreturn))
#else
#define PHYML_NORETURN
#endif

/* Ends the run with an error code (not PHYML_OK). Does not return. */
PHYML_NORETURN void PhyML_Raise(int code);

/* Ends the run on request, without error (e.g. after printing what an
   option asked for): PhyML_Run() returns PHYML_OK; without a recovery
   point the process exits with 'status'. Does not return. */
PHYML_NORETURN void PhyML_Stop(int status);

/* Sets the code that Exit() raises from now on and returns the previous
   one, to be restored at the end of the scope. */
int PhyML_Error_Scope(int code);
int PhyML_Error_Scope_Code(void);

/* Recovery point (used by PhyML_Run): errors raised while one is set
   longjmp to it with their code. */
int  PhyML_Error_Has_Recovery(void);
void PhyML_Error_Set_Recovery(jmp_buf *env);

/* Value passed to the recovery point by PhyML_Stop() */
#define PHYML_STOP_REQUESTED (-1)

/* A failed assert(): utilities.h defines assert() to call this. With a
   recovery point the run ends with PHYML_ERR_INTERNAL; without one the
   process aborts, as with the assert() of <assert.h>. */
PHYML_NORETURN void PhyML_Assert_Failed(const char *expr, const char *file,
                                        int line, const char *function);

/* Resources of a run. Between PhyML_Track_Begin() and PhyML_Track_End(),
   every block allocated and every file opened by the PhyML sources is
   recorded until it is freed or closed (utilities.h maps malloc, calloc,
   realloc, posix_memalign, free, fopen and fclose to the functions below).
   PhyML_Track_End() closes the files and frees the blocks that are left:
   those of a run stopped by an error, and what a complete run does not
   release. Outside these two calls the functions are the standard ones. */
void  PhyML_Track_Begin(void);
void  PhyML_Track_End(long *n_blocks, long *n_files);
void *PhyML_Malloc(size_t size);
void *PhyML_Calloc(size_t n, size_t size);
void *PhyML_Realloc(void *p, size_t size);
#ifndef WIN32
int   PhyML_Memalign(void **p, size_t alignment, size_t size);
#endif
void  PhyML_Free_Block(void *p);
FILE *PhyML_Fopen(const char *name, const char *mode);
int   PhyML_Fclose(FILE *fp);

/* Errors raised in OpenMP parallel regions. A thread cannot leave a
   parallel region with longjmp (nor an inactive one, which the runtime
   also tracks), so each thread of a region runs its work between
   PHYML_THREAD_TRY and PHYML_THREAD_END: an error raised in that work
   records its code and jumps to the end of the thread's work, the region
   ends normally, and PhyML_Thread_Errors_Raise(), called after the region,
   raises the recorded error. The guard is only set when PhyML is called as
   a library (PhyML_Run); in the phyml program an error in a thread exits
   the process as before, and the guard costs nothing. The work of a guarded
   thread must not contain an 'omp for' or a barrier, which a thread that
   jumped out would skip: give each thread its sites with Lk_Thread_Range.
   An error in a region without the guard still ends the process. */
#ifdef WIN32
typedef jmp_buf phyml_jmp_buf;
#define PHYML_SETJMP(env) setjmp(env)
#else
typedef sigjmp_buf phyml_jmp_buf; /* no signal mask: no system call */
#define PHYML_SETJMP(env) sigsetjmp(env, 0)
#endif

int  PhyML_Thread_Guard_Needed(void);
void PhyML_Thread_Guard(phyml_jmp_buf *env);
void PhyML_Thread_Errors_Raise_(void);

/* Without OpenMP there is no parallel region and the macros are empty: a
   function that calls setjmp is compiled with fewer optimisations (the
   likelihood loops of the default build were 3.5% slower). With OpenMP
   the region body is compiled as a separate function. */
#ifdef _OPENMP
#define PHYML_THREAD_TRY                                                   \
  {                                                                        \
    phyml_jmp_buf phyml_thread_env;                                        \
    const int     phyml_thread_guarded = PhyML_Thread_Guard_Needed();      \
    if (!phyml_thread_guarded || PHYML_SETJMP(phyml_thread_env) == 0)      \
    {                                                                      \
      if (phyml_thread_guarded) PhyML_Thread_Guard(&phyml_thread_env);

#define PHYML_THREAD_END                                                   \
    }                                                                      \
    if (phyml_thread_guarded) PhyML_Thread_Guard(NULL);                    \
  }
#define PhyML_Thread_Errors_Raise() PhyML_Thread_Errors_Raise_()
#else
#define PHYML_THREAD_TRY {
#define PHYML_THREAD_END }
#define PhyML_Thread_Errors_Raise() ((void)0)
#endif

/* For functions that return an error code instead of raising it: sets the
   code and the text of the last error (PhyML_Error_Code() and
   PhyML_Error_Message()), replacing the previous text, without printing
   anything, and returns the code. */
int PhyML_Error_Set(int code, const char *format, ...);

/* Prints the text of the last error on stderr and raises its code: for a
   caller that cannot go on after such a function failed. */
PHYML_NORETURN void PhyML_Error_Raise_Last(void);

/* Runs fn(arg) with its own recovery point: an error raised in it comes
   back as the returned code, with what it wrote on stderr as the error text
   (PhyML_Error_Message), instead of ending the run or the process. Exit()
   raises 'scope' meanwhile. A bridge for code that does not return codes
   yet: what fn allocated before the error is not released (PhyML_Run
   releases it at the end of the run). */
int PhyML_Try(void (*fn)(void *), void *arg, int scope);

/* Forgets the last error, and keeps the text written to stderr */
void PhyML_Error_Clear(void);
void PhyML_Error_Capture(const char *format, va_list ap);

#endif
