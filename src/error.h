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

/* Forgets the last error, and keeps the text written to stderr */
void PhyML_Error_Clear(void);
void PhyML_Error_Capture(const char *format, va_list ap);

#endif
