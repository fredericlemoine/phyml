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

/* Forgets the last error, and keeps the text written to stderr */
void PhyML_Error_Clear(void);
void PhyML_Error_Capture(const char *format, va_list ap);

#endif
