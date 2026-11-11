/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Error handling of libphyml: see error.h. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "error.h"

#define PHYML_ERROR_MSG_SIZE 4096

static jmp_buf *error_recovery = NULL;
static int      error_scope    = PHYML_ERR;
static int      error_code     = PHYML_OK;
static char     error_msg[PHYML_ERROR_MSG_SIZE];
static size_t   error_msg_len  = 0;

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* longjmp must not leave an OpenMP parallel region */

static int Error_Can_Recover(void)
{
  if (error_recovery == NULL) return 0;
#ifdef _OPENMP
  if (omp_in_parallel()) return 0;
#endif
  return 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PhyML_Raise(int code)
{
  if (code == PHYML_OK) code = PHYML_ERR;
  fflush(NULL);
  if (Error_Can_Recover())
  {
    error_code = code;
    longjmp(*error_recovery, code);
  }
  exit(1);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PhyML_Stop(int status)
{
  fflush(NULL);
  if (Error_Can_Recover()) longjmp(*error_recovery, PHYML_STOP_REQUESTED);
  exit(status);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int PhyML_Error_Scope(int code)
{
  int previous = error_scope;
  error_scope  = code;
  return previous;
}

int PhyML_Error_Scope_Code(void) { return error_scope; }

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int PhyML_Error_Has_Recovery(void) { return error_recovery != NULL; }

void PhyML_Error_Set_Recovery(jmp_buf *env) { error_recovery = env; }

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PhyML_Error_Clear(void)
{
  error_code    = PHYML_OK;
  error_scope   = PHYML_ERR;
  error_msg_len = 0;
  error_msg[0]  = '\0';
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Appends what is written to stderr to the error text, of which the last
   PHYML_ERROR_MSG_SIZE - 1 characters are kept. Only while a recovery point
   is set: the phyml program has no use for it. */

void PhyML_Error_Capture(const char *format, va_list ap)
{
  char   chunk[1024];
  int    n;
  size_t len, drop;

  if (error_recovery == NULL) return;

  n = vsnprintf(chunk, sizeof(chunk), format, ap);
  if (n <= 0) return;
  len = ((size_t)n < sizeof(chunk)) ? (size_t)n : sizeof(chunk) - 1;

#ifdef _OPENMP
#pragma omp critical(phyml_error_capture)
#endif
  {
    if (error_msg_len + len > PHYML_ERROR_MSG_SIZE - 1)
    {
      drop = error_msg_len + len - (PHYML_ERROR_MSG_SIZE - 1);
      if (drop > error_msg_len) drop = error_msg_len;
      memmove(error_msg, error_msg + drop, error_msg_len - drop);
      error_msg_len -= drop;
    }
    memcpy(error_msg + error_msg_len, chunk, len);
    error_msg_len += len;
    error_msg[error_msg_len] = '\0';
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

int PhyML_Error_Code(void) { return error_code; }

const char *PhyML_Error_Message(void) { return error_msg; }

const char *PhyML_Error_Name(int code)
{
  switch (code)
  {
  case PHYML_OK: return "PHYML_OK";
  case PHYML_ERR: return "PHYML_ERR";
  case PHYML_ERR_USAGE: return "PHYML_ERR_USAGE";
  case PHYML_ERR_INPUT: return "PHYML_ERR_INPUT";
  case PHYML_ERR_MEMORY: return "PHYML_ERR_MEMORY";
  case PHYML_ERR_INTERNAL: return "PHYML_ERR_INTERNAL";
  case PHYML_ERR_BUSY: return "PHYML_ERR_BUSY";
  }
  return "PHYML_ERR_UNKNOWN";
}
