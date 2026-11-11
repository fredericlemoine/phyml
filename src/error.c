/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Error handling of libphyml: see error.h. */

#include <stdint.h>
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

static void Error_Capture_Text(const char *text, size_t len);

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

void PhyML_Assert_Failed(const char *expr, const char *file, int line,
                         const char *function)
{
  char text[1024];

  // Same text as <assert.h>, on stderr and in the error text
  snprintf(text, sizeof(text),
           "\nAssertion failed: (%s), function %s, file %s, line %d.\n", expr,
           function, file, line);
  fputs(text, stderr);
  Error_Capture_Text(text, strlen(text));
  fflush(NULL);
  if (Error_Can_Recover())
  {
    error_code = PHYML_ERR_INTERNAL;
    longjmp(*error_recovery, PHYML_ERR_INTERNAL);
  }
  abort();
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

/* Appends to the error text, of which the last PHYML_ERROR_MSG_SIZE - 1
   characters are kept. Only while a recovery point is set: the phyml
   program has no use for it. */

static void Error_Capture_Text(const char *text, size_t len)
{
  size_t drop;

  if (error_recovery == NULL || len == 0) return;
  if (len > PHYML_ERROR_MSG_SIZE - 1)
  {
    text += len - (PHYML_ERROR_MSG_SIZE - 1);
    len = PHYML_ERROR_MSG_SIZE - 1;
  }

#ifdef _OPENMP
#pragma omp critical(phyml_error_capture)
#endif
  {
    if (error_msg_len + len > PHYML_ERROR_MSG_SIZE - 1)
    {
      drop = error_msg_len + len - (PHYML_ERROR_MSG_SIZE - 1);
      memmove(error_msg, error_msg + drop, error_msg_len - drop);
      error_msg_len -= drop;
    }
    memcpy(error_msg + error_msg_len, text, len);
    error_msg_len += len;
    error_msg[error_msg_len] = '\0';
  }
}

/* What PhyML_Fprintf writes to stderr */

void PhyML_Error_Capture(const char *format, va_list ap)
{
  char chunk[1024];
  int  n;

  if (error_recovery == NULL) return;

  n = vsnprintf(chunk, sizeof(chunk), format, ap);
  if (n <= 0) return;
  Error_Capture_Text(chunk, ((size_t)n < sizeof(chunk)) ? (size_t)n
                                                       : sizeof(chunk) - 1);
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Resources of a run (see error.h). A set of pointers: open addressing,
   linear probing, deletion by backward shift. It is only filled between
   PhyML_Track_Begin() and PhyML_Track_End(), so the phyml program, which
   never calls them, pays one test per allocation. */

typedef struct
{
  void **slot;
  size_t cap, n; /* cap is 0 or a power of two */
} ptr_set;

static int     track_on = 0;
static ptr_set track_blocks, track_files;

static size_t Set_Hash(const void *p, size_t cap)
{
  uint64_t x = (uint64_t)(uintptr_t)p;
  x ^= x >> 33;
  x *= 0xFF51AFD7ED558CCDULL;
  x ^= x >> 33;
  return (size_t)x & (cap - 1);
}

static int Set_Grow(ptr_set *s)
{
  size_t  i, j, cap = s->cap ? 2 * s->cap : 1024;
  void  **slot = (void **)calloc(cap, sizeof(void *));

  if (slot == NULL) return 0;
  for (i = 0; i < s->cap; i++)
    if (s->slot[i])
    {
      j = Set_Hash(s->slot[i], cap);
      while (slot[j]) j = (j + 1) & (cap - 1);
      slot[j] = s->slot[i];
    }
  free(s->slot);
  s->slot = slot;
  s->cap  = cap;
  return 1;
}

/* A pointer that cannot be recorded (no memory for the set) is simply not
   released by PhyML_Track_End(). */
static void Set_Add(ptr_set *s, void *p)
{
  size_t j;

  if (p == NULL) return;
  if (2 * (s->n + 1) > s->cap && !Set_Grow(s)) return;
  j = Set_Hash(p, s->cap);
  while (s->slot[j])
  {
    if (s->slot[j] == p) return;
    j = (j + 1) & (s->cap - 1);
  }
  s->slot[j] = p;
  s->n++;
}

static void Set_Remove(ptr_set *s, void *p)
{
  size_t i, j, k;

  if (p == NULL || s->cap == 0) return;
  i = Set_Hash(p, s->cap);
  while (s->slot[i] != p)
  {
    if (s->slot[i] == NULL) return; // not recorded
    i = (i + 1) & (s->cap - 1);
  }
  // Backward shift: move up the entries that the hole would cut off
  j = i;
  for (;;)
  {
    j = (j + 1) & (s->cap - 1);
    if (s->slot[j] == NULL) break;
    k = Set_Hash(s->slot[j], s->cap);
    if ((i <= j) ? (k <= i || k > j) : (k <= i && k > j))
    {
      s->slot[i] = s->slot[j];
      i          = j;
    }
  }
  s->slot[i] = NULL;
  s->n--;
}

static void Track_Add(ptr_set *s, void *p)
{
  if (!track_on) return;
#ifdef _OPENMP
#pragma omp critical(phyml_track)
#endif
  Set_Add(s, p);
}

static void Track_Remove(ptr_set *s, void *p)
{
  if (!track_on) return;
#ifdef _OPENMP
#pragma omp critical(phyml_track)
#endif
  Set_Remove(s, p);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void PhyML_Track_Begin(void)
{
  track_blocks.n = track_files.n = 0;
  if (track_blocks.slot)
    memset(track_blocks.slot, 0, track_blocks.cap * sizeof(void *));
  if (track_files.slot)
    memset(track_files.slot, 0, track_files.cap * sizeof(void *));
  track_on = 1;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Closes the files and frees the blocks still recorded, and gives their
   numbers (either pointer may be NULL). */

void PhyML_Track_End(long *n_blocks, long *n_files)
{
  size_t i;

  track_on = 0;
  if (n_blocks) *n_blocks = (long)track_blocks.n;
  if (n_files) *n_files = (long)track_files.n;

  for (i = 0; i < track_files.cap; i++)
    if (track_files.slot[i]) fclose((FILE *)track_files.slot[i]);
  for (i = 0; i < track_blocks.cap; i++)
    if (track_blocks.slot[i]) free(track_blocks.slot[i]);

  free(track_files.slot);
  free(track_blocks.slot);
  memset(&track_files, 0, sizeof(track_files));
  memset(&track_blocks, 0, sizeof(track_blocks));
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void *PhyML_Malloc(size_t size)
{
  void *p = malloc(size);
  Track_Add(&track_blocks, p);
  return p;
}

void *PhyML_Calloc(size_t n, size_t size)
{
  void *p = calloc(n, size);
  Track_Add(&track_blocks, p);
  return p;
}

void *PhyML_Realloc(void *p, size_t size)
{
  void *q = realloc(p, size);
  if (q != NULL && q != p)
  {
    Track_Remove(&track_blocks, p);
    Track_Add(&track_blocks, q);
  }
  return q;
}

#ifndef WIN32
int PhyML_Memalign(void **p, size_t alignment, size_t size)
{
  int rc = posix_memalign(p, alignment, size);
  if (rc == 0) Track_Add(&track_blocks, *p);
  return rc;
}
#endif

void PhyML_Free_Block(void *p)
{
  Track_Remove(&track_blocks, p);
  free(p);
}

FILE *PhyML_Fopen(const char *name, const char *mode)
{
  FILE *fp = fopen(name, mode);
  Track_Add(&track_files, fp);
  return fp;
}

int PhyML_Fclose(FILE *fp)
{
  Track_Remove(&track_files, fp);
  return fclose(fp);
}
