/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

/* Public interface of libphyml.

   The library holds everything but main(): data structures, models,
   likelihood, tree search and branch supports. This header gives the
   entry points that do not end the calling process on an error. The
   internal headers (utilities.h and the others) are available to programs
   that work on trees and models directly; errors raised from their
   functions also come back to the enclosing PhyML_Run().

   Example:

     #include <phyml.h>

     char *args[] = {"phyml", "-i", "seq.phy", "-m", "GTR", "--r_seed", "1"};
     int   rc     = PhyML_Run(7, args);
     if (rc != PHYML_OK)
       fprintf(stderr, "%s: %s", PhyML_Error_Name(rc), PhyML_Error_Message());
*/

#ifndef PHYML_H
#define PHYML_H

#ifdef __cplusplus
extern "C"
{
#endif

/* Error codes returned by PhyML_Run() */
enum
{
  PHYML_OK           = 0,
  PHYML_ERR          = 1, /* error with no more precise code */
  PHYML_ERR_USAGE    = 2, /* command line or menu: unknown or incompatible options */
  PHYML_ERR_INPUT    = 3, /* alignment, tree or XML input that cannot be read */
  PHYML_ERR_MEMORY   = 4, /* allocation failure */
  PHYML_ERR_INTERNAL = 5, /* internal check failed ("Err. in file ... line ...") */
  PHYML_ERR_BUSY     = 6  /* PhyML_Run() called while another run is in progress */
};

/* Runs PhyML as the command line does: argv[0] is the program name, the
   other arguments are PhyML's options. Returns PHYML_OK, or the code of the
   error that stopped the run; the process is not terminated. Output files
   are written as by the phyml program. Not reentrant: one run at a time.

   Limits of this first version: the memory and the files held by a run
   that fails are not released; an error raised inside an OpenMP parallel
   region, and a failed assert(), still end the process; some errors wait
   for a key on stdin ("Type any key to exit") before they come back; the
   MPI build initialises and finalises MPI in each run, so it can run once. */
int PhyML_Run(int argc, char **argv);

/* The same analysis without the recovery point: an error ends the process
   with status 1. This is the main() of the phyml program. */
int PhyML_Main(int argc, char **argv);

/* Code and text of the error that stopped the last PhyML_Run(): the text is
   what PhyML wrote to stderr during that run (its last 4 KB), empty if
   nothing was written. */
int         PhyML_Error_Code(void);
const char *PhyML_Error_Message(void);

/* Name of an error code, e.g. "PHYML_ERR_INPUT" */
const char *PhyML_Error_Name(int code);

#ifdef __cplusplus
}
#endif

#endif
