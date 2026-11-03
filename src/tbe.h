/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#ifndef TBE_H
#define TBE_H

#include "utilities.h"

#include <stdint.h>

/* Adds to the tdist_score of each branch of ref_tree its transfer index in
   boot_tree (Match_Tip_Numbers first). O(n) memory. */
void TBE_Add_Transfer_Index(t_tree *ref_tree, t_tree *boot_tree);
#endif
