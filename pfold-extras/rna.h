/*********************************************************************

  rna.h

  Contains general functions for handling RNA.

*********************************************************************/

#ifndef __rna_h__
#define __rna_h__

#include <ctype.h>

int IsRNAnuc(char nuc);
int IsRNAnucExt(char nuc);
int StdPair6(char nuc1, char nuc2);
int StdPair4(char nuc1, char nuc2);

#endif
