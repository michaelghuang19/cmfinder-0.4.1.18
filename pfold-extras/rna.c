/*********************************************************************

  rna.c

  Contains general functions for handling RNA.

  001016 Bjarne Knudsen (bk@daimi.au.dk)

  Copyright (C) 2000 Bjarne Knudsen

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.

*********************************************************************/

#include "rna.h"

/*
   Returns true if nuc is a knwon nucleotide or '-'. Otherwise, false.
*/
int IsRNAnuc(char nuc)
{
  nuc = toupper(nuc);
  if (nuc == 'A' || nuc == 'C' ||
      nuc == 'G' || nuc == 'U' || nuc == '-')
    return 1;
  return 0;
}

/*
   Returns true if nuc is a knwon nucleotide, from the extended set,
   or '-'. Otherwise, false.

   Extended nucleotide set:

   A, C, G, U, Y, R, K, M, B, D, H, V, N, X
*/
int IsRNAnucExt(char nuc)
{
  nuc = toupper(nuc);
  if (nuc == 'A' || nuc == 'C' || nuc == 'G' || nuc == 'U' ||
      nuc == 'Y' || nuc == 'R' || nuc == 'K' || nuc == 'M' ||
      nuc == 'A' || nuc == 'C' ||  /* correct these */
      nuc == 'B' || nuc == 'D' || nuc == 'H' || nuc == 'V' ||
      nuc == 'N' || nuc == 'X' || nuc == '-')
    return 1;
  return 0;
}

/* Returns true if nuc1 and nuc 2 forms a standard
   base-pair. Otherwise, false */
int StdPair6(char nuc1, char nuc2)
{
  nuc1 = toupper(nuc1);
  if (nuc1 == 'T') nuc1 = 'U';
  nuc2 = toupper(nuc2);
  if (nuc2 == 'T') nuc2 = 'U';

  if ((nuc1 == 'A' && nuc2 == 'U') ||
      (nuc1 == 'U' && nuc2 == 'A') ||
      (nuc1 == 'C' && nuc2 == 'G') ||
      (nuc1 == 'G' && nuc2 == 'C') ||
      (nuc1 == 'G' && nuc2 == 'U') ||
      (nuc1 == 'U' && nuc2 == 'G'))
    return 1;
  return 0;
}

/* Returns true if nuc1 and nuc 2 form one of the 4 strongest standard
   base-pairs.  Otherwise, false */
int StdPair4(char nuc1, char nuc2)
{
  nuc1 = toupper(nuc1);
  if (nuc1 == 'T') nuc1 = 'U';
  nuc2 = toupper(nuc2);
  if (nuc2 == 'T') nuc2 = 'U';

  if ((nuc1 == 'A' && nuc2 == 'U') ||
      (nuc1 == 'U' && nuc2 == 'A') ||
      (nuc1 == 'C' && nuc2 == 'G') ||
      (nuc1 == 'G' && nuc2 == 'C'))
    return 1;
  return 0;
}
