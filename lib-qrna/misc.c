/********************************************************************************************************
 * QRNA - Comparative analysis of biological sequences 
 *         with pair hidden Markov models, pair stochastic context-free
 *        grammars, and probabilistic evolutionary  models.
 *       
 * Version 2.0.0 (JUN 2003)
 *
 * Copyright (C) 2000-2003 Howard Hughes Medical Institute/Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ***********************************************************************************************************/

/* matrix.c
 *
 * ER, Fri Nov 22 11:51:28 CST 2002 [St. Louis]
 * 
 * Purpose:  Includes all the functions for matrix manipulations
 *           used in constructing the evolutionary models.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>

#include "misc.h"
#include "matrix.h"
#include "evolve.h"
#include "mat_globals.h"
#include "squid.h"
#include "mat_structs.h"

/* Function: PrintProbs()
 * Date:     ER, Tue Jun 13 13:42:27 CDT 2000 [St. Louis]
 *
 * Purpose:  Print probabilities
 *
 * Args:     othmodel -- the othmodel prob's, in log2 form
 *
 * Returns:  void. prints transition and emission probs, in [0,1] form.
 */
void
PrintProbs(FILE *ofp, double *prob, int L) 
{
  int x, y;

  for (x = 0; x < L; x++) {
    for (y = 0; y < L; y++) {
      fprintf(ofp, "%.4f ", prob[x*L+y]);
    }
    fprintf(ofp, "\n");
  }
}

void
PrintProbs2D(FILE *ofp, double **prob, int Lx, int Ly) 
{
  int x, y;

  for (x = 0; x < Lx; x++) {
    for (y = 0; y < Ly; y++) {
      fprintf(ofp, "%.4f ", prob[x][y]);
    }
    fprintf(ofp, "\n");
  }
}

void
PrintVectorProbs(FILE *ofp, double *prob, int L) 
{
  int x;

  for (x = 0; x < L; x++) 
    fprintf(ofp, "%.4f ", prob[x]);
  
  fprintf(ofp, "\n");
  
}

