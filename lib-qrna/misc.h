/************************************************************
 * QRNA - Comparative analysis of biological sequences 
 *        with pair hidden Markov models and pair stochastic context-free grammars
 * Copyright (C) 2000-Howard Hughes Medical Institute
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* evolve.h
 * ANSI prototypes for all external functions.
 * 
 */

#include <stdio.h>
#include "mat_structs.h"

extern void     PrintProbs(FILE *ofp, double *prob, int L);
extern void     PrintProbs2D(FILE *ofp, double **prob, int Lx, int Ly);
extern void     PrintVectorProbs(FILE *ofp, double *prob, int L);

