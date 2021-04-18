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
/* from er_math.c:
 */
extern void   CheckDoubleProb(double **pdouble, int dimx, int dimy);
extern void   CheckSingleProb(double *psingle, int size);
extern void   CheckSingleLog2Prob(double *psingle, int size);
extern void   CheckSYMLog2Prob(double *mat, int edge_size);
extern void   DistributionMeanVar(double *p, int dim, double *ret_mean, double *ret_var);
extern void   DistributionMean(double *p, int dim, double *ret_mean);
extern void   DistributionLogMean(double *p, int dim, double *ret_mean);
extern void   DistributionLogMeanVar(double *p, int dim, double *ret_mean, double *ret_var);
extern int    DENChoose(double *p, int N); 
extern void   DExp2(double *vec, int n);
extern void   DLog2(double *vec, int n);
extern int    DLog2Choose(double *p, int N);
extern void   DLog2Norm(double *log2p, int n);
extern void   DLog2SYMNorm(double *log2p, int n);
extern double DLog2Sum(double *log2p, int n);
extern double DLog2SYMSum(double *log2p, int n);
extern int    ILogSum(int *logp, int n);
extern int    ILog2Sum(int *log2p, int n);
extern void   FLog(float *vec, int n);
extern void   FLog2(float *vec, int n);
extern void   FExp(float *vec, int n);
extern void   FExp2(float *vec, int n);
extern void   FAddConst(float *vec, float plus, int n);
extern void   IAddConst(int *vec, int plus, int n);

