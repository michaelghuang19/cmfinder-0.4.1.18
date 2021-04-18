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

/* er_math.c
 * 
 * should go to a library
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#include <time.h>
#include <limits.h>
#include <float.h>

#include "er_math.h"
#include "misc.h"
#include "matrix.h"
#include "evolve.h"
#include "squid.h"
#include "mat_globals.h"
#include "mat_structs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: CheckDoubleProb()
 * Date:     ER, Tue May 25 13:18:02 CDT 1999  [St. Louis] 
 * 
 * Purpose:  Verify that \sum_{x,y} P(x,y) = 1 
 * 
 * Args:     pdouble - double-argument probability distribution to check 
 * 
 * Returns:  void.  
 */ 
void
CheckDoubleProb(double **pdouble, int dimx, int dimy)
{
  int    x;
  int    y;
  double sum = 0.0;

  for (x = 0; x < dimx; x++)
    for (y = 0; y < dimy; y++){
      if (pdouble[x][y] < -MARGIN) Die ("CheckDoubleProb(dimx=%d,dimy=%d): probabilities are getting too small here P = %f", dimx, dimy, pdouble[x][y]);
      sum += pdouble[x][y];
    }
  
  if (sum > 1.0+MARGIN || sum < 1.0-MARGIN) Die("CheckDoubleProb(dimx=%d,dimy=%d): sum_{x,y} P(x,y) is %f\n", dimx, dimy, sum);
}

/* Function: CheckSingleProb()
 * Date:     ER, Tue May 25 13:18:02 CDT 1999  [St. Louis]
 *
 * Purpose:  Verify that \sum_x P(x) = 1
 *
 * Args:     psingle - single-argument probability distribution to check
 *
 * Returns:  void. 
 */
void
CheckSingleProb(double *psingle, int size)
{
  int x;
  double sum = 0.0;

  for (x = 0; x < size; x++) {
    if (psingle[x] > 1.0+MARGIN) Die ("CheckSingleProb(L=%d): probabilities are getting too large here. P[%d] = %f", size, x, psingle[x]);
    if (psingle[x] < -MARGIN)    Die ("CheckSingleProb(L=%d): probabilities are getting too small here, P[%d] = %f", size, x, psingle[x]);
    sum += psingle[x];
  }
  
  if (sum > 1.0+MARGIN || sum < 1.0-MARGIN) Die ("CheckSingleProb(L=%d): sum_x P(x) is %f\n", size, sum);
}
void
CheckSingleLog2Prob(double *psingle, int size)
{
  int    x;
  double prob;
  double sum = 0.0;

  for (x = 0; x < size; x++) {
    if (psingle[x] > - 50.) {
      prob = EXP2(psingle[x]);
      if (prob > 1.0+MARGIN) Die ("CheckSingleLog2Prob(L=%d): probabilities are getting too large here. P[%d] = %f", size, x, prob);
      if (prob <    -MARGIN) Die ("CheckSingleLog2Prob(L=%d): probabilities are getting too small here. P[%d] = %f", size, x, prob);
      sum += prob;
    }
  }

  if (sum > 1.0+MARGIN || sum < 1.0-MARGIN) Die ("CheckSingleLog2Prob(L=%d): sum_x P(x) is %f\n", size, sum);
}

void
CheckSYMLog2Prob(double *mat, int edge_size)
{
  int    L;
  int    i;
  int    x, y;
  double prob;
  double sum = 0.0;

  L = edge_size*edge_size;

  for (i = 0; i < L; i++) {

      x = i/edge_size;
      y = i%edge_size;

      if (y <= x) {
	
	if (mat[matrix_index(y,x)] > -50.0) {
	  prob = EXP2(mat[matrix_index(y,x)]);
	  
	  if (prob > 1.0+MARGIN2) Die ("CheckSingleLog2Prob(): probabilities are getting too large here. P[%d] = %f", i, prob);
	  if (prob <    -MARGIN2) Die ("CheckSingleLog2Prob(): probabilities are getting too small here. P[%d] = %f", i, prob);
	  
	  sum += (x == y)? prob : 2.0*prob;
	}
	
      }
  }
  
  if (sum > 1.0+MARGIN2 || sum < 1.0-MARGIN2) Die ("CheckSYMLog2Prob(): sum_x P(x) is %f\n", sum);
}

/* Function: DistributionMeanVar()
 * Date:     ER, Wed Jul 16 15:19:40 CDT 2003   [St. Louis] 
 * 
 * Purpose:  Calculates \sum_{x} x P(x)  
 * 
 * Args:     p - distributio
 * 
 * Returns:  void.  
 */
void
DistributionMeanVar(double *p, int dim, double *ret_mean, double *ret_var)
{
  double mean = 0.0;
  double var  = 0.0;
  int    x;
  
  for (x = 0; x < dim; x++) {
    if (p[x] < -MARGIN) Die ("DistributionAverage(): probabilities are getting too small here P = %f", p[x]);
    mean += x*p[x];
    var  += x*x*p[x];
  }
  
  var -= mean * mean;
  
  var = sqrt(var);

  *ret_mean = mean;
  *ret_var  = var;
}
void
DistributionMean(double *p, int dim, double *ret_mean)
{
  double mean = 0.0;
  int    x;
  
  for (x = 0; x < dim; x++) {
    if (p[x] < -MARGIN) Die ("DistributionAverage(): probabilities are getting too small here P = %f", p[x]);
    mean += x*p[x];
  }
  
  *ret_mean = mean;
}

void
DistributionLogMean(double *p, int dim, double *ret_mean)
{
  double mean = 0.0;
  int    x;
  
  for (x = 0; x < dim; x++) 
    mean +=     x * EXP2(p[x]);

  *ret_mean = mean;
}

void
DistributionLogMeanVar(double *p, int dim, double *ret_mean, double *ret_var)
{
  double mean = 0.0;
  double var  = 0.0;
  int    x;
  
  for (x = 0; x < dim; x++) {
    mean +=     x * EXP2(p[x]);
    var  += x * x * EXP2(p[x]);
  }
  
  var -= mean * mean;
  
  var = sqrt(var);

  *ret_mean = mean;
  *ret_var  = var;
}

int
ILogSum(int *logp, int n)
{
  int x;
  float max = -1.e30;
  float sum;
  
  for (x = 0; x < n; x++)
    if ((float)logp[x] > max) max = (float)logp[x];

  sum = 0.0;
  for (x = 0; x < n; x++)
    if ((float)logp[x] > max - 50.)
      sum += exp((float)logp[x] - max);
  sum = log(sum) + max;

  return (int)sum;
}

int
ILog2Sum(int *log2p, int n)
{
  int x;
  float max = -1.e30;
  float sum;
  
  for (x = 0; x < n; x++)
    if ((float)log2p[x] > max) max = (float)log2p[x];

  sum = 0.0;
  for (x = 0; x < n; x++)
    if ((float)log2p[x] > max - 50.)
      sum += EXP2((float)log2p[x] - max);

  sum = LOG2(sum) + max;

  return (int)sum;
}

/* Functions: DChoose(), FChoose()
 *
 * Purpose:   Make a random choice from a normalized distribution.
 *            DChoose() is for double-precision vectors;
 *            FChoose() is for single-precision float vectors.
 *            Returns the number of the choice.
 */
int
DENChoose(double *p, int N)
{
  double *vec;                   /* vector to store exponetiated an normalized probs */
  double  sum = 0.;              /* integrated prob */
  int     i;                     /* counter over the probs */


  vec = (double *)MallocOrDie(sizeof(double) * N);

  /* exponentiate */
  for (i = 0; i < N; i++)
    {
      vec[i] = EXP2(p[i]);
      sum += vec[i];
    }

  /* normalize */
  for (i = 0; i < N; i++) 
    vec[i] /= sum; 

  return DChoose(vec, N);
}


void
DExp2(double *vec, int n)
{
  int i;
  
  for (i = 0; i < n; i++)
    if (vec[i] > -50.) vec[i] = EXP2(vec[i]);
    else               vec[i] = 0.0;
}

void
DLog2(double *vec, int n)
{
  int x;
  for (x = 0; x < n; x++)
    if (vec[x] > 0.) vec[x] = LOG2(vec[x]);
    else             vec[x] = -BIGFLOAT;
}

int
DLog2Choose(double *p, int N)
{
  double roll;                  /* random fraction */
  double sum;                   /* integrated prob */
  int    i;                     /* counter over the probs */

  roll    = sre_random();
  sum     = 0.0;
  for (i = 0; i < N; i++)
    {
      if ((sum += EXP2(p[i])) > 1.00001) Die ("DLog2Choose() sum larger than one");
      if (roll < sum) return i;
    }
  return (int) (sre_random() * N);   /* bulletproof */
}

void
DLog2Norm(double *log2p, int n)
{
  double dlogsum;
  int    x;
 
  dlogsum = DLog2Sum(log2p, n);  

  for (x = 0; x < n; x++) log2p[x] -= dlogsum;

}

void
DLog2SYMNorm(double *log2p, int edge_size)
{
  double dlogsum;
  int    L;
  int    i;
  int    x, y;
 
  L = edge_size*edge_size;

  dlogsum = DLog2SYMSum(log2p, edge_size);  

  for (i = 0; i < L; i++) {
    x = i/edge_size;
    y = i%edge_size;
    if (y <= x) 
      log2p[matrix_index(y,x)] -= dlogsum;
  }
  
}

double
DLog2Sum(double *log2p, int n)
{
  int x;
  double max = -1.e30;
  double sum;
  
  for (x = 0; x < n; x++)
    if (log2p[x] > max) max = log2p[x];

  sum = 0.0;
  for (x = 0; x < n; x++)
    if (log2p[x] > max - 50.)
      sum += EXP2(log2p[x] - max);

  sum = LOG2(sum) + max;

  return sum;
}
double
DLog2SYMSum(double *log2p, int edge_size)
{
  int L;
  int i;
  int x, y;
  double max = -1.e30;
  double sum;
  
  L = edge_size*edge_size;

  for (i = 0; i < L; i++) {
    x = i/edge_size;
    y = i%edge_size;

    if (y <= x) 
      if (log2p[matrix_index(y,x)] > max) max = log2p[matrix_index(y,x)];
  }
  
  sum = 0.0;
  for (i = 0; i < L; i++) {
    
    x = i/edge_size;
    y = i%edge_size;
    
    if (log2p[matrix_index(y,x)] > max - 50.)
      sum += EXP2(log2p[matrix_index(y,x)] - max);
  }
  
  sum = LOG2(sum) + max;
  
  return sum;
}



void
FLog2(float *vec, int n)
{
  int i;
  
  for (i = 0; i < n; i++)
    if (vec[i] > 0.0) vec[i] = LOG2(vec[i]);
    else Die("invalid data for FLog2");
}






void
FExp2(float *vec, int n)
{
  int i;
  
  for (i = 0; i < n; i++)
    vec[i] = exp(LN2*vec[i]);
}

void
FAddCOnst(float *vec, float plus, int n)
{
  int i;
  
  for (i = 0; i < n; i++)
    vec[i] = vec[i] + plus;
}

void
IAddConst(int *vec, int plus, int n)
{
  int i;
  
  for (i = 0; i < n; i++)
    vec[i] = vec[i] + plus;
}






