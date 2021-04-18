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

/* evolve.c
 *
 * ER, Fri Nov 22 10:53:20 CST 2002 [St. Louis]
 * 
 * Purpose:  Given:
 *                       - a set of joint probabilities P(a,b| t) at a given evolutionary time t,
 *
 *                       - set of conditional probabilities at time zero Q_0
 *         
 *           Calculate:
 *                       - P(a,b| r*t) using a Markov chain model for evolution
 *  
 *
 * Method:   (1) conditionals: P(a,b | t) --> Q(a | b,t) 
 *            
 *           (2) Q(a | b, r*t) = Q_0 * exp{r * log[Q_0^{-1} * Q}
 *
 *           (3) Q(a | b, r*t) --> P(a,b | r*t)
 *
 *
 *
 * ConditionalsEvolved()  creates:    Q(a | b, r*t)   from    P(a,b | t)
 *
 * Joint2Joint()          converts:   P(a | b, t)     to      P(a,b | r*t)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>

#include "er_math.h"
#include "misc.h"
#include "matrix.h"
#include "evolve.h"
#include "mat_globals.h"
#include "squid.h"
#include "mat_structs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static void    adjust_prob(double *psingle, int size);
static double  calfunc(double x, double *S, double *pml, double *pmr, int L);
static double  caldfunc(double x, double *S, double *pml, double *pmr, int L);
static void    check_three_probs (double *T, double *T_zero, double *T_infty, int dim);
static void    check_reversibility(double *QL, double *QR, double *ml, double *mr, int L);
static void    check_Q_0_reversibility (double *QL, double *QR, double *ml, double *mr, int L, int hasindel);

double
Cal_lambda(FILE *ofp, double *S, double *pml, double *pmr, int L, int verbose)
{
  double lambda;
  double func;
  double dfunc;
  double x = 0;
  double xnew = 1.0;
  
  while (fabs(xnew-x) > 1.0-accuracy1) {
    
    x = xnew;

    func  = calfunc (x, S, pml, pmr, L);
    dfunc = caldfunc(x, S, pml, pmr, L);
    
    if (dfunc != 0.0) xnew = x - func/dfunc;
    
  }

  lambda = xnew;
  
  return lambda;
}

double
calfunc(double x, double *S, double *pml, double *pmr, int L)
{
  double y = 0.0;
  int    i,j;

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) 
      y += exp(x*S[i*L+j]) * pml[i] * pmr[j];

  y -= 1.0;

  return y;
}
double
caldfunc(double x, double *S, double *pml, double *pmr, int L)
{
  double y = 0.0;
  int    i,j;

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) 
      y += S[i*L+j] * exp(x*S[i*L+j]) * pml[i] * pmr[j];

  return y;
}

/* Function: CalculateConditionalsAndMarginals()
 *
 *           ER, Wed Sep 17 10:36:04 CDT 2003 [St. Louis]
 *
 * Purpose:  Given       a LxL matrix P(i,j)  of joint probabilities, 
 *           calculate the LxL matrix Ql(i,j) of conditional probabilities.
 *           calculate the L   vecotr ml(i)   of marginal probabilities.
 *
 *           ml[i]     = \sum_k  P(i,k|t_o) 
 *
 *           Ql[i*L+j] = P(j|i, t_o) = P(i,j|t_o) / ml[i]
 *
 *           where P(i,j|t_o) are calculated according to 
 *           using a given pammodel t_o.
 *
 *           Notice that if P(i,j) are not real joint probs (ie non symmetric)
 *           this function calculates the LEFT conditionals and LEFT marginals.
 *
 * Args:     P  - LxL joint prob matrix (prealloc)
 *
 * Returns:  Ql(LxL), conditional probabilities.
 *           ml(L),   marginal probabilities.
 *           Q and ml are alocated here, freed by caller.
 */
void
CalculateConditionalsAndMarginals(FILE *ofp, double *P, double **ret_Ql, double **ret_Qr, 
				  double **ret_ml, double **ret_mr, int L, int hasindel, int verbose)
{
  double *ml;
  double *mr;
  double *Ql;
  double *Qr;
  double  sum;
  int     dim;
  int     i, j;
 
  if (hasindel) dim = L-1;
  else          dim = L;

  if (verbose) {
    fprintf(stdout, "Joint Probs\n");
    PrintProbs(stdout, P, L);
  }
 /* paranoia */
  CheckSingleProb(P, L*L);

  /* allocate memory
   */
  ml = (double *) MallocOrDie (sizeof(double) * L);
  mr = (double *) MallocOrDie (sizeof(double) * L);
  Ql = (double *) MallocOrDie (sizeof(double) * L * L);
  Qr = (double *) MallocOrDie (sizeof(double) * L * L);

  /* calculate LEFT/RIGHT marginals
   */
  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) {
      Ql[i*L+j] = 0.0;
      Qr[i*L+j] = 0.0;
    }

  for (i = 0; i < L; i++) {
    ml[i] = 0.0;
    mr[i] = 0.0;
  }

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) {
      ml[i] += P[i*L+j];
      mr[i] += P[j*L+i];
    }
  
  CheckSingleProb(ml, L);
  CheckSingleProb(mr, L);
  
  /* calculate LEFT conditionals
   */
  for (i = 0; i < L; i++) {
    sum = 0.;
    for (j = 0; j < L; j++) {
      if (ml[i] > 0.0) Ql[i*L+j] = P[i*L+j] / ml[i];
      sum += Ql[i*L+j];
    }

   /* normalize */
    for (j = 0; j < L; j++) 
      if (sum > 0.0) Ql[i*L+j] /= sum;
  }
   
  /* calculate RIGHT conditionals
   */
  for (i = 0; i < L; i++) {
    sum = 0.;
    for (j = 0; j < L; j++) {
      if (mr[i] > 0.0) Qr[i*L+j] = P[j*L+i] / mr[i];
      sum += Qr[i*L+j];
    }
    
    /* normalize */
    for (j = 0; j < L; j++) 
      if (sum > 0.0) Qr[i*L+j] /= sum;
  }

  /* consistency check, rows should add up to one
   */
  for (i = 0; i < dim; i++)
    CheckSingleProb(Ql+i*L, L);
  for (i = 0; i < dim; i++)
    CheckSingleProb(Qr+i*L, L);
  
  if (verbose) {
    printf("RIGHT  marginalized probabilities\n");
    PrintVectorProbs(stdout, mr, L);
    printf("LEFT marginalized probabilities\n");
    PrintVectorProbs(stdout, ml, L);
    fprintf(stdout, "RIGHT Conditionals\n");
    PrintProbs(stdout, Qr, L);
    fprintf(stdout, "LEFT Conditionals\n");
    PrintProbs(stdout, Ql, L);
  }
  
  *ret_ml = ml;
  *ret_mr = mr;
  *ret_Ql = Ql;
  *ret_Qr = Qr;
}

void
ComputeLRConditionalsAndMarginals(FILE *ofp, double *P, double *QL, double *QR, double *pml, double *pmr, int L, int hasindel, int verbose)
{
  double *Ql;
  double *Qr;
  double *ml;
  double *mr;
 
  CalculateConditionalsAndMarginals(ofp, P, &Ql, &Qr, &ml, &mr, L, hasindel, verbose);

   CopyMatrix (QL, Ql, L, L);
   CopyMatrix (QR, Qr, L, L);
   CopyVector (pml, ml, L);
   CopyVector (pmr, mr, L);

   free (Ql);
   free (Qr);
   free (ml);
   free (mr);
}

void
ComputeConditionalsAndMarginals(FILE *ofp, double *P, double *Q, double *pm, int L, int hasindel, int verbose)
{
  double *Ql;
  double *Qr;
  double *ml;
  double *mr;
 
  CalculateConditionalsAndMarginals(ofp, P, &Ql, &Qr, &ml, &mr, L, hasindel, verbose);

   CopyMatrix (Q, Ql, L, L);
   CopyVector (pm, ml, L);

   free (Ql);
   free (Qr);
   free (ml);
   free (mr);
}

int
Check_Accuracy(double *vec, int L)
{
  int    i;
  int    issig;
  double m;

  issig = FALSE;
  
  for (i = 0; i < L; i ++) 
    {
      m = (vec[i] > 0)? vec[i]:-vec[i];
      if (m > DBL_EPSILON) {
	issig = TRUE;
	break;
      }
    }
  
  return issig;
}

/* Function: ChangeFrequencies()
 * Date:     ER, Sat Sep 27 10:05:51 CDT 2003 [St. Louis]
 *
 * Purpose:  Given a set of marginals pml, pmr change them
 *           to targetfreq
 *
 *           If this is the case of having indels as extra characters, we do that
 *           with only the subspace (L-1)x(L-1).
 *
 *  
 * Args:              
 *
 * Returns:  (void)
 *           
 */
void
ChangeFrequencies(FILE *ofp, double *Q_0, double *R, double *pm, int L, double *targetfreq, int hasindel, int verbose)
{
  double *psat;
  double *target;
  int     dim;
  int     i;
  
  if (hasindel) dim = L-1;
  else          dim = L;
  
  if (verbose) {
    fprintf(ofp, "old marginals\n");
    PrintVectorProbs(ofp, pm, L);
  }

  target = (double *) MallocOrDie (sizeof(double) * dim);

  psat = SaturationProbs(R, Q_0, L, FALSE, verbose);
  
 if (verbose) {
    fprintf(ofp, "saturation probs\n");
    PrintVectorProbs(ofp, psat, L);
  }

 if (hasindel) {
   if (fabs(1.0-psat[dim]) > 1.0-accuracy1) {
     for (i = 0; i < dim; i++)  
       target[i] = psat[i] / (1.0 - psat[dim]);
   }
   else { /* singular case  we cannot infer the background freqs from the saturation probs */
     for (i = 0; i < dim; i++)  
       target[i] = targetfreq[i];
   }
 }
 else {
   for (i = 0; i < dim; i++)  
     target[i] = psat[i];
 }
 if (verbose) {
   fprintf(ofp, "ChangeFrequencies(): New Target marginals\n");
   PrintVectorProbs(ofp, target, dim);
 }
 CheckSingleProb(target, dim);
 
 for (i = 0; i < dim; i++) 
   if (hasindel)  pm[i] = target[i] * (1.0 - pm[dim]);
   else           pm[i] = target[i];
 
 if (verbose) {
   fprintf(ofp, "ChangeFrequencies(): New marginals\n");
   PrintVectorProbs(ofp, pm, L);
 }
 CheckSingleProb(pm, L);
 
 if (verbose) {
   fprintf(ofp, "ChangeFrequencies(): New marginals\n");
   PrintVectorProbs(ofp, pm, L);
   fprintf(ofp, "ChangeFrequencies(): New Target marginals\n");
   PrintVectorProbs(ofp, target, dim);
 }
 
 free(target);
 free (psat);
}

/* Function: ChangeConditionalMatrix()
 *
 * Date:     ER, Fri Aug 27 13:03:54 CDT 2004[ St. Louis at work, Coro with Maribel, her last day! ]
 *
 * Purpose:  Given a matrix of conditional probabilities at a given time instance,  modify the
 *           matrix so it stays reversible respect to the new frequecies.
 *
 *           If this is the case of having indels as extra characters, we do that
 *           with only the subspace (L-1)x(L-1).
 *  
 * Method:   
 *    
 *                    Q_ij = S_ij (p_i*p_j)^{alpja_{ij}} * p_j  for i \neq j
 *
 *          to  
 *                    ~Q_ij = S_ij (~p_i*~p_j/p_i*p_j)^{alpja_{ij}} * ~p_j/p_j  for i \neq j
 *
 *                    set diagonal by normalization.
 *
 * Args:              
 *
 * Returns:  (void)
 *           
 */
void
ChangeConditionalMatrix(FILE *ofp, double *Q_star, double *pml_star, double *pmr_star, int L, double *targetfreq, 
			int hasindel, int isleft, int verbose)
{
  double           *pstatl;
  double           *pstatr;
  double           *alpha;
  double            comesfrom;
  double            goesto;
  int               dim;
  int               i,j;

  if (hasindel) dim = L-1;
  else          dim = L;
  
 if (verbose) {
    fprintf(ofp, "Old Q Matrix\n");
    PrintProbs(ofp, Q_star, L);
 }

 /* paranoia */
  CheckSingleProb(pml_star, L);
  CheckSingleProb(pmr_star, L);
  CheckSingleProb(targetfreq, dim);
  
  pstatl = (double *) MallocOrDie (sizeof(double) * dim);
  pstatr = (double *) MallocOrDie (sizeof(double) * dim);

  for (i = 0; i < dim; i++) {
    if (hasindel) pstatl[i] = pml_star[i]/(1.0-pml_star[dim]);
    else          pstatl[i] = pml_star[i];

    if (hasindel) pstatr[i] = pmr_star[i]/(1.0-pmr_star[dim]);
    else          pstatr[i] = pmr_star[i];
  }
  
  alpha = (double *) MallocOrDie (sizeof(double) * dim * dim);

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) 
      alpha[i*dim+j] = 0.0;

  for (i = 0; i < dim; i++)
    for (j = i+1; j < dim; j++) {
      alpha[i*dim+j] = -0.5;      
      alpha[j*dim+i] = alpha[i*dim+j]; 
    }

  if (verbose) {
    fprintf(ofp, "alpha Matrix\n");
    PrintProbs(ofp, alpha, dim);
    fprintf(ofp, "stat L marginals\n");
    PrintVectorProbs(ofp, pstatl, L);
    fprintf(ofp, "stat R marginals\n");
    PrintVectorProbs(ofp, pstatr, L);
    fprintf(ofp, "target marginals\n");
    PrintVectorProbs(ofp, targetfreq, L);
  }
  
  /* Non diagonal "reduced" components of (Q_0*R) */
  for (i = 0; i < dim; i++) 
    for (j = 0; j < dim; j++) {
      
      if (isleft) {
	comesfrom = LOG2(targetfreq[i]) - LOG2(pstatl[i]);
	goesto    = LOG2(targetfreq[j]) - LOG2(pstatr[j]);
      }
      else {
	comesfrom = LOG2(targetfreq[i]) - LOG2(pstatr[i]);
	goesto    = LOG2(targetfreq[j]) - LOG2(pstatl[j]);
      }
      
     if (i != j ) Q_star[i*L+j] = EXP2( LOG2(Q_star[i*L+j]) + alpha[i*dim+j]*comesfrom + (1.0+alpha[i*dim+j])*goesto );
    }

   if (hasindel) {
    for (i = 0; i < dim; i++) 
     Q_star[dim*L+i] = targetfreq[i] * (1.0 - Q_star[dim*L+dim]);
 
  }

 /* Diagonal components of Q set by normalization*/
  for (i = 0; i < L; i++) { 
    Q_star[i*L+i] = 1.0;
    
    for (j = 0; j < L; j++) 
      if (j != i) Q_star[i*L+i] -= Q_star[i*L+j];
  }

  if (verbose) {
    fprintf(ofp, "New Q Matrix\n");
    PrintMatrix(ofp, Q_star, L);
  }

  /* consistency check, rows should add up to one
   */
  for (i = 0; i < L; i++)
    CheckSingleProb(Q_star+i*L, L);

  free(pstatl);
  free(pstatr);
  free(alpha);
}

/* Function: ChangeQ_0Matrix()
 * Date:     ER, Sat Sep 27 10:34:17 CDT 2003 [St. Louis]
 *           fixed: Thu Aug 26 20:12:31 CDT 2004 [st. Louis, at home, Coro cocinando con Sean ]
 *
 * Purpose:  Given a matrix of conditionals at time zero Q_0,  modify the
 *           rmatrix so it stays reversible respect to the new frequecies.
 *
 *           If this is the case of having indels as extra characters, we do that
 *           with only the subspace (L-1)x(L-1).
 *  
 * Method:   
 *    
 *                    (Q_0)_ij = S_ij (p_i*p_j)^{alpja_{ij}} * p_j  
 *
 *          to  
 *                    (Q_0) _ij = S_ij (p'_i*p'_j/p_i*p_j)^{alpja_{ij}} * p'_j/p_j  
 *
 *
 * Args:              
 *
 * Returns:  (void)
 *           
 */
void
ChangeQ_0Matrix(FILE *ofp, double *Q_0, double *ml, double *mr, int L, double *targetfreq, 
		int hasindel, int isleft, int verbose)
{
  double           *pstatl;
  double           *pstatr;
  double           *alpha;
  double            comesfrom;
  double            goesto;
  int               dim;
  int               i,j;

  if (hasindel) dim = L-1;
  else          dim = L;
  
 if (verbose) {
    fprintf(ofp, "Old Q_0 Matrix\n");
    PrintProbs(ofp, Q_0, L);
 }

 /* paranoia */
  CheckSingleProb(ml, L);
  CheckSingleProb(mr, L);
  CheckSingleProb(targetfreq, dim);
  
  pstatl = (double *) MallocOrDie (sizeof(double) * dim);
  pstatr = (double *) MallocOrDie (sizeof(double) * dim);

  for (i = 0; i < dim; i++) {
    if (hasindel) pstatl[i] = ml[i]/(1.0-ml[L-1]);
    else          pstatl[i] = ml[i];

    if (hasindel) pstatr[i] = mr[i]/(1.0-mr[L-1]);
    else          pstatr[i] = mr[i];
  }
  
  alpha = (double *) MallocOrDie (sizeof(double) * dim * dim);

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) 
      alpha[i*dim+j] = 0.0;

  for (i = 0; i < dim; i++)
    for (j = i+1; j < dim; j++) {
      alpha[i*dim+j] = -0.5;      
      alpha[j*dim+i] = alpha[i*dim+j]; 
    }

  if (verbose) {
    fprintf(ofp, "alpha Matrix\n");
    PrintProbs(ofp, alpha, dim);
    fprintf(ofp, "stat L marginals\n");
    PrintVectorProbs(ofp, pstatl, L);
    fprintf(ofp, "stat R marginals\n");
    PrintVectorProbs(ofp, pstatr, L);
    fprintf(ofp, "target marginals\n");
    PrintVectorProbs(ofp, targetfreq, L);
  }
  
 /* Non diagonal "reduced" components of (Q_0*R) */
  for (i = 0; i < dim; i++) 
    for (j = 0; j < dim; j++) {

      if (isleft) {
	comesfrom = LOG2(targetfreq[i]) - LOG2(pstatl[i]);
	goesto    = LOG2(targetfreq[j]) - LOG2(pstatr[j]);
      }
      else {
	comesfrom = LOG2(targetfreq[i]) - LOG2(pstatr[i]);
	goesto    = LOG2(targetfreq[j]) - LOG2(pstatl[j]);
      }
      
     if (i != j ) Q_0[i*L+j] = EXP2( LOG2(Q_0[i*L+j]) + alpha[i*dim+j]*comesfrom + (1.0+alpha[i*dim+j])*goesto );
    }

  if (hasindel) {
    for (i = 0; i < dim; i++) {
      Q_0[i*L+dim] = 0.0;
      Q_0[dim*L+i] = targetfreq[i] * (1.0 - Q_0[dim*L+dim]);
    }
  }

  /* Diagonal components of (Q_0*R) set by normalization*/
  for (i = 0; i < L; i++) { 
    Q_0[i*L+i] = 1.0;
    
    for (j = 0; j < L; j++) 
      if (j != i) Q_0[i*L+i] -= Q_0[i*L+j];
  }

  if (verbose) {
    fprintf(ofp, "New Q_0 Matrix\n");
    PrintMatrix(ofp, Q_0, L);
  }

  /* consistency check, rows should add up to one
   */
  for (i = 0; i < L; i++)
    CheckSingleProb(Q_0+i*L, L);
 
  free(pstatl);
  free(pstatr);
  free(alpha);
}

/* Function: ChangeQ_0MatrixNaive()
 * Date:    on Sep 29 18:25:19 CDT 2003  [St. Louis]
 *
 * Purpose:  Given a matrix of conditionals at time zero Q_0,  modify the
 *           rmatrix so it stays reversible respect to the new frequecies.
 *
 *           If this is the case of having indels as extra characters, we do that
 *           with only the subspace (L-1)x(L-1).
 *  
 * Method:   
 *    
 *
 *
 * Args:              
 *
 * Returns:  (void)
 *           
 */
void
ChangeQ_0MatrixNaive(FILE *ofp, double *QL_0, double *QR_0, int L, double *targetfreq, int hasindel, int isleft, int verbose)
{
  double           *P;
  double           *Id;
  double           *Qnewl;
  double           *Qnewr;
  double           *pml;
  double           *pmr;
  int               flag;
  int               dim;
  int               i,j;

  Id = (double *)Cal_Id(L);

  if (hasindel) dim = L-1;
  else          dim = L;
  
 if (verbose) {
    fprintf(ofp, "Old QL_0 Matrix\n");
    PrintProbs(ofp, QL_0, L);
    fprintf(ofp, "Old QR_0 Matrix\n");
    PrintProbs(ofp, QR_0, L);
 }

 flag = TRUE;
 for (i = 0; i < dim; i++)
   for (j = 0; j < dim; j++) 
     if (fabs(QL_0[i*L+j] - Id[i*L+j])  > MARGIN) flag = FALSE; 
     
 if (flag) return;

 /* paranoia */
  CheckSingleProb(targetfreq, dim);
  
  P = (double *) MallocOrDie (sizeof(double) * dim * dim);

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) 
      P[i*dim+j] = targetfreq[i] * 0.5 * (QL_0[i*L+j] + QR_0[j*L+i]);

 if (verbose) {
    fprintf(ofp, "new P Matrix\n");
    PrintProbs(ofp, P, dim);
 }
  
  CalculateConditionalsAndMarginals(ofp, P, &Qnewl, &Qnewr, &pml, &pmr, dim, FALSE, verbose);
  fprintf(ofp, "Naive newQl(i|j,t) probabilities\n");
  PrintProbs(ofp, Qnewl, dim);
  fprintf(ofp, "Naive newQr(i|j,t) probabilities\n");
  PrintProbs(ofp, Qnewr, dim);
  fprintf(ofp, "Naive marginals probabilities\n");
  PrintVectorProbs(ofp, pml, dim);
  PrintVectorProbs(ofp, pmr, dim);
  
  free(Qnewl);
  free(Qnewr);
  free(pml);
  free(pmr);
  
  free(P);
  free(Id);
}

void
ChangeQ_0MatrixIterate(FILE *ofp, double *Q_0, double *Q_0_inv, double *pml, double *pmr, int L, double *targetfreq, 
		       int hasindel, int isleft, int verbose)
{
  double           *P;
  double           *S;
  double           *Ql;
  double           *Qr;
  double           *Id;
  double           *pstatl;
  double           *pstatr;
  double            lambda;
  int               flag;
  int               dim;
  int               iterations = 0;
  int               i,j;

  if (hasindel) dim = L-1;
  else          dim = L;
  
 Id = Cal_Id(L);

 /* If the process is "quasi-reversible,
  * the only implication of reversability in the off diagonals
  * of Q_0 is that Q_0(i,5)= 0 which is always the case
  */
  flag = TRUE;
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) 
      if (fabs(Q_0[i*L+j] - Id[i*L+j])  > MARGIN) flag = FALSE; 
  if (flag) { free(Id); return; }
  
  /* paranoia */
  CheckSingleProb(pml, L);
  CheckSingleProb(pmr, L);
  CheckSingleProb(targetfreq, dim);
  
  P      = (double *) MallocOrDie (sizeof(double) * dim * dim); 
  Ql     = (double *) MallocOrDie (sizeof(double) * dim * dim); 
  Qr     = (double *) MallocOrDie (sizeof(double) * dim * dim); 
  S      = (double *) MallocOrDie (sizeof(double) * dim * dim); 
  pstatl = (double *) MallocOrDie (sizeof(double) * dim);
  pstatr = (double *) MallocOrDie (sizeof(double) * dim);


  for (i = 0; i < dim; i++) {
    if (hasindel) pstatl[i] = pml[i]/(1.0-pml[L-1]);
    else          pstatl[i] = pml[i];

    if (hasindel) pstatr[i] = pmr[i]/(1.0-pmr[L-1]);
    else          pstatr[i] = pmr[i];
  }

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) 
      if (isleft) { Ql[i*dim+j] = Q_0[i*L+j]; Qr[i*dim+j] = 0.0; }
      else        { Qr[i*dim+j] = Q_0[i*L+j]; Ql[i*dim+j] = 0.0; }

  if (isleft)
    Joint_From_Condi(ofp, P, Ql, pstatl, dim, verbose);
  else
    Joint_From_Condi(ofp, P, Qr, pstatr, dim, verbose);

  if (verbose) {
    fprintf(ofp, "Old Ql_0(i|j,t) probabilities\n");
    PrintProbs(ofp, Ql, dim);
    fprintf(ofp, "Old Qr_0(i|j,t) probabilities\n");
    PrintProbs(ofp, Qr, dim);
    fprintf(ofp, "Old marginals probabilities\n");
    PrintVectorProbs(ofp, pml, L);
    PrintVectorProbs(ofp, pmr, L);
    fprintf(ofp, "Old P_0 Matrix\n");
    PrintProbs(ofp, P, dim);
  }

  while (CompareFreqs(pstatl, pstatr, targetfreq, dim)) {
    iterations ++;

    for (i = 0; i < dim; i++) 
      for (j = 0; j < dim; j++) 
	S[i*dim+j] = log(P[i*dim+j]) - log(pstatl[i]) - log(pstatr[j]);
    
    lambda = Cal_lambda(ofp, S, targetfreq, targetfreq, dim, verbose);
    
    for (i = 0; i < dim; i++) 
      for (j = 0; j < dim; j++) 
	P[i*dim+j] = exp(lambda*S[i*dim+j])*targetfreq[i]*targetfreq[j];
    
    /* check they are probabilities */  
    CheckSingleProb(P, dim*dim);
    
    ComputeLRConditionalsAndMarginals(ofp, P, Ql, Qr, pstatl, pstatr, dim, FALSE, verbose);
  }
  
  if (verbose) {
    fprintf(ofp, "new P_0 (%d it) Matrix\n", iterations);
    PrintProbs(ofp, P, dim);
    fprintf(ofp, "new Ql_0(i|j,t) probabilities\n");
    PrintProbs(ofp, Ql, dim);
    fprintf(ofp, "new Qr_0(i|j,t) probabilities\n");
    PrintProbs(ofp, Qr, dim);
    fprintf(ofp, "new marginals probabilities\n");
    PrintVectorProbs(ofp, pstatl, dim);
    PrintVectorProbs(ofp, pstatr, dim);  
  }
  
    for (i = 0; i < dim; i++) 
      for (j = 0; j < dim; j++) 
	if (isleft) Q_0[i*L+j] = Ql[i*dim+j];
	else        Q_0[i*L+j] = Qr[i*dim+j];

  /* consistency check, rows should add up to one
   */
  for (i = 0; i < L; i++)
    CheckSingleProb(Q_0+i*L, L);
  
  Comp_M_Inv(ofp, Q_0, Q_0_inv, L, verbose);

  free(Ql);
  free(Qr);
  free(pstatl);
  free(pstatr);
  
  free(Id);
  free(S);
  free(P);
}

/* Function: ChangePairProbs()
 * Date:     ER, Wed Mar  5 10:04:54 CST 2003 [St. Louis]
 *
 * Purpose: To change the overall base composition, [we do this before dealing with gaps]
 *
 * calculate:         P(i|j) = P(i,j) / P(j)
 *
 * then use:         ^P(j) = f(j)
 *
 * to recalculate:   ^P(i,j) = P(i|j) * ^P(j)
 *
 *
 *
 */
void
ChangePairProbs (int L, double *pairprob, double *psingle_target, int hasindel, int verbose)
{
  double *Ql, *Qr;              /* [5x5] base-pair conditional probabilities +*/
  double *pml, *pmr;            /* [5] marginalization of the pair5probs     +*/
  int     i,j;                  /* Symbols for emission prob's               +*/
  int     LSQ = L*L;

  if (verbose) {
    fprintf(stdout, "Target Marginal probabilities \n");
    PrintVectorProbs(stdout, psingle_target, L);
  }

  /* paranoia */
  CheckSingleProb(psingle_target, L);

  CalculateConditionalsAndMarginals(stdout, pairprob, &Ql, &Qr, &pml, &pmr, L, hasindel, verbose);

  if (verbose) {
    fprintf(stdout, "old PAIR Matrix\n");
    PrintProbs(stdout, pairprob, L);
    fprintf(stdout, "old Ql(i|j,t) probabilities\n");
    PrintProbs(stdout, Ql, L);
    fprintf(stdout, "old Qr(i|j,t) probabilities\n");
    PrintProbs(stdout, Qr, L);
    fprintf(stdout, "old marginals probabilities\n");
    PrintVectorProbs(stdout, pml, L);
    PrintVectorProbs(stdout, pmr, L);  
  }

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++)
      pairprob[i*L+j] = Ql[i*L+j] * psingle_target[i];

  DNorm(pairprob, LSQ);

  if (verbose) {
    ComputeLRConditionalsAndMarginals(stdout, pairprob, Ql, Qr, pml, pmr, L, hasindel, verbose);
    
    fprintf(stdout, "new PAIR Matrix\n");
    PrintProbs(stdout, pairprob, L);
    fprintf(stdout, "new Ql(i|j,t) probabilities\n");
    PrintProbs(stdout, Ql, L);
    fprintf(stdout, "new Qr(i|j,t) probabilities\n");
    PrintProbs(stdout, Qr, L);
    fprintf(stdout, "new marginals probabilities\n");
    PrintVectorProbs(stdout, pml, L);
    PrintVectorProbs(stdout, pmr, L);  
    fprintf(stdout, "Target Marginal probabilities \n");
    PrintVectorProbs(stdout, psingle_target, L);
  }

  free(pml);
  free(pmr);
  free(Ql);
  free(Qr);
}

/* Function: ChangeStationaryProbs()
 *
 * Date:     ER,  Wed May  5 09:09:19 CDT 2004 [St. Louis, with Coro at work]
 *
 * Purpose: To change the overall base composition using my method.
 *
 *          Notice that this method is limited to targetfreq[i] >= 1/2P(i)
 *
 * calculate:         P(i|j) = P(i,j) / P(j)
 *
 * then use:         ^P(j) = 2*targetfreq(j) - P(j)
 *
 * to recalculate:   ^P(i,j) = 0.5 * [ P(i|j) * ^P(j) + P(j|i) * ^P(i) ]
 *
 *   then targetfreq are guaranteed to be the new marginal (stationary) probs. We check for that
 *
 */
void
ChangeStationaryProbs(int L, double *jointprob, double *targetfreq, int hasindel, int verbose)
{
  double *pcond;                 /* [4x4] base-pair conditional probabilities                                +*/
  double *pmar;                  /* [4] marginalization of the joint5probs                                   +*/
  double *pmix;                  /* [4]                                                                      +*/
  int     i,j;                   /* Symbols for emission prob's                                              +*/
  int     LSQ = L*L;
  double  sum1, sum2;

  /* paranoia */
  CheckSingleProb(targetfreq, L);
  CheckSingleProb(jointprob, LSQ);
  
  if (verbose) {
    fprintf(stdout, "Target Marginal probabilities \n");
    PrintVectorProbs(stdout, targetfreq, L);
    fprintf(stdout, "PJoint(i,j) Joint probabilities \n");
    PrintProbs(stdout, jointprob, L);
  }

  /* First, check that the given probabilies are joint probs ie symmetric otherwise make them symmetric and they add up to one
   */
  for (i = 0; i < L; i++)
    for (j = i+1; j < L; j++)
      if (fabs(jointprob[i*L+j] - jointprob[j*L+i]) > 0.0+MARGIN) 
	Die ("ChangeStationaryProbs(): given probabilities have to be symmetric\n");
  
  /* Calculate the marginals
   */
  pmar = (double *) MallocOrDie (sizeof(double) * L);
  MarginalizeJointProbs(jointprob, pmar, L, 2);

  /* Calculate conditionals
   */
  pcond = (double *) MallocOrDie (sizeof(double) * LSQ);
  for (i = 0; i < LSQ; i++)
    pcond[i] = 0.0;

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++)
      if (pmar[i] > 0.0) pcond[i*L+j] = jointprob[i*L+j] / pmar[i];
	      else               pcond[i*L+j] = 0.0;
  for (i = 0; i < L; i++)
    DNorm(pcond+i*L, L);

  if (verbose) {
    fprintf(stdout, "PJoint(i,j) After-symmetrization Joint probabilities \n");
    PrintProbs(stdout, jointprob, L);
    fprintf(stdout, "Conditionals probabilities \n");
    PrintProbs(stdout, pcond, L);
    fprintf(stdout, "Marginal probabilities \n");
    PrintVectorProbs(stdout, pmar, L);
  }
  
  /* Now the recalculation of the joint probabilities
   *
   *
   */
  pmix = (double *) MallocOrDie (sizeof(double) * L);
  for (i = 0; i < L; i++)
    if (2.0*targetfreq[i] - pmar[i] >= 0.0)  pmix[i] = 2.0*targetfreq[i] - pmar[i];
    else                                    
      Die ("ChangeStationaryProbs(): method does not apply for this case. i=%d p_old=%f p_new=%f\n", i, pmar[i], targetfreq[i]); 
  CheckSingleProb(pmix, L);

  if (verbose) {
    fprintf(stdout, "mix probabilities \n");
    PrintVectorProbs(stdout, pmix, L);
  }

  for (i = 0; i < L; i++) {
    sum1 = 0.0;
    sum2 = 0.0;

    for (j = 0; j < L; j++) {
      jointprob[i*L+j] = 0.5 * (pcond[i*L+j]*pmix[i] + pcond[j*L+i]*pmix[j]);
      sum1 += pcond[i*L+j]*pmix[i];
      sum2 += pcond[j*L+i]*pmix[j];
    }
    printf("i %d sum %f %f \n", i, sum1, sum2);
  }
  DNorm(jointprob, LSQ);

  /* check that the new marginals are the targetfreqs
   */
  MarginalizeJointProbs(jointprob, pmar, L, 2);
  if (verbose) {
    fprintf(stdout, "Pjoint(i,j) After--Joint probabilities \n");
    PrintProbs(stdout, jointprob, L);

    fprintf(stdout, "Pjoint Marginal probabilities \n");
    PrintVectorProbs(stdout, pmar, L);
  }
  for (i = 0; i < L; i++)
    if (targetfreq[i] - pmar[i] > MARGIN || targetfreq[i] - pmar[i] < -MARGIN) 
      Die ("ChangeStationaryProbs(): bad new marginals. pmar = %f targetfreq = %f\n", pmar[i], targetfreq[i]);

  free(pmar);
  free(pmix);
  free(pcond);

}

/* Function: ChangeStationaryProbsGW()
 *
 * Date:     ER, Mon May  3 11:50:16 CDT 2004 [St. Louis, with Coro at home]
 *
 * Purpose: To change the overall base composition,using the method of Goldman and Whelan
 *
 * calculate:         P(i|j) = P(i,j) / P(j)
 *
 * then use:         ^P(j) = 2*targetfreq(j) - P(j)
 *
 * to recalculate:   ^P(i,j) = 0.5 * [ P(i|j) * ^P(j) + P(j|i) * ^P(i) ]
 *
 *   then targetfreq are guaranteed to be the new marginal (stationary) probs. We check for that
 *
 */
void
ChangeStationaryProbsGW(int L, double *jointprob, double *targetfreq, int hasindel, int verbose)
{
  double *pcond;                 /* [4x4] base-pair conditional probabilities                                +*/
  double *Q;                     /* [4x4] modified conditional probabilities                                 +*/
  double *R;                     /* [4x4] Rate matrix                                                        +*/
  double *pmar;                  /* [4] marginalization of the joint5probs                                   +*/
  double *Id;
  int     dim;
  int     i,j;                   /* Symbols for emission prob's                                              +*/
  int     LSQ = L*L;

  if (hasindel) dim = L-1;
  else          dim = L;

  /* paranoia */
  CheckSingleProb(targetfreq, L);
  CheckSingleProb(jointprob, LSQ);
  
  if (verbose) {
    fprintf(stdout, "Target Marginal probabilities \n");
    PrintVectorProbs(stdout, targetfreq, L);
    fprintf(stdout, "PJoint(i,j) Joint probabilities \n");
    PrintProbs(stdout, jointprob, L);
  }

  /* First, check that the given probabilies are joint probs ie symmetric otherwise make them symmetric and they add up to one
   */
  for (i = 0; i < L; i++)
    for (j = i+1; j < L; j++)
      if (fabs(jointprob[i*L+j] - jointprob[j*L+i]) > 0.0+MARGIN) 
	Die ("ChangeStationaryProbsGW(): given probabilities have to be symmetric\n");
  
  /* Calculate the marginals
   */
  pmar = (double *) MallocOrDie (sizeof(double) * L);
  MarginalizeJointProbs(jointprob, pmar, L, 2);

  /* Calculate conditionals
   */
  pcond = (double *) MallocOrDie (sizeof(double) * LSQ);
  for (i = 0; i < LSQ; i++)
    pcond[i] = 0.0;

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++)
      if (pmar[i] > 0.0) pcond[i*L+j] = jointprob[i*L+j] / pmar[i];

  for (i = 0; i < dim; i++)
    DNorm(pcond+i*L, L);

  if (verbose) {
    fprintf(stdout, "Conditionals probabilities \n");
    PrintProbs(stdout, pcond, L);
    fprintf(stdout, "Marginal probabilities \n");
    PrintVectorProbs(stdout, pmar, L);
  }

  /* Calculate Rate Matrix
   */
  Id = Cal_Id(L);
  R = (double *) MallocOrDie (sizeof(double) * LSQ);
  for (i = 0; i < LSQ; i++)
    R[i] = 0.0;
  RateFromConditionals(stdout, R, pcond, Id, L, 1.0, hasindel, FALSE, verbose);

  if (verbose) {
    fprintf(stdout, "Rate Matrix \n");
    PrintProbs(stdout, R, L);
  }

  /* Modify  Rate Matrix
   */
  ChangeRateMatrix(stdout, R, Id, R, pmar, pmar, L, targetfreq, hasindel, FALSE, verbose);
  if (verbose) {
    fprintf(stdout, "New Rate Matrix \n");
    PrintProbs(stdout, R, L);
  }
    
  /* Now the recalculation of the joint probabilities
   *
   *
   */
  Q = (double *) MallocOrDie (sizeof(double) * LSQ);
  for (i = 0; i < LSQ; i++)
    Q[i] = 0.0;

  Condi_From_Rate(&Q, R, Id, 1.0, L, FALSE, verbose);
  Joint_From_Condi(stdout, jointprob, Q, targetfreq, L, verbose);
  DNorm(jointprob, LSQ);

  /* check that the new marginals are the targetfreqs
   */
  MarginalizeJointProbs(jointprob, pmar, L, 2);
  if (verbose) {
    fprintf(stdout, "Pjoint(i,j) After--Joint probabilities \n");
    PrintProbs(stdout, jointprob, L);

    fprintf(stdout, "Pjoint Marginal probabilities \n");
    PrintVectorProbs(stdout, pmar, L);
  }
  for (i = 0; i < L; i++)
    if (fabs(targetfreq[i] - pmar[i]) > MARGIN) 
      Die ("ChangeStationaryProbsGW(): bad new marginals. pmar = %f targetfreq = %f\n", pmar[i], targetfreq[i]);

  free(pmar);
  free(pcond);
  free(R);
  free(Q);
  free(Id);

}


void
ChangePairProbsIterate(FILE *ofp, double *pairprobs, int L, double *targetfreq, int hasindel, int verbose)
{
  double           *S;
  double           *Ql;
  double           *Qr;
  double           *pml;
  double           *pmr;
  double            lambda;
  int               dim;
  int               iterations = 0;
  int               i,j;

  if (hasindel) dim = L-1;
  else          dim = L;

  CalculateConditionalsAndMarginals(ofp, pairprobs, &Ql, &Qr, &pml, &pmr, L, hasindel, verbose);

  if (verbose) {
    fprintf(ofp, "Target marginals probabilities\n");
    PrintVectorProbs(ofp, targetfreq, L);
    fprintf(ofp, "Old PAIR Matrix \n");
    PrintProbs(ofp, pairprobs, L);
    fprintf(ofp, "Old Ql(i|j,t) probabilities\n");
    PrintProbs(ofp, Ql, L);
    fprintf(ofp, "Old Qr(i|j,t) probabilities\n");
    PrintProbs(ofp, Qr, L);
    fprintf(ofp, "Old marginals probabilities\n");
    PrintVectorProbs(ofp, pml, L);
    PrintVectorProbs(ofp, pmr, L);
  }
  
  S = (double *) MallocOrDie (sizeof(double) * L * L); 
  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) 
      S[i*L+j] = 0.0;
  
  while (CompareFreqs(pml, pmr, targetfreq, L)) {
    iterations ++;
    
    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++) 
	if (pml[i] > 0.0 && pmr[j] > 0.0 && pairprobs[i*L+j] > 0.0)
	  S[i*L+j] = log(pairprobs[i*L+j]) - log(pml[i]) - log(pmr[j]);
    
    lambda = Cal_lambda(ofp, S, targetfreq, targetfreq, L, verbose);
    
    for (i = 0; i < L; i++)
      for (j = 0; j < L; j++)
        pairprobs[i*L+j] = exp(lambda*S[i*L+j])*targetfreq[i]*targetfreq[j];

    /* check they are probabilities */  
    CheckSingleProb(pairprobs, L*L);
    
    ComputeLRConditionalsAndMarginals(ofp, pairprobs, Ql, Qr, pml, pmr, L, hasindel, verbose);
  }
  
  if (fabs(lambda) < MARGIN && iterations > 0) 
    Warn ("ChangePairProbsIterate(): lambda = 0 iterations = %d. trivial", iterations);
  
  if (verbose) {
    fprintf(ofp, "new PAIR (%d it) Matrix\n", iterations);
    PrintProbs(ofp, pairprobs, L);
    fprintf(ofp, "new Ql(i|j,t) probabilities\n");
    PrintProbs(ofp, Ql, L);
    fprintf(ofp, "new Qr(i|j,t) probabilities\n");
    PrintProbs(ofp, Qr, L);
    fprintf(ofp, "new marginals probabilities\n");
    PrintVectorProbs(ofp, pml, L);
    PrintVectorProbs(ofp, pmr, L);  
  }
  
  free(Ql);
  free(Qr);
  free(pml);
  free(pmr);
  
  free(S);
}

/* Function: ChangeRateMatrix()
 * Date:     ER, Wed Sep 17 15:38:25 CDT 2003 [St. Louis]
 *
 * Purpose:  Given a rate matrix R,  modify the
 *           rate matrix so it stays reversible respect to the new frequecies.
 *
 *           If this is the case of having indels as extra characters, we do that
 *           with only the subspace (L-1)x(L-1).
 *  
 * Method:   
 *    
 *                    (Q_0 K)_ij = S_ij (p_j/p_i)^{1/2}  i \neq j 
 *
 *          to  
 *                    (Q_0 K') _ij = S_ij (p'_j/p'_i)^{1/2}
 * *          so         K' = Q_0_inv * (Q_0 K')
 *
 * Args:              
 *
 * Returns:  (void)
 *           
 */
void
ChangeRateMatrix(FILE *ofp, double *Q_0R, double *Q_0_inv, double *R, double *ml, double *mr, int L, double *targetfreq, 
		 int hasindel, int isleft, int verbose)
{
  double *pstatl;
  double *pstatr;
  double *tdep_target;
  double *old_psat;
  double *new_psat;
  double *new_target;
  double *alpha;
  double *Rev;
  double *Id;
  double  comesfrom;
  double  goesto;
  double  sum;
  int     dim;
  int     i,j;

  Id = Cal_Id(L);

  if (hasindel) dim = L-1;
  else          dim = L;
  
  /* paranoia */
  CheckSingleProb(ml, L);
  CheckSingleProb(mr, L);
  CheckSingleProb(targetfreq, dim); 
  
  /* allocate */
  pstatl     = (double *) MallocOrDie (sizeof(double) * dim);
  pstatr     = (double *) MallocOrDie (sizeof(double) * dim);
  new_target = (double *) MallocOrDie (sizeof(double) * dim);

  old_psat = SaturationProbs(R, Id, L, FALSE, verbose);

 if (verbose) {
    fprintf(ofp, "Old Rate Matrix\n");
    PrintProbs(ofp, R, L);
    fprintf(ofp, "Old Q_0R Matrix\n");
    PrintProbs(ofp, Q_0R, L);
    fprintf(ofp, "Old Saturations probs\n");
    PrintVectorProbs(ofp, old_psat, L);
 }

  for (i = 0; i < dim; i++) {
    if (hasindel) pstatl[i] = ml[i]/(1.0-ml[L-1]);
    else          pstatl[i] = ml[i];

    if (hasindel) pstatr[i] = mr[i]/(1.0-mr[L-1]);
    else          pstatr[i] = mr[i];
  }
  
  alpha = (double *) MallocOrDie (sizeof(double) * dim * dim);

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) 
      alpha[i*dim+j] = 0.0;

  for (i = 0; i < dim; i++)
    for (j = i+1; j < dim; j++) {
      alpha[i*dim+j] = -0.5;      
      alpha[j*dim+i] = alpha[i*dim+j]; 
    }

 if (verbose) {
    fprintf(ofp, "alpha Matrix\n");
    PrintProbs(ofp, alpha, dim);
    fprintf(ofp, "stat L marginals\n");
    PrintVectorProbs(ofp, pstatl, L);
    fprintf(ofp, "stat R marginals\n");
    PrintVectorProbs(ofp, pstatr, L);
    fprintf(ofp, "target marginals\n");
    PrintVectorProbs(ofp, targetfreq, L);
  }
 
 /* Non diagonal "reduced" components of (Q_0*R) */
  for (i = 0; i < dim; i++) 
    for (j = 0; j < dim; j++) {

      if (isleft) {
	comesfrom = LOG2(targetfreq[i]) - LOG2(pstatl[i]);
	goesto    = LOG2(targetfreq[j]) - LOG2(pstatr[j]);
      }
      else {
	comesfrom = LOG2(targetfreq[i]) - LOG2(pstatr[i]);
	goesto    = LOG2(targetfreq[j]) - LOG2(pstatl[j]);
      }
      
     if (i != j ) Q_0R[i*L+j] = EXP2( LOG2(Q_0R[i*L+j]) + alpha[i*dim+j]*comesfrom + (1.0+alpha[i*dim+j])*goesto );
    }

  if (hasindel) {
    for (i = 0; i < dim; i++) {

      if (isleft) {
	comesfrom = LOG2(targetfreq[i]) - LOG2(pstatl[i]);
	goesto    = LOG2(targetfreq[i]) - LOG2(pstatr[i]);
      }
      else {
	comesfrom = LOG2(targetfreq[i]) - LOG2(pstatr[i]);
	goesto    = LOG2(targetfreq[i]) - LOG2(pstatl[i]);
      }

      Q_0R[dim*L+i] = -targetfreq[i] * Q_0R[dim*L+dim];
    }

  }

  /* Diagonal components of (Q_0*R) */
  for (i = 0; i < L; i++) { 
    Q_0R[i*L+i] = 0.0;
    
    for (j = 0; j < L; j++) 
      if (j != i) Q_0R[i*L+i] -= Q_0R[i*L+j];
  }
  /* consistency check, rows should add up to zero
   */	
  for (i = 0; i < L; i++) {
    sum = 0.0;
    for (j = 0; j < L; j++)
      sum += Q_0R[i*L+j];
    if (sum > MARGIN || sum < -MARGIN) Die("ChangeRateMatrix(): column %d bad Q_0R matrix (sum = %f)\n", i, sum);
  }
  
  if (verbose) {
    fprintf(ofp, "New Q_0R Matrix(L=%d)\n", L);
    PrintProbs(ofp, Q_0R, L);
  }
  
  /* Now calculate R = Q_0_inv * (Q_0*R)
   */
  Comp_M_N_Prod(ofp, Q_0_inv, Q_0R, L, verbose); /*  multiply Q_0_inv*M it in M */
  CopyMatrix(R, Q_0R, L, L);
  /* consistency check, rows should add up to zero
   */	
  for (i = 0; i < L; i++) {
    sum = 0.0;
    for (j = 0; j < L; j++)
      sum += R[i*L+j];
    if (sum > MARGIN || sum < -MARGIN) Die("ChangeRateMatrix(): column %d bad rate matrix (sum = %f)\n", i, sum);
  }

  /* Calculate the new saturation probabilities
   */
  new_psat = SaturationProbs(R, Id, L, FALSE, verbose);

  if (verbose) {
    fprintf(ofp, "new_sat  probs\n");
    PrintVectorProbs(ofp, new_psat, L);
  }

  for (i = 0; i < dim; i++) 
    if (hasindel) {
      if (fabs(1.0-new_psat[dim]) > 1.0-accuracy1) new_target[i] = new_psat[i] / (1.0 - new_psat[dim]);
      else                                         new_target[i] = targetfreq[i];
    }
    else 
      new_target[i] = new_psat[i];
  
  if (verbose) {
    fprintf(ofp, "New Rate Matrix(L=%d)\n", L);
    PrintProbs(ofp, R, L);
    fprintf(ofp, "New Saturations probs\n");
    PrintVectorProbs(ofp, new_psat, L);
    fprintf(ofp, "Target  probs\n");
    PrintVectorProbs(ofp, targetfreq, dim);
    fprintf(ofp, "Found Target  probs\n");
    PrintVectorProbs(ofp, new_target, dim);
  }
  
  /* paranoia */
  if (hasindel && fabs(new_psat[dim]-old_psat[dim]) > MARGIN) 
    Die ("ChangeRateMatrix():\nbad change of background frequencies, with indels, last component should not change");
  
  if (hasindel) {
    if (fabs(1.0-new_psat[dim]) > MARGIN) 
      for (i = 0; i < dim; i++) 
	if (fabs(new_target[i]-targetfreq[i]) > MARGIN) 
	  Warn ("ChangeRateMatrix():\nbad change of background frequencies -- indels, i=%d new_target=%f targetfreq=%f", i, new_target[i], targetfreq[i]);
  }
  else {
    for (i = 0; i < dim; i++) 
      if (fabs(new_target[i]-targetfreq[i]) > MARGIN) 
	Die ("ChangeRateMatrix():\nbad change of background frequencies -- no indels, i=%d new_target=%f targetfreq=%f", i, new_target[i], targetfreq[i]);
  }

  if (verbose) {
    Rev         = (double *) MallocOrDie (sizeof(double) * L * L);
    tdep_target = (double *) MallocOrDie (sizeof(double) * L);
    
    for (i = 0; i < L; i++) 
      if (hasindel) {
	if (i < dim) 
	  if (isleft) tdep_target[i] = targetfreq[i] * (1.0 - ml[L-1]); 
	  else        tdep_target[i] = targetfreq[i] * (1.0 - mr[L-1]);
	else         
	  if (isleft) tdep_target[i] = ml[L-1]; 
	  else        tdep_target[i] = mr[L-1];
      }
      else 
	tdep_target[i] = targetfreq[i];
    
    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++)
	Rev[i*L+j] = Q_0R[i*L+j] * tdep_target[i];
	
    fprintf(ofp, "(Q_0*R * tdep_target) reversibility\n");
    PrintProbs(ofp, Rev, L);
    fprintf(ofp, "pi(t)\n");
    PrintVectorProbs(ofp, tdep_target, L);

    for (i = 0; i < L; i++) 
      if (hasindel) {
      if (i < dim) 
	if (isleft) tdep_target[i] = new_target[i] * (1.0 - ml[L-1]); 
	else        tdep_target[i] = new_target[i] * (1.0 - mr[L-1]);
      else         
	if (isleft) tdep_target[i] = ml[L-1]; 
	else        tdep_target[i] = mr[L-1];
      }
      else 
	tdep_target[i] = new_target[i];
    
    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++)
	Rev[i*L+j] = Q_0R[i*L+j] * tdep_target[i];
	
    fprintf(ofp, "(Q_0*R * newtarget) reversibility\n");
    PrintProbs(ofp, Rev, L);
    fprintf(ofp, "pi(t)\n");
    PrintVectorProbs(ofp, tdep_target, L);

    free(Rev);
    free(tdep_target);
  }
  
  free(pstatl);
  free(pstatr);
  free(old_psat);
  free(new_psat);
  free(new_target);
  free(alpha);
  free(Id);
}


/* Function: Condi_From_Joint()
 * Date:     ER, Wed Mar  1 11:04:24 CST 2000 [St. Louis]
 *
 * Purpose:  Given       a LL matrix M(i,j) of joint probabilities, 
 *           calculate the LL matrix Q(i,j) of conditional probabilities.
 *
 *           Q[i*L+j] = Q(j|i) = P(i,j|t_o) / \sum_k  P(i,k|t_o) 
 *
 *           where P(i,j|t_o) are calculated according to 
 *           using a given pammodel t_o.
 *
 *           Notice that if P(i,j) are not real joint probs (ie non symmetric)
 *           this function calculates the RIGHT conditionals and LEFT marginals.
 *
 * Args:     P  - LxL joint prob matrix (prealloc)
 *
 * Returns:  Q(LxL), conditional probabilities.
 *           Q is alocated here, freed by caller.
 */
double *
Condi_From_Joint(FILE *ofp, double *P, int L, int verbose)
{
  double *Q;
  double *pm;
  double  sum;
  int     i, j;
  
  if (verbose) {
    fprintf(ofp, "P(i,j, t_o) probabilities\n");
    PrintProbs(ofp, P, L);
  }
  
  /* allocate Q[L*L]
   */
  Q  = (double *) MallocOrDie (sizeof(double) * L * L);
  pm = (double *) MallocOrDie (sizeof(double) * L);

  for (i = 0; i < L; i++)
    pm[i] = 0.0;
  
  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) 
      pm[i] += P[i*L+j];
  
  for (i = 0; i < L; i++) {
    sum = 0.;
    for (j = 0; j < L; j++) {
      Q[i*L+j] = P[i*L+j] / pm[i];
      sum += Q[i*L+j];
    }
    /* normalize */
    for (j = 0; j < L; j++) 
      Q[i*L+j] /= sum;
  }
 
  if (verbose) { 
    fprintf(ofp, "Q(j|i, t_o) probabilities\n");
    PrintProbs(ofp, Q, L);
  }
  
  /* consistency check, rows should add up to one
   */
  for (i = 0; i < L; i++)
    CheckSingleProb(Q+i*L, L);
  
  free (pm);

  return Q;
}

double *
CondiR_From_Joint(FILE *ofp, double *P, int L, int verbose)
{
  double *QR;
  double *pmr;
  double  sum;
  int     i, j;
  
  if (verbose) {
    fprintf(ofp, "P(i,j, t_o) probabilities\n");
    PrintProbs(ofp, P, L);
  }
  
  /* allocate Q[L*L]
   */
  QR  = (double *) MallocOrDie (sizeof(double) * L * L);
  pmr = (double *) MallocOrDie (sizeof(double) * L);

  for (i = 0; i < L; i++)
    pmr[i] = 0.0;
  
  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) 
      pmr[i] += P[j*L+i];
  
  for (i = 0; i < L; i++) {
    sum = 0.;
    for (j = 0; j < L; j++) {
      QR[i*L+j] = P[j*L+i] / pmr[i];
      sum += QR[i*L+j];
    }
    /* normalize */
    for (j = 0; j < L; j++) 
      QR[i*L+j] /= sum;
  }
 
  if (verbose) { 
    fprintf(ofp, "QR(j|i, t_o) probabilities\n");
    PrintProbs(ofp, QR, L);
  }
  
  /* consistency check, rows should add up to one
   */
  for (i = 0; i < L; i++)
    CheckSingleProb(QR+i*L, L);
  
  free (pmr);

  return QR;
}

/* Function: Cond_From_Rate()
 * Date:     ER, Thu Sep 25 09:38:36 CDT 2003 [St. Louis]
 *
 * Purpose:  Given a Rate R 
 *           and set of conditional probabilities at time zero Q_0
 *         
 *           calculate coditionals; Q(t) = Q_0 * exp(tR)
 *  
 * Args:              
 *
 * Returns:  void
 *           
 */
void 
Condi_From_Rate(double **ret_Q, double *R, double *Q_0, double tfactor, int L, int pedantic, int verbose)
{
  double *Qr;
  double *Q;
  double  tmax;
  double  xtime;
  int     n_times;
  int     i;
  
  if (verbose) {
    fprintf(stdout, "Rate(L=%d)\n", L);
    PrintProbs(stdout, R, L);
  }

  tmax = TMAX;

  n_times = (int)tfactor / (int)tmax;
  xtime   = tfactor - tmax*(double)n_times;

  /* calculate Q(t) = Q_0 * exp{tfactor*K}
   */
  Q = (double *)Cal_Id(L);
  while (n_times > 0) {
    Qr = Cal_M_Exp(stdout, R, L, tmax, verbose);
    Comp_M_N_Prod(stdout, Qr, Q, L, FALSE); /*  multiply Qr^(n_times*tmax), dump it in Qp */
    n_times --;
    free (Qr);
  }
  Qr = Cal_M_Exp(stdout, R, L, xtime, verbose);
  Comp_M_N_Prod(stdout, Qr, Q, L, verbose); /*  multiply Qr^xtime, dump it in Qp */
 
  Comp_M_N_Prod(stdout, Q_0, Q, L, verbose); /* multiply Q_0*Qp, dump it in Qp */
  
  if (verbose) {
    fprintf(stdout, "Q(L=%d, t=%f)\n", L, tfactor);
    PrintProbs(stdout, Q, L);
  }

  for (i = 0; i < L; i++) 
    QOMRegularizationAlgorithm(Q+i*L, L, verbose);
  
  if (verbose) {
    fprintf(stdout, "Q(L=%d, t=%f)\n", L, tfactor);
    PrintProbs(stdout, Q, L);
  }

  /* consistency check, rows should add up to one
   */	
  for (i = 0; i < L; i++)
    CheckSingleProb(Q+i*L, L);
 
  *ret_Q = Q;

  free (Qr);
}

void
ConditionalsFromRate(double *Q, double *R, double *Q_0, double tfactor, int L, int pedantic, int verbose)
{
  double *pcond;

  Condi_From_Rate(&pcond, R, Q_0, tfactor, L, pedantic, verbose);

  CopyMatrix(Q, pcond, L, L);

  free(pcond);
}

void
ZCondi_From_Rate(double **ret_Q, struct zmatrix_s *R, struct zmatrix_s *Q_0, double tfactor, int L, int pedantic, int verbose)
{
  struct zmatrix_s *Qr;
  struct zmatrix_s *Q;
  double           *Qreal;
  double            tmax;
  double            xtime;
  int               n_times;
  int               i,j;
  
  tmax = TMAX;

  n_times = (int)tfactor / (int)tmax;
  xtime   = tfactor - tmax*(double)n_times;

  /* calculate Q(t) = Q_0 * exp{tfactor*K}
   */
  Q = (struct zmatrix_s *)ZCal_Id(L);
  while (n_times > 0) {
    Qr = ZCal_M_Exp(stdout, R, L, tmax, verbose);
    ZComp_M_N_Prod(stdout, Qr, Q, L, FALSE); /*  multiply Qr^(n_times*tmax), dump it in Qp */
    n_times --;
    FreeZmatrix (Qr);
  }
  Qr = ZCal_M_Exp(stdout, R, L, xtime, verbose);
  ZComp_M_N_Prod(stdout, Qr, Q, L, verbose); /*  multiply Qr^xtime, dump it in Qp */
 
  ZComp_M_N_Prod(stdout, Q_0, Q, L, verbose); /* multiply Q_0*Qp, dump it in Qp */
  
  if (verbose) {
    fprintf(stdout, "Conditional Zmatrix\n");
    PrintZmatrix(stdout, Q, L);
  }

  /* Is the conditional matrix real ?
   */
  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++)
      if (fabs(Q->imag[i*L+j]) > MARGIN2) 
	Die ("ZCondi_From_Rate(): your conditional matrix is not real!");
 
  Qreal = (double *) MallocOrDie (sizeof(double) * L * L); 
  CopyMatrix (Qreal, Q->real, L, L);

  for (i = 0; i < L; i++)
    QOMRegularizationAlgorithm(Qreal+i*L, L, verbose);

   /* consistency check, rows should add up to one
   */	
  for (i = 0; i < L; i++)
    CheckSingleProb(Qreal+i*L, L);
  
 
  *ret_Q = Qreal;

  FreeZmatrix(Q);
  FreeZmatrix(Qr);
}


/* Function: ConditionalsEvolved()
 * Date:     ER, Wed Apr 12 18:11:23 CDT 2000 [St. Louis]
 *
 * Purpose:  Given a set of conditionals probabilities Q(a|b, t) 
 *           at a given evolutionary time t, 
 *           and set of conditional probabilities at time zero Q_0
 *         
 *           calculate Q(a,b| r*t) using a Markov chain model for evolution
 *  
 * Method:   Q(a | b, r*t) = Q_0 * exp{r * log[Q_0^{-1} * Q}
 *
 *
 *
 *
 * Args:              
 *
 * Returns: Qr(LxL). 
 *          Qr is allocated here, freed by caller.
 *           
 */
void
ConditionalsEvolved(FILE *ofp, double *Q, double *QR, double *Q_0, double *QR_0, double *ml, double *mr, 
		    int L, double tfactor, double *targetfreq, int changefreq, int hasindel, int pedantic, int verbose)
{
  double  *Q_0_inv;
  double  *QR_0_inv;
  double  *Q_0K;
  double  *QR_0KR;
  double  *Qp;
  double  *QRp;
  double  *K;
  double  *KR;
  double  *S;
  double   tmax;
  double   xtime;
  int      n_times;
  int      i;

  tmax = TMAX;

  n_times = (int)tfactor / (int)tmax;
  xtime   = tfactor - tmax*(double)n_times;

  /* consistency check, rows should add up to one
   */	
  for (i = 0; i < L; i++)  CheckSingleProb(Q+i*L, L);
  
 if (verbose) {
    fprintf(ofp, "QL(j|i, t_*)  probabilities\n");
    PrintProbs(ofp, Q, L);
    fprintf(ofp, "QR(j|i, t_*)  probabilities\n");
    PrintProbs(ofp, QR, L);
    fprintf(ofp, "QL_0(j|i)  probabilities\n");
    PrintProbs(ofp, Q_0, L);
    fprintf(ofp, "QR_0(j|i)  probabilities\n");
    PrintProbs(ofp, QR_0, L);
    fprintf(ofp, "LEFT Marginal  probabilities\n");
    PrintVectorProbs(ofp, ml, L);
    fprintf(ofp, "RIGHT  Marginal  probabilities\n");
    PrintVectorProbs(ofp, mr, L);  
  }
 
 if (verbose) {
   check_Q_0_reversibility(Q_0, QR_0, ml, mr, L, hasindel);
   check_reversibility    (Q,   QR,   ml, mr, L);
 }
 
 /* calculate Q_0_inv
  */
  Q_0_inv  = Cal_M_Inv(ofp, Q_0,  L, verbose);
  QR_0_inv = Cal_M_Inv(ofp, QR_0, L, verbose);

  /* construct K = log [Q_0_inv * Q]
   */
  Comp_M_N_Prod(ofp, Q_0_inv,  Q,  L, verbose); /*  multiply Q_0_inv*Q, dump it in Q */
  Comp_M_N_Prod(ofp, QR_0_inv, QR, L, verbose); /*  multiply Q_0_inv*Q, dump it in Q */

  S = Cal_Id(L);
  Comp_M_N_Sum(ofp, Q, S, L, L, FALSE, verbose);   /* S = Q - Id  */
  
  /* check that Q is consistent with a markovian stationary 
   * model of evolution. 
   * In plain words it checks whether logQ is going to exist.
   *
   * Basically, it checks whether  all eigenvalues of Q - I
   * are between -1 and 1. 
   *
   * remember log (1+x) = x - x^2/2 + ...   converges for -1 < x <= 1
   *
   */
  
  if (verbose) IsQConsistent(ofp, S, L, pedantic, verbose);
  K  = Cal_M_Taylor_Log(ofp, Q,  L, verbose);
  KR = Cal_M_Taylor_Log(ofp, QR, L, verbose);
  
  /* if (IsQConsistent(ofp, S, L, pedantic, verbose)) { 
     K  = Cal_M_Taylor_Log(ofp, Q,  L, verbose);
     KR = Cal_M_Taylor_Log(ofp, QR, L, verbose);
     }
     else {
     K  = Cal_M_Vandermonde_Log(ofp, Q,  L, verbose);
     KR = Cal_M_Vandermonde_Log(ofp, QR, L, verbose);
     }*/
  
  if (verbose) {
    fprintf(ofp, "R  rate matrix\n");
    PrintProbs(ofp, K, L);
  }
 
  if (pedantic)
    IslogQConsistent(ofp, Q_0, K, L, L, L, verbose); /*  check that K is consistent with a 
						      *  markovian stationary model of evolution. 
						      *
						      *  Basically, it checks whether all entries of  
						      *  Q_0 + epsilon*Q_0*K are none negatives.
						      */

  /* If we are told to change the stationary frequencies we do it here.
   * This will also modify the rate matrix K, while keeping the exchangeability parameters S
   * defined as
   *            (Q_0 K)_ij  = S_ij (p_j/p_i)^{1/2}  i \neq j 
   *
   *   to  
   *            (Q_0 K')_ij = S_ij (p'_j/p'_i)^{1/2}
   *
   *    so       K' = Q_0_inv * (Q_0 K')
   *
   */
  if (changefreq) {
    
    Q_0K   = Cal_M_N_Prod(ofp, Q_0,  K,  L, L, L, verbose);
    QR_0KR = Cal_M_N_Prod(ofp, QR_0, KR, L, L, L, verbose);
    
    ChangeQ_0MatrixIterate(ofp, Q_0,  Q_0_inv,  ml, mr, L, targetfreq, hasindel, TRUE,  verbose);
    ChangeQ_0MatrixIterate(ofp, QR_0, QR_0_inv, ml, mr, L, targetfreq, hasindel, FALSE, verbose);
    
    ChangeRateMatrix(ofp, Q_0K,   Q_0_inv,  K,  ml, mr, L, targetfreq, hasindel, TRUE,  verbose);
    ChangeRateMatrix(ofp, QR_0KR, QR_0_inv, KR, ml, mr, L, targetfreq, hasindel, FALSE, verbose);
 
    free(Q_0K);
    free(QR_0KR);
  }
  
   /* if this process involves indels, we have to evolve the marginal probabilities too
    */
  if (hasindel) {
    EvolveIndelMarginals(stdout, K,  Q_0,  mr, tfactor, L, FALSE, verbose);
    EvolveIndelMarginals(stdout, KR, QR_0, ml, tfactor, L, FALSE, verbose);
  }

  if (changefreq) {
    ChangeFrequencies(ofp, Q_0,  K,  mr, L, targetfreq, hasindel, verbose); /* Left  conditionals saturate to right marginals */
    ChangeFrequencies(ofp, QR_0, KR, ml, L, targetfreq, hasindel, verbose); /* Right conditionals saturate to left  marginals */
  }
  
  /* calculate Qp = Q_0 * exp{tfactor*K}
   */
  Condi_From_Rate(&Qp,  K,  Q_0,  tfactor, L, pedantic, verbose);
  Condi_From_Rate(&QRp, KR, QR_0, tfactor, L, pedantic, verbose);
  CopyMatrix(Q,  Qp,  L, L);
  CopyMatrix(QR, QRp, L, L);
  
  /* consistency check, rows should add up to one
   */	
  for (i = 0; i < L; i++)
    CheckSingleProb(Q+i*L, L);
  for (i = 0; i < L; i++)
    CheckSingleProb(QR+i*L, L);
  
  if (verbose) {
    fprintf(ofp, "Q(j|i, r*t_o=%f)  probabilities\n", tfactor);
    PrintProbs(ofp, Q, L);
    fprintf(ofp, "QR(j|i, r*t_o=%f)  probabilities\n", tfactor);
    PrintProbs(ofp, QR, L);
  }
  
  if (verbose) {
   check_Q_0_reversibility(Q_0, QR_0, ml, mr, L, hasindel);
   check_reversibility    (Q,   QR,   ml, mr, L);
  }
  
  
  free(S);
  free(K);
  free(Q_0_inv);
  free(Qp);
  free(QRp);
}



/* Function: ConditionalsEvolved_2()
 * Date:     ER, Mon May  8 10:33:18 CDT 2000 [St. Louis]
 *
 * Purpose:  Given a set of conditionals probabilities Q(a|b, t) 
 *           at a given evolutionary time t, 
 *           and set of conditional probabilities at time zero Q_0
 *         
 *           calculate Q(a,b| r*t) 
 *           this is a more general model than the previous one:
 *  
 * Method:   Q(a | b, r*t) = Q_0 +  exp{ r * log[I + R^{-1}(Q-Q_0)] }
 *
 *           for once, the property Q(t+s) = Q(t)*Q(s) is not satisfied any more.
 *
 * Args:              
 *
 * Returns: Qr(LxL).  *          Qr is allocated here, freed by caller.
 *           
 */
void
ConditionalsEvolved_2(FILE *ofp, double *Q, double *Q_0, double *R, int L, double tfactor, double *targetfreq, 
		      int changefreq, int pedantic, int verbose)
{
  double  *R_inv;
  double  *Qr;
  double  *Qp;
  double  *K;
  double  *Id;
  double   tmax;
  double   xtime;
  int      n_times;
  int      i;

  tmax = TMAX;

  n_times = (int)tfactor / (int)tmax;
  xtime   = tfactor - tmax*(double)n_times;

  /* consistency check, rows should add up to one
   */	
  for (i = 0; i < L; i++)
    CheckSingleProb(Q+i*L, L);
  
  Id  = Cal_Id(L);
  
  /* calculate R_inv 
   */
  R_inv = Cal_M_Inv(ofp, R, L, verbose);
  
  /* calculate I + R_inv*(Q-Q_0)  dump it in Q 
   */
  Comp_M_N_Sum(ofp, Q, Q_0, L, L, FALSE, verbose);    /*  add:      Q - Q_0,   dump it in Q */
  CopyMatrix(Q, Q_0, L, L);
  Comp_M_N_Prod(ofp, R_inv, Q, L, verbose);           /*  multiply: R_inv * Q, dump it in Q */

  if (pedantic)
    IsQConsistent(ofp, Q, L, pedantic, verbose);  /* check that Q is consistent with a markovian stationary 
						   * modelof evolution. 
						   * In plain words it checks whether log(I+Q) is going to exist.
						   *
						   * Basically, it checks whether  all eigenvalues of Q
						   * are between -1 and 1. 
						   *
						   * remember log (1+x) = x - x^2/2 + ...   converges for -1 < x <= 1
						   *
						   */

  Comp_M_N_Sum(ofp, Id, Q, L, L, TRUE, verbose);   /*  add:      Id + Q,    dump it in Q */

  /* construct K = log [I + R_inv*(Q-Q_0)]   */
  K = Cal_M_Taylor_Log(ofp, Q, L, FALSE);
  
  if (pedantic)
    IslogQConsistent_2(ofp, Q_0, R, K, L,  L,  L, verbose); /*  check that K is consistent with a 
							     *  markovian stationary model of evolution. 
							     *
							     *  Basically, it checks whether all entries of  
							     *  Q_0 + epsilon*R*K are none negatives.
							     */
  
  /* calculate Qr = Q_0 + R * [ exp{tfactor*K} - Id ]
   */
  Qp = (double *)Cal_Id(L);
  while (n_times > 0) {
    Qr = Cal_M_Exp(ofp, K, L, tmax, verbose);
    Comp_M_N_Prod(ofp, Qr, Qp, L, FALSE); /*  multiply Qr^(n_times*tmax), dump it in Qp */
    n_times --;
    free (Qr);
  }
  Qr = Cal_M_Exp(ofp, K, L, xtime, verbose);
  Comp_M_N_Prod(ofp, Qr, Qp, L, FALSE); /*  multiply Qr^xtime, dump it in Qp */
  
  Comp_M_N_Prod(ofp, R, Qr, L, FALSE);               /* multiply: R*Qp,   dump it in Qp */
  Comp_M_N_Sum(ofp, Q_0, Qr, L, L, TRUE, FALSE);     /* add:      Q_0+Qp, dump it in Qp */
  
  if (verbose) {
    fprintf(ofp, "Q(j|i, r*t_o)  probabilities\n"); 
    PrintProbs(ofp, Qp, L);
  }
  
  CopyMatrix(Q, Qp, L, L);

  /* consistency check, rows should add up to one
   */	
  for (i = 0; i < L; i++)
    CheckSingleProb(Q+i*L, L);
  
  free(K);
  free(Id);
  free(R_inv);
  free(Qr);
  free(Qp);

}

/* Function: EvolveIndelMarginals()
 * Date:     ER, Wed Sep 17 14:35:07 CDT 2003 [St. Louis]
 *
 * Purpose:  Given       a LxL matrix R(i,j)   of rates              , 
 *                       a LxL matrix Q_0(i,j) of conditional probabilities at t=0
 *
 *                       a L vector of marginal probs of the form   pmar    = [p_i(1-lambda*), lamda*]
 *
 *           Use R and Q_0 to calculate the marginals at time infty p_infty = [p_i(1-lambda_infty), lamda_infty] 
 *
 *           evolve lambda(t) = lambda_infty[1-exp{t*log(1-lambda*\/lamda_infty)]
 *
 *           construct the evolved marginals as pmar(t) = [p_i(1-lambda(t)), lamda(t)]
 *
 *
 * Returns:  void.
 */
void
EvolveIndelMarginals(FILE *ofp, double *R, double *Q_0, double *pm, double tfactor, int L, int pedantic, int verbose)
{
  double *pinfty;
  double  lambda_star;
  double  lambda_infty;
  int     i;

  pinfty = (double *) MallocOrDie (sizeof(double) * L); 

   if (verbose) {
    fprintf(ofp, "RATE\n");
    PrintProbs(ofp, R, L);
    fprintf(ofp, "t^* marginals\n");
    PrintVectorProbs(ofp, pm, L);
  }

  /* saturation probabilitie
   */
  pinfty = SaturationProbs(R, Q_0, L, FALSE, verbose);
 
  /* evolve the sigle probs -- only evolve the last component
   */
  lambda_star  = pm[L-1];
  lambda_infty = pinfty[L-1];
  
  if (verbose) {
    fprintf(ofp, "EvolveIndelMarginals()\n");
    fprintf(ofp, "Lambda_star = %f\n", lambda_star);
    fprintf(ofp, "Lambda_infty = %f\n", lambda_infty);
    fprintf(ofp, "INFTY-marginals\n");
    PrintVectorProbs(ofp, pinfty, L);
  }

  if (lambda_star > lambda_infty) 
    Die ("EvolveIndelMarginals(): bad marginals: lambda_star = %.10f lambda_infty = %.10f", lambda_star, lambda_infty);

  pm[L-1] = pinfty[L-1]*(1.0-EXP2(tfactor*LOG2(1.0-lambda_star/lambda_infty)));
  
  for (i = 0; i < L-1; i++) 
      pm[i] *= (1.0-pm[L-1]) / (1.0-lambda_star);
 
  CheckSingleProb(pm, L);
  
  if (verbose) {
    fprintf(ofp, "marginals evolved %f time\n", tfactor);
    PrintVectorProbs(ofp, pm, L);
    fprintf(ofp, "INFTY-marginals\n");
    PrintVectorProbs(ofp, pinfty, L);
  }

  free(pinfty);
}


/* Function: Joint_From_Condi()
 * Date:     ER, Mon Jun 12 14:41:47 CDT 2000 [St. Louis]
 *
 * Purpose:  Given       a LxL matrix Q(i,j) of conditional probabilities, 
 *           calculate the LxL matrix P(i,j) of joint       probabilities.
 *
 *           P(i,j) = Q(j|i) * pm[j] 
 *
 *           where P(i,j|t_o) are calculated according to 
 *           using a given pammodel t_o.
 *
 * Args:     Q  - LxL conditional prob matrix (prealloc)
 *           pm - L   single emision probabilities (prealloc)
 *
 * Returns:  P(LxL), conditional probabilities.
 *           P is alocated here, freed by caller.
 */
void
Joint_From_Condi(FILE *ofp, double *P, double *Q, double *pm, int L, int verbose)
{
  double *Qnewl;
  double *Qnewr;
  double *pml;
  double *pmr;
  double  norm = 0.0;
  int     i, j;
  int     hasindel;
  
  if (verbose) {
    fprintf(ofp, "Q(i|j,t) probabilities\n");
    PrintProbs(ofp, Q, L);
    fprintf(ofp, "marginals probabilities\n");
    PrintVectorProbs(ofp, pm, L);
  }
  
  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) {
      P[i*L+j] = Q[i*L+j] * pm[i]; 
      norm += P[i*L+j];
    }

  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) 
      P[i*L+j] /= norm;
  
  if (verbose) { 
    fprintf(ofp, "P(i,j,t) probabilities\n");
    PrintProbs(ofp, P, L);

    if (L%4) hasindel = TRUE;
    else     hasindel = FALSE;

    CalculateConditionalsAndMarginals(ofp, P, &Qnewl, &Qnewr, &pml, &pmr, L, hasindel, verbose);
    fprintf(ofp, "newQl(i|j,t) probabilities\n");
    PrintProbs(ofp, Qnewl, L);
    fprintf(ofp, "newQr(i|j,t) probabilities\n");
    PrintProbs(ofp, Qnewr, L);
    fprintf(ofp, "marginals probabilities\n");
    PrintVectorProbs(ofp, pml, L);
    PrintVectorProbs(ofp, pmr, L);
    
    free(Qnewl);
    free(Qnewr);
    free(pml);
    free(pmr);
  }
  
  /* consistency check
   */
  CheckSingleProb(P, L*L);
  
}


 
void
Joint_From_CondiR(FILE *ofp, double *P, double *QR, double *pr, int L, int verbose)
{
  double *Qnewl;
  double *Qnewr;
  double *pml;
  double *pmr;
  int     i, j;
  
  if (verbose) {
    fprintf(ofp, "QR(i|j,t) probabilities\n");
    PrintProbs(ofp, QR, L);
    fprintf(ofp, "Rmarginals probabilities\n");
    PrintVectorProbs(ofp, pr, L);
  }
  
  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++)
      P[i*L+j] = QR[j*L+i] * pr[j]; 
  
  if (verbose) { 
    fprintf(ofp, "P(i,j,t) probabilities\n");
    PrintProbs(ofp, P, L);

    CalculateConditionalsAndMarginals(ofp, P, &Qnewl, &Qnewr, &pml, &pmr, L, FALSE, verbose);
    fprintf(ofp, "newQl(i|j,t) probabilities\n");
    PrintProbs(ofp, Qnewl, L);
    fprintf(ofp, "newQr(i|j,t) probabilities\n");
    PrintProbs(ofp, Qnewr, L);
    fprintf(ofp, "marginals probabilities\n");
    PrintVectorProbs(ofp, pml, L);
    PrintVectorProbs(ofp, pmr, L);
    
    free(Qnewl);
    free(Qnewr);
    free(pml);
    free(pmr);
  }
  
  /* consistency check
   */
  CheckSingleProb(P, L*L);
  
}
 
void
Joint_From_Condi_Symmetrize(FILE *ofp, double *P, double *Q, double *p, int L, int verbose)
{
  int  i, j;

  Joint_From_Condi (ofp, P, Q, p, L, verbose);
  
  for (i = 0; i < L; i ++) 
    for (j = i+1; j < L; j ++) 
      {
	P[i*L+j] = 0.5 * (P[i*L+j] + P[j*L+i]);
	P[j*L+i] = P[i*L+j];
      }

  if (verbose) {
    fprintf(ofp, "\nPSYM(i,j) Joint probabilities\n");
    PrintProbs(stdout, P, L);
  }

}

/* Function: Joint2Joint()
 * Date:     ER, Wed Apr 12 18:11:23 CDT 2000 [St. Louis]
 *
 * Purpose:  Given a set of joint probabilities P(a,b| t_o) 
 *           at a given evolutionary time t_o, 
 *           and set of conditional probabilities at time zer0 P^o
 *         
 *           calculate P(a,b| tfactor*t_o) using a Markov chain model for evolution
 *  
 * Method:   (1) calculate conditionals: P(a,b | t_o) --> Q(a | b,t_o) 
 *            
 *           (2) calculate marginals: P(a,b | t_o) --> pmar(a | b,t_o)  
 *            
 *           (3) evolve conditionals Q(a | b, r*t_o) = Q_0 * exp{ r * log[Q_0^{-1} * Q] }
 *
 *           (4) construct joints using same marginals Q(a | b, r*t_o), pmar(a | b,t_o)  --> P(a,b | r*t_o)
 *
 * Args:              
 *
 * Returns:  (void)
 *           
 */
void
Joint2Joint(double *p_joint, double *ql_nought, double *qr_nought, int L, double tfactor, double *targetfreq, 
	    int changefreq, int hasindel, int ispairprob, int pedantic, int verbose)
{
  double *p_margl;
  double *p_margr;
  double *p_condl;
  double *p_condr;
  
  if (verbose) {
    fprintf(stdout, "p_joint(i,j,t) probabilities\n");
    PrintProbs(stdout, p_joint, L);
  }

  /* calculate the matrix of conditional probabilities (Q) and marginals from joint probabilities (P)
   */
  CalculateConditionalsAndMarginals(stdout, p_joint, &p_condl, &p_condr, &p_margl, &p_margr, L, FALSE, verbose);

  /* evolve conditionals, change the marginals if asked to do so
   */
  ConditionalsEvolved(stdout, p_condl, p_condr, ql_nought, qr_nought, p_margl, p_margr, L, tfactor, 
		      targetfreq, changefreq, hasindel, pedantic, verbose);
 
  /* Reconstruct the new evolved Joint probabilities
   */
  Joint_From_Condi(stdout, p_joint, p_condl, p_margl, L, verbose);
  
  free(p_margl);
  free(p_margr);
  free(p_condl);
  free(p_condr);
}


/* Function: Joint2JointGap()
 * Date:     ER, Tue Jan 14 11:24:14 CST 2003 [St. Louis]
 *
 * Purpose:  Given a set of joint probabilities P(a,b| t_*)
 *           at a given evolutionary time t_*,
 *           the set of conditional probabilities at time zero P^o,
 *           and the set of single probs at the time of interests pv
 *
 *           calculate P(a,b| tfactor*t_o) using a Markov chain model for evolution.
 *           For probabilities that include gaps.
 *
 * Method:   (1) calculate conditionals: P(a,b | t_*) --> Q(a | b,t_*)
 *
 *           (4) evolve conditionals Q(a | b, r*t_*) = Q_0 * exp{ r * log[Q_0^{-1} * Q] }
 *
 *           (6) construct joints with new probs Q(a | b, r*t_*), pv(a | b, r*t_*)  --> P(a,b | r*t_*)
 *
 * Args:
 *
 * Returns:  (void)
 *
 */
void
Joint2JointGap(double *p_joint, double *q_nought, double *pv, int L, double tfactor, double *targetfreq, 
	       int changefreq,  int pedantic, int verbose)
{
  double *p_cond;
  double *p_condr;
  double *q_noughtr;
  double *pvr;

  /* calculate the matrix of conditional probabilities (Q) from joint probabilities (P)
   */
  p_cond  = Condi_From_Joint (stdout, p_joint, L, verbose);
  p_condr = CondiR_From_Joint(stdout, p_joint, L, verbose);

  q_noughtr = (double *) MallocOrDie (sizeof(double) * L * L);
  pvr       = (double *) MallocOrDie (sizeof(double) * L);
  CopyMatrix(q_noughtr, q_nought, L, L);
  CopyVector(pvr, pv, L);

  /* evolve conditionals
   */
  ConditionalsEvolved(stdout, p_cond, p_condr, q_nought, q_noughtr, pv, pvr, L, tfactor, targetfreq, changefreq, FALSE, pedantic, verbose);

  Joint_From_Condi(stdout, p_joint, p_cond, pv, L, verbose);

  CheckSingleProb(p_joint, L*L);

  if (verbose) {
    fprintf(stdout, "P_COND(i,j| r*t_o=%f)  Cond probabilities\n", tfactor);
    PrintProbs(stdout, p_cond, L);
    fprintf(stdout, "Psingle(i| r*t_o=%f) single nt probabilities\n", tfactor);
    PrintVectorProbs(stdout, pv, L);
    fprintf(stdout, "P_JOINT(i,j| r*t_o=%f)  Joint probabilities\n", tfactor);
    PrintProbs(stdout, p_joint, L);
 }

  free(pvr);
  free(p_cond);
  free(p_condr);
  free(q_noughtr);
}


/* Function: MarginalizeJointProbs()
 * Date:     ER, Wed Aug 21 13:19:04 CDT 2002 [St. Louis]
 *
 * Purpose:  Given a jp[idx(dim,dim)] joint probability distrubution, 
 *           marginalize to calculate sp[dim].
 *
 * Args:    jp - dimxdim A..UxA..U joint prob matrix (prealloc)
 *          sp - dim seqX emission prob (prealloc)
 *
 * Returns:  (void)
 *           Fills in sp (already allocated)
 */
void
MarginalizeJointProbs(double *jp, double *sp, int dim, int which)
{
  int x;
  int y;

  /* Zero null model
   */
  for (x = 0; x < dim; x++)
    sp[x] = 0.0;

  /* Marginalize and average over Y positions
   */
  for (x = 0; x < dim; x++)
    for (y = 0; y < dim; y++)
      if      (which == 0) sp[x] += jp[x*dim+y];
      else if (which == 1) sp[x] += jp[y*dim+x];
      else if (which == 2) sp[x] += 0.5*(jp[x*dim+y]+jp[y*dim+x]);
      else Die ("You have to marginalize X(=0) or Y (=1), or averaging (2) nothing else allowed\n");

  CheckSingleProb(sp, dim);
}

/* Function: MatrixMethodTransitionsFromRate()
 * Date:     ER, Fri Feb 18 11:14:23 CST 2005 [St. Louis]
 *
 * Purpose:  Given a rate R a Q_zero and a Q_infty calculate Q(t)
 *  
 *  
 * Method:   Q(t) = Q_infty + (q)zero-Q_infty)*exp{tR}
 *
 *
 * Args:              
 *
 * Returns: Q(t), allocated here. Freed by caller
 *           
 */
void
MatrixMethodTransitionsFromRate (double time, int Lx, int Ly, double **ret_Q, 
				 double *R, double *Q_zero, double *Q_infty, int verbose)
{
  double *Q;    /* the conditional matrix */
  double *expR; /* the exp(tR) matrix */
  double *M;    /* M = Q_0-Q_infty */
  double *aux;  /* an auxiliary matrix */
  double  tmax;
  double  xtime;
  int     ntimes;
  int     i, j;

  tmax = 1.0;

  ntimes = (int)time / (int)tmax;
  xtime  = time - tmax*(double)ntimes;

  if (verbose) 
    printf("\ntime = %f tmax = %f ntimes = %d xtime = %f\n", time, tmax, ntimes, xtime);

  /* M = (Q_0  - Q_infty)
   */
  M = (double *)Cal_M_N_Sum(stdout, Q_zero, Q_infty, Lx, Ly, FALSE, FALSE);
  if (verbose) {
    printf("M = (Q_zero-Q_infty)\n");
    for (i = 0; i < Lx; i++) 
      PrintVectorProbs(stdout, M+i*Ly, Ly);
  }  

  expR = (double *)Cal_Id(Ly);
  if (ntimes > 0) {
    aux = (double *)Cal_M_Exp(stdout, R, Ly, tmax, FALSE);
    while (ntimes > 0) {
      Comp_M_N_Prod(stdout, aux, expR, Ly, FALSE); /*  calculate R^(ntimes*tmax), dump it in  expR*/
      ntimes --;
    }
    free (aux);
  }

  aux = (double *)Cal_M_Exp(stdout, R, Ly, xtime, FALSE);
  Comp_M_N_Prod(stdout, aux, expR, Ly, FALSE); /*  multiply R^xtime *R^(ntimes*tmax) , dump it in expR */

  if (verbose) {
    printf("\nexp(tR)(t=%f)\n", time);
    for (i = 0; i < Ly; i++) 
      PrintVectorProbs(stdout, expR+i*Ly, Ly);
  }

  /* (Q_0-Q_infty)*e^{tR} */
  Q = (double *)Cal_M_N_Prod(stdout, M, expR, Lx, Ly, Ly, FALSE);

  /* Q_infty + (Q_0-Q_infty)*e^{tR} */
  Comp_M_N_Sum(stdout, Q_infty, Q, Lx, Ly, TRUE, FALSE);
  
  if (verbose) {
    printf("\nBEFORE REG\nMatrix Q(t=%f) = Q_infty + (Q_0-Q_infty)*e^{tR}\n", time);
   for (i = 0; i < Lx; i++) {
      for (j = 0; j < Ly; j++) 
	printf(" %f ", Q[i*Ly+j]);
      printf("\n");
    }
  }

  /* Regularize in case there are negative entries
   */
  for (i = 0; i < Lx; i++) 
    if (verbose && QOMRegularizationAlgorithm(Q+i*Ly, Ly, FALSE)) 
      printf("regularization at time = %f\n",  time);

  if (verbose) {
    printf("\nMatrix Q(t=%f) = Q_infty + (Q_0-Q_infty)*e^{tR}\n", time);
    for (i = 0; i < Lx; i++) 
      PrintVectorProbs(stdout, Q+i*Ly, Ly);
  }

 /* check it's a matrix of conditional probabilites */
  for (i = 0; i < Lx; i++) 
    CheckSingleProb(Q+i*Ly, Ly);

  free(aux);
  free(expR);

  *ret_Q = Q;
}

/* Function: RateFromConditionals()
 * Date:     ER, Fri Sep 12 10:09:25 CDT 2003 [St. Louis]
 *
 * Purpose:  Given a set of conditionals probabilities Q(a|b, t) 
 *           calculate the rate matrix
 *  
 *  
 * Method:   R = log[Q_0^{-1} * Q}
 *
 *
 * Args:              
 *
 * Returns: Qr(LxL). 
 *          Qr is allocated here, freed by caller.
 *           
 */
void
RateFromConditionals(FILE *ofp, double *R, double *Q, double *Q_0, int L, double tfactor, int hasindel, int pedantic, int verbose)
{
  double  *AUX;
  double  *Q_0_inv;
  double  *K;
  double  *S;
  int      dim;
  int      i;

  if (hasindel) dim = L-1;
  else          dim = L;

  /* consistency check, rows should add up to one
   */	
  fprintf(ofp, "Q(j|i, %.2f t_*)  probabilities\n", tfactor);
  PrintProbs(ofp, Q, L);
  
  for (i = 0; i < dim; i++)
    CheckSingleProb(Q+i*L, L);
  
  /* calculate Q_0_inv
   */
  Q_0_inv = Cal_M_Inv(ofp, Q_0, L, verbose);
  if (verbose) {
    fprintf(ofp, "Q_0_inv(j|i)  probabilities\n");
    PrintProbs(ofp, Q_0_inv, L);
  }
  
  /* construct K = log [Q_0_inv * Q]
   */
  AUX = Cal_M_N_Prod(ofp, Q_0_inv, Q, L, L, L, verbose); /*  multiply Q_0_inv*Q, dump it in AUX */
  
  S = Cal_Id(L);
  Comp_M_N_Sum(ofp, AUX, S, L, L, FALSE, verbose);   /* S = AUX - Id  */
 
  if (verbose) 
    IsQConsistent(ofp, S, L, pedantic, verbose);  /* check that Q is consistent with a markovian stationary 
						   * model of evolution. 
						   * In plain words it checks whether logQ is going to exist.
						   *
						   * Basically, it checks whether  all eigenvalues of Q - I
						   * are between -1 and 1. 
						   *
						   * remember log (1+x) = x - x^2/2 + ...   converges for -1 < x <= 1
						   *
						   */
  K = Cal_M_Taylor_Log(ofp, AUX, L, verbose);
  MultiplyMatrix ((1.0/tfactor), L, K);

  fprintf(ofp, "R  rate matrix\n");
  PrintProbs(ofp, K, L);
  
 
  if (pedantic)
    IslogQConsistent(ofp, Q_0, K, L, L, L, verbose); /*  check that K is consistent with a 
						      *  markovian stationary model of evolution. 
						      *
						      *  Basically, it checks whether all entries of  
						      *  Q_0 + epsilon*Q_0*K are none negatives.
						      */
  
  CopyMatrix(R, K, L, L);

  free(K);
  free(S);
  free(AUX);
  free(Q_0_inv);
}

/* Function: SaturationProbs()
 * Date:     ER, Thu Sep 11 16:39:06 CDT 2003 [St. Louis]
 *
 * Purpose:  Given a set of joint P(a,b| t_o) 
 *           at a given evolutionary time t_o, 
 *           and set of conditional probabilities at time zer0 P^o
 *         
 *           calculate P(a,b| tfactor*t_o) using a Markov chain model for evolution
 *  
 * Method:   (1) calculate conditionals: P(a,b | t_o) --> Q(a | b,t_o) 
 *            
 *           (3) evolve conditionals at time infity Q(a | b, r*t_o) = Q_0 * exp{ r * log[Q_0^{-1} * Q] }
 *
 *           (4) extract the saturation probs.
 *
 * Args:              
 *
 * Returns:  saturation probabilities
 *           
 */
double *
SaturationProbs(double *R, double *Q_0, int L, int pedantic, int verbose)
{
  double *psat;
  double *Q;
  double  tfactor = 10000.0;
  double  margin = 0.0005;
  int     i, j;
  
  psat = (double *) MallocOrDie (sizeof(double) * L);

  /* calculate Q^\infty = Q_0 * exp{tfactor*K}
   */
  Condi_From_Rate(&Q, R, Q_0, tfactor, L, pedantic, verbose);

  if (verbose) {
    fprintf(stdout, "R \n");
    PrintProbs(stdout, R, L);
    fprintf(stdout, "Q_0 \n");
    PrintProbs(stdout, Q_0, L);
    fprintf(stdout, "Q(%f) probs\n", tfactor);
    PrintProbs(stdout, Q, L);
  }
  
  for (j = 0; j < L; j++) {
    psat[j] = Q[j*L+j];
    for (i = 1; i < L; i++)
      if (psat[j]-Q[i*L+j]>margin) Die ("SaturationProbs(): you have not reached saturation yet\n");
  }
  
  /* normalize */
  DNorm(psat, L);
 
  free (Q);

  return psat;
}

void
ComputeSaturationProbs(double *R, double *Q_0, double *pm, int L, int pedantic, int verbose) 
{
  double *psat;
  
  psat = SaturationProbs(R, Q_0, L, pedantic, verbose); 

  CopyVector(pm, psat, L);

  free (psat);
}

/* Function: TimeIdCorrelation()
 *
 * Date:     ER, Thu Nov 14 14:16:12 CST 2002 [St. Louis]
 *
 *
 *     id = div.id_zero - div.a * (time-div.t_zero) ^ div.b ==>  log (time-div.t_zero) = 1/div.b log[(div.id_zero-id)/div.a]
 *
 *     time = div.t_zero + exp{ 1/div.b log[(div.id_zero-id)/div.a] }
 * 
 *       using 3 different fitting curves
 *
 */
double 
TimeIdCorrelation(struct three_divergence_s div, double id)
{
  double time;

  if (id == 100.) { time = 0.0; return time; }

  /* check that the curve parameters make sense
   */
  if (div.prxmal.a == 0.0 || div.medial.a == 0 || div.distal.a == 0) 
    Die ("TimeIdCorrelation() error. No time dependence. a parameter is zero\n");
  if (div.prxmal.b == 0.0 || div.medial.b == 0 || div.distal.b == 0) 
    Die ("TimeIdCorrelation() error. No time dependence. b parameter is zero\n");
  if (div.prxmal.id_cutoff >= div.prxmal.id_zero || div.medial.id_cutoff >= div.medial.id_zero || div.distal.id_cutoff >= div.distal.id_zero) 
    Die ("TimeIdCorrelation() error. No ID range\n");
  
  if      (id >= div.prxmal.id_cutoff)
    time = div.prxmal.t_zero + EXP2( LOG2( (div.prxmal.id_zero-id) / div.prxmal.a) / div.prxmal.b );
  else if (id >= div.medial.id_cutoff)
    time = div.medial.t_zero + EXP2( LOG2( (div.medial.id_zero-id) / div.medial.a) / div.medial.b );
  else if (id >= div.distal.id_cutoff)
    time = div.distal.t_zero + EXP2( LOG2( (div.distal.id_zero-id) / div.distal.a) / div.distal.b );
  else
    Die ("TimeIdCorrelation() error. ID (%d) is our of range\n", id);
    
  if (id < 100.0 && time < TMIN) time = TMIN;

  return time;
}

/* Function: TimeIdCorrelation3Models()
 *
 * Date:     ER, Thu Apr 29 13:35:59 CDT 2004 [St. Louis/ home with Coro and Rocio]
 *
 *     time/id correlations for the three models using 3 different fitting curves
 *
 */
struct three_times_s 
TimeIdCorrelation3Models(struct three_divergence_s div_o, struct three_divergence_s div_c, struct three_divergence_s div_r, double id)
{

  struct three_times_s time;
  
  time.oth = TimeIdCorrelation(div_o, id);
  time.cod = TimeIdCorrelation(div_c, id);
  time.rna = TimeIdCorrelation(div_r, id);

  return time;

}

/* Function: TransitionsEvolvedLinear()
 * Date:     ER, Fri Jul 26 15:20:09 CDT 2002 [St. Louis]
 *
 * Purpose:  Given a set of conditionals probabilities Q(a|b, t) 
 *           at a given evolutionary time t, 
 *           and set of conditional probabilities at time zero Q_0
 *         
 *           calculate Q(a,b| r*t) using a Markov chain model for evolution
 *  
 * Method:   
 *
 *                    q^i(t) = q^i_0 + r^i [ exp{tA^i} - I ]  ---   A(nxn), q(1xn), q_0(1xn), r(1xn)
 *
 *
 *              | q^1 |                           | q^1_o |                     |  0  |                           | 1 |
 *              |  .  |                           |   .   |                     |  0  |                           | 1 |
 *              |  .  |                           |   .   |                     |  0  |                           | 1 |
 *              |  .  |                           |   .   |                     |  0  |                           | 1 |
 *        Q^* = | q^i |                     Q_0 = | q^i_o |               R^i = |  0  |          A^i = 1/t* log [ | 1 | (q^i - q^I_o) + I ]
 *              |  .  |                           |   .   |                     | r^i |                           | 1 |
 *              |  .  |                           |   .   |                     |  0  |                           | 1 |
 *              |  .  |                           |   .   |                     |  0  |                           | 1 |
 *              | q^n |                           | q^n_o |                     |  0  |                           | 1 |
 *
 *
 *
 *     THEN,      Q(t)  =  Q_0  + Sum_i R^i [ exp(tA^i) - I ]
 *
 * Args:              
 *
 * Returns: Qr(LxL). 
 *          Qr is allocated here, freed by caller.
 *           
 */
void
TransitionsEvolvedLinear(FILE *ofp, double *q, double *q_zero, double *r, double *r_diag, double *r_fix, 
			 int L, double tfactor, int pedantic, int verbose)
{
  double  *Id;        /* LxL array */
  double  *A;         /* LxL array */
  double  *K;         /* LxL array */
  double  *Qr;        /* 1xL array */
  int      i, j;

  Id = Cal_Id(L);

  /* Allocate 
   */
  A = (double *) MallocOrDie (sizeof(double) * L * L);

  /* 
   *
   *                       | 1 |           | K1                |
   *                       | 1 |           |    K2             |
   *                       | 1 |           |       K3          |
   *    Calculate  A = I + | 1 | r_fix +   |          .        | 
   *                       | 1 |           |            .      | 
   *                       | 1 |           |              .    | 
   *                       | 1 |           |                Kn | 
   *
   *
   *          with Ki = r_diag^i =  ( q^i - q^i_zero - beta * r^i_fix) / r^i
   *
   */

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) {
      A[i*L+j] = r_fix[j];  

      if (i==j) A[i*L+j] += r_diag[i];
    }
  
  if (pedantic)
    IsQConsistent(ofp, A, L, pedantic, verbose);  /* check that A is consistent with a markovian stationary 
						   * model of evolution. 
						   * In plain words it checks whether log(I+A) is going to exist.
						   *
						   * Basically, it checks whether  all eigenvalues of A
						   * are between -1 and 1. 
						   *
						   * remember log (1+x) = x - x^2/2 + ...   converges for -1 < x <= 1
						   *
						   */

  /* Now add the identity to complete the 
   * definition of A
   */
  for (i = 0; i < L; i++) A[i*L+i] += 1.0;

  if (FALSE) {
    fprintf(ofp, "\n EXP[ t^* A ] \n");

    for (i = 0; i < L; i++) {
      for (j = 0; j < L; j++) 
        fprintf(ofp, "%.4f ", A[i*L+j]);
      fprintf(ofp, "\n");
    }
    fprintf(ofp, "\n\n");
  }

  /* construct K = log [A]
   */
  K = Cal_M_Taylor_Log(ofp, A, L, verbose);
  
  if (pedantic)
    IslogQConsistent_2(ofp, q_zero, r, K, 1, L, L, verbose); /*  check that K is consistent with a 
							      *  markovian stationary model of evolution. 
							      *
							      *  Basically, it checks whether all entries of  
							      *  q_zero + epsilon*q_zero*K are none negatives.
							      */
  
  /* calculate Qr = Q_0 + R * [ exp{tfactor*K} - Id]
   */    
  Comp_M_Exp(ofp, K, L, tfactor, verbose);                /* exp:      exp(t*K),  dump it in K      */
  Comp_M_N_Sum(ofp, K, Id, L, L, FALSE, verbose);             /* add:      K-Id,      dump it in K      */
  CopyMatrix(K, Id, L, L);
  Qr = Cal_M_N_Prod(ofp, r, K, 1, L, L, verbose);             /* multiply: r*K,       dump it in Qr     */
  Comp_M_N_Sum(ofp, q_zero, Qr, 1, L, TRUE, verbose);         /* add:      q_zero+Qr, dump it in Qr     */
  
  
  if (verbose) {
    fprintf(ofp, "\nT(i->j | t)  Transition probabilities at time %.3f\n", tfactor);
    for (i = 0; i < L; i++) 
        fprintf(ofp, "%.4f ", Qr[i]);
    fprintf(ofp, "\n----------------------------------------------------------\n");
  }
  
  CopyMatrix(q, Qr, 1, L);
  /* consistency check, rows should add up to one
   */	
  CheckSingleProb(q, L);
  
  free(K);
  free(A);
  free(Id);
  free(Qr);

}

/* Function: TransitionsEvolved()
 * Date:     ER, Tue Aug 13 16:17:27 CDT 2002 [St. Louis]
 *
 * Purpose:  Given a set of conditionals probabilities Q(a|b, t) 
 *           at a given evolutionary time t, 
 *           and set of conditional probabilities at time zero Q_0
 *         
 *           calculate Q(a,b| r*t) using a Markov chain model for evolution
 *  
 * Method:   
 *
 *                    q^i(t) = q^i_0 + r^i [ exp{tA^i} - I ]  ---   A(nxn), q(1xn), q_0(1xn), r(1xn)
 *
 *
 *              | q^1 |                           | q^1_o |                     |  0  |                           | 1 |
 *              |  .  |                           |   .   |                     |  0  |                           | 1 |
 *              |  .  |                           |   .   |                     |  0  |                           | 1 |
 *              |  .  |                           |   .   |                     |  0  |                           | 1 |
 *        Q^* = | q^i |                     Q_0 = | q^i_o |               R^i = |  0  |          A^i = 1/t* log [ | 1 | r^fix + Diag(e^(-K^1)) ]
 *              |  .  |                           |   .   |                     | r^i |                           | 1 |
 *              |  .  |                           |   .   |                     |  0  |                           | 1 |
 *              |  .  |                           |   .   |                     |  0  |                           | 1 |
 *              | q^n |                           | q^n_o |                     |  0  |                           | 1 |
 *
 *
 *
 *     THEN,      Q(t)  =  Q_0  + Sum_i R^i [ exp(tA^i) - I ]
 *
 * Args:              
 *
 * Returns: void
 *           
 */
void
TransitionsEvolved(FILE *ofp, double *q, double *q_zero, double *q_infty, double *r_diag, 
		   int L, double tfactor, int pedantic, int verbose)
{
  double  *Id;        /* LxL array */
  double  *A;         /* LxL array */
  double  *K;         /* LxL array */
  double  *Q;         /* 1xL array */
  double  *r;         /* 1xL array */
  double   w = 0.0;
  int      isdiagonal = 1;
  int      i, j;

  Id = Cal_Id(L);

  if (verbose) {
    fprintf(ofp, "\nT(i->j | t)  Transition probabilities at time zero\n");
    for (i = 0; i < L; i++) 
      fprintf(ofp, "%.4f ", q_zero[i]);
    fprintf(ofp, "\nT(i->j | t)  Transition probabilities at time star\n");
    for (i = 0; i < L; i++) 
      fprintf(ofp, "%.4f ", q[i]);
    fprintf(ofp, "\nT(i->j | t)  Transition probabilities at time infty\n");
    for (i = 0; i < L; i++) 
      fprintf(ofp, "%.4f ", q_infty[i]);
    printf("\n");
  }
  
  /* Allocate 
   */
  A = (double *) MallocOrDie (sizeof(double) * L * L);

  /* 
   *
   *                      | e^(-K^1)                            |
   *                      |         e^(-K2)                     |
   *                      |                e^(-K3)              |
   *    Calculate  A =    |                       .             | 
   *                      |                          .          | 
   *                      |                             .       | 
   *                      |                             e^(-Kn) | 
   *
   *
   *          with   e^(-K^i) = 1 + r_diag^i =  1 + ( q^i - q^i_zero) / r^i
   *
   */

  for (i = 0; i < L; i++) A[i*L+i] = 1.0 + r_diag[i];

  if (verbose) {
    fprintf(ofp, "\n EXP[ t^* A ] \n");

    for (i = 0; i < L; i++) {
      for (j = 0; j < L; j++) 
        fprintf(ofp, "%.4f ", A[i*L+j]);
      fprintf(ofp, "\n");
    }
    fprintf(ofp, "\n\n");
  }


  /* Calculate exp (t A), 
   *
   *   If A is diagonal it is very simple:  exp(tA) =  exp[ t^* A ] * t
   *
   *   Otherwise have to calculate the log and then exponentiate.
   *
   */
  if (isdiagonal) { 
    K = (double *) MallocOrDie (sizeof(double) * L * L);
    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++) 
	if (i==j) {
	  if (tfactor > 0.0) K[i*L+j] = exp(log(A[i*L+j]) * tfactor); 
	  else               K[i*L+j] = 1.0;
	}
	else
	  K[i*L+j] = 0.0;
  } 
  
  else {
    /* construct K = log [A]
     */
    K = Cal_M_Taylor_Log(ofp, A, L, verbose);
    
    if (pedantic)
      IslogQConsistent_2(ofp, q_zero, r, K, 1, L, L, verbose); /*  check that K is consistent with a 
								*  markovian stationary model of evolution. 
								*
								*  Basically, it checks whether all entries of  
							        *  q_zero + epsilon*q_zero*K are none negatives.
							        */
    
    /* exp{tfactor*K} 
     */    
    Comp_M_Exp(ofp, K, L, tfactor, verbose);                /* exp:      exp(t*K),  dump it in K      */
  }
  

  /* Calculate w(t) = u^T exp(tA) r
   */
  r = Cal_M_N_Sum(ofp, q_zero, q_infty, 1, L, FALSE, verbose);         /* add:      r = q_0 - q_infty,           */  

  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++)  
      w += K[i*L+j] * r[j];
  if (verbose) {
    printf("normalization factor %f\n", w);
  }

  /* Calculate 
   *
   *           q(t) = q_0 + [ exp(tA) - I ] r
   *
   */
  Comp_M_N_Sum(ofp, K, Id, L, L, FALSE, verbose);                   /* add:      K-Id,      dump it in K      */
  CopyMatrix (K, Id, L, L);
  Q = Cal_M_N_Prod(ofp, r, K, 1, L, L, verbose);                    /* multiply: r*K,       dump it in q      */
  Comp_M_N_Sum(ofp, q_zero, Q, 1, L, TRUE, verbose);                /* add:      q_zero+q,  dump it in q      */
  
  
  /* Normalize the probabilities
   */
  for (i = 0; i < L; i++) 
    Q[i] /= (1. + w);

  if (verbose) {
    for (i = 0; i < L; i++) 
      fprintf(ofp, "%.4f ", Q[i]);
    fprintf(ofp, "\n");
  }
  
  CopyMatrix(q, Q, 1, L);

  /* consistency check, rows should add up to one
   */	
  CheckSingleProb(q, L);
  
  free(r);
  free(K);
  free(A);
  free(Id);
  free(Q);

 }


/* Function: TransitionsDirectoryCosines()
 * Date:     ER, Tue Aug 13 17:01:49 CDT 2002 [St. Louis]
 *
 * Purpose:  Fill R_diag ant T_rate of the transfer matrix for the OTH model
 *
 * Args:     T
 *           T_zero
 *
 *
 * Returns:  void. Calculates R_diag and T_rate. Both allocated here, freed by caller.
 */
void
TransitionsDirectoryCosines(double *q, double *q_zero, double *q_rate, double *R_diag, int dim)
{
  double *q_infty;
  double  sum = 0.0;
  double  angle;
  double  shift;
  int     num = 0;
  int     col;

  /* Allocate q^infty */
  q_infty = (double  *) MallocOrDie (sizeof(double) * dim);

  /* fill R_diag.  
   * 
   *  R_diag_i = cos(theta)
   * 
   *  cos (theta) = (q_i - q^0_i) / (q^0_i - q^infty_i) 
   * 
   * 
   *  -1 < cos(theta) <= 0 ===> pi/2 <= theta < pi --- 
   * 
   *  also because 0 <= q_0 - r <= 1  ===>  cos(theta) <= (q^* / q^0) - 1
   * 
   * 
   */
  angle = 3.0*PI/4.0;
  shift = 0.05;

  for (col = 0; col < dim; col++) {
    if (q[col] - q_zero[col] < MARGIN && q[col] - q_zero[col] > -MARGIN) { 
      R_diag[col] = 0.0;           
      sum += q_zero[col]; 
      num ++; 
    }
    else {
      if (cos(angle) <= q[col]/q_zero[col] - 1.0) R_diag[col] = cos(angle);
      else                                        R_diag[col] = q[col]/q_zero[col] - 1.0;
    }
  }
  
  /* fill T_rate.  
   * 
   * q^rate = q^infty - q^0
   * 
   * q^rate = (q_i - q^0_i) / cos (theta)
   * 
   *  If  cos (theta) = 0, q^rate is free with the condition    0 <= q^0 - q^rate <= 1
   * 
   */
  sum /= num;
  for (col = 0; col < dim; col++) {
    if (q[col] - q_zero[col] < MARGIN && q[col] - q_zero[col] > -MARGIN) 
      q_rate[col] = q_zero[col] - sum; 
    else {                                                                
      if (cos(angle) <= q[col]/q_zero[col] - 1.0) q_rate[col] = (q[col] - q_zero[col]) / R_diag[col];
      else                                        q_rate[col] = q_zero[col];
    }

    q_infty[col] = q_zero[col] - q_rate[col];
  }
 
  if (FALSE) {
    printf("RA_diag - OTH transfer matrix\n");
    for (col = 0; col < dim; col++) 
      printf("%f ", R_diag[col]);
    
    printf("\n");
    
    printf("TA_rate - OTH transfer matrix\n");
    for (col = 0; col < dim; col++) 
      printf("%f ", q_rate[col]);
    
    printf("\n");
    
    
    printf("TA_infty - OTH transfer matrix\n");
    for (col = 0; col < dim; col++) 
      printf("%f ", q_infty[col]);
    
    printf("\n");
  }

  /* paranoia - check that q^infty = q^0 + q_rate add up to one
   */
  CheckSingleProb(q_infty, dim);
    
  free(q_infty);
}

/* Function: TransitionsExpVector()
 * Date:     ER, Wed Aug 14 17:48:05 CDT 2002 [St. Louis]
 *
 * Purpose:  Fill R_diag vector and r vector
 *                 
 *               from the expressions:        r_i = q^0_i - q^infty_i
 *
 *                                       R^diag_i = ( q_i - q^0_i ) / r_i
 *                 
 * Args:     
 *
 * Returns:  void
 */
double *
TransitionsExpVector(double *T, double *T_zero, double *T_infty, int dim)
{
  double  *R_exp;
  double   sum   = 0.;
  int      col;
 
  /* check the three input probability vectors are in good order
   *
   * to have things well behaving one of these two conditions have to occur
   *
   *
   *            T^0_i < T^*_i < T^infty_i    or     T^infty_i < T^*_i < T^0_i  
   *
   */
  /*check_three_probs (T, T_zero, T_infty, dim);*/

  /* allocate R_diag
   */
  R_exp = (double *) MallocOrDie (sizeof(double) * dim);

  /* initialize
   */
  for (col = 0; col < dim; col++)
    R_exp[col] = 0.0;
  

  /* check consistency
   */
  for (col = 0; col < dim; col++) {
    if (T_zero[col] >= T[col] && T[col] < T_infty[col]) 
      Die ("TransitionsExpVector(): unconsistent probabilities. dim=%d col=%d P(t=0)=%f P(t=t^*)=%f P(t=infty)=%f\n",
	   dim, col, T_zero[col], T[col], T_infty[col]);
    if (T_zero[col] <= T[col] && T[col] > T_infty[col]) 
      Die ("TransitionsExpVector(): unconsistent probabilities. dim=%d col=%d P(t=0)=%f P(t=t^*)=%f P(t=infty)=%f\n",
	   dim, col, T_zero[col], T[col], T_infty[col]);
  }

  /* Calculate 
   *
   *    
   *
   */
  for (col = 0; col < dim; col++)  
    if (T_zero[col] - T_infty[col] < MARGIN && T_zero[col] - T_infty[col] > -MARGIN) 
      R_exp[col] = 0.0;     
    else
      R_exp[col] = (T[col] - T_zero[col]) / (T_zero[col] - T_infty[col]); 
  
  if (FALSE) {
    printf("R^exp - matrix\n");
    for (col = 0; col < dim; col++) 
      printf("%.4f ", R_exp[col]);
    printf("\n");
  }
  
  /* paranoia: 
   *
   *           Define the scalar product K.R = sum_k K[k] * R[k]
   *        
   *           IS R^diag.R = 0 ????
   *
   *           Vectors R and R^diag being orthogonals  guarantees that q(t) add up to one.
   *
   */
  for (col = 0; col < dim; col++)
     sum += (T_zero[col] - T_infty[col]) * R_exp[col];
  if (sum < -MARGIN && sum > MARGIN) Die ("Bad Diag Matrix, sum_col K[col] R[col] = %f", sum);
    
  return R_exp;
}

/* Function: TransitionsDiagVector()
 * Date:     ER, Fri Aug  9 15:43:08 CDT 2002 [St. Louis]
 *
 * Purpose:  Fill R_diag matrix 
 *                 
 *               from the expression: R^diag_i * r_i = q_i - q^0_i - (sum_l r_l) R^fix_i  
 *                 
 *                 
 *           R^diag_k = (q_k - q^0_k - R^fix_k ) / r_k
 *
 * Args:     R_diag -- transfer matrix for the trnasitions  of the model
 *
 * Returns:  R_diag is allocated and filled here. R_diag freed by caller.
 */
double *
TransitionsDiagVector(double *T, double *T_zero, double *T_rate, double *R_fix, int dim)
{
  double *R_diag;
  double  beta  = 0.;
  double  sum   = 0.;
  int     col;

  /* allocate R_diag
   */
  R_diag = (double *) MallocOrDie (sizeof(double) * dim);

  /* initialize
   */
  for (col = 0; col < dim; col++)
    R_diag[col] = 0.0;
  
  /* Calculate beta = sum_k R[k]
   *
   */
  for (col = 0; col < dim; col++) beta += T_rate[col];
  if (beta < -MARGIN || beta > MARGIN) Die ("Bad Rate Matrix in TransitionsDiagVector(), sum_col R[col] = %f", beta);
  
  /* Calculate R^diag
   *
   *                       R^diag_k = (q_k - q^0_k - R^fix_k ) / r_k
   *
   */
  for (col = 0; col < dim; col++) 
    R_diag[col] = (T[col] - T_zero[col] - R_fix[col]) / T_rate[col]; 
  
  if (FALSE) {
    printf("R^diag - matrix\n");
    for (col = 0; col < dim; col++) 
      printf("%.4f ", R_diag[col]);
    printf("\n");
  }
  
  /* paranoia: 
   *
   *           Define the scalar product K.R = sum_k K[k] * R[k]
   *        
   *           IS R^diag.R = 0 ????
   *
   *           Vectors R and R^diag being orthogonals  guarantees that q(t) add up to one.
   *
   */
  for (col = 0; col < dim; col++)
     sum += T_rate[col] * R_diag[col];
  if (sum < -MARGIN && sum > MARGIN) Die ("Bad Diag Matrix, sum_col K[col] R[col] = %f", sum);
    
  return R_diag;
}


/* Function: TransitionsFixVector()
 * Date:     ER, Fri Aug  9 15:32:29 CDT 2002 [St. Louis]
 *
 * Purpose:  Fill R_fix matrix 
 *                 
 *               from the expression: R^diag_i * r_i = q_i - q^0_i - (sum_l r_l) R^fix_i  
 *                 
 *                 
 *           R^fix_k = (q_k - q^0_k) / beta - R^diag_k * r_k
 *
 *           where  beta = (sum_l r_l) 
 *
 * Args:     R_fix -- transfer matrix for the trnasitions  of the model
 *
 * Returns:  R_fix is allocated and filled here. R_fix freed by caller.
 */
double *
TransitionsFixVector(double *T, double *T_zero, double *T_rate, double *R_diag, int dim)
{
  double *R_fix;
  double  beta  = 0.;
  double  sum   = 0.;
  int     col;
  
  /* allocate R_fix
   */
  R_fix = (double *) MallocOrDie (sizeof(double) * dim);
  
  /* initialize
   */
  for (col = 0; col < dim; col++)
    R_fix[col] = 0.0;
  
  /* Calculate beta = sum_k R[k]
   *
   */
  for (col = 0; col < dim; col++) beta += T_rate[col];
  if (beta < -MARGIN && beta > MARGIN) Die ("Bad Rate Matrix in TransitionsFixVector(), sum_col R[col] = %f", beta);
  
  /* Calculate R^fix
   *
   *        R^fix_k = (q_k - q^0_k) / beta - R^diag_k * r_k
   *  
   */
  for (col = 0; col < dim; col++) 
    R_fix[col] = (T[col] - T_zero[col]) / beta - R_diag[col] * T_rate[col]; 
  
  if (FALSE) {
    printf("R^fix - matrix\n");
    for (col = 0; col < dim; col++) 
      printf("%.4f ", R_fix[col]);
    printf("\n");
  }
  
  /* paranoia: 
   *            elements of R^fix should add up to zero
   *
   */
  for (col = 0; col < dim; col++)
    sum += R_fix[col];
  if (sum < -MARGIN || sum > MARGIN) Die ("Bad Fix Matrix, sum_col R_fix[col] = %f", sum);
  
  return R_fix;
}


/* Function: TransitionsRateVector()
 * Date:     ER, Wed Aug  7 09:52:32 CDT 2002 [St. Louis]
 *
 * Purpose:  Fill R the rate matrix 
 *                 
 *               from the expression: R^diag_i * r_i = q_i - q^0_i - (sum_l r_l) R^fix_i  
 *                 
 *                 
 *           r_k = (q_k - q^0_k) / R^diag_k - beta * R^fix_k / R^diag_k
 *
 *           where  beta = (sum_l r_l) = [ sum_l (q_l - q^0_l) / R^diag_l ] / 1 + sum_l R^fix_k / R^diag_k
 *
 * Args:     R -- transfer matrix for the trnasitions  of the model
 *
 * Returns:  R is allocated and filled here. R freed by caller.
 */
double *
TransitionsRateVector(double *T, double *T_zero, double *R_fix, double *R_diag, int dim)
{
  double *R;
  double  beta  = 0.;
  double  sum_q = 0.;
  double  sum_i = 0.;
  double  sum   = 0.;
  int     col;

  /* allocate R 
   */
  R = (double *) MallocOrDie (sizeof(double) * dim);
  
  /* initialize
   */
  for (col = 0; col < dim; col++)
    R[col] = 0.0;
  
  /* Calculate beta = sum_k R[k]
   *
   */
  for (col = 0; col < dim; col++) 
    if (R_diag[col] > MARGIN || R_diag[col] < -MARGIN) {
      sum_q += (T[col] -  T_zero[col]) / R_diag[col]; 
      sum_i += R_fix[col]              / R_diag[col];
    }
  
  if (sum_q != 0.) {
    if (1.0 + sum_i < MARGIN && 1.0 + sum_i > -MARGIN) 
      Die ("There is an inconsistency in your Rate matrix. Revise your directory cosines. sum_q = %f   1+sum_i = %f", sum_q, 1.0+sum_i);
    else
      beta += sum_q / (1.0 + sum_i);
  }
  
  /* fill R 
   *           R_k = ( q_k - q^0_k - beta R^fix_k ) / R^diag_k
   *  
   */
  for (col = 0; col < dim; col++) 
    if (R_diag[col] > MARGIN || R_diag[col] < -MARGIN) 
      R[col] += ( T[col] -  T_zero[col] - beta * R_fix[col] ) / R_diag[col]; 
  
  if (FALSE) {
    printf("R - rate matrix -  beta = %f\n", beta);
    for (col = 0; col < dim; col++) 
      printf("%.4f ", R[col]);
    printf("\n");
  }
  
  /* paranoia: 
   *
   *           Define the scalar product K.R = sum_k K[k] * R[k]
   *        
   *           IS R^diag.R = 0 ????
   *
   *           Vectors R and R^diag being orthogonals  guarantees that q(t) add up to one.
   *
   */
  for (col = 0; col < dim; col++)
    sum += R[col] * R_diag[col];
  if (sum < -MARGIN || sum > MARGIN) Die ("Bad Rate Matrix, sum_col K[col] R[col] = %f", sum);
  
  return R;
}

/* Function: check_three_probs()
 * Date:     ER, Wed Aug 14 18:12:30 CDT 2002 [St. Louis]
 *
 * Purpose:  check the three input probability vectors are in good order
 *
 *             to have things well behaving one of these two conditions have to occur
 *
 *
 *              T^0_i < T^*_i < T^infty_i    or     T^infty_i < T^*_i < T^0_i  
 *
 */
void
check_three_probs (double *T, double *T_zero, double *T_infty, int dim)
{
  int i;

  /* check for one of the two conditions
   *
   *   q_0 < q_* < q_infty    or  q_0 > q_* > q_infty
   */
  for (i = 0; i < dim; i++) 
    if      (T_infty[i] <= T[i] && T[i] <= T_zero[i]) ;
    else if (T_infty[i] >= T[i] && T[i] >= T_zero[i]) ;
    else Die ("check your parameters for state %d. Time zero:(%f), star:(%f) and infty:(%f) \n", i, T_zero[i], T[i], T_infty[i]);
  
  
  /*
   * Also check that if two of the probs are equal, the third it is also equal.
   */
  for (i = 0; i < dim; i++) {
    if ( T_zero[i] - T[i]  < MARGIN*MARGIN && T_zero[i]  - T[i] > -MARGIN*MARGIN  &&
	 (T_infty[i] - T[i] > MARGIN*MARGIN || T_infty[i] - T[i] < -MARGIN*MARGIN)  )
      Die ("check your parameters for state %d. Time zero:(%f), star:(%f) and infty:(%f) \n", i, T_zero[i], T[i], T_infty[i]);
  }
}


/* Function: IsQConsistent()
 * Date:     ER, Thu Apr 20 14:19:42 CDT 2000 [St. Louis]
 *
 * Purpose:  Q(t) ---a matrix of conditional probabilities---
 *           is consistent with being derived from a markovian stationary
 *           model of evolution from:
 *                    Q_0 - a matrix of conditional probabilities at time zero
 *                    K - a matrix of rate of evolution
 *
 *           as
 *                     Q(t) = Q_0 * exp {t K}
 *
 *
 *           IF: Q_0 + epsilon * Q_0*K  can play the role of an "instant matrix of evolution"
 *
 *               that is, if all the entries are positive.
 *
 * Args:      Q_0 - LxL matrix of conditional prob matrix at time zero (prealloc)
 *           logQ - LxL log matrix of conditional prob matrix (prealloc)
 *
 * Returns:  1 if consistent, 0 otherwise.
 */
int
IsQConsistent(FILE *ofp, double *Q, int L, int pedantic, int verbose)
{
  struct eigenvalue_s *eigen;
  double              *H;
  double               norm;
  double               isconsistent = TRUE;
  int                  i;
  
  H = HessenbergForm(ofp, Q, L, verbose);
  Hessenberg2Eigenvalues(ofp, H, L, &eigen, pedantic, verbose);

  if (verbose) {
      fprintf(ofp, "Q matrix\n"); 
      PrintProbs(ofp, Q, L);

      fprintf(ofp, "Hessenberg form\n"); 
      PrintProbs(ofp, H, L);
      
      fprintf(ofp, "Eigenvalues for Q\n"); 
      for (i = 0; i < L; i++) 
	fprintf(ofp, "%.4f %.4f \n", eigen->real[i], eigen->imag[i]);
  }
  
  /* The condition is that ALL the eigenvalues of Q are of norm smaller
   * than one. 
   */
  for (i = 0; i < L; i++) {
    norm = sqrt(eigen->real[i]*eigen->real[i]+eigen->imag[i]*eigen->imag[i]);
    if (norm > 1.0+MARGIN2) {
	Warn ("Q not consistent with a markovian stationary model of evolution.\n Norm of eigenvalue (%f, %f) is = %f", 
	      eigen->real[i], eigen->imag[i], norm);
	isconsistent = FALSE;
    }
}

  free(H);
  free(eigen->real);
  free(eigen->imag);
  free(eigen);

  return isconsistent;
}
  
/* Function: IslogQConsistent()
 * Date:     ER, Fri Apr 14 12:06:24 CDT 2000 [St. Louis]
 *
 * Purpose:  Q(t) ---a matrix of conditional probabilities---
 *           is consistent with being derived from a markovian stationary
 *           model of evolution from:
 *                    Q_0 - a matrix of conditional probabilities at time zero
 *                    K - a matrix of rate of evolution
 *
 *           as
 *                     Q(t) = Q_0 * exp {t K}
 *
 *
 *           IF: Q_0 + epsilon * Q_0*K  can play the role of an "instant matrix of evolution"
 *
 *               that is, if all the entries are positive.
 *
 * Args:      Q_0 - LxL matrix of conditional prob matrix at time zero (prealloc)
 *           logQ - LxL log matrix of conditional prob matrix (prealloc)
 *
 * Returns:  void.
 */
void
IslogQConsistent(FILE *ofp, double *Q_0, double *K, int Lx, int Ly, int Lz, int verbose)
{
  double *R;
  int     i, j;

  R = Cal_M_N_Prod(ofp, Q_0, K, Lx, Ly, Lz, FALSE); /*  multiply Q_0*K, dump it in R */

  for (i = 0; i < Lx; i++)
    for (j = 0; j < Lz; j++) 
      if ( (R[i*Lz+j] < 0.0 && Q_0[i*Lz+j] == 0.0)  )
	Warn ("Ah! Q_0log[Q_0^{-1}*Q] not consistent with a markovian stationary model of evolution[L=%d] (%d,%d: %f %f).", 
	      Lz, i, j, Q_0[i*Lz+j], R[i*Lz+j]);     
  
  if (verbose) { 
    fprintf(ofp, "RA Instant evolution probabilities\n");
    for (i = 0; i < Lx; i++) {
      for (j = 0; j < Lz; j++) {
	fprintf(ofp, "%.4f ", R[i*Lz+j]);
      }
      fprintf(ofp, "\n");
    }

    fprintf(ofp, "Q_0\n");
    for (i = 0; i < Lx; i++) {
      for (j = 0; j < Lz; j++) {
	fprintf(ofp, "%.4f ", Q_0[i*Lz+j]);
      }
      fprintf(ofp, "\n");
    }

  }
  
  free(R);
}

/* Function: IslogQConsistent_2()
 * Date:     ER, Mon May  8 12:55:29 CDT 2000 [St. Louis]
 *
 * Purpose:  Q(t) ---a matrix of conditional probabilities---
 *           is consistent with being derived from a markovian stationary
 *           model of evolution from:
 *                    Q_0 - a matrix of conditional probabilities at time zero
 *                    A - a matrix of rate of evolution
 *                    K - a matrix of rate of evolution
 *
 *           as
 *                     Q(t) = Q_0 + A[ exp{t K} - Id ]
 *
 *
 *           IF: Q_0 + epsilon * A*K  can play the role of an "instant matrix of evolution"
 *
 *               that is, if all the entries are positive.
 *
 * Args:     Q_0  - LxL matrix of conditional prob matrix at time zero (prealloc)
 *           A    - LxL matrix  (prealloc)
 *           logQ - LxL log matrix of conditional prob matrix (prealloc)
 *
 * Returns:  void.
 */
void
IslogQConsistent_2(FILE *ofp, double *Q_0, double *A, double *K, int Lx, int Ly, int Lz, int verbose)
{
  double *R;
  int     i, j;

  R = Cal_M_N_Prod(ofp, A, K, Lx, Ly, Lz, FALSE); /*  multiply A*K, dump it in R */

  for (i = 0; i < Lx; i++)
    for (j = 0; j < Lz; j++)
      if ( (R[i*Lz+j] < 0.0 && Q_0[i*Lz+j] < 0.0)  )
	Warn ("Ah! Q_0log[Q_0^{-1}*Q] not consistent with a markovian stationary model of evolution (%d,%d: %f %f).", i, j, Q_0[i*Lz+j], R[i*Lz+j]);     
  
  if (verbose) { 
    fprintf(ofp, "RA Instant evolution probabilities\n");
    for (i = 0; i < Lx; i++) {
      for (j = 0; j < Lz; j++) {
	fprintf(ofp, "%.4f ", R[i*Lz+j]);
      }
      fprintf(ofp, "\n");
    }
  
    fprintf(ofp, "Q_0\n");
    for (i = 0; i < Lx; i++) {
      for (j = 0; j < Lz; j++) {
	fprintf(ofp, "%.4f ", Q_0[i*Lz+j]);
      }
      fprintf(ofp, "\n");
    }

  }
  
  free(R);
}


void
adjust_prob(double *psingle, int size)
{
  double sum     = 0.0;
  double min_neg = 0.0;
  int    x;
  int    flag     = FALSE;
  int    pedantic = FALSE;

  for (x = 0; x < size; x++) {
    if (psingle[x] > 1.0+MARGIN) Die ("adjust_prob(): probabilities are getting too large here. P[%d] = %f", x, psingle[x]);
    if (psingle[x] < -0.0) {
      flag = TRUE;

      if (psingle[x] < min_neg) min_neg = psingle[x];

      if (psingle[x] < -MARGIN) Warn ("adjust_prob(): probabilities are getting too small here, P[%d] = %.20f", x, psingle[x]);
    }

    if (psingle[x] < 0.0) 
      psingle[x] = 0.0;
    
   sum += psingle[x];
  }

  if (flag && pedantic)  Warn ("adjust_prob(): there are negative probabilities. size = %d, abs_min = %.40f\n", size, min_neg);

  for (x = 0; x < size; x++) psingle[x] /= sum;
 
}

void 
OrderVector(int **ret_perm, double *b, int L, int isdescending, int verbose)
{
  int *perm;
  int  x;
  int  i, j;

  perm = (int *) MallocOrDie (sizeof(int) * L); 
  for (i = 0; i < L; i++) perm[i] = i;

  for (i = 0; i < L-1; i++) {
    for (j = i+1; j < L; j++) 
      if (isdescending) {
	if (b[perm[j]] > b[perm[i]]) { 
	  x = perm[i];
	  perm[i] = perm[j]; 
	  perm[j] = x;
	}
      }
      else {
	if (b[perm[j]] < b[perm[i]]) { 
	  x = perm[i];
	  perm[i] = perm[j]; 
	  perm[j] = x;
	}
      }
  }

  for (i = 0; i < L-1; i++) {
    if (isdescending  && b[perm[i]] < b[perm[i+1]]) Die("bad permutatiom");
    if (!isdescending && b[perm[i]] > b[perm[i+1]]) Die("bad permutatiom");
  }

  if (verbose) {
    printf("\noriginal vector\n");
    PrintVectorProbs(stdout, b, L);
    printf("permutation\n");
    for (i = 0; i < L; i++) printf("%d  ", perm[i]);
    printf("\n");
    printf("permutated vector\n");
    for (i = 0; i < L; i++) printf("%f ", b[perm[i]]);
    printf("\n");
  }

  *ret_perm = perm;
}

/* Function: QOMRegularizationAlgorithm()
 *
 * Date:     ER, Tue Nov 30 15:33:36 CST 2004 [St. Louis]
 *
 * Purpose:  regularize a conditonal matrix
 *           taken from Kreinin and Sidelnikova, Algo Research Quarterly 4 (2001) 23-40.
 *
 *
 * (1) Construct the vector w, such that w(i) = r(i)-\lambda where  \lambda = \frac{1}{n}(\sum_i r(i) - 1\).
 *
 * (2) If all w(i) are non negative, r <---- w, is the new regularized row.
 *
 * (3) Otherwise, calculate the permutation w_p= P(w) such that w_p(i)\geq w_p(i+1).
 *
 * (4) Construct C(k) = \sum{i=1}^{k} w_p(i) - k*w_p(k),, for k=1,\ldots,n.
 *
 * (5) Calculate k_{max} = \max [ k; k\geq 1 such that  C(k) \leq 1 ]
 *
 * (6) Construct vector
 *
 *                |  w_p(i) + 1/k_max [1 - \sum_{j=1}^{k_max} b(j) ]    if 1\leq i\leq k_{max}
 *                |
 *  \hat{r}(i) =  |
 *                |
 *                |  0                                                  otherwise
 *
 * (7) The regularized row is given by r <--- P^{-1}(\hat{r}).
 * 
 * Args:      b - a row of a conditional matrix
 *
 * Returns:  void.
 */
int 
QOMRegularizationAlgorithm(double *b, int size, int verbose)
{
  double *c;
  double  corr;
  double  min_neg = 0.0;
  double  lambda;
  int    *perm;
  int     hasneg = FALSE;
  int     pedantic;
  int     i, j;
  int     k_max;
  int     isdescending = TRUE;

  pedantic = FALSE;

  /* Check for negative probabilities
   */
  for (i = 0; i < size; i++) {

    if (b[i] < 0.0) {
      hasneg = TRUE; 

      if (b[i] < min_neg) min_neg = b[i];
      /*if (b[i] < -MARGIN2) 
	Die ("QOMRegularizationAlgorithm(): probabilities are getting too small here, P[%d] = %.20f", i, b[i]);*/
    }
  }

  if (hasneg && pedantic)  
    Warn("QOMRegularizationAlgorithm(): there are negative probabilities. size = %d, abs_min = %.40f\n", size, min_neg);

  if (verbose && hasneg) {
    printf("\nOriginal conditional vector\n");
    PrintVectorProbs(stdout, b, size);
  }

  /* Construct b(i) = b(i) - lambda ,, for lambda = 1/n { sum_i b(i) - 1}
   */
  lambda = 0.0;
  for (i = 0; i < size; i++) lambda += b[i];
  lambda /= size;
  lambda -= 1.0;
  if (verbose) printf("lambda = %f\n", lambda);

  for (i = 0; i < size; i++) b[i] -= lambda;		       
  if (verbose) {
    printf("\nVector -lambda\n");
    PrintVectorProbs(stdout, b, size);
  }

  /* calculate perm[i]
   */
  c = (double *) MallocOrDie (sizeof(double) * size); 
  for (i = 0; i < size; i++) c[i] = 0.0;
  
  OrderVector(&perm, b, size, isdescending, verbose);

  /* c[i] = sum_{0}^{i} a_i - (i+1) a_i */
  for (i = 0; i < size; i++) 
    {
      c[i] = - (i+1) * b[perm[i]];

      for (j = 0; j <= i; j++) c[i] += b[perm[j]];
    }

  /* calculate k_max = \max [ k; k\geq 1 such that  C(k) \leq 1 ]
   */
  k_max = 0;
  for (i = size-1; i >= 0; i--) {
    if (c[i] <= 1.0) { k_max = i; break; }
  }

  /* regularize:
   *
   *  if (1<=i<=k_max) b(i) <--- b(i) + corr,,   where: corr = 1/k_max [1 - \sum_{j=1}^{k_max} b(j) ]
   *  else             b(i) <--- 0,,  
   *
   */
  corr = 1.0;
  for (i = 0; i <= k_max ; i++) corr -= b[perm[i]];
  corr /= (k_max+1);
  
  if (verbose) printf("k_max %d size %d corr %f\n", k_max, size, corr);

  for (i = 0; i < size; i++) 
    if (i <= k_max) b[perm[i]] += corr;
    else            b[perm[i]]  = 0.0;

  if (verbose && hasneg) {
    printf("\nQOM-regularized vector\n");
    PrintVectorProbs(stdout, b, size);
  }

  free(c);
  free(perm);

  return hasneg;
}

/* Function: QOGRegularizationAlgorithm()
 *
 * Date:     ER, Tue Nov 30 15:52:12 CST 2004 [St. Louis]
 *
 * Purpose:  regularize a rate matrix.
 *           taken from Kreinin and Sidelnikova, Algo Research Quarterly 4 (2001) 23-40.
 *
 *            this algorithm had a typo in the paper in step (3)
 *
 * (1) Permutate row vector so that r(1) = R(whichrow, whichrow)
 *
 * (2) Construct the vector w, such that w(i) = r(i)-\lambd$ where  \lambda = \frac{1}{n}\sum_i r(i).
 *
 * (3) Calculate the permutation w_p= P(w) such that w_p(i) < w_p(i+1). # wrong in the paper
 *
 * (4) Construct C(k) = w_p(k) + \sum_{i=0}^{n-k-1}w_p(n-i) - (n-k+1)w_p(k+1) for k=2,\ldots,n-1.
 *
 * (5) Calculate k_{min} = \min [ k;  2\leq k\leq n-1 such that C(k) <= 0 ]
 * 
 * (6) Construct 
 *
 *                |  0                                                    if 2\leq i\leq k_{min}
 *                |
 *  \hat{r}(i) =  |
 *                |
 *                |  w_p(i) - 1/(n-k_{min}+1) \sum_{k_{min}} w_p(j)        otherwise
 *
 * (7) The regularized row is given by r\leftarrow P^{-1}(\hat{r}).
 * 
 * Args:      b - a row of a rate matrix
 *
 * Returns:  void.
 */
void
QOGRegularizationAlgorithm(double *b, int size, int whichrow, int verbose)
{
  double *c;
  double *r_1;
  double  corr;
  double  lambda;
  double  sum = 0.0;
  int    *perm;
  int     hasneg = FALSE;
  int     pedantic;
  int     i, j;
  int     k_min;
  int     isdescending = FALSE;

  pedantic = FALSE;

  /* Check for off-diagonal negative entries
   */
  for (i = 0; i < size; i++) 
    if (i != whichrow && b[i] < 0.0) hasneg = TRUE; 

  if (hasneg && pedantic)  
    Warn("QOGRegularizationAlgorithm(): there are negative off--diagonals for row (%d). size = %d\n", whichrow, size);

  if (verbose && hasneg) {
    printf("\nOriginal Rate vector\n");
    PrintVectorProbs(stdout, b, size);
  }
  if (!hasneg) return;

  /* Reorder row so that r(1) = R(whichrow, whichrow)
   */
  r_1 = (double *) MallocOrDie (sizeof(double) * size); 
  for (i = 0; i < size; i++) r_1[i] = b[i];
    
  if (whichrow > 0) {
    b[0] = r_1[whichrow];

    for (i = 1; i < size; i++) if (i <= whichrow) b[i] = r_1[i-1];
  }
  if (verbose) {
    printf("\nOriginal reordered Rate vector\n");
    PrintVectorProbs(stdout, b, size);
  }

  /* Construct b(i) = b(i) - lambda ,, for lambda = 1/n sum_i b(i)
   */
  lambda = 0.0;
  for (i = 0; i < size; i++) lambda += b[i];
  lambda /= size;
  if (verbose) printf("lambda = %f\n", lambda);

  for (i = 0; i < size; i++) b[i] -= lambda;		       
  if (verbose) {
    printf("\nVector -lambda\n");
    PrintVectorProbs(stdout, b, size);
  }

  /* calculate perm[i]
   */
  c = (double *) MallocOrDie (sizeof(double) * size); 
  for (i = 0; i < size; i++) c[i] = 100.0;
  
  OrderVector(&perm, b, size, isdescending, verbose);

  /* c[i] = w_p(1) + \sum_{j=0}^{n-i-1}w_p(n-j) - (n-i+1)w_p(i+1) for i=2,\ldots,n-1
   */
  for (i = 1; i < size-1; i++) 
    {
      c[i] = b[perm[0]] - (size-i+1)*b[perm[i+1]];

      for (j = i+1; j < size; j++) c[i] += b[perm[j]];
    }
  if (verbose) {
    printf("\nC vector \n");
    PrintVectorProbs(stdout, c, size);
  }

  /* calculate k_min = \min [ k;  2\leq k\leq n-1 such that C(k) <= 0 ]
   */
  k_min = size;
  for (i = 1; i < size-1; i++) {
    if (c[i] <= 0.0) { k_min = i; break; }
  }

  /* regularize:
   *
   *  if (2<=i<=k_min)  b(i) <--- 0
   *  else              b(i) <--- b(i) - corr,,   where: corr = 1/(n-k_{min}+1) \sum_{j=k_min+1} w_p(j)
   *
   */
  corr = b[perm[0]];
  for (i = k_min+1; i < size; i++) corr += b[perm[i]];
  corr /= (size-k_min+1);
  
  if (verbose) printf("k_min %d size %d corr %f\n", k_min, size, corr);

  for (i = 1; i < size; i++) 
    if (i <= k_min) b[perm[i]]  = 0.0;
    else            b[perm[i]] -= corr;

  if (verbose) {
    printf("\nregularized vector\n");
    PrintVectorProbs(stdout, b, size);
  }

  /* Reorder row to original form
   */
  for (i = 0; i < size; i++) r_1[i] = b[i];
    
  if (whichrow > 0) {
    b[whichrow] = r_1[0];

    for (i = 0; i < size; i++) if (i < whichrow) b[i] = r_1[i+1];
  }
  if (verbose) {
    printf("\nregularized reordered vector\n");
    PrintVectorProbs(stdout, b, size);
  }

  /* consistency test
   */
  sum = 0.0;
  for (i = 0; i < size; i++) {
    sum += b[i];

    if (i == whichrow) { 
      if (b[i] > 0.0) Die("QOGRegularizationAlgorithm(): positive diagnonal i=%d val=%f", i, b[i]); 
    }
    else {
      if (b[i] < -MARGIN) Die("QOGRegularizationAlgorithm(): negative off-diagnonal i=%d val=%f", i, b[i]); 
    }
  }
  if (fabs(sum) > MARGIN) Die("QOGRegularizationAlgorithm(): sum is different from zero sum=%f", sum); 

  if (verbose) {
    printf("\nQOG_regularized vector\n");
    PrintVectorProbs(stdout, b, size);
  }

  free(c);
  free(r_1);
  free(perm);
}

void
check_reversibility (double *QL, double *QR, double *ml, double *mr, int L)
{
  double *slr;
  double *srl;
  int     i,j;

  slr = (double *) MallocOrDie (sizeof(double) * L * L); 
  srl = (double *) MallocOrDie (sizeof(double) * L * L); 

  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) {
      slr[i*L+j] = QL[i*L+j]/mr[j];
      srl[i*L+j] = QR[i*L+j]/ml[j];
    }

    fprintf(stdout, "Slr probabilities\n");
    PrintProbs(stdout, slr, L);
    fprintf(stdout, "Srl probabilities\n");
    PrintProbs(stdout, srl, L);

    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++) 
	if (slr[i*L+j] - srl[j*L+i] > MARGIN || slr[i*L+j] - srl[j*L+i] < -MARGIN) 
	  Warn("reversibility check failed\n");
    
  free(slr);
  free(srl);
}

void
check_Q_0_reversibility (double *QL, double *QR, double *ml, double *mr, int L, int hasindel)
{
  double *slr;
  double *srl;
  double *pstatl;
  double *pstatr;
  int     dim;
  int     i,j;

  fprintf(stdout, "QL_0 probabilities\n");
  PrintProbs(stdout, QL, L);
  fprintf(stdout, "QR_0 probabilities\n");
  PrintProbs(stdout, QR, L);

  if (hasindel) dim = L-1;
  else          dim = L;
  
  slr    = (double *) MallocOrDie (sizeof(double) * dim * dim); 
  srl    = (double *) MallocOrDie (sizeof(double) * dim * dim); 
  pstatl = (double *) MallocOrDie (sizeof(double) * dim); 
  pstatr = (double *) MallocOrDie (sizeof(double) * dim); 
  
  for (i = 0; i < dim; i++) {
    if (hasindel) pstatl[i] = ml[i]/(1.0-ml[L-1]);
    else          pstatl[i] = ml[i];
    
    if (hasindel) pstatr[i] = mr[i]/(1.0-mr[L-1]);
    else          pstatr[i] = mr[i];
  }
  for (i = 0; i < dim; i++) 
    for (j = 0; j < dim; j++) {
      slr[i*dim+j] = QL[i*L+j]/pstatr[j];
      srl[i*dim+j] = QR[i*L+j]/pstatl[j];
    }
  
  fprintf(stdout, "Slr probabilities\n");
  PrintProbs(stdout, slr, dim);
  fprintf(stdout, "Srl probabilities\n");
  PrintProbs(stdout, srl, dim);
  
  for (i = 0; i < dim; i++) 
    for (j = 0; j < dim; j++) 
      if (slr[i*dim+j] - srl[j*dim+i] > MARGIN || slr[i*dim+j] - srl[j*dim+i] < -MARGIN) 
	Warn("reversibility check failed\n");
  
  free(slr);
  free(srl);
  free(pstatl);
  free(pstatr);
}

int
CompareFreqs(double *pml, double *pmr, double *targetfreq, int L)
{
  int flag = FALSE;
  int i;
  
  for (i = 0; i < L; i ++) 
    if (fabs(pml[i]-targetfreq[i]) > 1.0-accuracy1 || fabs(pmr[i]-targetfreq[i]) > 1.0-accuracy1) 
      flag = TRUE;

  return flag;
}
