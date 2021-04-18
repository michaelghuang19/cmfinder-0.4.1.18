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

#include "er_math.h"
#include "misc.h"
#include "matrix.h"
#include "evolve.h"
#include "mat_globals.h"
#include "squid.h"
#include "mat_structs.h"


struct eigenvalue_s *
AllocEigenvalues (int L)
{
  struct eigenvalue_s *eigen;
  int                  i;

  eigen       = (struct eigenvalue_s *) MallocOrDie (sizeof(struct eigenvalue_s));
  eigen->real = (double              *) MallocOrDie (sizeof(double) * L);
  eigen->imag = (double              *) MallocOrDie (sizeof(double) * L);

  /* initialize eigenvalues */
  for (i = 0; i < L; i++) {
    eigen->real[i] = 0.;
    eigen->imag[i] = 0.;
  }
  
  return eigen;
}

struct zmatrix_s *
AllocZmatrix (int L)
{
  struct zmatrix_s *zmatrix;
  int               i;

  zmatrix       = (struct zmatrix_s *) MallocOrDie (sizeof(struct zmatrix_s));
  zmatrix->real = (double           *) MallocOrDie (sizeof(double) * L);
  zmatrix->imag = (double           *) MallocOrDie (sizeof(double) * L);

  /* initialize eigenvalues */
  for (i = 0; i < L; i++) {
    zmatrix->real[i] = 0.;
    zmatrix->imag[i] = 0.;
  }
  
  return zmatrix;
}

/* Function: Cal_Id()
 * Date:     ER, Thu Apr 13 10:28:16 CDT 2000 [St. Louis]
 *
 * Purpose:  Creates a Id(LxL) matrix
 *         
 * Args:    L - dimension          
 *
 * Returns: Id(LxL) = \delta_{ij}
 *          Id is allocated here, freed by caller.
 *           
 */
double *
Cal_Id(int L)
{
  double *Id;
  int     i, j;

  Id = (double *) MallocOrDie (sizeof(double) * L * L);

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) 
      Id[i*L+j] = ((i == j) ? 1. : 0.);

  return Id;
}

struct zmatrix_s *
ZCal_Id(int L)
{
  struct zmatrix_s *Id;
  int     i, j;

  Id = AllocZmatrix(L*L);

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) {
      Id->real[i*L+j] = ((i == j) ? 1.0 : 0.0);
      Id->imag[i*L+j] = 0.0;
    }

  return Id;
}

void
Comp_Id(double *Id, int L)
{
  int i, j;

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) 
      Id[i*L+j] = ((i == j) ? 1. : 0.);
}

/* Function: Cal_LUP_Decomp()
 * Date:     ER, Mon Apr 17 13:43:07 CDT 2000 [St. Louis]
 *
 * Purpose:  given M calculate a LUP decomposition, stored in LU
 *  
 *           (Implemented after reading chapter 31 of "introduction to algorithms"
 *            by Cormen et al.)
 *            Find matrices P, L(lower triangular), U(uppper triangular)
 *            such that P*M = L*U
 *
 *           if l_ij with i<j  are the elements of L  (l_ii = 1)
 *              u_ij with i>=j are the elements of U  (u_ii = m_ii)
 *
 *                 | l_ij for i  > j (we forget about the one's in the diagonal of L)
 *           LU  = |
 *                 | u_ij for i <= j
 *
 * Method: 
 *           Inspect by column k  
 *           Inspect rows from k on, until u_ik is a max value.
 *           Store the permutation in P.
 *           Permute rows k with i both in m and LU matrices.
 *
 *           u_ik =            m_ik - \sum_{j=1}^{j=i-1} [ l_ij * u_jk ]
 *
 *           l_ik = 1/u_kk * { m_ik - \sum_{j=1}^{j=k-1} [ l_ij * u_jk ]}
 *
 *
 *
 * Args:     M - LxL non-singular matrix.
 *
 * Returns:  LUP(LxL) = M^{-1}(LxL)
 *           LUP is alocated here, freed by caller.
 */
double *
Cal_LUP_Decomp(FILE *ofp, double *M, int *perm, int L, int verbose)
{
  double *LU;
  double *m;
  double *Per;
  double *Pro;
  double  val;
  double  exchange;
  int     new_row, new_perm;
  int     L2;
  int     i, j, k;
 
  L2 = L * L;

  /* allocate LUP[L2]
   */
  LU = (double *) MallocOrDie (sizeof(double) * L2);
  m  = (double *) MallocOrDie (sizeof(double) * L2);


  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) {
      m[i*L+j]  = M[i*L+j];
      LU[i*L+j] = 0.;
    }

  for (i = 0; i < L; i++) 
    perm[i] = i; /* initializa to no permutations */

  /* start moving by columns
   */
  for (k = 0; k < L; k++) {

    val = 0.;
    /* For a given column move rows to have the largest possible value
     * in the diagonal.
     */
    for (i = k; i < L; i++) {
      if (fabs(m[i*L+k]) > val) {
	val = fabs(m[i*L+k]);
	new_row = i;
      }
    }
    if (val == 0.) Die( "this matrix is singular. No LUP decomposition possible");

    /* register in the P matrix that you did a permutation
     */
    new_perm      = perm[k];
    perm[k]       = perm[new_row];
    perm[new_row] = new_perm;
    
    /* exchange values of the two rows both for LU and m
     */
    for (i = 0; i < L; i++) {
      exchange       = m[k*L+i];
      m[k*L+i]       = m[new_row*L+i];
      m[new_row*L+i] = exchange;

      exchange        = LU[k*L+i];
      LU[k*L+i]       = LU[new_row*L+i];
      LU[new_row*L+i] = exchange;
    }
    
    LU[k*L+k] = m[k*L+k];
    
    for (i = k+1; i < L; i++) {
      LU[i*L+k] = m[i*L+k]/LU[k*L+k];  /* L part of the LUP matrix */
      LU[k*L+i] = m[k*L+i];            /* U part of the LUP matrix */
    }
    
    /* redefine rest of matrix as m_ij -= l_ik * u_kj
     */
    for (i = k+1; i < L; i++) 
      for (j = k+1; j < L; j++) 
	m[i*L+j] -= LU[i*L+k]*LU[k*L+j];
	
  }
  if (verbose) {
    Per = (double *) MallocOrDie (sizeof(double) * L2);
    Pro = (double *) MallocOrDie (sizeof(double) * L2);

    fprintf(ofp, "LU matrix\n");
    PrintProbs(ofp, LU, L);

    fprintf(ofp, "L * U matrix\n");
    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++) 
	m[i*L+j] = 0.;

    for (i = 0; i < L; i++) 
      for (k = 0; k < L; k++) {
	if (i <= k) m[i*L+k] += LU[i*L+k];

	for (j = 0; j <= ((i-1 < k)?i-1:k); j++) 
	  m[i*L+k] += LU[i*L+j] * LU[j*L+k];
      }
    
    PrintProbs(ofp, m, L);
    
    fprintf(ofp, "Permutation matrix\n");
    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++) 
	Per[i*L+j] = (perm[i] == j)? 1. : 0.;
    PrintProbs(ofp, Per, L);
    
    /* test Per*M = L*U ? 
     */
    fprintf(ofp, "Per*M matrix\n");
    Pro = (double *)Cal_M_N_Prod(ofp, Per, M, L, L, L, verbose);
    PrintProbs(ofp, Pro, L);
    
    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++) {
	if (fabs(Pro[i*L+j]-m[i*L+j]) > MARGIN )
	    Die("Cal_LUP_Decomp(): bad LUP decomposion");
      }
    
    free (Per);
    free (Pro);
  }
  
  free(m);

  return LU;
}

struct zmatrix_s *
ZCal_LUP_Decomp(FILE *ofp, struct zmatrix_s *M, int *perm, int L, int verbose)
{
  struct zmatrix_s *LU;
  struct zmatrix_s *m;
  struct zmatrix_s *Per;
  struct zmatrix_s *Pro;
  double  val;
  double  nval;
  double  exchange_real;
  double  exchange_imag;
  int     new_row, new_perm;
  int     L2;
  int     i, j, k;

  L2 = L * L;

 /* allocate LUP[L2]
   */
  LU = AllocZmatrix(L2);
  m  = AllocZmatrix(L2);

  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) {
      m->real[i*L+j]  = M->real[i*L+j];
      m->imag[i*L+j]  = M->imag[i*L+j];
      LU->real[i*L+j] = 0.;
      LU->imag[i*L+j] = 0.;
    }

  for (i = 0; i < L; i++) 
    perm[i] = i; /* initializa to no permutations */

  /* start moving by columns
   */
  for (k = 0; k < L; k++) {

    val = 0.;

    /* For a given column move rows to have the largest possible value
     * in the diagonal.
     */
    for (i = k; i < L; i++) {
      nval = m->real[i*L+k]*m->real[i*L+k]+ m->imag[i*L+k]*m->imag[i*L+k];

      if (nval > val) {
	val = nval;
	new_row = i;
      }
    }
    if (val == 0.) Die( "this matrix is singular. No LUP decomposition possible");

    /* register in the P matrix that you did a permutation
     */
    new_perm      = perm[k];
    perm[k]       = perm[new_row];
    perm[new_row] = new_perm;
    
    /* exchange values of the two rows both for LU and m
     */
    for (i = 0; i < L; i++) {
      exchange_real        = m->real[k*L+i];
      exchange_imag        = m->imag[k*L+i];
      m->real[k*L+i]       = m->real[new_row*L+i];
      m->imag[k*L+i]       = m->imag[new_row*L+i];
      m->real[new_row*L+i] = exchange_real;
      m->imag[new_row*L+i] = exchange_imag;

      exchange_real         = LU->real[k*L+i];
      exchange_imag         = LU->imag[k*L+i];
      LU->real[k*L+i]       = LU->real[new_row*L+i];
      LU->imag[k*L+i]       = LU->imag[new_row*L+i];
      LU->real[new_row*L+i] = exchange_real;
      LU->imag[new_row*L+i] = exchange_imag;
    }
    
    LU->real[k*L+k] = m->real[k*L+k];
    LU->imag[k*L+k] = m->imag[k*L+k];
    
    for (i = k+1; i < L; i++) {
      /* L part of the LUP matrix */
      LU->real[i*L+k] = (m->real[i*L+k] * LU->real[k*L+k] + m->imag[i*L+k] * LU->imag[k*L+k]) / 
	(LU->real[k*L+k] * LU->real[k*L+k] + LU->imag[k*L+k] * LU->imag[k*L+k]);  

      LU->imag[i*L+k] = (m->imag[i*L+k] * LU->real[k*L+k] - m->real[i*L+k]  * LU->imag[k*L+k]) / 
	(LU->real[k*L+k] * LU->real[k*L+k] + LU->imag[k*L+k] * LU->imag[k*L+k]);

      /* U part of the LUP matrix */
      LU->real[k*L+i] = m->real[k*L+i];           
      LU->imag[k*L+i] = m->imag[k*L+i];         
    }
    
    /* redefine rest of matrix as m_ij -= l_ik * u_kj
     */
    for (i = k+1; i < L; i++) 
      for (j = k+1; j < L; j++) {
	m->real[i*L+j] -= LU->real[i*L+k]*LU->real[k*L+j] -  LU->imag[i*L+k]*LU->imag[k*L+j];
	m->imag[i*L+j] -= LU->real[i*L+k]*LU->imag[k*L+j] +  LU->imag[i*L+k]*LU->real[k*L+j];
      }	

  }

  if (TRUE) {
    Per = AllocZmatrix(L2);
    Pro = AllocZmatrix(L2);
    
    fprintf(ofp, "LU matrix\n");
    PrintZmatrix(ofp, LU, L);

    fprintf(ofp, "L * U matrix\n");
    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++) {
	m->real[i*L+j] = 0.;
	m->imag[i*L+j] = 0.;
      }

    for (i = 0; i < L; i++) 
      for (k = 0; k < L; k++) {
	if (i <= k) {
	  m->real[i*L+k] += LU->real[i*L+k];
	  m->imag[i*L+k] += LU->imag[i*L+k];
	}

	for (j = 0; j <= ((i-1 < k)?i-1:k); j++) {
	  m->real[i*L+k] += LU->real[i*L+j] * LU->real[j*L+k] - LU->imag[i*L+j] * LU->imag[j*L+k];
	  m->imag[i*L+k] += LU->imag[i*L+j] * LU->real[j*L+k] + LU->real[i*L+j] * LU->imag[j*L+k];
	}
      }
    
    PrintZmatrix(ofp, m, L);
    
    fprintf(ofp, "Per matrix\n");
    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++) {
	Per->real[i*L+j] = (perm[i] == j)? 1. : 0.;
	Per->imag[i*L+j] = 0.;
      }
    PrintZmatrix(ofp, Per, L);


    /* test Per*M = L*U ? 
     */
    fprintf(ofp, "Per*M matrix\n");
    Pro = (struct zmatrix_s *)ZCal_M_N_Prod(ofp, Per, M, L, L, L, verbose);
    PrintZmatrix(ofp, Pro, L);
     
    for (i = 0; i < L; i++) 
      for (j = 0; j < L; j++) {
	if (fabs(Pro->real[i*L+j]-m->real[i*L+j]) > MARGIN || fabs(Pro->imag[i*L+j]-m->imag[i*L+j]) > MARGIN )
	    Die ("ZCal_LUP_Decomp(): bad LUP decompositon %d %d (%f %f) (%f %f)", 
		i, j, Pro->real[i*L+j], Pro->imag[i*L+j], m->real[i*L+j], m->imag[i*L+j]);
      }
    
    FreeZmatrix(Per);
    FreeZmatrix(Pro);
  }
  
  FreeZmatrix(m);

  return LU;
}


/* Function: Cal_M_Exp()
 * Date:     ER, Wed Mar  1 11:04:24 CST 2000 [St. Louis]
 *
 * Purpose:  Given a matrix M, calculate exp{r*M}
 *
 *            exp{rM} = \sum_{n=0}^{\infty} [ r^n * M^n / n!   ] 
 *
 * Args:     M  - LL joint prob matrix (prealloc)
 *
 * Returns:  Q(LxL) = exp{rM(LxL)}
 *           Q is alocated here, freed by caller.
 */
double *
Cal_M_Exp(FILE *ofp, double *M, int L, double r, int verbose)
{
  double *Qr;        /* Qr = exp{r*M} matrix */
  double *taylorQr;  /* next term for Qr in the taylor expansion */
  double *Mr;        /* holds r*M */
  double *Mpower;    /* holds at a given n (Mr)^n */
  double  coeff;
  int     L2;
  int     stop = FALSE;
  int     i;
  int     n = 0;
  
  L2 = L * L;

  /* allocate Qr[L2]
   */
  Mr       = (double *) MallocOrDie (sizeof(double) * L2);
  taylorQr = (double *) MallocOrDie (sizeof(double) * L2);

  /* Mr - r*M */
  CopyMatrix (Mr, M, L, L);
  MultiplyMatrix (r, L, Mr);

  /* inititalize Mpower to I and Qr to I
   */
  Mpower = Cal_Id(L);
  Qr     = Cal_Id(L);

  coeff  = 0.0; 
  while (!stop) {
    n++;

    coeff += log(1.0/(double)n);
    
    /* calculate Mr^{n-1}*Mr
     */
    Comp_M_N_Prod(ofp, Mr, Mpower, L, FALSE); /* multipL Mpower*M */
    
    for (i = 0; i < L2; i++) {

      if (verbose) {
	fprintf(ofp, "Mpower  coeff = %f\n", coeff);
	PrintProbs(ofp, Mpower, L);
      }

      taylorQr[i] = exp(coeff) * Mpower[i];
      if (abs(taylorQr[i]) > DBL_MAX || abs(Qr[i]-taylorQr[i]) > DBL_MAX) 
	Die ("sorry I did not reach convergence before double-float limit in Cal_M_Exp()\n");
    }
    
    if (Check_Accuracy(taylorQr, L*L)) {
      for (i = 0; i < L2; i++) { 
	Qr[i] += taylorQr[i];
      }
    }
    else {
      stop = TRUE; 
      if (verbose) fprintf(ofp, "exp(r*K) convergence in %d terms\n", n);
    }   
  }
  
  if (verbose) {
    fprintf(ofp, "exp{r*M} matrix\n");
    PrintProbs(ofp, Qr, L);
  }

  free(Mr);
  free(Mpower);
  free(taylorQr);

  return Qr;
}

struct zmatrix_s *
ZCal_M_Exp(FILE *ofp, struct zmatrix_s *M, int L, double r, int verbose)
{
  struct zmatrix_s *Qr;        /* Qr = exp{r*M} matrix */
  struct zmatrix_s *taylorQr;  /* next term for Qr in the taylor expansion */
  struct zmatrix_s *Mpower;    /* holds at a given n (M_I)^n */ double  coeff;
  int     L2;
  int     stop = FALSE;
  int     i;
  int     n = 0;
  
  L2 = L * L;

  /* allocate Qr[L2]
   */
  taylorQr = AllocZmatrix(L2);

  /* inititalize Mpower to I and Qr to I
   */
  Mpower = ZCal_Id(L);
  Qr     = ZCal_Id(L);

  coeff  = 0.0; 
  while (!stop) {
    n++;

    coeff += log(r/n);
    
    /* calculate M^{n-1}*M
     */
    ZComp_M_N_Prod(ofp, M, Mpower, L, FALSE); /* multipL Mpower*M */
    
    for (i = 0; i < L2; i++) {
      taylorQr->real[i] = exp(coeff) * Mpower->real[i];
      taylorQr->imag[i] = exp(coeff) * Mpower->imag[i];

      if (abs(taylorQr->real[i]) > DBL_MAX || abs(Qr->real[i]-taylorQr->real[i]) > DBL_MAX) 
	Die ("sorry I did not reach convergence before double-float limit in ZCal_M_Exp()\n");
      if (abs(taylorQr->imag[i]) > DBL_MAX || abs(Qr->imag[i]-taylorQr->imag[i]) > DBL_MAX) 
	Die ("sorry I did not reach convergence before double-float limit in ZCal_M_Exp()\n");
    }
    
    if (Check_Accuracy(taylorQr->real, L*L)) {
      for (i = 0; i < L2; i++) { 
	Qr->real[i] += taylorQr->real[i];
	Qr->imag[i] += taylorQr->imag[i];
      }
    }
    else if (Check_Accuracy(taylorQr->imag, L*L)) {
      for (i = 0; i < L2; i++) { 
	Qr->real[i] += taylorQr->real[i];
	Qr->imag[i] += taylorQr->imag[i];
      }
    }
    else {
      stop = TRUE; 
      if (verbose) fprintf(ofp, "exp(r*K) convergence in %d terms\n", n);
    }   
  }
  
  if (verbose) {
    fprintf(ofp, "exp{r*M} matrix\n");
    PrintZmatrix(ofp, Qr, L);
  }

  FreeZmatrix(Mpower);
  FreeZmatrix(taylorQr);

  return Qr;
}

void
Comp_M_Exp(FILE *ofp, double *M, int L, double r, int verbose)
{
  double *Qr;        /* Qr = exp{r*M} matrix */

  Qr = Cal_M_Exp(ofp, M, L, r, verbose);

  CopyMatrix (M, Qr, L, L);

  if (verbose) {
    fprintf(ofp, "exp{r*M} matrix\n");
    PrintProbs(ofp, M, L);
  }

  free(Qr);
}

/* Function: Cal_M_Inv()
 * Date:     ER, Mon Apr 17 13:40:50 CDT 2000 [St. Louis]
 *
 * Purpose:  given M calculate M^{-1} using a LUP decomposition
 *
 * Method:   PA = LU => LUA^{-1} = P
 *
 * Args:     M^{-1} - LxL conditional prob matrix (prealloc)
 *
 * Returns:  Minv(LxL) = M^{-1}(LxL)
 *           Minv is alocated here, freed by caller.
 */
double *
Cal_M_Inv(FILE *ofp, double *M, int L, int verbose)
{
  double *Minv;     /* Minv  = M^{-1} matrix        */
  double *LU;       /* LU = L*U = PM matrix         */
  double *b;        /* LU = L*U = M matrix          */
  int    *perm;     /* P matrix                     */
  int     L2;
  int     i, j, k;
  double  sum;

  L2 = L * L;

  if (verbose) {
    fprintf(ofp, "M matrix\n");
    PrintProbs(ofp, M, L);
 }

  /* allocate Minv[L2], LU
   */
  Minv  = (double *) MallocOrDie (sizeof(double) * L2);
  b     = (double *) MallocOrDie (sizeof(double) * L);
  perm  = (int    *) MallocOrDie (sizeof(int)    * L);

  LU = Cal_LUP_Decomp(ofp, M, perm, L, verbose);

  for (j = 0; j < L; j++) {
    for (i = 0; i < L; i++) 
      if (j == perm[i]) b[i] = 1.;
      else              b[i] = 0.;
    
    Solve_LU_Eq(ofp, LU, b, L, verbose);

    for (i = 0; i < L; i++) 
      Minv[i*L+j] = b[i];
  }
  
  if (verbose) {
     fprintf(ofp, "M^{-1} matrix\n");
     PrintProbs(ofp, Minv, L);

    fprintf(ofp, "M * M^{-1} matrix\n");
    for (i = 0; i < L; i++) {
      for (k = 0; k < L; k++) {
	sum = 0.;
	for (j = 0; j < L; j++) 
	  sum += M[i*L+j] * Minv[j*L+k];
	
	fprintf(ofp, "%.4f ", sum);

	if (i == k && (sum > 1.0+MARGIN || sum < 1.0-MARGIN)) Die ("Cal_M_Inv(): bad inverse matrix");
	if (i != k && (fabs(sum) >MARGIN))                    Die ("Cal_M_Inv(): bad inverse matrix");
      }
      fprintf(ofp, "\n");
    }
  }
  
  free(LU);
  free(perm);
  free(b);

  return Minv;
}

struct zmatrix_s *  
ZCal_M_Inv(FILE *ofp, struct zmatrix_s *M, int L, int verbose)
{
  struct zmatrix_s *Minv;     /* Minv  = M^{-1} matrix        */
  struct zmatrix_s *LU;       /* LU = L*U = PM matrix         */
  struct zmatrix_s *b;        /* LU = L*U = M matrix          */
  int    *perm;               /* P matrix                     */
  int     L2;
  int     i, j, k;
  double  sum_real;
  double  sum_imag;

  L2 = L * L;

  if (verbose) {
    fprintf(ofp, "M matrix\n");
    PrintZmatrix(ofp, M, L);
 }

  /* allocate Minv[L2], LU
   */
  Minv  = AllocZmatrix(L2);
  b     = AllocZmatrix(L);
  perm  = (int    *) MallocOrDie (sizeof(int)    * L);

  LU = ZCal_LUP_Decomp(ofp, M, perm, L, verbose);

  for (j = 0; j < L; j++) {
    for (i = 0; i < L; i++) 
      if (j == perm[i]) { b->real[i] = 1.; b->imag[i] = 0.; }
      else              { b->real[i] = 0.; b->imag[i] = 0.; }
    
    ZSolve_LU_Eq(ofp, LU, b, L, TRUE);

    for (i = 0; i < L; i++) {
      Minv->real[i*L+j] = b->real[i];
      Minv->imag[i*L+j] = b->imag[i];
    }
  }
  
  if (TRUE) {
     fprintf(ofp, "M^{-1} matrix\n");
     PrintZmatrix(ofp, Minv, L);

    fprintf(ofp, "M * M^{-1} matrix\n");
    for (i = 0; i < L; i++) {
      for (k = 0; k < L; k++) {
	sum_real = 0.;
	sum_imag = 0.;
	for (j = 0; j < L; j++) {
	  sum_real += M->real[i*L+j] * Minv->real[j*L+k] - M->imag[i*L+j] * Minv->imag[j*L+k];
	  sum_imag += M->real[i*L+j] * Minv->imag[j*L+k] + M->imag[i*L+j] * Minv->real[j*L+k];
	}
	
	fprintf(ofp, "(%.4f %.4f) ", sum_real, sum_imag);

	if (i == k && (sum_real > 1.0+MARGIN2 || sum_real < 1.0-MARGIN2)) Die ("Cal_M_Inv(): bad inverse matrix");
	if (i != k && (fabs(sum_real) > MARGIN2))                         Die ("Cal_M_Inv(): bad inverse matrix");
	if (fabs(sum_imag) > MARGIN2)                                     Die ("Cal_M_Inv(): bad inverse matrix");

      }
      fprintf(ofp, "\n");
    }
  }
  
  FreeZmatrix(LU);
  FreeZmatrix(b);
  free(perm);

  return Minv;

}

void
Comp_M_Inv(FILE *ofp, double *M, double *MI, int L, int verbose)
{
  double *Minv;     /* Minv  = M^{-1} matrix        */

  Minv = Cal_M_Inv(ofp, M, L, verbose);

  CopyMatrix (MI, Minv, L, L);

  if (verbose) {
    fprintf(ofp, "M^{-1} matrix\n");
    PrintProbs(ofp, MI, L);
  }

  free(Minv);
}


/* Function: Cal_M_Eigenvalues()
 * Date:     ER, Mon Sep 29 08:55:28 CDT 2003 [St. Louis]
 *
 * Purpose:  given M calculate its eigenvalues using a QR_decomposition
 *
 * Args:     M  - LxL matrix (prealloc)
 *
 * Returns:  K(LxL) = logM(LxL)
 *           K is alocated here, freed by caller.
 */
struct eigenvalue_s *
Cal_M_Eigenvalues(FILE *ofp, double *M, int L, int verbose)
{
  struct eigenvalue_s *eigen;
  double *H;
 
  H  = (double *)HessenbergForm(ofp, M, L, verbose);
  Hessenberg2Eigenvalues(ofp, H, L, &eigen, FALSE, verbose);

  free(H);

  return eigen;
}

/* Function: Cal_M_Vandemonde_Log()
 * Date:     ER, Mon Sep 29 08:50:49 CDT 2003 [St. Louis]
 *
 * Purpose:  given M calculate K = log{M} using eigenvalues and the Vandemonde matrix.
 *
 * Args:     M  - LxL conditional prob matrix (prealloc)
 *
 * Returns:  K(LxL) = logM(LxL)
 *           K is alocated here, freed by caller.
 */
struct zmatrix_s *
Cal_M_Vandermonde_Log(FILE *ofp, double *M, int L, int verbose)
{
  struct eigenvalue_s *eigen;
  struct zmatrix_s    *ZM;
  struct zmatrix_s    *Id;
  struct zmatrix_s    *Vandemonde;
  struct zmatrix_s    *Vandemonde_inv;
  struct zmatrix_s    *K;               /* K = logM matrix */
  struct zmatrix_s    *log_eval;
  struct zmatrix_s    *v;               
  struct zmatrix_s    *m;
  struct zmatrix_s    *mv;
  double              *Q;
  double               sum;
  double               rad;
  double               angle;
  int                  i,j;
  
  /* allocate memory */
  ZM         = AllocZmatrix(L*L);
  Vandemonde = AllocZmatrix(L*L);
  v          = AllocZmatrix(L);
  log_eval   = AllocZmatrix(L);
  mv         = AllocZmatrix(L*L);

  /* calculate eigenvalues */
  eigen = Cal_M_Eigenvalues(ofp, M, L, verbose);

  for (i = 0; i < L; i++) 
    for (j = i+1; j < L; j++)     
      if (fabs(eigen->real[i]-eigen->real[j]) < 1.0-accuracy1 && 
	  fabs(eigen->imag[i]-eigen->imag[j]) < 1.0-accuracy1 ) 
	Die ("Cal_M_Vandemonde_Log(): only applies to non-degenearte eigenvalues: diff = %f ", eigen->real[i]-eigen->real[j]);
  
  for (i = 0; i < L; i++) {
    rad   = sqrt(eigen->real[i]*eigen->real[i] + eigen->imag[i]*eigen->imag[i]);
    angle = atan(eigen->imag[i]/eigen->real[i]);
    
    if (eigen->real[i] > 0.0 ) {
      log_eval->real[i] = log(rad);
      log_eval->imag[i] = angle;
    }
    
    if (eigen->real[i] < 0.0 ) {
      log_eval->real[i] = log(rad);
      log_eval->imag[i] = PI+angle;
    }
    
  }
  if (verbose) {
    fprintf(ofp, "Log evalues\n");
    PrintZvector(ofp, log_eval, L);
  }
    
  /* The Vandermonde matrix */
  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) {
      Vandemonde->real[i*L+j] = exp(j*log_eval->real[i]) * cos(j*log_eval->imag[i]);
      Vandemonde->imag[i*L+j] = exp(j*log_eval->real[i]) * sin(j*log_eval->imag[i]);
    }
  
  if (TRUE) {
    fprintf(ofp, "Vandermonde Matrix\n");
    PrintZmatrix(ofp, Vandemonde, L);
  }
  
  Vandemonde_inv = ZCal_M_Inv(ofp, Vandemonde,  L, verbose);
  
  v = (struct zmatrix_s *)ZCal_M_N_Prod(ofp, Vandemonde_inv, log_eval, L, L, 1, verbose); 
  if (TRUE) {
    fprintf(ofp, "coeff\n");
    PrintZvector(ofp, v, L);
  }
  sum = 0.0;
  for (i = 0; i < L; i++) 
    sum += v->real[i];
  if(fabs(sum) > MARGIN2) Die ("real_coeff should add up to zero, sum = %f", sum);
  sum = 0.0;
  for (i = 0; i < L; i++) 
    sum += v->imag[i];
  if(fabs(sum) > MARGIN2) Die ("inag_coeff should add up to zero, sum = %f", sum);

  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) 
      ZM->real[i*L+j] = M[i*L+j];

  m = (struct zmatrix_s *)ZCal_Id(L);
  K = (struct zmatrix_s *)ZCal_M_num_Prod(ofp, m, v->real[0], v->imag[0], L, L, verbose);

  for (i = 1; i < L; i++) {
    ZComp_M_N_Prod(ofp, ZM, m, L, verbose); /* m *= M */

    CopyZmatrix(mv, m, L, L);          

    MultiplyZmatrix(v->real[i], v->imag[i], L, mv);  /* mv = m*v_i */

    ZComp_M_N_Sum(ofp, mv, K, L, L, TRUE, verbose); /* K += mv */
  }

  if (TRUE) { 
    fprintf(ofp, "Vandermonde Rate Matrix\n");
    PrintZmatrix(ofp, K, L);
  }

  /* reconstruct the exponential matrix
   */
  if (TRUE) {
    Id = ZCal_Id(L);
    ZCondi_From_Rate(&Q, K, Id, 1.0, L, FALSE, verbose);
   
    fprintf(ofp, "EXP Matrix\n");
    PrintMatrix(ofp, Q, L);
    
    Comp_M_N_Sum(ofp, M, Q, L, L, FALSE, verbose); /* M-Q */

    free (Q);
    FreeZmatrix(Id);
  }
  
  FreeEigenvalues(eigen);
  FreeZmatrix(log_eval);
  FreeZmatrix(ZM);
  FreeZmatrix(Vandemonde);
  FreeZmatrix(Vandemonde_inv);
  FreeZmatrix(m);
  FreeZmatrix(mv);
  FreeZmatrix(v);

  return K;
}


/* Function: Cal_M_Taylor_Log()
 * Date:     ER, Wed Mar  1 11:04:24 CST 2000 [St. Louis]
 *
 * Purpose:  given M calculate K = log{M} using a taylor series.
 *
 *           K = \sum_{n=1}^{\infty} [ (-)^{n+1}(M-I)^n / n]
 *
 * Args:     Q  - LxL conditional prob matrix (prealloc)
 *
 * Returns:  K(LxL) = logM(LxL)
 *           K is alocated here, freed by caller.
 */
double *
Cal_M_Taylor_Log(FILE *ofp, double *M, int L, int verbose)
{
  double *K;        /* K = logM matrix                         */
  double *taylorK;  /* next term for K in the taylor expansion */
  double *Mmin;     /* holds M - I                             */
  double *Mpower;   /* holds at a given n (M_I)^n              */
  double  coeff;
  double  epsilon;
  int     L2;
  int     stop = FALSE;
  int     i;
  int     n = 0;

  L2 = L * L;

  epsilon = 1.0;

  /* allocate K[L2]
   */  
  taylorK = (double *) MallocOrDie (sizeof(double) * L2);

  /* Multiply epsilon*M
   */
  MultiplyMatrix(epsilon, L, M);

  /* substract I from M
   */
  Mmin = (double *)Cal_M_MinId(M, L);  /* M --> Mmin = M-I */

  if (verbose) {
    fprintf(ofp, "Mmin matrix\n");
    PrintProbs(ofp, Mmin, L);
  }

  /* inititalize Mpower to I and K to -log(epsilon)I
   */
  Mpower = Cal_Id(L);
  K      = Cal_Id(L); 

  for (i = 0; i < L; i++) 
    K[i*L+i] = -log(epsilon);
  
  while (!stop) {
    n ++;
    
    Comp_M_N_Prod(ofp, Mmin, Mpower, L, FALSE); /* multiply Mpower*Mmin */
    
    coeff = (n%2==1)? 1./n : -1./n;

    for (i = 0; i < L2; i++) {
      taylorK[i] = coeff * Mpower[i];
      if (taylorK[i] > DBL_MAX) 
	Die ("sorry I did not reach convergence before double-float limit in Cal_M_Taylor_Log() \n");
    }
     
  if (FALSE) {
    fprintf(ofp, "logQ taylor matrix\n");
    PrintProbs(ofp, taylorK, L);
  }

    if (Check_Accuracy(taylorK, L*L)) {
      for (i = 0; i < L2; i++) 
	K[i] += taylorK[i];
    }
    else {
      stop = TRUE; 
      if (verbose) fprintf(ofp, "logM convergence in %d terms\n", n);
    }   
  }
  
  if (verbose) {
    fprintf(ofp, "logM matrix\n");
    PrintProbs(ofp, K, L);
  }
  
  free(taylorK);
  free(Mmin);
  free(Mpower);

  return K;
}

void
Comp_M_Taylor_Log(FILE *ofp, double *M, int L, int verbose)
{
  double *K;        /* K = logM matrix                         */

  K = Cal_M_Taylor_Log(ofp, M, L, verbose);

   CopyMatrix (M, K, L, L);

  if (verbose) {
    fprintf(ofp, "log{M} matrix\n");
    PrintProbs(ofp, M, L);
  }

 free(K);
}

/* Function: Cal_M_MinId()
 * Date:     ER, Thu Apr 13 10:30:38 CDT 2000 [St. Louis]
 *
 * Purpose:  given M calculates M - Id.
 *
 * Args:     M  - LxL (prealloc)
 *
 * Returns:  Mmin(LixLj) = M(LxL) - ID
 *           Mmin is alocated here, freed by caller.
 */
double *
Cal_M_MinId(double *M, int L)
{
  double *Mmin;
  int     i, j;

  Mmin = (double *) MallocOrDie (sizeof(double) * L * L);
  
  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) 
      Mmin[i*L+j] =  M[i*L+j] - ((i==j)? 1. : 0);
  
  return Mmin;
}

void
Comp_M_MinId(double *M, int L)
{
  int i, j;

  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) 
      M[i*L+j] -= ((i==j)? 1. : 0);
  
}

/* Function: Cal_M_num_Prod()
 * Date:     ER,  [St. Louis]
 *
 * Purpose:  given M (LixLj) and N (LixLj) calculate Msum(LixLj) = M + N.
 *
 * Args:     M    - LixLj (prealloc)
 *           num  - a double
 *
 * Returns:  Mprod(LixLj) = M(LixLj) * num
 *           Mprod is alocated here, freed by caller.
 */
double *
Cal_M_num_Prod(FILE *ofp, double *M, double num, int Li, int Lj, int verbose)
{
  int     i, j;
  double *Mprod;
  
  Mprod = (double *) MallocOrDie (sizeof(double) * Li * Lj);
  
  for (i = 0; i < Li; i++)
    for (j = 0; j < Lj; j++)
	Mprod[i*Lj+j] = M[i*Lj+j] * num;
  
  if (verbose) {
    fprintf(ofp, "num * M matrix\n");
    
    for (i = 0; i < Li; i++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "%.4f ", Mprod[i*Lj+j]);
      fprintf(ofp, "\n");
    }
  }
  return Mprod;
}

struct zmatrix_s *
ZCal_M_num_Prod(FILE *ofp, struct zmatrix_s *M, double real, double imag, int Li, int Lj, int verbose) 
{
  struct zmatrix_s *Mprod;
  int             i, j;
  
  Mprod = AllocZmatrix(Li*Lj);

  for (i = 0; i < Li; i++)
    for (j = 0; j < Lj; j++) {
	Mprod->real[i*Lj+j] = M->real[i*Lj+j] * real - M->imag[i*Lj+j] * imag;
	Mprod->imag[i*Lj+j] = M->real[i*Lj+j] * imag + M->imag[i*Lj+j] * real;
    }
  
  if (verbose) {
    fprintf(ofp, "num * M zmatrix\n");
    
    for (i = 0; i < Li; i++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "(%.4f %.4f) ", Mprod->real[i*Lj+j], Mprod->imag[i*Lj+j]);
      fprintf(ofp, "\n");
    }
  }

  return Mprod;
}

void
Comp_M_num_Prod(FILE *ofp, double *M, double num, int Li, int Lj, int verbose)
{
  int  i, j;
  
  for (i = 0; i < Li; i++)
    for (j = 0; j < Lj; j++)
	M[i*Lj+j] *= num;
  
  if (verbose) {
    fprintf(ofp, "num * M matrix\n");
    
    for (i = 0; i < Li; i++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "%.4f ", M[i*Lj+j]);
      fprintf(ofp, "\n");
    }
  }
}


/* Function: Cal_M_N_Prod()
 * Date:     ER, Thu Apr 13 10:30:38 CDT 2000 [St. Louis]
 *
 * Purpose:  given M (LixLk) and N (LkxLj) calculate Mprod(LixLk) = M * N.
 *
 * Args:     M  - LixLk (prealloc)
 *           N  - LkxLj (prealloc)
 *
 * Returns:  Mprod(LixLj) = M(LixLk) * N(LkxLj)
 *           Mprod is alocated here, freed by caller.
 */
double *
Cal_M_N_Prod(FILE *ofp, double *M, double *N, int Li, int Lk, int Lj, int verbose)
{
  int     i, j, k;
  double  prod;
  double *Mprod;

  Mprod = (double *) MallocOrDie (sizeof(double) * Li * Lj);
  
  if (verbose) {
    fprintf(ofp, "M matrix\n");
    
    for (i = 0; i < Li; i++) {
      for (k = 0; k < Lk; k++) 
	fprintf(ofp, "%.4f ", M[i*Lk+k]);
      fprintf(ofp, "\n");
    }

    fprintf(ofp, "N matrix\n");
    
    for (k = 0; k < Lk; k++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "%.4f ", N[k*Lj+j]);
      fprintf(ofp, "\n");
    }
  }

  for (i = 0; i < Li; i++)
    for (j = 0; j < Lj; j++) {
      prod = 0.;
      for (k = 0; k < Lk; k++) {
	prod += M[i*Lk+k] * N[k*Lj+j];
	if (prod > DBL_MAX) 
	  Die("sorry, I am getting out of bounds in Cal_M_N_Prod() prod[%d][%d] = %f\n", i, j, prod);
      }
      Mprod[i*Lj+j] = prod;
    }

  if (verbose) {
    fprintf(ofp, "M*N matrix\n");
    
    for (i = 0; i < Li; i++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "%.4f ", Mprod[i*Lj+j]);
      fprintf(ofp, "\n");
    }
  }
  return Mprod;
}

struct zmatrix_s *
ZCal_M_N_Prod(FILE *ofp, struct zmatrix_s *M, struct zmatrix_s *N, int Li, int Lk, int Lj, int verbose)
{
  struct zmatrix_s *Mprod;
  double            prod_real;
  double            prod_imag;
  int               i, j, k;

  Mprod = AllocZmatrix(Li * Lj);
  
  if (verbose) {
    fprintf(ofp, "M matrix\n");
    
    for (i = 0; i < Li; i++) {
      for (k = 0; k < Lk; k++) 
	fprintf(ofp, "(%.4f %.4f) ", M->real[i*Lk+k], M->imag[i*Lk+k]);
      fprintf(ofp, "\n");
    }

    fprintf(ofp, "N matrix\n");
    
    for (k = 0; k < Lk; k++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "(%.4f %.4f) ", N->real[k*Lj+j], N->imag[k*Lj+j]);
      fprintf(ofp, "\n");
    }
  }

  for (i = 0; i < Li; i++)
    for (j = 0; j < Lj; j++) {
      prod_real = 0.;
      prod_imag = 0.;
      for (k = 0; k < Lk; k++) {
	prod_real += M->real[i*Lk+k] * N->real[k*Lj+j] - M->imag[i*Lk+k] * N->imag[k*Lj+j] ;
	prod_imag += M->real[i*Lk+k] * N->imag[k*Lj+j] + M->imag[i*Lk+k] * N->real[k*Lj+j] ;

	if (prod_real > DBL_MAX || prod_imag > DBL_MAX) 
	  Die("sorry, I am getting out of bounds in Cal_M_N_Prod() prod[%d][%d] = (%f %f)\n", i, j, prod_real, prod_imag);
      }
      Mprod->real[i*Lj+j] = prod_real;
      Mprod->imag[i*Lj+j] = prod_imag;
    }

  if (verbose) {
    fprintf(ofp, "M*N matrix\n");
    
    for (i = 0; i < Li; i++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "(%.4f %.4f) ", Mprod->real[i*Lj+j], Mprod->imag[i*Lj+j]);
      fprintf(ofp, "\n");
    }
  }
  return Mprod;
}

void
Comp_M_N_Prod(FILE *ofp, double *M, double *N, int L, int verbose)
{
  double *Mprod;

  Mprod = Cal_M_N_Prod(ofp, M, N, L, L, L, verbose);

  CopyMatrix (N, Mprod, L, L);

  if (verbose) {
    fprintf(ofp, "M*N matrix\n");
    PrintMatrix(ofp, N, L);
  }

  free(Mprod);
}

void
Comp_M_N_ProdGen(FILE *ofp, double *M, double *N, int Li, int Lk, int Lj, int verbose)
{
  double *Mprod;
  int     i;

  Mprod = Cal_M_N_Prod(ofp, M, N, Li, Lk, Lj, verbose);

  CopyMatrix (N, Mprod, Li, Lj);

  if (verbose) {
    fprintf(ofp, "M*N matrix\n");
    for (i = 0; i < Li; i++) 
      PrintVectorProbs(stdout, Mprod+i*Lj, Lj);
  }

  free(Mprod);
}

void
ZComp_M_N_Prod(FILE *ofp, struct zmatrix_s *M, struct zmatrix_s *N, int L, int verbose)
{
  struct zmatrix_s *Mprod;

  Mprod = ZCal_M_N_Prod(ofp, M, N, L, L, L, verbose);

  CopyZmatrix (N, Mprod, L, L);

  if (verbose) {
    fprintf(ofp, "M*N matrix\n");
    PrintZmatrix(ofp, N, L);
  }

  FreeZmatrix(Mprod);
}

/* Function: Cal_V_M_Prod()
 * Date:     ER, Mon Sep 15 15:19:36 CDT 2003 [St. Louis]
 *
 * Purpose:  given V (L) and M (LxL) calculate v_i M_ij
 *
 * Args:     V  - L (prealloc)
 *           M  - LxL (prealloc)
 *
 * Returns:  Vprod(L) = V(L) * M(LxL)
 *           Msum is alocated here, freed by caller.
 */
double *
Cal_V_M_Prod(FILE *ofp, double *vec, double *M, int L, int verbose)
{
  int     i, j;
  double *Vprod;
  
  Vprod = (double *) MallocOrDie (sizeof(double) * L);
  
  for (i = 0; i < L; i++)
    for (j = 0; j < L; j++) 
	Vprod[j] = vec[i] * M[i*L+j];
  
  if (verbose) {
    fprintf(ofp, "vM matrix\n");
    
    for (i = 0; i < L; i++) 
 	fprintf(ofp, "%.4f ", Vprod[i]);
      fprintf(ofp, "\n");    
  }

  return Vprod;
}

/* Function: Cal_M_N_Sum()
 * Date:     ER, Mon May  8 10:44:33 CDT 2000 [St. Louis]
 *
 * Purpose:  given M (LixLj) and N (LixLj) calculate Msum(LixLj) = M + N.
 *
 * Args:     M  - LixLj (prealloc)
 *           N  - LixLj (prealloc)
 *
 * Returns:  Msum(LixLj) = M(LixLj) * N(LixLj)
 *           Msum is alocated here, freed by caller.
 */
double *
Cal_M_N_Sum(FILE *ofp, double *M, double *N, int Li, int Lj, int add, int verbose)
{
  int     i, j;
  double *Msum;
  
  Msum = (double *) MallocOrDie (sizeof(double) * Li * Lj);
  
  for (i = 0; i < Li; i++)
    for (j = 0; j < Lj; j++) 
      if (add)
	Msum[i*Lj+j] = M[i*Lj+j] + N[i*Lj+j];
      else
	Msum[i*Lj+j] = M[i*Lj+j] - N[i*Lj+j];
  
  if (verbose) {
    fprintf(ofp, "M+N matrix\n");
    
    for (i = 0; i < Li; i++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "%.4f ", Msum[i*Lj+j]);
      fprintf(ofp, "\n");
    }
  }
  return Msum;
}

struct zmatrix_s *
ZCal_M_N_Sum(FILE *ofp, struct zmatrix_s *M, struct zmatrix_s *N, int Li, int Lj, int add, int verbose)
{
  int               i, j;
  struct zmatrix_s *Msum;
  
  Msum = AllocZmatrix(Li*Lj);
  
  for (i = 0; i < Li; i++)
    for (j = 0; j < Lj; j++) 
      if (add) {
	Msum->real[i*Lj+j] = M->real[i*Lj+j] + N->real[i*Lj+j];
	Msum->imag[i*Lj+j] = M->imag[i*Lj+j] + N->imag[i*Lj+j];
      }
      else {
	Msum->real[i*Lj+j] = M->real[i*Lj+j] - N->real[i*Lj+j];
	Msum->imag[i*Lj+j] = M->imag[i*Lj+j] - N->imag[i*Lj+j];
      }
  
  if (verbose) {
    fprintf(ofp, "M+N matrix\n");
    
    for (i = 0; i < Li; i++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "(%.4f %.4f) ", Msum->real[i*Lj+j], Msum->imag[i*Lj+j]);
      fprintf(ofp, "\n");
    }
  }

  return Msum;
}

void
Comp_M_N_Sum(FILE *ofp, double *M, double *N, int Li, int Lj, int add, int verbose)
{
  int     i, j;
  double *Msum;
  
  Msum = Cal_M_N_Sum(ofp, M, N, Li, Lj, add, verbose);
  
  CopyMatrix (N, Msum, Li, Lj);

  if (verbose) {
    fprintf(ofp, "M + N matrix\n");
    for (i = 0; i < Li; i++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "%.4f ", N[i*Lj+j]);
      fprintf(ofp, "\n");
    }
  }

  free(Msum);
}

void
ZComp_M_N_Sum(FILE *ofp, struct zmatrix_s *M, struct zmatrix_s *N, int Li, int Lj, int add, int verbose)
{
  int               i, j;
  struct zmatrix_s *Msum;
  
  Msum = ZCal_M_N_Sum(ofp, M, N, Li, Lj, add, verbose);

  CopyZmatrix (N, Msum, Li, Lj);

  if (verbose) {
    fprintf(ofp, "M+N matrix\n");
    
    for (i = 0; i < Li; i++) {
      for (j = 0; j < Lj; j++) 
	fprintf(ofp, "(%.4f %.4f) ", Msum->real[i*Lj+j], Msum->imag[i*Lj+j]);
      fprintf(ofp, "\n");
    }
  }

  FreeZmatrix(Msum);
}

double *
Col_From_M (double *A, int Lrow, int Lcol, int idx)
{
  double *col;
  int     row_i;

  if (idx >= Lcol) Die ("cannot give you that col");

  col = (double *) MallocOrDie (sizeof(double) * Lcol);
 
  for (row_i=0; row_i< Lrow; row_i++) col[row_i] = A[row_i*Lcol+idx];

  return col;
}

void
CopyMatrix (double *copyQ, double *Q, int N, int M)
{
  int col, row;

  for (row = 0; row < N; row++) 
    for (col = 0; col < M; col++) 
      copyQ[row*M+col] = Q[row*M+col];

}

void
CopyZmatrix (struct zmatrix_s *copyQ, struct zmatrix_s *Q, int N, int M)
{
  int col, row;

  for (row = 0; row < N; row++) 
    for (col = 0; col < M; col++) {
      copyQ->real[row*M+col] = Q->real[row*M+col];
      copyQ->imag[row*M+col] = Q->imag[row*M+col];
    }

}

void
CopyVector (double *copyVec, double *Vec, int N)
{
  int row;

  for (row = 0; row < N; row++) 
    copyVec[row] = Vec[row];

}

void
FreeEigenvalues(struct eigenvalue_s *eigen)
{
  free(eigen->real);
  free(eigen->imag);
  free(eigen);
}

void
FreeZmatrix(struct zmatrix_s *zmatrix)
{
  free(zmatrix->real);
  free(zmatrix->imag);
  free(zmatrix);
}

int
IsDiagonal (double *A, int dim)
{
  int flag = 1;
  int row, col;
  
  for (row = 0; row < dim; row++) 
    for (col = 0; col < dim; col++) {
      if (col != row && (A[col*dim+row] > MARGIN || A[col*dim+row] < - MARGIN)) flag = 0; break;
    }

  return flag;
}

/* Function: HessenbergForm()
 *
 * Date:    ER, Fri Apr 21 16:15:23 CDT 2000 [St. Louis]
 *
 * Purpose: Given a real matrix M (LxL) obtain its Hessenberg form.
 *          This is an intermediate step to calculate the eigenvalues
 *          using the QL algorithm for real non-symmetric matrices.
 *
 *          The Hessenberg form has zeros everywhere below the diagonal,
 *          expect for the subdiagonal.
 *
 *
 *          Implemented after reading chapter 11 of "numerical recipies in C"
 *
 * Method:  -- pick a column k (k=0, k<L-1)
 *
 *          -- look at rows  i (i>k, i<L)
 *
 *                 find row i_o with larger value M(i_o,k)
 *
 *          -- exchange rows i_o <--> k+1
 *          -- exchange cols i_o <--> k+1
 *
 *          -- for rows  i (i=k+2, i<L)
 *
 *                 row_i <-- row_i - row_{k+1}*M(i,k)/M(k+1,k)
 *
 *             notice that for the new row_i:
 *                 M(i,k) = 0 so we have put zeros in column k, for all rows under the subdiagonal (i > k+1)
 *
 *         -- to make the elimination a similarity transformation, also reassign:
 *
 *               for rows  i (i=k+2, i<L)
 *                 col_{k+1} <-- col_{k+1} - col_{i}*M(i,k)/M(k+1,k)
 *
 * Args:              
 *
 * Returns:  H hessenberg form.
 *           H is allocated here, and freed by caller.
 *
 */
double *
HessenbergForm(FILE *ofp, double *M, int L, int verbose)
{
  double *H;
  double  val;
  double  exchange;
  int     i, j, k, k1;
  int     new_row;


  H  = (double *) MallocOrDie (sizeof(double) * L * L);

  for (i = 0; i < L; i++) 
    for (j = 0; j < L; j++) 
      H[i*L+j] = M[i*L+j];

  /* start moving by columns
   */
  for (k = 0; k < L-1; k++) {

    k1 = k + 1;
    val = 0.;
    /* For a given column move rows to have the largest possible value
     * in the diagonal.
     */
    new_row = k1; /* initialize */
    for (i = k1; i < L; i++) {
      if (fabs(H[i*L+k]) > val) {
	val = fabs(H[i*L+k]);
	new_row = i;
      }
    }
    
    if (k1 != new_row) {
      for (i = 0; i < L; i++) {
	/* exchange values of the two rows (k+1 and new_row) 
	 */
	exchange       = H[k1*L+i];
	H[k1*L+i]      = H[new_row*L+i];
	H[new_row*L+i] = exchange;
      }
      for (i = 0; i < L; i++) {
	/* also exchange values for the columns (k+1 and new_row)
	 */
	exchange       = H[i*L+k1];
	H[i*L+k1]      = H[i*L+new_row];
	H[i*L+new_row] = exchange;
      }
    }
    
    if (val != 0.) {
      for (i = k1+1; i < L; i++) {
	for (j = 0; j < L; j++) 
	  H[i*L+j]  -= H[i*L+k] * H[k1*L+j] / val;
	for (j = 0; j < L; j++) 
	  H[j*L+k1] += H[i*L+k] * H[j*L+i] / val;
     }
    }
  } /* for every column k */
  
 if (verbose) {
    fprintf(ofp, "Hessenberg\n");
    PrintProbs(ofp, H, L);
 }

  return H;
}

/* Function: Hessenberg2Eigenvalues()
 *
 * Date:    ER, Fri Apr 21 21:07:27 CDT 2000 [St. Louis]
 *
 * Purpose: Given a real non-symmetric matrix M (LxL) in its Hessenberg form
 *          calculate its eigenvalues using the QR algorithm.
 *
 *          Implemented after reading chapter 6 of Matrix Methods (R. Bronson).
 *
 *          A_o = H
 *
 *          Use QR decomposition to calculate    A_o - (A_o)_nn I = Q_o * R_o
 *
 *          then,                                A_1              = R_o * Q_o + (A_o)_nn I
 *
 *          The QR decomposition preserves the Hessenberg form (which makes it faster).
 *
 *          At the end A_k is of the form       |   S      T |
 *                                              |            |
 *                                              | 0...0    a |    ===> a is an eigenvalue. Continue QR with submatrix S
 *
 *          OR
 *                                              |   S       T  |
 *                                              |              |
 *                                              | 0...0    b c |
 *                                              | 0...0    d e | ===> 2 complex eigenvalues. Continue QR with submatrix S
 *
 *
 * Args:    H -- (LxL) Hessenberg matrix          
 *
 * Returns: eigen(2L) vector with real and imaginary part of eigenvalues
 *
 */
void
Hessenberg2Eigenvalues(FILE *ofp, double *H, int L, struct eigenvalue_s **ret_eigen, int pedantic, int verbose)
{
  struct eigenvalue_s *eigen;
  double *A;         /* the series matrix that is equivalent to the original matrix */
  double *Q;
  double *R;
  double *I;
  double *last_row;
  double *nxtl_row;
  double  Ann;       /* nn element of the A matrix                                  */
  double  a, b, c;   /* coefficients of second degree equation ax^2+bx+c=0          */
  int     dim;       /* dimension of the square matrix A                            */
  int     it = 0;    /* number of iterations                                        */
  int     idx;       /* order of the eigenvalues being calculated                   */
  int     i,j;
  int     flag = 0;
  
  /* memory allocation */
  A =  (double *) MallocOrDie (sizeof(double) * L * L);
  Q =  (double *) MallocOrDie (sizeof(double) * L * L);
  R =  (double *) MallocOrDie (sizeof(double) * L * L);
  eigen = AllocEigenvalues (L);

  /* initialize A_o matrix Q and R*/
  for (i = 0; i < L*L; i++) {
    A[i] = H[i];
    Q[i] = 0.0;
    R[i] = 0.0;
  }
  
  /* initialize dimension of space */
  dim = L;
  
  /* initialize number of eigenvalues */
  idx = 0;
  
  /* do the iteration */
  
  while (dim > 2) {
    flag = 0;
    it ++;
    if (verbose) 
      printf("Hessenberg2Eigenvalues()\nnumber of iterations = %d\ndimension = %d\n", it, dim);
    
    last_row = (double *)Row_From_M (A, dim, dim, dim-1);
    nxtl_row = (double *)Row_From_M (A, dim, dim, dim-2);
    
    for (i=0; i<dim-1; i++) 
      if (last_row[i] > MARGIN || last_row[i] < -MARGIN) { flag = 1; break; }
    
    if (flag == 0) { /* one real eigenvalue */    
      eigen->real[idx] = last_row[dim-1];
      eigen->imag[idx] = 0.0;

      if (verbose) printf("  [%d]->real eigenvalue %f\n", idx, eigen->real[idx]);

      idx ++; /* one more eigenvalue */    

      /* reduce matrix */
      for (i = 0; i < dim-1; i++) 
	for (j = 0; j < dim-1; j++) 
	  A[i*(dim-1)+j] = A[i*dim+j];
      
      dim --; /* reduce the dimension of the matrix by one */
    }
    else 
      {
	flag = 0;
	for (i=0; i<dim-2; i++) 
	  if (last_row[i] > MARGIN || last_row[i] < -MARGIN || 
	      nxtl_row[i] > MARGIN || nxtl_row[i] < -MARGIN  ) { flag = 1; break; }

	if (flag == 0) { /* two (posibly complex)  eigenvalues */   
	  a = 1.0;
	  b = - nxtl_row[dim-2] - last_row[dim-1];
	  c = nxtl_row[dim-2]*last_row[dim-1] - nxtl_row[dim-1]*last_row[dim-2];
	  
	  SecondDegree_Solutions (a, b, c, eigen, idx);
	  
	  if (verbose) {
	    printf("  [%d]->pair eigenvalue (%f, %f) \n", idx,   eigen->real[idx],   eigen->imag[idx]);
	    printf("  [%d]->pair eigenvalue (%f, %f) \n", idx+1, eigen->real[idx+1], eigen->imag[idx+1]);
	  }

	  idx += 2; /* two more eigenvalues */

	  /* reduce matrix */
	  for (i = 0; i < dim-2; i++) 
	    for (j = 0; j < dim-2; j++) 
	      A[i*(dim-2)+j] = A[i*dim+j];
	  
	  dim -= 2; /* reduce the dimension of the matrix by 2 */
	  
	}
	else { /* ok, do the actual QR decomposition */ 

	  /* shift matrix */
	  Ann = -A[(dim-1)*dim+(dim-1)];
	  I   = Cal_Id (dim);
	  Comp_M_num_Prod (ofp, I, Ann, dim, dim, verbose);
	  Comp_M_N_Sum (ofp, I, A, dim, dim, TRUE, verbose);    /* substract the shift              */
	  
	  QR_Decomposition (ofp, A, Q, R, dim, pedantic, verbose);  /* QR decomposition of A      */
	  
	  Comp_M_N_Prod (ofp, R, Q, dim, verbose);               /* calculate new A            */
	  CopyMatrix(A, Q, dim, dim);
	  Comp_M_N_Sum (ofp, A, I, dim, dim, FALSE, verbose);    /* add the shift              */
	  CopyMatrix(A, I, dim, dim);

	}
      }
    
    free(last_row);
    free(nxtl_row);
    
  }

  if (dim == 2) {

    a = 1.0;
    b = - A[0] - A[3];
    c = A[0]*A[3] - A[1]*A[2];

    SecondDegree_Solutions (a, b, c, eigen, idx);

    if (verbose) {
      printf("  [%d]->pair eigenvalue (%f, %f) \n", idx,   eigen->real[idx],   eigen->imag[idx]);
      printf("  [%d]->pair eigenvalue (%f, %f) \n", idx+1, eigen->real[idx+1], eigen->imag[idx+1]);
    }
    idx += 2;  /* two eigenvalues */
  }
  else if (dim == 1) {
    eigen->real[idx] = A[0];
    eigen->imag[idx] = 0.0;

    if (verbose) printf("  [%d]->real eigenvalue %f\n", idx, eigen->real[idx]);
    idx ++; /* one eigenvalue */
  }

  /* paranoia */
  if (idx != L) Die ("You have not calculated all the eigenvalues\n");

  /* printout eigenvalues */
  if (verbose) 
    {
      fprintf(ofp, "Eigenvalues for Q\n"); 
      for (i = 0; i < L; i++)  
	fprintf(ofp, "%.4f %.4f \n", eigen->real[i], eigen->imag[i]);
    }
 
  free (A);
  free (Q);
  free (R);
  free (I);

  *ret_eigen = eigen;
}

void
MultiplyMatrix (double tfactor, int L, double *K)
{
  int i,j;
  
  for (i = 0; i < L; i ++)
    for (j = 0; j < L; j ++)
      K[i*L+j] *= tfactor;

}

void
MultiplyZmatrix (double real, double imag, int L, struct zmatrix_s *K)
{
  double Kreal;
  double Kimag;
  int    i,j;
  
  for (i = 0; i < L; i ++)
    for (j = 0; j < L; j ++) {
      Kreal = K->real[i*L+j] * real - K->imag[i*L+j] * imag;
      Kimag = K->real[i*L+j] * imag + K->imag[i*L+j] * real;

      K->real[i*L+j] = Kreal;
      K->imag[i*L+j] = Kimag;
    }

}

void
PrintMatrix(FILE *ofp, double *m, int L) 
{
  int x, y;

  for (x = 0; x < L; x++) {
    for (y = 0; y < L; y++) {
      fprintf(ofp, "%.4f ", m[x*L+y]);
    }
    fprintf(ofp, "\n");
  }
}

void
PrintZmatrix(FILE *ofp, struct zmatrix_s *m, int L) 
{
  int x, y;

  for (x = 0; x < L; x++) {
    for (y = 0; y < L; y++) {
      fprintf(ofp, "(%.4f %.4f) ", m->real[x*L+y], m->imag[x*L+y]);
    }
    fprintf(ofp, "\n");
  }
}

void
PrintZvector(FILE *ofp, struct zmatrix_s *v, int L) 
{
  int x;

  for (x = 0; x < L; x++) 
    fprintf(ofp, "(%.4f %.4f) ", v->real[x], v->imag[x]);
  
  fprintf(ofp, "\n");
  
}

/* Function: QR_Decomposition()
 *
 * Date:    ER, Thu Jul 25 14:35:13 CDT 2002  [St. Louis]
 *
 * Purpose: Given a real non-symmetric matrix M (nxn) 
 *          calculate the QR decomposition.
 *
 *          Implemented after reading chapter 6 of Matrix Methods (R. Bronson).
 *
 *          X = [x_1, ..., x_n]
 *
 *          for each i [1,...,n]
 *                          
 *              - for each j [i,...,n]
 *                          r_ij = <x_i,x_j>  ---->   r_i = (0,...,0,r_ii,...,r_in)
 *                          
 *              - q_i = 1/r_ii x_i
 *                          
 *              - for each j [i+1,...,n]
 *                         x_j = x_j - r_ij q_i
 *
 *         Then define Q = [q_n,...,q_n]
 *           
 *                         | r_1 |
 *                         |  .  |
 *                     R = |  .  |                      and X = QR
 *                         |  .  |
 *                         | r_n |
 *
 *         Then define new X = RQ
 *
 *
 * Args:    X -- (LxL) Hessenberg matrix     
 *          Q
 *          R
 *          L - dimension
 *
 * Returns: void. Calculates the new X matrix
 *
 */
void 
QR_Decomposition (FILE *ofp, double *X, double *Q, double *R, int n, int pedantic, int verbose)
{
  double *Check;
  double *Xdup;
  int     i, j, k;
  
  Xdup = (double *) MallocOrDie (sizeof(double) * n * n);

  /* initialize matrices Q and R*/
  for (i = 0; i < n*n; i++) {
    Q[i] = 0.0;
    R[i] = 0.0;
  }

  for (i = 0; i < n*n; i++) Xdup[i] = X[i];
  
  for (i=0; i<n; i++)
    {
      /* 1. calculate Rii = sqrt <x_i,x_i>
       */
      for (k=0; k<n; k++) R[i*n+i] += X[k*n+i] * X[k*n+i];
      R[i*n+i] = sqrt(R[i*n+i]);
      
      /* 2. calculate q_i = 1/Rii x_i
       */
      for (k=0; k<n; k++)   if (R[i*n+i] != 0.0) Q[k*n+i] = X[k*n+i] / R[i*n+i];

      /* 3. calculate Rij = <x_j,q_i>
       */
      for (j=i+1; j<n; j++) 
	for (k=0; k<n; k++) R[i*n+j] += X[k*n+j] * Q[k*n+i];

      /* 4. redefinition vector x_j by x_j - Rij * q_i
       */
      for (j=i+1; j<n; j++) 
	for (k=0; k<n; k++) X[k*n+j] -= R[i*n+j] * Q[k*n+i];
    }
  
  /* check is X = QR ? */
  Check = Cal_M_N_Prod(ofp, Q, R, n, n, n, FALSE);
  
  if (verbose) {    
    fprintf(ofp, "QR decomposition\n");
    fprintf(ofp, "Q Matrix\n");
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	fprintf(ofp, "%.4f ", Q[i*n+j]);
      }
      fprintf(ofp, "\n");
    }
    
    fprintf(ofp, "R Matrix\n");
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	fprintf(ofp, "%.4f ", R[i*n+j]);
      }
      fprintf(ofp, "\n");
    }

    fprintf(ofp, "X Matrix\n");
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	fprintf(ofp, "%.4f ", Xdup[i*n+j]);
      }
      fprintf(ofp, "\n");
    }

    fprintf(ofp, "\nQR Matrix\n");
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	fprintf(ofp, "%.4f ", Check[i*n+j]);
      }
      fprintf(ofp, "\n");
    }
    fprintf(ofp, "\n");
    
  }
    
  for (i = 0; i < n*n; i++) 
    if (abs(Check[i]) > 0.0) {
      if (Check[i]/Xdup[i] < accuracy || Check[i]/Xdup[i] > 2.0-accuracy) Die ("Bad QR decomposition %d %f %f", i, Xdup[i], Check[i]);
    }
    else
      if (Check[i] - Xdup[i] > MARGIN || Check[i] - Xdup[i] < -MARGIN) Die ("Bad QR decomposition %d %f %f", i, Xdup[i], Check[i]);
  

  free (Check);
  free (Xdup);
  
}

double *
Row_From_M (double *A, int Lrow, int Lcol, int idx)
{
  double *row;
  int     col_i;

  if (idx >= Lcol) Die ("cannot give you that row");

  row = (double *) MallocOrDie (sizeof(double) * Lrow);
 
  for (col_i=0; col_i< Lcol; col_i++) row[col_i] = A[idx*Lcol+col_i];

  return row;
}

  
  
/* Function: SecondDegree_Solutions()
 * Date:     ER, Thu Jul 25 13:22:05 CDT 2002 [St. Louis]
 *
 * Purpose:  calculate the 2 solutions of equation ax^2 + bx + c = 0
 *         
 * Args:     a
 *           b
 *           c
 *           eigen
 *
 * Returns: void
 *           
 */
void
SecondDegree_Solutions(double a, double b, double c, struct eigenvalue_s *eigen, int idx)
{
  double discriminant;
  double real;
  double imag;

  real = 0.0;
  imag = 0.0;
  discriminant = b*b - 4.0*a*c;

  if      (discriminant >= 0) real += sqrt(+discriminant);
  else if (discriminant <  0) imag += sqrt(-discriminant);

  b    /= 2.0 * a;
  real /= 2.0 * a;
  imag /= 2.0 * a;

  eigen->real[idx]   = -b + real;
  eigen->real[idx+1] = -b - real;

  eigen->imag[idx]   = +imag;
  eigen->imag[idx+1] = -imag;

}



/* Function: Solve_LUP_Eq()
 * Date:     ER, Mon Apr 17 17:22:31 CDT 2000 [St. Louis]
 *
 * Purpose:  Given LU x = b solve for x
 *  
 *           (Implemented after reading chapter 31 of "introduction to algorithms"
 *            by Cormen et al.)
 *
 * Method:   (1) solve equation  Ly = b for y
 *                     trivial because L is a lower triangular matrix.
 *
 *                     y_0 = b_0/l_00 = b_0
 *
 *                     y_i = 1/l_ii [ b_i - \sum_{j=1}^{i-1} l_ij*y_j ] for i=2,...,L-1
 *
 *                         = b_i - \sum_{j=1}^{i-1} l_ij*y_j 
 *
 *           (2) solve equation  Ux = y for x
 *                     trivial because U is an upper triangular matrix.
 *
 *                     x_{L-1} = y_{L-1}/u_{L-1,L-1}
 *
 *                     x_i = 1/u_ii [ y_i - \sum_{j=i+1}^{L-1} u_ij*x_j ] for i=L-2,...,0
 *
 *
 * Args:     LU - LxL LU  matrix.
 *
 * Returns:  (void)
 *           values for x are dumped in b.
 */
void
Solve_LU_Eq(FILE *ofp, double *LU, double *b, int L, int verbose)
{
  int     i, j;
  double  val, add;

  if (verbose) 
    PrintVectorProbs(ofp, b, L);

  /* calculate y (dump values into b)
   */
  for (i = 1; i < L; i++) {
    val = b[i];
    add = 0.0;
      for (j = 0; j < i; j++)
	{
	  add -= LU[i*L+j] * b[j];
	}
    val += add;
    b[i] = val; 
  }

  /* calculate x (dump values into b)
   */
 b[L-1] /= LU[(L-1)*L+L-1];

  for (i = L-2; i >= 0; i--) {
    val = b[i];
    add = 0.0;
      for (j = i+1; j < L; j++)
	{
	  add -= LU[i*L+j] * b[j];
	}
    val += add;
    b[i] = val / LU[i*L+i];
  }

  if (verbose) 
    PrintVectorProbs(ofp, b, L);

}

void
ZSolve_LU_Eq(FILE *ofp, struct zmatrix_s *LU, struct zmatrix_s *b, int L, int verbose)
{
  int     i, j;
  double  val_real;
  double  val_imag;
  double  add_real;
  double  add_imag;

  if (verbose) {
    fprintf(ofp, "\nInput b (LUx=b)\n");
    PrintZvector(ofp, b, L);
  }

  /* calculate y (dump values into b)
   */
  for (i = 1; i < L; i++) {
    val_real = b->real[i];
    val_imag = b->imag[i];
    add_real = 0.0;
    add_imag = 0.0;
      for (j = 0; j < i; j++)
	{
	  add_real -= LU->real[i*L+j] * b->real[j] - LU->imag[i*L+j] * b->imag[j];
	  add_imag -= LU->real[i*L+j] * b->imag[j] + LU->imag[i*L+j] * b->real[j];
	}
    val_real += add_real;
    val_imag += add_imag;
    b->real[i] = val_real; 
    b->imag[i] = val_imag; 
  }
  if (verbose) {
    fprintf(ofp, "\nPartial sol y (Ly=b)\n");
    PrintZvector(ofp, b, L);
  }


  /* calculate x (dump values into b)
   */
  printf("LU %f %f\n", LU->real[(L-1)*L+L-1], LU->imag[(L-1)*L+L-1]);
  val_real = (b->real[L-1] * LU->real[(L-1)*L+L-1] + b->imag[L-1] * LU->imag[(L-1)*L+L-1]) / 
    (LU->real[(L-1)*L+L-1] * LU->real[(L-1)*L+L-1] + LU->imag[(L-1)*L+L-1] * LU->imag[(L-1)*L+L-1]);
  val_imag = (b->imag[L-1] * LU->real[(L-1)*L+L-1] - b->real[L-1] * LU->imag[(L-1)*L+L-1]) / 
    (LU->real[(L-1)*L+L-1] * LU->real[(L-1)*L+L-1] + LU->imag[(L-1)*L+L-1] * LU->imag[(L-1)*L+L-1]);

  b->real[L-1] = val_real; 
  b->imag[L-1] = val_imag; 
    
  for (i = L-2; i >= 0; i--) {
    val_real = b->real[i];
    val_imag = b->imag[i];
    add_real = 0.0;
    add_imag = 0.0;
      for (j = i+1; j < L; j++)
	{
	  add_real -= LU->real[i*L+j] * b->real[j] - LU->imag[i*L+j] * b->imag[j];
	  add_imag -= LU->real[i*L+j] * b->imag[j] + LU->imag[i*L+j] * b->real[j];
	}
    val_real += add_real;
    val_imag += add_imag;

    b->real[i] = (val_real *  LU->real[i*L+i] + val_imag *  LU->imag[i*L+i]) / 
      (LU->real[i*L+i] * LU->real[i*L+i] + LU->imag[i*L+i] * LU->imag[i*L+i]);
    b->imag[i] = (val_imag *  LU->real[i*L+i] - val_real *  LU->imag[i*L+i]) / 
      (LU->real[i*L+i] * LU->real[i*L+i] + LU->imag[i*L+i] * LU->imag[i*L+i]);
  }

  if (verbose) {
    fprintf(ofp, "\nSol x (Ux=y)\n");
    PrintZvector(ofp, b, L);
  }

}

