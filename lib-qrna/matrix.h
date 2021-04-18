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

/* matrix.h
 * ANSI prototypes for all external functions.
 * 
 */

#include <stdio.h>
#include "mat_structs.h"

/* from matrix.c: functions for matrix manipulations
 */
extern struct eigenvalue_s *AllocEigenvalues(int L);
extern struct zmatrix_s    *AllocZmatrix(int L);
extern double              *Cal_Id(int L);
extern double              *Cal_LUP_Decomp(FILE *ofp, double *M, int *perm, int L, int verbose);
extern struct eigenvalue_s *Cal_M_Eigenvalues(FILE *ofp, double *M, int L, int verbose);
extern double              *Cal_M_Exp(FILE *ofp, double *M, int L, double power, int verbose);
extern double              *Cal_M_Inv(FILE *ofp, double *M, int L, int verbose);
extern struct zmatrix_s    *Cal_M_Vandermonde_Log(FILE *ofp, double *M, int L, int verbose);
extern double              *Cal_M_Taylor_Log(FILE *ofp, double *M, int L, int verbose);
extern double              *Cal_M_MinId(double *M, int L);
extern double              *Cal_M_num_Prod(FILE *ofp, double *M, double num, int Li, int Lj, int verbose);
extern double              *Cal_V_M_Prod(FILE *ofp, double *vec, double *M, int L, int verbose);
extern double              *Cal_M_N_Prod(FILE *ofp, double *M, double *N, int Li, int Lj, int Lk, int verbose);
extern double              *Cal_M_N_Sum(FILE *ofp, double *M, double *N, int Li, int Lj, int add, int verbose);
extern double              *Col_From_M(double *A, int Lrow, int Lcol, int col_i);
extern void                 Comp_Id(double *Id, int L);
extern void                 Comp_M_Exp(FILE *ofp, double *M, int L, double power, int verbose);
extern void                 Comp_M_Inv(FILE *ofp, double *M, double *Minv, int L, int verbose);
extern void                 Comp_M_Taylor_Log(FILE *ofp, double *M, int L, int verbose);
extern void                 Comp_M_MinId(double *M, int L);
extern void                 Comp_M_num_Prod(FILE *ofp, double *M, double num, int Li, int Lj, int verbose);
extern void                 Comp_M_N_Prod(FILE *ofp, double *M, double *N, int L, int verbose);
extern void                 Comp_M_N_ProdGen(FILE *ofp, double *M, double *N, int Li, int Lj, int Lk, int verbose);
extern void                 Comp_M_N_Sum(FILE *ofp, double *M, double *N, int Li, int Lj, int add, int verbose);
extern void                 CopyMatrix (double *copyQ, double *Q, int N, int M);
extern void                 CopyVector(double *copyVec, double *Vec, int N);
extern int                  IsDiagonal (double *A, int dim);
extern void                 FreeEigenvalues(struct eigenvalue_s *eigen);
extern void                 FreeZmatrix(struct zmatrix_s *zmatrix);
extern double              *HessenbergForm(FILE *ofp, double *A, int L, int verbose);
extern void                 Hessenberg2Eigenvalues(FILE *ofp, double *H, int L, struct eigenvalue_s **ret_eigen, int pedantic, int verbose);
extern void                 MultiplyMatrix (double tfactor, int L, double *K);
extern void                 PrintMatrix(FILE *ofp, double *m, int L);
extern void                 PrintZmatrix(FILE *ofp, struct zmatrix_s *m, int L);
extern void                 PrintZvector(FILE *ofp, struct zmatrix_s *v, int L) ;
extern void                 QR_Decomposition (FILE *ofp, double *X, double *Q, double *R, int n, int pedantic, int verbose);
extern double              *Row_From_M(double *A, int Lrow, int Lcol, int row_i);
extern void                 SecondDegree_Solutions(double a, double b, double c, struct eigenvalue_s *eigen, int idx);
extern void                 Solve_LU_Eq(FILE *ofp, double *LU, double *b, int L, int verbose);
extern struct zmatrix_s    *ZCal_Id(int L);
extern struct zmatrix_s    *ZCal_LUP_Decomp(FILE *ofp, struct zmatrix_s *M, int *perm, int L, int verbose);
extern struct zmatrix_s    *ZCal_M_Exp(FILE *ofp, struct zmatrix_s *M, int L, double r, int verbose);
extern struct zmatrix_s    *ZCal_M_Inv(FILE *ofp, struct zmatrix_s  *M, int L, int verbose);
extern struct zmatrix_s    *ZCal_M_num_Prod(FILE *ofp, struct zmatrix_s *M, double real, double imag, int Li, int Lj, int verbose);
extern struct zmatrix_s    *ZCal_M_N_Prod(FILE *ofp, struct zmatrix_s *M, struct zmatrix_s *N, int Li, int Lk, int Lj, int verbose);
extern struct zmatrix_s    *ZCal_M_N_Sum(FILE *ofp, struct zmatrix_s *M, struct zmatrix_s *N, int Li, int Lj, int add, int verbose);
extern void                 ZComp_M_N_Prod(FILE *ofp, struct zmatrix_s *M, struct zmatrix_s *N, int L, int verbose);
extern void                 ZComp_M_N_Sum(FILE *ofp, struct zmatrix_s *M, struct zmatrix_s *N, int Li, int Lj, int add, int verbose);
extern void                 CopyZmatrix (struct zmatrix_s *copyQ, struct zmatrix_s *Q, int N, int M);
extern void                 MultiplyZmatrix (double real, double imag, int L, struct zmatrix_s *K);
extern void                 ZSolve_LU_Eq(FILE *ofp, struct zmatrix_s *LU, struct zmatrix_s *b, int L, int verbose);
