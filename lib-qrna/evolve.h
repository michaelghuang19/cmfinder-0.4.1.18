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

extern double                Cal_lambda(FILE *ofp, double *S, double *pml, double *pmr, int L, int verbose);
extern void                  CalculateConditionalsAndMarginals(FILE *ofp, double *P, double **ret_QL, double **ret_QR,
                                                               double **ret_ml, double **ret_mr, int L, int hasindel, int verbose);
extern void                  ComputeConditionalsAndMarginals(FILE *ofp, double *P, double *Q, double *pm, int L, int hasindel, int verbose);
extern void                  ComputeLRConditionalsAndMarginals(FILE *ofp, double *P, double *Ql, double *Qr, double *ml, double *mr,
                                                               int L, int hasindel, int verbose);
extern void                  ChangeFrequencies(FILE *ofp, double *Q_0, double *R, double *pm, int L, double *targetfreq, int hasindel,int verbose);
extern void                  ChangePairProbs(int L5, double *pair5prob, double *psingle_target, int hasindel, int verbose);
extern void                  ChangeStationaryProbs(int L, double *jointprob, double *targetfreq, int hasindel, int verbose);
extern void                  ChangeStationaryProbsGW(int L, double *jointprob, double *targetfreq, int hasindel, int verbose);
extern void                  ChangeStationaryPAMProbsGW(double **pammodel, double *targetfreq, int verbose);
extern void                  ChangePairProbsIterate(FILE *ofp, double *pairprobs, int L, double *targetfreq, int hasindel, int verbose);
extern void                  ChangeConditionalMatrix(FILE *ofp, double *Qstar, double *ml, double *mr, int L, double *targetfreq,
                                                     int hasindel, int isleft, int verbose);
extern void                  ChangeQ_0Matrix(FILE *ofp, double *Q_0, double *ml, double *mr, int L, double *targetfreq,
                                             int hasindel, int isleft, int verbose);
extern void                  ChangeQ_0MatrixNaive(FILE *ofp, double *Q_0, double *QR_0, int L, double *targetfreq,
                                                  int hasindel, int isleft, int verbose);
extern void                  ChangeQ_0MatrixIterate(FILE *ofp, double *Q_0, double *Q_0_inv, double *pml, double *pmr, 
						    int L, double *targetfreq, int hasindel, int isleft, int verbose);
extern void                  ChangeRateMatrix(FILE *ofp, double *Q_0R, double *Q_0_inv, double *R,
                                              double *ml, double *mr, int L, double *targetfreq, int hasindel, int isleft, int verbose);
extern int                   Check_Accuracy(double *vec, int L);
extern int                   CompareFreqs(double *pml, double *pmr, double *targetfreq, int L);
extern double               *Condi_From_Joint(FILE *ofp, double *P, int L, int verbose);
extern double               *CondiR_From_Joint(FILE *ofp, double *P, int L, int verbose);
extern void                  Condi_From_Rate(double **ret_Q, double *R, double *Q_0, double tfactor, int L, int pedantic, int verbose);
extern void                  ConditionalsFromRate(double *Q, double *R, double *Q_0, double tfactor, int L, int pedantic, int verbose);
extern void                  ZCondi_From_Rate(double **ret_Q, struct zmatrix_s *R, struct zmatrix_s *Q_0, double tfactor,
                                              int L, int pedantic, int verbose);
extern void                  ConditionalsEvolved(FILE *ofp, double *Ql, double *Qr, double *Ql_0, double *Qr_0, double *ml, double *mr,
                                                 int L, double tfactor, double *targetfreq,
                                                 int changefreq, int hasindel, int pedantic, int verbose);
extern void                  ConditionalsEvolved_2(FILE *ofp, double *Q, double *Q_0, double *R, int L, double tfactor, double *targetfreq,
                                                   int changefreq, int pedantic, int verbose);
extern void                  EvolveIndelMarginals(FILE *ofp, double *R, double *Q_0, double *pm, double tfactor,
                                                  int L, int pedantic, int verbose);
extern int                   IsQConsistent(FILE *ofp, double *Q, int L, int pedantic, int verbose);
extern void                  IslogQConsistent(FILE *ofp, double *Q_0, double *K, int Lx, int Ly, int Lz, int verbose);
extern void                  IslogQConsistent_2(FILE *ofp, double *Q_0, double *A, double *K, int Lx, int Ly, int Lz, int verbose);
extern void                  Joint_From_Condi(FILE *ofp, double *P, double *Q, double *pm, int L, int verbose);
extern void                  Joint_From_CondiR(FILE *ofp, double *P, double *Q, double *pm, int L, int verbose);
extern void                  Joint_From_Condi_Symmetrize(FILE *ofp, double *P, double *Q, double *p, int L, int verbose);
extern void                  Joint2Joint(double *p_joint, double *ql_nought, double *qr_nought, int L, double tfactor, double *targetfreq,
                                         int changefreq, int hasindel, int ispairprob, int pedantic, int verbose);
extern void                  Joint2JointGap(double *p_joint, double *q_nought, double *pv, int L, double tfactor, double *targetfreq,
                                            int changefreq,  int pedantic, int verbose);
extern void                  MatrixMethodTransitionsFromRate (double time, int Lx, int Ly, double **ret_Q, 
							      double *R, double *Q_zero, double *Q_infty, int verbose);
extern void                  MarginalizeJointProbs(double *jp, double *sp, int dim, int which);
extern void                  OrderVector(int **ret_perm, double *b, int L, int isdescending, int verbose);
extern int                   QOMRegularizationAlgorithm(double *psingle, int size, int verbose);
extern void                  QOGRegularizationAlgorithm(double *b, int size, int whichrow, int verbose);
extern void                  RateFromConditionals(FILE *ofp, double *R, double *Q, double *Q_0, int L, double tfactor,
                                                  int hasindel, int pedantic, int verbose);
extern double               *SaturationProbs(double *R, double *q_nought, int L, int pedantic, int verbose);
extern void                  ComputeSaturationProbs(double *R, double *Q_0, double *pm, int L, int pedantic, int verbose);
extern double                TimeIdCorrelation (struct three_divergence_s div, double id);
extern struct three_times_s  TimeIdCorrelation3Models(struct three_divergence_s div_o, struct three_divergence_s div_c,
                                                      struct three_divergence_s div_r, double id);
extern void                  TransitionsEvolvedLinear(FILE *ofp, double *q, double *q_0, double *r, double *r_diag, double *r_fix,
                                                      int L, double tfactor,
                                                      int pedantic, int verbose);
extern void                  TransitionsEvolved(FILE *ofp, double *q, double *q_0, double *q_infty, double *r_diag, int L, double tfactor,
                                                int pedantic, int verbose);
extern void                  TransitionsDirectoryCosines(double *q, double *q_zero, double *q_infty, double *R_diag, int dim);
extern double               *TransitionsExpVector(double *T, double *T_zero, double *T_infty, int dim);
extern double               *TransitionsDiagVector(double *T, double *T_zero, double *T_rate, double *A_fix,  int dim);
extern double               *TransitionsFixVector (double *T, double *T_zero, double *T_rate, double *A_diag, int dim);
extern double               *TransitionsRateVector(double *T, double *T_zero, double *A_fix,  double *A_diag, int dim);

