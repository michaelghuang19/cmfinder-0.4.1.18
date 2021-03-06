/* ratematrix.h
 * Header for evolutionary rate matrix routines in ratemx.c
 * 
 * SRE, Tue Jul 13 16:09:05 2004 [St. Louos]
 * SVN $Id: ratematrix.h 11 2005-01-06 11:44:17Z eddy $
 */
#ifndef ESL_RATEMATRIX_INCLUDED
#define ESL_RATEMATRIX_INCLUDED

extern int esl_ratemx_Symm2Q(ESL_DMATRIX *s, double *pi, ESL_DMATRIX *Q);
extern int esl_ratemx_Normalize(ESL_DMATRIX *Q, double *pi, double x);
extern int esl_ratemx_TaylorExp(ESL_DMATRIX *Q, double t, ESL_DMATRIX *P);

extern ESL_DMATRIX *esl_ratemx_CreateHKY(double *f, double alpha, double beta);

#endif /*ESL_RATEMATRIX_INCLUDED*/
