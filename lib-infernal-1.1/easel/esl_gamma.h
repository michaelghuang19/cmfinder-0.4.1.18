/* Gamma distributions.
 * 
 * SRE, Wed Nov 16 19:15:33 2005 [St. Louis]
 * SVN $Id: esl_gamma.h 664 2011-02-27 17:08:36Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/infernal/1.1/esl_gamma.h $
 */
#ifndef eslGAMMA_INCLUDED
#define eslGAMMA_INCLUDED

#ifdef eslAUGMENT_RANDOM
#include <esl_random.h>
#endif

double esl_gam_pdf    (double x, double mu, double lambda, double tau);
double esl_gam_logpdf (double x, double mu, double lambda, double tau);
double esl_gam_cdf    (double x, double mu, double lambda, double tau);
double esl_gam_logcdf (double x, double mu, double lambda, double tau);
double esl_gam_surv   (double x, double mu, double lambda, double tau);
double esl_gam_logsurv(double x, double mu, double lambda, double tau);
double esl_gam_invcdf (double p, double mu, double lambda, double tau);

double esl_gam_generic_pdf   (double x, void *params);
double esl_gam_generic_cdf   (double x, void *params);
double esl_gam_generic_surv  (double x, void *params);
double esl_gam_generic_invcdf(double x, void *params);

extern int esl_gam_Plot(FILE *fp, double mu, double lambda, double tau,
			double (*func)(double x, double mu, double lambda, double tau), 
			double xmin, double xmax, double xstep);

#ifdef eslAUGMENT_RANDOM
extern double esl_gam_Sample(ESL_RANDOMNESS *r, double mu, double lambda, double tau);
#endif

extern int esl_gam_FitComplete(double *x, int n, double mu, double *ret_lambda, double *ret_tau);

#endif /*eslGAMMA_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version i1.1rc2; December 2012
 * Copyright (C) 2012 HHMI Janelia Farm Research Campus
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
