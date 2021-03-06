/* esl_mixgev.c
 * Statistical routines for mixtures of generalized extreme value 
 * distributions.
 * 
 * SRE, Mon Aug 15 08:48:19 2005 [St. Louis]
 * xref STL9/139  
 * SVN $Id: esl_mixgev.c 82 2005-12-13 20:21:07Z eddy $
 */

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <easel.h>
#include <esl_vectorops.h>
#include <esl_stats.h>
#include <esl_gev.h>
#include <esl_mixgev.h>

#ifdef eslAUGMENT_RANDOM
#include <esl_dirichlet.h>  /* for uniform sampling of a probability vector */
#include <esl_random.h>
#endif 

#ifdef eslAUGMENT_MINIMIZER
#include <esl_minimizer.h>
#endif


/****************************************************************************
 * Routines for the ESL_MIXGEV object
 ****************************************************************************/ 

/* Function:  esl_mixgev_Create()
 * Incept:    SRE, Mon Aug 15 08:47:37 2005 [St. Louis]
 *
 * Purpose:   Creates an object to hold parameters for a <K>-component
 *            mixture of generalized extreme value distributions.
 *
 *            Parameters in the object are initialized ($q_k =
 *            \frac{1}{K}$, $\mu_k = 0$, $\lambda_k = 1$, $\alpha_k =
 *            0$), but the caller will want to set these according to
 *            its own purposes.
 *
 *            After an object is created, the caller can constrain any
 *            of the components to be a Gumbel (that is, constrain
 *            $\alpha_k = 0$ by calling <esl_mixgev_ForceGumbel(obj,
 *            k)>.
 *
 * Args:      K  - number of components in the mixture
 *
 * Returns:   ptr to newly allocated/initialized <ESL_MIXGEV> object.
 *
 * Throws:    NULL on allocation failure.
 */
ESL_MIXGEV *
esl_mixgev_Create(int K)
{
  ESL_MIXGEV *mg;
  int         k;

  if ((mg = malloc(sizeof(ESL_MIXGEV))) == NULL)           goto FAILURE;
  mg->q = mg->mu = mg->lambda = mg->alpha = mg->wrk = NULL;
  mg->isgumbel = NULL;
  mg->K = K;

  if ((mg->q        = malloc(sizeof(double) * K)) == NULL) goto FAILURE;
  if ((mg->mu       = malloc(sizeof(double) * K)) == NULL) goto FAILURE;
  if ((mg->lambda   = malloc(sizeof(double) * K)) == NULL) goto FAILURE;
  if ((mg->alpha    = malloc(sizeof(double) * K)) == NULL) goto FAILURE;
  if ((mg->wrk      = malloc(sizeof(double) * K)) == NULL) goto FAILURE;
  if ((mg->isgumbel = malloc(sizeof(int)    * K)) == NULL) goto FAILURE;

  for (k = 0; k < K; k++)
    {
      mg->q[k]        = 1. / (double) K;
      mg->mu[k]       = 0.;
      mg->lambda[k]   = 1.;
      mg->alpha[k]    = 0.;
      mg->isgumbel[k] = FALSE;
    }
  return mg;
  
 FAILURE:
  esl_mixgev_Destroy(mg);
  ESL_ERROR_NULL(eslEMEM, "malloc failed");
}

/* Function:  esl_mixgev_Destroy()
 * Incept:    SRE, Mon Aug 15 08:57:48 2005 [St. Louis]
 *
 * Purpose:   Deallocates the mixture GEV parameter object <mg>.
 *
 * Args:      mg  - ptr to the object to be deallocated.
 *
 * Returns:   (void)
 */
void
esl_mixgev_Destroy(ESL_MIXGEV *mg)
{
  if (mg == NULL) return;

  if (mg->q        != NULL) free(mg->q);
  if (mg->mu       != NULL) free(mg->mu);
  if (mg->lambda   != NULL) free(mg->lambda);
  if (mg->alpha    != NULL) free(mg->alpha);
  if (mg->wrk      != NULL) free(mg->wrk);
  if (mg->isgumbel != NULL) free(mg->isgumbel);

  free(mg);
}

/* Function:  esl_mixgev_Copy()
 * Incept:    SRE, Mon Aug 15 09:04:10 2005 [St. Louis]
 *
 * Purpose:   Makes a copy of the mixture GEV parameter object <src>
 *            in <dest>. Caller must have already allocated <dest> to have
 *            (at least) the same number of components as <src>.
 *
 * Args:      src   - object to be copied
 *            dest  - allocated object to copy <src> into
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINCOMPAT> if <dest> isn't allocated with enough
 *            components to hold a copy of <src>.
 */
int
esl_mixgev_Copy(ESL_MIXGEV *src, ESL_MIXGEV *dest)
{
  int k;

  if (dest->K < src->K) 
    ESL_ERROR(eslEINCOMPAT, "mixture GEV too small to copy into");

  for (k = 0; k < src->K; k++)
    {
      dest->q[k]        = src->q[k];
      dest->mu[k]       = src->mu[k];
      dest->lambda[k]   = src->lambda[k];
      dest->alpha[k]    = src->alpha[k];
      dest->isgumbel[k] = src->isgumbel[k];
    }
  dest->K = src->K;
  return eslOK;
}

/* Function:  esl_mixgev_ForceGumbel()
 * Incept:    SRE, Mon Aug 15 09:08:35 2005 [St. Louis]
 *
 * Purpose:   Constrain component <which> of the mixture GEV <mg>
 *            to be a Gumbel (that is, constrain $\alpha=0$ for
 *            that component. This constraint will be obeyed by
 *            any subsequent calls to parameter fitting routines.
 *
 *            Normally would be called just after creating the <mg>
 *            object, as part of its configuration before trying to 
 *            fit some observed data to a mixture GEV.
 *
 * Args:      mg    - mixture GEV object being configured
 *            which - which component to constrain to a Gumbel
 *
 * Returns:   <eslOK> on success.
 */
int
esl_mixgev_ForceGumbel(ESL_MIXGEV *mg, int which)
{
  mg->isgumbel[which] = TRUE;
  return eslOK;
}
/*----------------- end ESL_MIXGEV object maintenance ----------------------*/




/****************************************************************************
 * Routines for evaluating densities and distributions
 ****************************************************************************/ 

/* Function:  esl_mixgev_pdf()
 * Incept:    SRE, Mon Aug 15 09:23:03 2005 [St. Louis]
 *
 * Purpose:   Returns the probability density function $P(X=x)$ for
 *            quantile <x>, given mixture GEV parameters <mg>.
 */
double
esl_mixgev_pdf(double x, ESL_MIXGEV *mg)
{
  double pdf = 0.;
  int    k;

  for (k = 0; k < mg->K; k++)
    pdf += mg->q[k] * esl_gev_pdf(x, mg->mu[k], mg->lambda[k], mg->alpha[k]);
  return pdf;
}

/* Function:  esl_mixgev_logpdf()
 * Incept:    SRE, Mon Aug 15 09:30:38 2005 [St. Louis]
 *
 * Purpose:   Returns the log of the PDF ($\log P(X=x)$) for quantile <x>,
 *            given mixture GEV parameters <mg>.
 */
double
esl_mixgev_logpdf(double x, ESL_MIXGEV *mg)
{
  int k;
  for (k = 0; k < mg->K; k++)
    if (mg->q[k] == 0.0) 
      mg->wrk[k] = -eslINFINITY;
    else 
      mg->wrk[k] =  log(mg->q[k]) +
	esl_gev_logpdf(x, mg->mu[k], mg->lambda[k], mg->alpha[k]);

  return esl_vec_DLogSum(mg->wrk, mg->K);
}

/* Function:  esl_mixgev_cdf()
 * Incept:    SRE, Mon Aug 15 09:51:36 2005 [St. Louis]
 *
 * Purpose:   Returns the cumulative distribution function $P(X \leq x)$
 *            for quantile <x>, given mixture GEV parameters <mg>.
 */
double
esl_mixgev_cdf(double x, ESL_MIXGEV *mg)
{
  double cdf = 0.;
  int    k;

  for (k = 0; k < mg->K; k++)
    cdf += mg->q[k] * esl_gev_cdf(x, mg->mu[k], mg->lambda[k], mg->alpha[k]);
  return cdf;
}

/* Function:  esl_mixgev_logcdf()
 * Incept:    SRE, Mon Aug 15 09:56:10 2005 [St. Louis]
 *
 * Purpose:   Returns the log of the CDF $\log P(X \leq x)$
 *            for quantile <x>, given mixture GEV parameters <mg>.
 */
double
esl_mixgev_logcdf(double x, ESL_MIXGEV *mg)
{
  int k;

  for (k = 0; k < mg->K; k++)
    if (mg->q[k] == 0.0) 
      mg->wrk[k] = -eslINFINITY;
    else
      mg->wrk[k] = log(mg->q[k]) + 
	esl_gev_logcdf(x, mg->mu[k], mg->lambda[k], mg->alpha[k]);

  return esl_vec_DLogSum(mg->wrk, mg->K);
}

/* Function:  esl_mixgev_surv()
 * Incept:    SRE, Mon Aug 15 10:00:13 2005 [St. Louis]
 *
 * Purpose:   Returns the survivor function $P(X > x)$ (1-CDF)
 *            for quantile <x>, given mixture GEV parameters <mg>.
 */
double
esl_mixgev_surv(double x, ESL_MIXGEV *mg)
{
  double srv = 0.;
  int    k;

  for (k = 0; k < mg->K; k++)
    srv += mg->q[k] * esl_gev_surv(x, mg->mu[k], mg->lambda[k], mg->alpha[k]);
  return srv;
}

/* Function:  esl_mixgev_logsurv()
 * Incept:    SRE, Mon Aug 15 10:03:55 2005 [St. Louis]
 *
 * Purpose:   Returns the log survivor function $\log P(X > x)$ (log(1-CDF))
 *            for quantile <x>, given mixture GEV parameters <mg>.
 */
double
esl_mixgev_logsurv(double x, ESL_MIXGEV *mg)
{
  int k;
  for (k = 0; k < mg->K; k++)
    {
      mg->wrk[k] =  log(mg->q[k]);
      mg->wrk[k] += esl_gev_logsurv(x, mg->mu[k], mg->lambda[k], mg->alpha[k]);
    }
  return esl_vec_DLogSum(mg->wrk, mg->K);
}

/* Function:  esl_mixgev_invcdf()
 * Incept:    SRE, Sun Aug 21 14:32:53 2005 [St. Louis]
 *
 * Purpose:   Calculates the inverse CDF for a mixture GEV <mg>,
 *            returning the quantile <x> at which the CDF is <p>,
 *            where $0 < p < 1$.
 *            
 *            The inverse CDF of a mixture model has no analytical
 *            expression as far as I'm aware. The calculation here is
 *            a brute force bisection search in <x> using the CDF
 *            function. It will suffice for a small number of calls
 *            (for plotting applications, for example), but beware, it is not
 *            efficient.
 */
double
esl_mixgev_invcdf(double p, ESL_MIXGEV *mg)
{
  double x1, x2, xm;		/* low, high guesses at x */
  double f1, f2, fm;
  double tol = 1e-6;

  x2 = esl_vec_DMin(mg->mu, mg->K);
  x1 = x2 - 1.;
  do {				/* bracket, left side */
    x1 = x1 + 2.*(x2-x1);
    f1 = esl_mixgev_cdf(x1, mg);
  } while (f1 > p);
  do {				/* bracket, right side */
    x2 = x2 + 2.*(x2-x1);
    f2 = esl_mixgev_cdf(x2, mg);
  } while (f2 < p);		

  do {				/* bisection */
    xm = (x1+x2) / 2.;
    fm = esl_mixgev_cdf(xm, mg);
    
    if      (fm > p) x2 = xm;
    else if (fm < p) x1 = xm;
    else return xm;		/* unlikely case of fm==p */
  } while ( (x2-x1)/(x1+x2+1e-9) > tol);

  xm = (x1+x2) / 2.;
  return xm;
}
/*-------------------- end densities & distributions ------------------------*/




/****************************************************************************
 * Generic API routines: for general interface w/ histogram module
 ****************************************************************************/ 

/* Function:  esl_mixgev_generic_pdf()
 * Incept:    SRE, Thu Aug 25 08:03:24 2005 [St. Louis]
 *
 * Purpose:   Generic-API wrapper around <esl_mixgev_pdf()>, taking
 *            a void ptr to a <ESL_MIXGEV> parameter structure.
 */
double
esl_mixgev_generic_pdf(double x, void *params)
{
  ESL_MIXGEV *mg = (ESL_MIXGEV *) params;
  return esl_mixgev_pdf(x, mg);
}

/* Function:  esl_mixgev_generic_cdf()
 * Incept:    SRE, Sun Aug 21 14:44:06 2005 [St. Louis]
 *
 * Purpose:   Generic-API wrapper around <esl_mixgev_cdf()>, taking
 *            a void ptr to a <ESL_MIXGEV> parameter structure.
 */
double
esl_mixgev_generic_cdf(double x, void *params)
{
  ESL_MIXGEV *mg = (ESL_MIXGEV *) params;
  return esl_mixgev_cdf(x, mg);
}

/* Function:  esl_mixgev_generic_surv()
 * Incept:    SRE, Thu Aug 25 08:03:57 2005[St. Louis]
 *
 * Purpose:   Generic-API wrapper around <esl_mixgev_surv()>, taking
 *            a void ptr to a <ESL_MIXGEV> parameter structure.
 */
double
esl_mixgev_generic_surv(double x, void *params)
{
  ESL_MIXGEV *mg = (ESL_MIXGEV *) params;
  return esl_mixgev_surv(x, mg);
}

/* Function:  esl_mixgev_generic_invcdf()
 * Incept:    SRE, Sun Aug 21 14:44:59 2005 [St. Louis]
 *
 * Purpose:   Generic-API wrapper around <esl_mixgev_invcdf()>, taking
 *            a void ptr to a <ESL_MIXGEV> parameter structure.
 */
double
esl_mixgev_generic_invcdf(double p, void *params)
{
  ESL_MIXGEV *mg = (ESL_MIXGEV *) params;
  return esl_mixgev_invcdf(p, mg);
}
/*------------------------ end generic API ---------------------------------*/




/****************************************************************************
 * Routines for dumping plots to xmgrace XY format
 ****************************************************************************/ 

/* Function:  esl_mixgev_Plot()
 * Incept:    SRE, Mon Aug 15 10:06:35 2005 [St. Louis]
 *
 * Purpose:   Plot some function <func> (for instance, <esl_mixgev_pdf()>)
 *            for mixture GEV parameters <mg>, for a range of
 *            quantiles x from <xmin> to <xmax> in steps of <xstep>;
 *            output to an open stream <fp> in xmgrace XY input format.
 *
 * Returns:   <eslOK>.
 */
int
esl_mixgev_Plot(FILE *fp, ESL_MIXGEV *mg,
		double (*func)(double x, ESL_MIXGEV *mg), 
		double xmin, double xmax, double xstep)
{
  double x;
  for (x = xmin; x <= xmax; x += xstep)
    fprintf(fp, "%f\t%g\n", x, (*func)(x, mg));
  fprintf(fp, "&\n");
  return eslOK;
}
/*-------------------- end plot dumping routines ---------------------------*/





/****************************************************************************
 * Routines for sampling (requires augmentation w/ random module)
 ****************************************************************************/ 
#ifdef eslAUGMENT_RANDOM

/* Function:  esl_mixgev_Sample()
 * Incept:    SRE, Mon Aug 15 10:14:23 2005 [St. Louis]
 *
 * Purpose:   Sample a random variate x from a mixture GEV <mg>, 
 *            given random number source <r>.
 */
double
esl_mixgev_Sample(ESL_RANDOMNESS *r, ESL_MIXGEV *mg)
{
  int k;	
  k = esl_rnd_DChoose(r, mg->q, mg->K);
  return esl_gev_Sample(r, mg->mu[k], mg->lambda[k], mg->alpha[k]);
}

#endif /*eslAUGMENT_RANDOM*/
/*--------------------------- end sampling ---------------------------------*/






/****************************************************************************
 * Maximum likelihood fitting to mixture GEV distributions
 ****************************************************************************/ 
#ifdef eslAUGMENT_MINIMIZER

struct mixgev_data {
  double *x;
  int     n;
  double *wrk;	    /* workspace vector               [0..K-1]*/
  ESL_MIXGEV *mg;
};

/* Given mixture GEV parameters in <mg>;
 * do appropriate c.o.v.'s to unconstrained real parameters
 * and fill in the packed parameter vector <p>.
 * 
 * First K-1 are Q_1..Q_K-1 mixture coefficient parameters; Q_0 implicitly 0;
 *  cov is q_k = e^{Q_k} / \sum_j e^{Q_j},
 *  so     Q_k = log(q_k) - log(q_0).
 * Then K components: mu, lambda, optional alpha;
 * mu, alpha are already unconstrained real;
 * lambda cov is lambda = e^w, w = log(lambda).
 */
static void
mixgev_pack_paramvector(double *p, int np, ESL_MIXGEV *mg)
{
  int    i;			/* counter in parameter vector p */
  int    k;			/* counter in mixture components */
  double z;			/* tmp variable */

  /* mixture coefficients */
  z = log(mg->q[0]);
  i = 0;
  for (k = 1; k < mg->K; k++) 
    p[i++] = log(mg->q[k]) - z;
  
  /* gev parameters */
  for (k = 0; k < mg->K; k++)
    {
      p[i++] = mg->mu[k];
      p[i++] = log(mg->lambda[k]);
      if (! mg->isgumbel[k]) p[i++] = mg->alpha[k];
    }
  /* assert(i==np) in debugging, if you want */
}

/* Same as above but in reverse: given parameter vector <p>,
 * do appropriate c.o.v. back to desired parameter space, and
 * fill in the mixture GEV structure <mg>.
 */
static void
mixgev_unpack_paramvector(double *p, int np, ESL_MIXGEV *mg)
{
  int    i;			/* counter in parameter vector p */
  int    k;			/* counter in mixture components */
  double z;			/* tmp variable  */

  /* Fetch the params in their c.o.v. space first
   */
  i = 0;
  mg->q[0] = 0;	/* implicitly */
  for (k = 1; k < mg->K; k++) 
    mg->q[k] = p[i++]; 
  for (k = 0; k < mg->K; k++)
    {
      mg->mu[k]     = p[i++];
      mg->lambda[k] = p[i++];
      if (!mg->isgumbel[k]) mg->alpha[k]  = p[i++];
      else                  mg->alpha[k]  = 0.;
    }
  assert(i==np);
  
  /* Convert mix coefficients back to probabilities;
   * their  c.o.v. is q_k = e^{Q_k} / \sum_k e^{Q_k}
   * which rearranges to exp(Q_k - log[\sum_k e^Q_k]),
   * and we have the DLogSum() function to compute the log sum.
   */
  z = esl_vec_DLogSum(mg->q, mg->K);
  for (k = 0; k < mg->K; k++)
    mg->q[k] = exp(mg->q[k] - z);
  
  /* lambda c.o.v. is \lambda = e^w
   */
  for (k = 0; k < mg->K; k++)
    mg->lambda[k] = exp(mg->lambda[k]);
}

static double
mixgev_complete_func(double *p, int np, void *dptr)
{
  struct mixgev_data *data = (struct mixgev_data *) dptr;
  ESL_MIXGEV         *mg   = data->mg;
  int    i;
  double logL;

  /* Use the current parameter vector (in its unconstrained
   * real c.o.v. space) to deduce what the current mixture GEV
   * parameters are:
   */
  mixgev_unpack_paramvector(p, np, mg);

  /* Calculate the log likelihood:
   */
  logL = 0;
  for (i = 0; i < data->n; i++)
    logL += esl_mixgev_logpdf(data->x[i], mg);

  /* return the NLL
   */
  return -logL;
}


/* Make initial guesses at parameters.
 * We make "central" guesses based on method-of-moments for a single Gumbel.
 * Mixture coefficients are sampled uniformly.
 * mu, lambda are randomized +/- 10% around central
 * alpha is randomized +/- 0.1 around 0.
 */
int
esl_mixgev_FitGuess(ESL_RANDOMNESS *r, double *x, int n, ESL_MIXGEV *mg)
{
  double mean, variance;
  double mu, lambda;
  int    k;

  esl_stats_Mean(x, n, &mean, &variance);
  lambda = eslCONST_PI / sqrt(6.*variance);
  mu     = mean - 0.57722/lambda;

  esl_dirichlet_SampleUniform(r, mg->K, mg->q);
  for (k = 0; k < mg->K; k++)
    {
      mg->mu[k]     = mu     + 0.2 * mu     * (esl_random(r) - 0.5);
      mg->lambda[k] = lambda + 0.2 * lambda * (esl_random(r) - 0.5);
      if (mg->isgumbel[k]) mg->alpha[k] = 0.;
      else mg->alpha[k] = 0.2 * (esl_random(r) - 0.5);
    }
  return eslOK;
}


/* 
 *   <mg> is set to a starting point in parameter space;
 *   see <esl_mixgev_FitGuess()>, for example.
 */
int
esl_mixgev_FitComplete(double *x, int n, ESL_MIXGEV *mg)
{
  struct mixgev_data data;
  int     status;
  double *p;
  double *u;
  double *wrk;
  double  tol;
  int     np;
  double  fx;
  int     k;
  int     i;

  tol = 1e-6;

  /* Determine number of free parameters and allocate 
   */
  np = mg->K-1;			/* K-1 mix coefficients free */
  for (k = 0; k < mg->K; k++)
    np += (mg->isgumbel[k])? 2 : 3;
  p   = malloc(sizeof(double) * np);
  u   = malloc(sizeof(double) * np);
  wrk = malloc(sizeof(double) * np * 4);

  /* Copy shared info into the "data" structure
   */
  data.x   = x;
  data.n   = n;
  data.wrk = wrk;
  data.mg  = mg;

  /* From mg, create the parameter vector.
   */
  mixgev_pack_paramvector(p, np, mg);

  /* Define the step size vector u.
   */
  i = 0;
  for (k = 1; k < mg->K; k++) u[i++] = 1.0;
  for (k = 0; k < mg->K; k++)
    {
      u[i++] = 1.0;
      u[i++] = 1.0;
      if (! mg->isgumbel[k]) u[i++] = 0.02;
    }
  assert(np == i);

  /* Feed it all to the mighty optimizer.
   */

  status = esl_min_ConjugateGradientDescent(p, u, np, &mixgev_complete_func, NULL,
					    (void *) (&data), tol, wrk, &fx);

  /* Convert the final parameter vector back to a mixture GEV
   */
  mixgev_unpack_paramvector(p, np, mg);
  
  free(p);
  free(u);
  free(wrk);
  return status;
}
#endif /*eslAUGMENT_MINIMIZER*/
/*--------------------------- end fitting ----------------------------------*/




/****************************************************************************
 * Example main()
 ****************************************************************************/ 

#ifdef eslMIXGEV_EXAMPLE
/*::cexcerpt::mixgev_example::begin::*/
/* compile: 
   gcc -g -Wall -I. -I ~/src/easel -L ~/src/easel -o example -DeslMIXGEV_EXAMPLE\
     mixgev.c -leasel -lm 
 * run:     ./example
 */
#include <stdio.h>
#include <easel.h>
#include <mixgev.h>
#include <esl_random.h>

int
main(int argc, char **argv)
{
  FILE *fp;
  ESL_RANDOMNESS *r;		/* source of random numbers   */
  ESL_MIXGEV *mg;		/* mixture GEV to sample from */
  ESL_MIXGEV *emg;		/* estimated mixture GEV      */
  double     *x;		/* sampled dataset            */
  int         n = 100000;	/* number of samples          */
  int         i;
  int         k;
  double      nll;
  double      min, max;

  r  = esl_randomness_Create(42);
  mg = esl_mixgev_Create(2);
  mg->q[0]      = 0.85;   mg->q[1]      = 0.15; 
  mg->mu[0]     = -2.72;  mg->mu[1]     = -2.0; 
  mg->lambda[0] = 2.5;    mg->lambda[1] = 1.0;  
  mg->alpha[0]  = 0.;     mg->alpha[1]  = 0.09; 

  nll = 0.;
  min = 99999;
  max = -99999;
  x = malloc(sizeof(double) * n);
  for (i = 0; i < n; i++)
    {
      x[i] = esl_mixgev_Sample(r, mg);
      nll -= esl_mixgev_logpdf(x[i], mg);
      if (x[i] > max) max = x[i];
      if (x[i] < min) min = x[i];
    }
  printf("NLL of known mixGEV: %g\n", nll);

  /* Dump the raw data samples to an R file.
   */
  fp = fopen("data.out", "w");
  fprintf(fp, "     val\n");
  for (i = 0; i < n; i++)
    fprintf(fp, "%d   %f\n", i+1, x[i]);
  fclose(fp);

  emg = esl_mixgev_Create(2);
  esl_mixgev_FitGuess(r, x, n, emg); 
  /*  esl_mixgev_Copy(mg, emg); */
  esl_mixgev_SetGumbelConstraint(emg, 0); 
  esl_mixgev_FitComplete(x, n, emg);

  printf("Component   q      mu   lambda  alpha\n");
  for (k=0; k < 2; k++)
    printf("%d\t%7.4f\t%7.2f\t%7.4f\t%7.4f\n", 
	   k, emg->q[k], emg->mu[k], emg->lambda[k], emg->alpha[k]);

  nll = 0.;
  for (i = 0; i < n; i++)
    nll -= esl_mixgev_logpdf(x[i], emg);
  printf("NLL of fitted mixGEV: %g\n", nll);

  /* Dump some R commands for showing these distributions
   */
  printf("library(ismev)\n");
  printf("library(evd)\n");

  printf("d <- read.table(\"data.out\")$val\n");
  printf("plot(density(d,bw=0.2), log=\"y\")\n");
  printf("min <- %f\n", min);
  printf("max <- %f\n", max);
  printf("xax <- seq(min-2, max+5, by=0.1)\n");
  printf("cc <- xax - xax\n");
  printf("zc <- xax - xax\n");
  for (k = 0; k < mg->K; k++)
    {
      printf("c%d  <- %f * dgev(xax, %f, %f, %f)\n", 
	     k, mg->q[k], mg->mu[k], 1./mg->lambda[k], mg->alpha[k]);
      printf("cc   <- cc + c%d\n", k);
      printf("lines(xax, c%d, col=\"blue\")\n", k);
    }
  for (k = 0; k < emg->K; k++)
    {
      printf("z%d  <- %f * dgev(xax, %f, %f, %f)\n", 
	     k, emg->q[k], emg->mu[k], 1./emg->lambda[k], emg->alpha[k]);
      printf("zc   <- zc + z%d\n", k);
      printf("lines(xax, z%d, col=\"blue\")\n", k);
    }
  printf("lines(xax, cc, col=\"green\")\n");
  printf("lines(xax, zc, col=\"red\")\n");

  esl_mixgev_Destroy(mg);
  esl_mixgev_Destroy(emg);
  esl_randomness_Destroy(r);
  free(x);
  return 0;
}

/*::cexcerpt::mixgev_example::end::*/
#endif /*eslMIXGEV_EXAMPLE*/


/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/
