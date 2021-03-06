/* esl_histogram.h
 * Collection and display of score histograms.
 * 
 * SRE, Fri Jul  1 13:22:45 2005 [St. Louis]
 * SVN $Id: esl_histogram.h 82 2005-12-13 20:21:07Z eddy $
 */
#ifndef ESL_HISTOGRAM_INCLUDED
#define ESL_HISTOGRAM_INCLUDED

#include <math.h>   /* floor() is in one of the macros */


/* Structure: ESL_HISTOGRAM
 * 
 * Keeps a score histogram, in which scores are counted into bins of
 * size (width) w. 
 *   histogram starts at bmin <  floor(xmin/w) * w
 *   histogram ends at   bmax >= ceil(xmax/w)*w
 *   nb = (bmax-bmin)/w
 *   each score x is counted into bin b = nb - (int) (bmax-x)/w
 *   each bin b contains scores bw+bmin < x <= (b+1)w + bmin
 */  
typedef struct {
  /* The histogram is kept as counts in fixed-width bins.
   */
  double  xmin, xmax;	/* smallest, largest sample value observed          */
  int     n;            /* total number of raw data samples                 */
  int    *obs;		/* observed counts in bin b, 0..nb-1 (dynamic)      */
  double  bmin, bmax;	/* histogram bounds: all x satisfy bmin < x <= bmax */
  int     imin, imax;	/* smallest, largest bin that contain obs[i] > 0    */
  int     nb;           /* number of bins                                   */
  double  w;		/* fixed width of each bin                          */

  /* Optionally, in a "full" h, we can also keep all the raw samples in x.
   */
  double *x;		/* optional: raw sample values x[0..n-1]            */
  int     nalloc;	/* current allocated size of x                      */

  /* SetExpect() uses some theoretical distribution to set expect[0..nb-1].
   * If that dist was fitted to a censored tail, the caller uses SetCensoring()
   * to inform the histogram of that, setting phi and z.
   */
  double *expect;	/* expected counts in bin b, 0..nb-1 (not resized)  */
  double  phi;		/* censoring value; all x_i > phi                   */
  int     cmin;		/* smallest bin index that contains uncensored data */
  int     z;		/* # of censored values <= phi                      */
  int     Nc;	        /* # samples in complete data (including unobs)     */
  int     No;		/* # of samples in observed data                    */
  int     Nx;		/* # of samples modeled by a fitted distribution    */

  /* Some status flags
   */
  int is_full;		/* TRUE when we're keeping raw data in x           */
  int is_sorted;	/* TRUE if x is sorted smallest-to-largest         */
  int is_rounded;	/* TRUE to use/display only bins, not raw          */
  int is_tailfit;	/* TRUE if expected dist only describes tail       */
  enum { COMPLETE, VIRTUAL_CENSORED, TRUE_CENSORED } dataset_is; 

} ESL_HISTOGRAM;

#define esl_histogram_Bin2LBound(h,b)  ((h)->w*(b) + (h)->bmin)
#define esl_histogram_Bin2UBound(h,b)  ((h)->w*((b)+1) + (h)->bmin)
#define esl_histogram_Score2Bin(h,x)   ((int) ceil( ((x - (h)->bmin) / h->w) - 1.))


extern ESL_HISTOGRAM *esl_histogram_Create    (double bmin, double bmax, double w);
extern ESL_HISTOGRAM *esl_histogram_CreateFull(double bmin, double bmax, double w);
extern void           esl_histogram_Destroy(ESL_HISTOGRAM *h);

extern int esl_histogram_Add(ESL_HISTOGRAM *h, double x);

extern int esl_histogram_Sort(ESL_HISTOGRAM *h);
extern int esl_histogram_GetScoreAtRank(ESL_HISTOGRAM *h, int rank, double *ret_x);

extern int esl_histogram_TrueCensoring    (ESL_HISTOGRAM *h, int z, double phi);
extern int esl_histogram_VirtCensorByValue(ESL_HISTOGRAM *h, double phi);
extern int esl_histogram_VirtCensorByMass (ESL_HISTOGRAM *h, double tfrac);
extern int esl_histogram_DeclareTailfitting(ESL_HISTOGRAM *h);
extern int esl_histogram_DeclareRounding   (ESL_HISTOGRAM *h);
extern int esl_histogram_SetExpect(ESL_HISTOGRAM *h, 
				   double (*cdf)(double x, void *params),
				   void *params);


extern int esl_histogram_Print       (FILE *fp, ESL_HISTOGRAM *h);
extern int esl_histogram_Plot        (FILE *fp, ESL_HISTOGRAM *h);
extern int esl_histogram_PlotSurvival(FILE *fp, ESL_HISTOGRAM *h);
extern int esl_histogram_PlotTheory  (FILE *fp, ESL_HISTOGRAM *h, 
				      double (*fx)(double, void *), void *params);
extern int esl_histogram_PlotQQ      (FILE *fp, ESL_HISTOGRAM *h, 
				      double (*invcdf)(double, void *), void *params);

#ifdef eslAUGMENT_STATS
extern int esl_histogram_Goodness(ESL_HISTOGRAM *h, 
				  double (*cdf)(double x, void *params),
				  void *params, int nfitted, int *ret_nbins,
				  double *ret_G,  double *ret_Gp,
				  double *ret_X2, double *ret_X2p);
#endif



#endif /*!ESL_HISTOGRAM_INCLUDED*/
/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/
