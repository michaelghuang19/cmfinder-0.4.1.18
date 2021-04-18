/* esl_stats.h
 * Foundation for the statistics modules.
 * 
 * SRE, Tue Jul 19 11:35:28 2005
 * SVN $Id: esl_stats.h 82 2005-12-13 20:21:07Z eddy $
 */
#ifndef ESL_STATS_INCLUDED
#define ESL_STATS_INCLUDED

extern int esl_stats_Mean(double *x, int n, double *ret_mean, double *ret_var);
extern int esl_stats_LogGamma(double x, double *ret_answer);
extern int esl_stats_Psi(double x, double *ret_answer);
extern int esl_stats_IncompleteGamma(double a, double x, double *ret_pax, double *ret_qax);
extern int esl_stats_ChiSquaredTest(int v, double x, double *ret_answer);

#endif /*ESL_STATS_INCLUDED*/
/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/
