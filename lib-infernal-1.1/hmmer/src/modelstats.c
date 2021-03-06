/* Miscellaneous summary statistics calculated for HMMs and profiles.
 * 
 * SRE, Fri May  4 11:43:20 2007 [Janelia]
 * SVN $Id: modelstats.c 4103 2012-06-24 02:09:43Z wheelert $
 */

#include "p7_config.h"

#include "easel.h"
#include "esl_vectorops.h"

#include "hmmer.h"



/* Function:  p7_MeanMatchInfo()
 * Incept:    SRE, Fri May  4 11:43:56 2007 [Janelia]
 *
 * Purpose:   Calculate the mean information content per match state
 *            emission distribution, in bits:
 *            
 *            \[
 *              \frac{1}{M} \sum_{k=1}^{M}
 *                \left[ 
 *                     \sum_x p_k(x) \log_2 p_k(x) 
 *                   - \sum_x f(x) \log_2 f(x)
 *                \right] 
 *            \]
 *            
 *            where $p_k(x)$ is emission probability for symbol $x$
 *            from match state $k$, and $f(x)$ is the null model's
 *            background emission probability for $x$.
 */
double
p7_MeanMatchInfo(const P7_HMM *hmm, const P7_BG *bg)
{
  return esl_vec_FEntropy(bg->f, hmm->abc->K) - p7_MeanMatchEntropy(hmm);
}

/* Function:  p7_MeanMatchEntropy()
 * Incept:    SRE, Fri May  4 13:37:15 2007 [Janelia]
 *
 * Purpose:   Calculate the mean entropy per match state emission
 *            distribution, in bits:
 *            
 *            \[
 *              - \frac{1}{M} \sum_{k=1}^{M} \sum_x p_k(x) \log_2 p_k(x)
 *            \]
 *       
 *            where $p_k(x)$ is emission probability for symbol $x$
 *            from match state $k$.
 */
double
p7_MeanMatchEntropy(const P7_HMM *hmm)
{
  int    k;
  double H = 0.;

  for (k = 1; k <= hmm->M; k++)
    H += esl_vec_FEntropy(hmm->mat[k], hmm->abc->K);
  H /= (double) hmm->M;
  return H;
}


/* Function:  p7_MeanMatchRelativeEntropy()
 * Incept:    SRE, Fri May 11 09:25:01 2007 [Janelia]
 *
 * Purpose:   Calculate the mean relative entropy per match state emission
 *            distribution, in bits:
 *            
 *            \[
 *              \frac{1}{M} \sum_{k=1}^{M} \sum_x p_k(x) \log_2 \frac{p_k(x)}{f(x)}
 *            \]
 *       
 *            where $p_k(x)$ is emission probability for symbol $x$
 *            from match state $k$, and $f(x)$ is the null model's 
 *            background emission probability for $x$. 
 */
double
p7_MeanMatchRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg)
{
  int    k;
  double KL = 0.;

#if 0
  p7_bg_Dump(stdout, hmm->bg);
  for (k = 1; k <= hmm->M; k++)
    printf("Match %d : %.2f %.2f\n", k, 
	   esl_vec_FRelEntropy(hmm->mat[k], hmm->bg->f, hmm->abc->K),
	   esl_vec_FEntropy(bg->f, hmm->abc->K) - esl_vec_FEntropy(hmm->mat[k], hmm->abc->K));
#endif

  for (k = 1; k <= hmm->M; k++)
    KL += esl_vec_FRelEntropy(hmm->mat[k], bg->f, hmm->abc->K);
  KL /= (double) hmm->M;
  return KL;
}



double
p7_MeanForwardScore(const P7_HMM *hmm, const P7_BG *bg)
{
  int             L   = 350;
  int             N   = 100;
  P7_PROFILE     *gm  = p7_profile_Create(hmm->M, hmm->abc);
  P7_GMX         *gx  = p7_gmx_Create(gm->M, L);
  ESL_SQ         *sq  = esl_sq_CreateDigital(hmm->abc);
  ESL_RANDOMNESS *r   = esl_randomness_CreateFast(0);
  float           fsc;
  float           nullsc;
  double          bitscore;
  double          sum = 0.;
  int             i;

  if (p7_ProfileConfig (hmm, bg, gm, L, p7_LOCAL)        != eslOK) p7_Die("failed to configure profile");
  for (i = 0; i < N; i++)
    {
      if (p7_ReconfigLength(gm, L)                        != eslOK) p7_Die("failed to reconfig profile length");
      if (p7_ProfileEmit(r, hmm, gm, bg, sq, NULL)        != eslOK) p7_Die("failed to emit sequence");
      if (p7_ReconfigLength(gm, sq->n)                    != eslOK) p7_Die("failed to reconfig profile length");
      if (p7_gmx_GrowTo(gx, gm->M, sq->n)                 != eslOK) p7_Die("failed to grow the matrix");
      if (p7_GForward(sq->dsq, sq->n, gm, gx, &fsc)       != eslOK) p7_Die("failed to run Forward");
      if (p7_bg_NullOne(bg, sq->dsq, sq->n, &nullsc)      != eslOK) p7_Die("failed to run bg_NullOne()");
      bitscore = (fsc - nullsc) / eslCONST_LOG2;

      sum += bitscore;
    }
  
  esl_randomness_Destroy(r);
  esl_sq_Destroy(sq);
  p7_gmx_Destroy(gx);
  p7_profile_Destroy(gm);
  return (sum / (double) N);
}


/* Function:  p7_MeanPositionRelativeEntropy()
 * Synopsis:  Calculate the mean score per match position, including gap cost.
 * Incept:    SRE, Thu Sep  6 10:26:14 2007 [Janelia]
 *
 * Purpose:   Calculate the mean score (relative entropy) in bits per 
 *            match (consensus) position in model <hmm>, given background
 *            model <bg>.
 *            
 *            More specifically: the mean bitscore is weighted by
 *            match state occupancy (match states that aren't used
 *            much are downweighted), and the log transitions into
 *            that match state from the previous M, D, or I are
 *            counted against it, weighted by their probability.
 *            
 *            This isn't a complete accounting of the average score
 *            per model position nor per aligned residue; most
 *            notably, it doesn't include the contribution of
 *            entry/exit probabilities. So don't expect to approximate
 *            average scores by multiplying <*ret_entropy> by <M>.
 *
 * Returns:   <eslOK> on success, and <*ret_entropy> is the result.
 *
 * Throws:    <eslEMEM> on allocation failure, and <*ret_entropy> is 0.
 */
int
p7_MeanPositionRelativeEntropy(const P7_HMM *hmm, const P7_BG *bg, double *ret_entropy)
{
  int     status;
  float  *mocc = NULL;
  int     k;
  double  mre, tre;
  double  xm, xi, xd;
  
  ESL_ALLOC(mocc, sizeof(float) * (hmm->M+1));
  if ((status = p7_hmm_CalculateOccupancy(hmm, mocc, NULL)) != eslOK) goto ERROR;
  
  /* mre = the weighted relative entropy per match emission */
  for (mre = 0., k = 1; k <= hmm->M; k++)
    mre += mocc[k] * esl_vec_FRelEntropy(hmm->mat[k], bg->f, hmm->abc->K);
  mre /= esl_vec_FSum(mocc+1, hmm->M);

  /* The weighted relative entropy per match entry transition, 2..M 
   */
  for (tre = 0., k = 2; k <= hmm->M; k++)
    {
      xm = mocc[k-1]*hmm->t[k-1][p7H_MM] * log(hmm->t[k-1][p7H_MM] / bg->p1);
      xi = mocc[k-1]*hmm->t[k-1][p7H_MI] * (log(hmm->t[k-1][p7H_MM] / bg->p1) + log(hmm->t[k-1][p7H_IM] / bg->p1));
      xd = (1.-mocc[k-1])*hmm->t[k-1][p7H_DM] * log(hmm->t[k-1][p7H_DM] / bg->p1);
      tre += (xm+xi+xd) / eslCONST_LOG2;
    }
  tre /= esl_vec_FSum(mocc+2, hmm->M-1);

  free(mocc);
  *ret_entropy = mre+tre;
  return eslOK;

 ERROR:
  if (mocc != NULL) free(mocc);
  *ret_entropy = 0.;
  return status;
}


/* Function:  p7_hmm_CompositionKLDist()
 * Synopsis:  A statistic of model's composition bias.
 * Incept:    SRE, Mon Jul  2 08:40:12 2007 [Janelia]
 *
 * Purpose:   Calculates the K-L distance between the average match
 *            state residue composition in model <hmm> and a
 *            background frequency distribution in <bg>, and
 *            return it in <ret_KL>. 
 *            
 *            Optionally return the average match state residue
 *            composition in <opt_avp>. This vector, of length
 *            <hmm->abc->K> is allocated here and becomes the caller's
 *            responsibility if <opt_avp> is non-<NULL>. 
 *            
 *            The average match composition is an occupancy-weighted
 *            average (see <p7_hmm_CalculateOccupancy()>.
 *            
 *            The `K-L distance' <*ret_KL> is the symmetricized
 *            Kullback-Leibler distance in bits (log base 2).
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
p7_hmm_CompositionKLDist(P7_HMM *hmm, P7_BG *bg, float *ret_KL, float **opt_avp)
{
  int    K   = hmm->abc->K;
  float *occ = NULL;
  float *p   = NULL;
  int    status;
  int    k;

  ESL_ALLOC(occ, sizeof(float) * (hmm->M+1));
  ESL_ALLOC(p,   sizeof(float) * K);
  p7_hmm_CalculateOccupancy(hmm, occ, NULL);

  esl_vec_FSet(p, K, 0.);
  for (k = 1; k <= hmm->M; k++)
    esl_vec_FAddScaled(p, hmm->mat[k], occ[k], K);
  esl_vec_FNorm(p, K);

  *ret_KL = (esl_vec_FRelEntropy(p, bg->f, K) + esl_vec_FRelEntropy(bg->f, p, K)) / (2.0 * eslCONST_LOG2);
  if (opt_avp != NULL) *opt_avp = p;  else free(p); 
  free(occ);
  return eslOK;
  
 ERROR:
  if (occ != NULL) free(occ);
  if (p   != NULL) free(p);
  *ret_KL = 0.0;
  if (opt_avp != NULL) *opt_avp = NULL;
  return status;
}




/* Function:  p7_hmm_GetSimpleRepeats()
 * Synopsis:  Print locations of regions identified as containing
 *            simple tandem repeats
 *
 * Purpose:   Scan the model <hmm>, searching for regions with simple
 *            tandem repeats, and return their locations in the model
 *            using <ranges>.
 *
 *            Tests for repeats of up to length <maxK>, and requires
 *            at least <min_rep> consecutive (ungapped) instances,
 *            where the average relative entropy (depends on <hmm>
 *            and <bg>) of the length-K repeat must exceed <relent>
 *            bits.
 *
 *            Basic idea: build a circular list for each length
 *            1..maxK. Suppose we're dealing with list of length
 *            K=3. Then node 1 will count the number of times
 *            each character is observed in the model's consensus
 *            at positions n*3+1 over the most recent k*min_rep
 *            states. Once k*min_rep states have been added to the
 *            node counts, compute avg relative entropy, and compare
 *            to threshold. If passes, that's a repeat seed (see
 *            below), otherwise, remove the count associated with
 *            position 1, and add the next position. Repeat,
 *            removing count for position 2 and adding for
 *            k*min_rep+1, then removing 3 and adding k*min_rep+2,
 *            and so on.
 *
 *            Once complete, ignore long seeds that overlap short
 *            seeds, then for each seed, attempt to extend by lining
 *            up consecutive length-K blocks on either side of the
 *            seed, allowing up to length-2 indels.
 *
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation error.
 */

int
p7_hmm_GetSimpleRepeats(P7_HMM *hmm, int maxK, int min_rep, int min_length, float relent_thresh, P7_HMM_WINDOWLIST *ranges)
{
  int    K   = hmm->abc->K;
  int    i, j, k, n, x;
  int  **counts;
  int    pos;  // where we are in the circular list
  int    status;
  float  relent;
  P7_BG *bg;
  int    reps;


  // can't hope to find a repeat under this definition
  if (hmm->M < min_rep)
    return eslOK;

  if ((bg = p7_bg_Create(hmm->abc))              == NULL)  esl_fatal("failed to create null model");

  // for each circular-list length
  for (n=2; n<=maxK; n++) {

    reps = min_rep;
    //go to the next # of reps required to reach min-length with equal # chars from each node
    if (n*reps < min_length) {
      reps = ceil(min_length/n);
    }


    ESL_ALLOC(counts, n*sizeof(int*)) ;  // there are n nodes in the list

    //establish count vector for each node
    for (i=0; i<n; i++) {
      ESL_ALLOC(counts[i], K*sizeof(int));
      esl_vec_ISet(counts[i], K, 0);
    }

    pos = 0;
    //for each model position
    for (i = 1; i <= hmm->M; i++) {
      // add to count vector
      x = esl_vec_FArgMax(hmm->mat[i], hmm->abc->K);
      counts[pos][x]++;  //TODO: could use emission probabilities instead, using  hmm->mat[k]
      pos++;
      pos %= n;

      //if we've accrued enough counts, check relent
      if (i >= n*reps  ) {
        //that means go through each node, computing relent, add them up, then divide by # nodes
        relent = 0.0;
        for (j=0; j<n; j++) { //for each node

          for (k=0; k<K; k++) {
            float p       =  (float)(counts[j][k])/(float)reps;
            float logodds = log((p+.000000001) / bg->f[k]);  //natural log  add that tiny amount to avoid log-of-zero trouble ... just means the definition of 0*(log(0)) = 0
            relent += p * logodds / log(2);
          }
        }

        if (relent/n > relent_thresh) {//avg relent
          //printf ("%d : %.2f  (%d-mer)\n", i, rel_ent, n);
          p7_hmmwindow_new(ranges, 0, i, 0, 0, n*reps, relent/n, 0);
        }

        //remove the first character of the preceding (n*min_rep), since it'll be replaced with a new one
        x = esl_vec_FArgMax(hmm->mat[i - (n*reps) + 1], hmm->abc->K);
        counts[pos][x]--;
      }

    }

    for (i=0; i<n; i++)
      if (counts[i] != NULL) free(counts[i]) ;
    free(counts);

  }

  return eslOK;

 ERROR:

 if (counts != NULL) {
   for (i=0; i<n; i++)
     if (counts[i] != NULL) free(counts[i]) ;
   free(counts);
 }
  return status;
}


/*****************************************************************
 * Infernal - inference of RNA secondary structure alignments
 * Version 1.1rc2; December 2012
 * Copyright (C) 2012 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Infernal is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 *****************************************************************/
