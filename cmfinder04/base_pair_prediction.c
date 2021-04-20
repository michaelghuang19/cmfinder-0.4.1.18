#include "esl_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_distance.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msaweight.h"
#include "esl_msacluster.h"
#include "esl_stack.h"
#include "esl_stopwatch.h"
#include "esl_tree.h"
#include "esl_vectorops.h"

#include "hmmer.h"

#include "infernal.h"

#include "cand.h"
#include "cmfinder.h"

#include "part_func.h"
#include "utils.h"
#include "fold_vars.h"
#include "fold.h"

/* adapted from CalcCols in pscore, but I've decided not to keep all that complexity, especially since I think the fuzzy close-to-degen idea is probably less safe for CMfinder
 * if !use_fragmentary, then start/end will just be the first/last columns of the MSA
 */
void CalcFragmentStartEnds (ESL_MSA *digital_msa,int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col)
{
	int i,j;
	
	if ((digital_msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}

	if (use_fragmentary) {
		for (j=0; j < digital_msa->nseq; j++) {
			for (i=0; i < digital_msa->alen;i++) { /* from 5' end */
				if (!esl_abc_XIsGap(abc,digital_msa->ax[j][i+1])) {
					break;
				}
			}
			if (i==0) {
				/* no evidence of truncation */
				seq_fragment_start_col[j]=-1;
			}
			else {
				seq_fragment_start_col[j]=i;
			}
			for (i=digital_msa->alen-1; i>=0; i--) { /* from 3' end */
				if (!esl_abc_XIsGap(abc,digital_msa->ax[j][i+1])) {
					break;
				}
			}
			if (i==digital_msa->alen-1) {
				seq_fragment_end_col[j]=-1;
			}
			else {
				seq_fragment_end_col[j]=i;
			}
		}
	}
	else {
		for (j=0; j < digital_msa->nseq; j++) {
			seq_fragment_start_col[j]=-1;
			seq_fragment_end_col[j]=-1;
		}
	}
}

double standard_singlet_freq[]={0.238,0.235,0.261,0.266};  // from pscore's default files
double standard_pair_freq[]={0.0081,0.0160,0.0091,0.1554,0.0160,0.0084,0.2315,0.0098,0.0091,0.2315,0.0093,0.0596,0.1554,0.0098,0.0596,0.0114};
double get_standard_singlet_freq (int nuc) { assert(nuc>=0 && nuc<4); return standard_singlet_freq[nuc]; }
double get_standard_pair_freq (int leftNuc,int rightNuc) { assert(leftNuc>=0 && leftNuc<4); assert(rightNuc>=0 && rightNuc<4); return standard_pair_freq[leftNuc*4+rightNuc]; }

/* use column-only (independent) probabilities.  we just have a pair model and singlet model, and we calculate lod scores based on the
 * likelihoods for the two models within all sequences.  because we're calculating the likelihoods of all seqs within the alignment at
 * a given position pair, the individual probabilities are multiplied together, but we do it in log scale.  the weighting is multiplicative
 * relative to the log-scale likelihoods, because it's kind of saying how many times we observed that sequence pair.
 * 
 * we count gaps as degenerate 'N' nucleotides, since that's an easy solution, and should penalize them.
 * however, we ignore double gaps
 * we also ignore seqs with a truncated gap (if we're in truncation-aware mode)
 * 
 * It's potentially tricky to treat the pair model and singlet model separately because of the possibility of truncations that could
 * invalidate a pair but not only of the single nucleotides involved in the pair.  Therefore, for simplicity, we calculate both models
 * at the same time for each pair, so the singlet model gets calculated more times than might be necessary.  (The extra CPU time should
 * be trivial.)
 * 
 */
void prxy_column_only(ESL_MSA *msa,double **bp_pr,CmfinderVars *cmfinderVars,int *seq_fragment_start_col,int *seq_fragment_end_col,double ***ret_prxy)
{
	int      i,j,k;
	double **lod;
	ESL_STOPWATCH *timer_pf=esl_stopwatch_Create();
	int nucN=esl_abc_DigitizeSymbol(abc,'N');

	esl_stopwatch_Start(timer_pf);

	lod = DoubleAlloc2DArray(msa->alen+1);
	
	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}

	for (j = 1; j < msa->alen; j++) {
		for (i = 0; i < j; i++) {
			
			double sum_lod=0;
			double sum_weight=0;

			for (k = 0; k < msa->nseq; k++) {

				int nuci=msa->ax[k][i+1];
				int nucj=msa->ax[k][j+1];
				
				// check if we should process this seq
				int count_this_weight=1;
				if (cmfinderVars->use_fragmentary) {
					if (seq_fragment_start_col[k]!=-1 && i<seq_fragment_start_col[k]) {  /* i<j, so if j is left-truncated, then i will be also */
						assert(esl_abc_XIsGap(abc,nuci)); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
						count_this_weight=0;
					}
					if (seq_fragment_end_col[k]!=-1 && j>seq_fragment_end_col[k]) {
						assert(esl_abc_XIsGap(abc,nucj)); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
						count_this_weight=0;
					}
				}
				if (esl_abc_XIsGap(abc,nuci) && esl_abc_XIsGap(abc,nucj)) {
					// double gap, so we skip it
					count_this_weight=0;
				}
				
				if (count_this_weight) {
					int ni,nj;
					int numDegenI,numDegenJ,numDegenPair;
					double degenSumI,degenSumJ,degenSumPair;
					assert(esl_abc_XIsGap(abc,nuci) || esl_abc_XIsResidue(abc,nuci));
					assert(esl_abc_XIsGap(abc,nucj) || esl_abc_XIsResidue(abc,nucj));
					
					// treat gaps as the degenerate 'N' nuc
					if (esl_abc_XIsGap(abc,nuci)) {
						nuci=nucN;
					}
					if (esl_abc_XIsGap(abc,nucj)) {
						nucj=nucN;
					}

					numDegenI=0;
					numDegenJ=0;
					numDegenPair=0;
					degenSumI=0;
					degenSumJ=0;
					degenSumPair=0;
					for (ni=0; ni<NUM_NUCS; ni++) {
						if (abc->degen[nuci][ni]) {
							degenSumI += get_standard_singlet_freq(ni);
							numDegenI++;
							
							for (nj=0; nj<NUM_NUCS; nj++) {
								if (abc->degen[nucj][nj]) {
									numDegenPair++;
									degenSumPair += get_standard_pair_freq(ni,nj);
								}
							}
						}
					}
					for (nj=0; nj<NUM_NUCS; nj++) {
						if (abc->degen[nucj][nj]) {
							degenSumJ += get_standard_singlet_freq(nj);
							numDegenJ++;
						}
					}
					
					double probPair=degenSumPair/(double)(numDegenPair);
					double probI=degenSumI/(double)(numDegenI);
					double probJ=degenSumJ/(double)(numDegenJ);
					double thisOddsRatio=probPair/(probI*probJ);
					
					double thisLod=log(thisOddsRatio);
					double thisLodWeighted=thisLod*msa->wgt[k];

					sum_lod += thisLodWeighted;
					sum_weight += msa->wgt[k];
				}
			}

			if (sum_weight==0) {
				lod[j+1][i+1]=INT_MINSCORE;
			}
			else {
				
				double avgLod=sum_lod/sum_weight; // BTW, this is just because the weights don't necessarily sum to 1, so we're normalizing for that after the fact.  And they certainly don't sum to 1 if some seqs were skipped because of truncation, or double gaps
				
				lod[j+1][i+1] = 0.5 * 144.269504 * avgLod;
				if (lod[j+1][i+1] < INT_MINSCORE) {
					lod[j+1][i+1] = INT_MINSCORE;
				}
			}
			if (cmfinderVars->useIntsForMiAndPrLikeOldCmfinder) {
				// ignore this, since we're doing a function that wasn't in the old CMfinder anyway
			}
			if (0) {
				printf("prxy  , i=%d , j=%d , lod=%lg\n",i,j,lod[j+1][i+1]);
			}
		}
	}

	esl_stopwatch_Stop(timer_pf);
	esl_stopwatch_Include(cmfinderVars->timer_partition_function,timer_pf);
	esl_stopwatch_Destroy(timer_pf);

	*ret_prxy = lod;
}

double* bppr_seq(char *seq)
{
	int i,j;
	int size;
	double* bp_pr;
	char* structure;
	int len = (int)(strlen(seq));
	if (len == 0) return NULL;
	size = TriIndex(len,len-1) - 1;
	bp_pr = (double*) MallocOrDie(sizeof(double)*size);
	memset(bp_pr, 0, sizeof(double) * size);
	structure = (char *) malloc( sizeof(char) * (len + 1));
	init_pf_fold(len);
	pf_fold(seq, structure);
	for(j = 1; j < len; j++) {
		for(i = 0; i < j; i++){
			bp_pr[TriIndex(i,j)] += pr[iindx[i + 1] - (j + 1)];
		}
	}
	free(structure);
	free_pf_arrays();
	return (bp_pr);
}

void
Avg_bppr(ESL_MSA *text_msa,
	 double  **bp_pr,
	 int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col,
	 double ***ret_avg_bppr,
	 double **ret_avg_any_bppr)        
{
	double **avg_bppr;
	double *avg_any_bppr;
	double *sum_any_weight;
	double **sum_weight;
	char   *nogap_seq;
	int    *idx_map;
	int    *gotAnyBpprAtThisPosition;
	int    k,i,j, i1, j1;
	int count_this_weight;

	if ((text_msa->flags & eslMSA_DIGITAL)!=0) {
		esl_fatal("msa should be text");
	}

	//SaveMsa("text_msa.sto",text_msa);

	avg_any_bppr=(double *)MallocOrDie(sizeof(double)*(text_msa->alen+1));
	sum_any_weight=(double *)MallocOrDie(sizeof(double)*(text_msa->alen+1));
	gotAnyBpprAtThisPosition=(int *)MallocOrDie(sizeof(int)*(text_msa->alen+1));
	avg_bppr = DoubleAlloc2DArray(text_msa->alen+1);
	sum_weight=DoubleAlloc2DArray(text_msa->alen+1);
	
	if (bp_pr!=NULL) {
		esl_fatal("sorry, I'm not sure this works any more, and there's a risk of memory leak.  Must check more carefully.");
	}

	for(j=1; j < text_msa->alen+1; j++) {
		avg_any_bppr[j]=0;
		sum_any_weight[j]=0;
		for(i=0; i <j ; i++) {
			avg_bppr[j][i] = 0;
			sum_weight[j][i] = 0;
		}
	}

	if (bp_pr==NULL) {
		bp_pr = (double**)MallocOrDie(sizeof(double*) * text_msa->nseq);
		memset(bp_pr, 0, sizeof(double**) * text_msa->nseq);
	}
	for (k = 0; k < text_msa->nseq; k++) {
		nogap_seq = remove_gap(text_msa->aseq[k], &idx_map);
		if (bp_pr[k]==NULL) {
			bp_pr[k]= bppr_seq(nogap_seq);
		}
		if (bp_pr[k] == NULL) continue;
		for (j=0; j < text_msa->alen; j++) {
			gotAnyBpprAtThisPosition[j+1]=0;
		}
		for(j = 1; j < text_msa->alen; j++) {
			for(i = 0; i < j; i++) {
				i1 = idx_map[i];
				j1 = idx_map[j];
				
				/* in original CMfinder code, gaps are counted against pairing, even double-gaps
				 * the truncation-aware code preserves this property, except
				 * that we don't count base pairs that are truncated
				 */
				count_this_weight=1;
				if (use_fragmentary) {
					if (seq_fragment_start_col[k]!=-1 && i<seq_fragment_start_col[k]) {  /* i<j, so if j is left-truncated, then i will be also */
						assert(esl_abc_CIsGap(abc,text_msa->aseq[k][i])); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
						count_this_weight=0;
					}
					if (seq_fragment_end_col[k]!=-1 && j>seq_fragment_end_col[k]) {
						assert(esl_abc_CIsGap(abc,text_msa->aseq[k][j])); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
						count_this_weight=0;
					}
				}
				if (count_this_weight) {
					sum_weight[j+1][i+1] += text_msa->wgt[k];
					gotAnyBpprAtThisPosition[i+1]=1;
					gotAnyBpprAtThisPosition[j+1]=1;
				}
				if (i1 >=0 && j1 >= 0) {
					double add=bp_pr[k][TriIndex(i1,j1)] * text_msa->wgt[k];
					assert(count_this_weight); // if we're not counting the weight, why are we counting it in avg_bppr?
					avg_bppr[j+1][i+1] += add;
					
					avg_any_bppr[i+1] += add;
					avg_any_bppr[j+1] += add;
                                        if (0) {
                                          printf("seq=%d : i,j=%d,%d (i1,j1=%d,%d), add=%lg , weight=%lg\n",k,i,j,i1,j1,add,text_msa->wgt[k]);
                                        }
				}
			}
		}
		
		for (j=0; j < text_msa->alen; j++) {
			if (gotAnyBpprAtThisPosition[j+1]) {
				sum_any_weight[j+1] += text_msa->wgt[k];
			}
		}

		free(nogap_seq);
		free(idx_map);
	}
	
	/* divide avg_bppr by sum of weights */
	for(j = 1; j < text_msa->alen; j++) {
		for(i = 0; i < j; i++) {
			if (sum_weight[j+1][i+1]==0) {
				avg_bppr[j+1][i+1]=0; /* bad evidence for pairing, doesn't matter.  I believe this is implicit from above */
			}
			else {
				avg_bppr[j+1][i+1] /= sum_weight[j+1][i+1];
				if (0) {
					printf("frag=%d , i=%d , j=%d , avg_bppr=%lg\n",use_fragmentary,i,j,avg_bppr[j+1][i+1]);
				}
			}
		}
	}

	if (use_fragmentary) {
		for(i=0; i < text_msa->alen; i++) {
			if (avg_any_bppr[i+1]==0) {
				/* okay, just leave it at zero */
			}
			else {
				avg_any_bppr[i+1] /= sum_any_weight[i+1];
			}
		}
	}
	else {
		for(i=0; i < text_msa->alen; i++) {
			avg_any_bppr[i+1]=0;
			for(j = 0; j < text_msa->alen; j++) {
				if (i != j) {
					double add=i<j ? avg_bppr[j+1][i+1] : avg_bppr[i+1][j+1];
					avg_any_bppr[i+1] += add;
					if (0 && (i==11 || i==79)) {
						printf("add, i=%d , j=%d , add=%lg\n",i,j,add);
					}
				}
			}
		}
	}

	/* free bp_pr, which for now I assume is allocated here.  the strategy for saving partition function calculations will make this trickier */
	for (k = 0; k < text_msa->nseq; k++) {
		if (bp_pr[k]!=NULL) {
			free(bp_pr[k]);
		}
	}
	free(bp_pr);
	
	Free2DArray((void **)sum_weight,text_msa->alen+1);
	free(sum_any_weight);
	free(gotAnyBpprAtThisPosition);

	*ret_avg_bppr = avg_bppr;
	*ret_avg_any_bppr = avg_any_bppr;
}

/*
Calculate base-pairing probabilities according to the average partition function
(which is calculated by the Avg_bppr function).
Then calculate LOD scores for favoring base pairing, based on the partition function
evidence.  Converted to the same scale as the mixy function.
*/
void
prxy (ESL_MSA *text_msa,
	  double **bp_pr,
	  CmfinderVars *cmfinderVars,
	  int *seq_fragment_start_col,int *seq_fragment_end_col,
	  double  ***ret_prxy)        /* RETURN: prxy array           */
{
	int      i,j;
	double **lod;
	double  *px;
	double **pxy;

	if ((text_msa->flags & eslMSA_DIGITAL)!=0) {
		esl_fatal("msa should be text");
	}

	ESL_STOPWATCH *timer_pf=esl_stopwatch_Create();

	esl_stopwatch_Start(timer_pf);

	Avg_bppr(text_msa, bp_pr, cmfinderVars->use_fragmentary && !cmfinderVars->nonFragmentaryAvgBppr, seq_fragment_start_col, seq_fragment_end_col, &pxy, &px);

	lod = DoubleAlloc2DArray(text_msa->alen+1);

	for(j = 1; j < text_msa->alen; j++) {
		for(i = 0; i < j; i++) {
			assert(px[i+1]>=0 && px[i+1]<=1);
			assert(px[j+1]>=0 && px[j+1]<=1);
			assert(pxy[j+1][i+1]>=0 && pxy[j+1][i+1]<=1);
			if (pxy[j+1][i+1]==0) { /* doesn't really matter what px[i+1] or px[j+1] are, we shouldn't select this as a pair anyway */
				lod[j+1][i+1]=INT_MINSCORE;
			}
			else {
				lod[j+1][i+1] = 0.5 * 144.269504 * log(pxy[j+1][i+1] /( (1-px[i+1]) * (1- px[j+1])));
				if (lod[j+1][i+1] < INT_MINSCORE) {
					lod[j+1][i+1] = INT_MINSCORE;
				}
			}
			if (cmfinderVars->useIntsForMiAndPrLikeOldCmfinder) {
				lod[j+1][i+1]=(int)(lod[j+1][i+1]);
			}
			if (0) {
				printf("prxy  , i=%d , j=%d , lod=%lg , pxy=%lg , px=%lg,%lg\n",i,j,lod[j+1][i+1],pxy[j+1][i+1],px[i+1],px[j+1]);
			}
		}
	}

	free(px);
	Free2DArray((void **)pxy, text_msa->alen+1);
	*ret_prxy = lod;
	
	esl_stopwatch_Stop(timer_pf);
	esl_stopwatch_Include(cmfinderVars->timer_partition_function,timer_pf);
	esl_stopwatch_Destroy(timer_pf);
}

/* Function: mixy()
*
* Originally from COVE package, and copied from CMfinder 0.3.
*
* Purpose:  given a set of N aligned sequences aseq, calculate
*           pairwise covariances (mutual information). ret_mixy
*           is allocated, filled, and returned, as a diagonal 2D
*           (NxN) matrix of values. It must be freed by
*           the caller. It is a lower diagonal matrix mxy[j][i],
*           j > i, 0..alen-1 by 0..j-1.
*
*           The values are the average
*           secondary structure information content (i.e. weighted for
*           the number of pairs actually occurring in columns i,j)
*           in bits, to two decimal places (i.e. info*100).
*           (In COVE and CMfinder 0.3, they were integers.)
* 
*           For truncation/fragmentary: gaps are already ignored, so this deals with truncation properly already (we want to ignore truncated nucleotides, and truncated nucleotides would be gaps, so they're already ignored)
*
* Returns:  mxy, which must be free'd by caller with free_mixy().
*/
void
mixy(ESL_MSA  *digital_msa,
     CmfinderVars *cmfinderVars,
	 double ***ret_mxy)
{
	double **mxy;                 /* RETURN: diagonal covariance matrix  */
	double   fx[NUM_NUCS];        /* singlet frequency vector        */
	double   fy[NUM_NUCS];	/* another singlet frequency vector    */
	double   fxy[NUM_NUCS][NUM_NUCS]; /* pairwise frequency 2D array  */
	int     idx;			/* counter for sequences               */
	int     i, j;			/* counters for columns x,y            */
	int     symi, symj;		/* counters for symbols                */
	double  pairs;		/* counter for pairs in which there are no gaps */
	
	if ((digital_msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}

	mxy = DoubleAlloc2DArray(digital_msa->alen + 1);

	/* calculate mxy
	*/
	for (j = 2; j <= digital_msa->alen; j++) {
		for (i = 1; i < j; i++) {
			/* zero counter array */
			for (symj = 0; symj < NUM_NUCS; symj++) {
				fx[symj] = fy[symj] = 0.0;
				for (symi = 0; symi < NUM_NUCS; symi++) {
					fxy[symj][symi] = 0.0;
				}
			}
			/* count symbols in a column */
			pairs = 0;
			for (idx = 0; idx < digital_msa->nseq; idx++) {
				symi = digital_msa->ax[idx][i];
				symj = digital_msa->ax[idx][j];
				/* Gaps are ignored */
				if (esl_abc_XIsGap(abc,symi) || esl_abc_XIsGap(abc,symj)) {/* if either nuc is a gap character */
					continue;
				}
				if (!esl_abc_XIsCanonical(abc,symi) || !esl_abc_XIsCanonical(abc,symj)) {
					continue; /* skip non-canonicals too */
				}
				fx[symi] += digital_msa->wgt[idx];
				fy[symj] += digital_msa->wgt[idx];
				fxy[symi][symj] += digital_msa->wgt[idx];
				pairs += digital_msa->wgt[idx];
			}

			/* convert to frequencies */
			if (pairs > 0) {
				for (symi = 0; symi < NUM_NUCS; symi++) {
					fx[symi] /=  pairs;
					fy[symi] /=  pairs;
					for (symj = 0; symj < NUM_NUCS; symj++) {
						fxy[symi][symj] /=  pairs;
					}
				}
			}
			else {
				mxy[j][i] =  INT_NEGINFINITY;
				continue;
			}
			/* calculate mxy. 144.269504 is a conversion of ln's into
			* bits * 100: i.e. 100 * (1/log(2))
			*/
			mxy[j][i] = 0;
			for (symi = 0; symi < NUM_NUCS; symi++) {
				for (symj = 0; symj < NUM_NUCS; symj++) {
					if (fxy[symi][symj] > 0.0) {
						double f=144.269504 * fxy[symi][symj] * log((fxy[symi][symj] / (fx[symi] * fy[symj])));
						if (cmfinderVars->useIntsForMiAndPrLikeOldCmfinder) {
							mxy[j][i] += (int)f;
						}
						else {
							mxy[j][i] += f;
						}
					}
				}
			}

			/* Sat Jul 17 22:17:17 1993:  We weight by pairs to get an expected score
			* over all the sequences. Fixes a problem that columns with few symbols
			* could dominate the calculation just because of noise.
			*/
			mxy[j][i] =  (mxy[j][i] * pairs) / (double)(digital_msa->nseq);
		}
	}

	*ret_mxy = mxy;
}

/* This function is copied from summarize.c in CMfinder 0.3
 * It has some differences with the above 'mixy' function (which was from autommaker.c in CMfinder 0.3)
 * and I've decided not to try to merge the functions.
 * 
 * The differences are:
 * (1) summarize_mxy loops only over pairs in the SS_cons line, while mixy loops over all pairs of columns.
 * (2) summarize_mxy marginalizes DIGITAL_GAP symbols over all possible base pairs, while mixy ignores
 * base pairs that have a gap.
 * (3) summarize_mxy has a nullModel for the prior distribution of the 4 nucs, but mixy does not.  However,
 * the nullModel is only used for marginalizing DIGITAL_GAP.  (I now also use it for marginalizing other degenerate
 * nucs.)
 * (4) mixy does a correction (see end of function) by multiplying mxy by weight-of-pairs-without-gap / num-seqs.
 * summarize_mxy does not do this.
 * (5) mixy uses int's to store the mxy table, while summarize_mxy uses floats.
 * (6) I converted summarize_mxy to just using the 'msa', while I've left mixy in its older form.
 * (7) summarize_mxy just returns the total mutual information, while mixy returns the mutual info for each pair of columns.
 */
double summarize_mxy(ESL_MSA *msa, double *nullModel, double tot_weight,int *pt)
{
	int     i, j, idx;
	double  fx[NUM_NUCS];        /* singlet frequency vector            */
	double  fy[NUM_NUCS];	/* another singlet frequency vector    */
	double  fxy[NUM_NUCS][NUM_NUCS]; /* pairwise frequency 2D array    */
	int     symi, symj;		/* counters for symbols                */
	double  total_mxy=0;
	double  mxy;

	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}

	for (i = 0; i < msa->alen; i++) {
		if (pt[i] <= i) continue;
		j = pt[i];
		
		/* zero counter array */
		for (symj = 0; symj < NUM_NUCS; symj++) {
			fx[symj] = fy[symj] = 0.0;
			for (symi = 0; symi < NUM_NUCS; symi++) {
				fxy[symj][symi] = 0.0;
			}
		}
		
		/* count symbols in a column */
		for (idx = 0; idx < msa->nseq; idx++) {
			double fx_norm,fy_norm,fxy_norm;
			int l,r;
			int nuci=msa->ax[idx][i+1];
			int nucj=msa->ax[idx][j+1];
			if (esl_abc_XIsGap(abc,nuci)) {
				nuci=esl_abc_DigitizeSymbol(abc,'N'); /* this is really what the original summarize.c code was doing: treating a gap as an 'N', and marginalizing over the 4 nucs */
			}
			if (esl_abc_XIsGap(abc,nucj)) {
				nucj=esl_abc_DigitizeSymbol(abc,'N');
			}
			
			fx_norm=0;
			fy_norm=0;
			fxy_norm=0;
			for (l=0; l<NUM_NUCS; l++) {
				if (!(abc->degen[nuci][l])) {
					continue;
				}
				fx_norm += nullModel[l];
			}
			for (r=0; r<NUM_NUCS; r++) {
				if (!(abc->degen[nucj][r])) {
					continue;
				}
				fy_norm += nullModel[r];
			}
			/* fxy_norm=fx_norm * fy_norm, but rather than making sure, I'd rather just do the loop again */
			for (l=0; l<NUM_NUCS; l++) {
				if (!(abc->degen[nuci][l])) {
					continue;
				}
				for (r=0; r<NUM_NUCS; r++) {
					if (!(abc->degen[nucj][r])) {
						continue;
					}
					fxy_norm += nullModel[l]*nullModel[r];
				}
			}
			
			for (l=0; l<NUM_NUCS; l++) {
				if (!(abc->degen[nuci][l])) {
					continue;
				}
				fx[l] += msa->wgt[idx]*nullModel[l]/fx_norm;
			}
			for (r=0; r<NUM_NUCS; r++) {
				if (!(abc->degen[nucj][r])) {
					continue;
				}
				fy[r] += msa->wgt[idx]*nullModel[r]/fy_norm;
			}
			for (l=0; l<NUM_NUCS; l++) {
				if (!(abc->degen[nuci][l])) {
					continue;
				}
				for (r=0; r<NUM_NUCS; r++) {
					if (!(abc->degen[nucj][r])) {
						continue;
					}
					fxy[l][r] += msa->wgt[idx] * nullModel[l]*nullModel[r]/fxy_norm;
				}
			}
		}
		
		/* convert to frequencies */
		for (symi = 0; symi < NUM_NUCS; symi++) {
			fx[symi] /=  tot_weight;
			fy[symi] /=  tot_weight;
			for (symj = 0; symj < NUM_NUCS; symj++) {
				fxy[symj][symi] /=  tot_weight;
			}
		}
		
		/* calculate mxy. 144.269504 is a conversion of ln's into
		* bits * 100: i.e. 100 * (1/log(2)) 
		*/
		mxy = 0;
		for (symi = 0; symi < NUM_NUCS; symi++) {
			for (symj = 0; symj < NUM_NUCS; symj++) {
				if (fxy[symi][symj] > 0.0) {
					mxy +=  1.44269504 * fxy[symi][symj] * log((fxy[symi][symj] / (fx[symi] * fy[symj])));
				}
			}
		}
		if (mxy < -0.00001) {
			esl_fatal("Error ! Column %d  %d mxy = %f", i, j, mxy);
		}    
		total_mxy += mxy;    
	}  
	return total_mxy;  
}

/* combine mutual information and the prior using partition function 
 * xy=combined info, mxy=mutual info, pxy=partition-function-based info
 * xy=mxy+pxy, except xy=0 if pxy=0 (i.e., if the partition function doesn't favor the pair, we never favor it)
 */
void merge(double **mxy, double **pxy, int alen, double ***ret_xy)
{
	double **xy;
	int i,j;

	xy = DoubleAlloc2DArray(alen+1);

	for(j=2; j <= alen; j++) {
		for(i=1; i < j; i++) {
			xy[j][i] = mxy[j][i] + pxy[j][i] ;
			if (pxy[j][i] < 0 ) {
				xy[j][i] = 0;
			}
			if (0) {
				printf("merge , i=%d , j=%d , xy=%lg\n",i-1,j-1,xy[j][i]);
			}
		}
	}
	
	*ret_xy = xy;
}

/* stuff from summarize.c in CMfinder 0.3 */

double entropy(ESL_MSA *msa, double *nullModel, double tot_weight)
{
	int    i, idx,sym;
	double total_e, e;
	double f[NUM_NUCS];        /* singlet frequency vector            */
	double one_over_ln_2=1.44269504;

	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}

	total_e = 0;  
	for (i = 0; i < msa->alen; i++) {    
		for (sym = 0; sym <NUM_NUCS; sym++) {
			f[sym] = 0;
		}

		for (idx = 0; idx < msa->nseq; idx++) {
			int thisNuc=msa->ax[idx][i+1];
			if (esl_abc_XIsGap(abc,thisNuc)) {
				for (sym = 0; sym < NUM_NUCS; sym++) {
					f[sym] += nullModel[sym] * msa->wgt[idx];
				}
			}
			else {
				/* handles both canonical and degenerate cases */
				if (esl_abc_DCount(abc,f,thisNuc,msa->wgt[idx])!=eslOK) { /* importantly, we already know it's not a gap */
					esl_fatal("esl_abc_DCount failed");
				}
			}
		}
		e = 0;
		for (sym = 0; sym < NUM_NUCS; sym++) {
			f[sym] /= tot_weight;
			e -=  nullModel[sym]*log(nullModel[sym]);
			if (f[sym] > 0) {
				e +=  f[sym] * log(f[sym]);
			}
		}
		total_e += e *  one_over_ln_2;
	}
	return total_e;
}

double weighted_base_pair (ESL_MSA* msa,int *pt,double **bp_pr)
{
	int    i;
	double total_bp = 0;
	for (i=0; i < msa->alen; i++) {
		if (pt[i] > i) {
			total_bp += bp_pr[pt[i]+1][i+1];
		}
	}
	return total_bp;
}

double average_base_pair(ESL_MSA* msa,int *pt)
{
	int    i;
	double total_bp = 0;
	for (i=0; i < msa->alen; i++) {
		if (pt[i] > i) {
			total_bp ++;
		}
	}
	return total_bp;
}

double average_score(ESL_MSA *msa, float tot_weight)
{
	float total_score = 0;
	float score;
	int i;
	int start, end;
	for (i=0; i < msa->nseq; i++) {
		if (msa->sqdesc && msa->sqdesc[i]) {
			if (sscanf(msa->sqdesc[i], "%d..%d\t%f", &start, &end, &score)<3) {
				esl_fatal("#=GR ... DE tag is not in the correct format for sequence#%d",i);
			}
			total_score += score * msa->wgt[i];
		}
		else {
			esl_fatal("no #=GR ... DE tag for sequence #%d",i);
		}
	}
	return total_score / tot_weight;
}

double average_seq_id (ESL_MSA *msa, float tot_weight,int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col)
{
	int i, j, k;
	double total_id = 0;
	double id = 0;
	int    len;
	double pair = 0;
	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}
	if (msa->nseq==1) {
		// this is a weird case, because this is avg _pairwise_ seq_id, so if there's only 1 seq, then there's no pairwise comparisons.  but I think 1.0 is as reasonable an answer as any, and anyway the MSA should be rejected by other things in the pipeline
		return 1.0;
	}
	for(i=1; i < msa->nseq; i++) {
		for(j=0; j < i; j++) {
			id = 0;
			len=0;
			for (k=0; k < msa->alen; k++) {
				int nuci=msa->ax[i][k+1];
				int nucj=msa->ax[j][k+1];
				if (esl_abc_XIsGap(abc,nuci) && esl_abc_XIsGap(abc,nucj)) {
					continue;
				}
				if (esl_abc_XIsDegenerate(abc,nuci) || esl_abc_XIsDegenerate(abc,nucj)) { // if either is degen, then just skip it
					continue;
				}
				if (use_fragmentary) {
					/* ignore if either nuci or nucj is truncated by fragmentary */
					if (seq_fragment_start_col[i]!=-1 && k<seq_fragment_start_col[i]) {
						assert(esl_abc_XIsGap(abc,nuci)); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
						continue;
					}
					if (seq_fragment_start_col[j]!=-1 && k<seq_fragment_start_col[j]) {
						assert(esl_abc_XIsGap(abc,nucj)); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
						continue;
					}
					if (seq_fragment_end_col[i]!=-1 && k>seq_fragment_end_col[i]) {
						assert(esl_abc_XIsGap(abc,nuci)); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
						continue;
					}
					if (seq_fragment_end_col[j]!=-1 && k>seq_fragment_end_col[j]) {
						assert(esl_abc_XIsGap(abc,nucj)); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
						continue;
					}
				}
				len++;
				if (nuci==nucj) {
					id++;
				}
			}
			if (len<10) {
				// on the off chance there's very little overlap, skip it
			}
			else {
				total_id += (id / len) * msa->wgt[i] * msa->wgt[j];
				pair += msa->wgt[i] * msa->wgt[j];
			}
		}
	}
	assert(pair>0);
	return total_id / pair;
}

double average_seq_len(ESL_MSA *msa, float tot_weight)
{
	double total_len=0;
	int    i,j;
	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}
	for (i=0; i < msa->nseq; i++) {
		int length=0;
		for (j=0; j < msa->alen; j++) {
			int nuc=msa->ax[i][j+1];
			if (!esl_abc_XIsGap(abc,nuc)) {
				length++;
			}
		}
		total_len += length * msa->wgt[i];
	}
	return total_len / tot_weight;
}

double average_energy(ESL_MSA  *msa, float tot_weight)
{
	char* ss;
	char* seq;
	char  *sp1, *sp2;
	int i,j;
	double energy;
	double tot_energy=0;
	if ((msa->flags & eslMSA_DIGITAL)!=0) {
		esl_fatal("msa should be text");
	}
	seq = (char *)malloc(sizeof(char) * (msa->alen + 1));
	ss = (char *)malloc(sizeof(char) * (msa->alen + 1));
	for(i=0; i < msa->nseq; i++) {
		sp1 = seq;
		sp2 = ss;
		for(j=0; j < msa->alen; j++) {
			if (!esl_abc_CIsGap(abc,msa->aseq[i][j])) {
				*(sp1++)= msa->aseq[i][j];
				*sp2= msa->ss[i][j];
				if (pair_left(*sp2)) *sp2= '('; /* Vienna likes parentheses only */
				if (pair_right(*sp2)) *sp2= ')';
				sp2++;
			}
		}
		*sp1=*sp2='\0';
		initialize_fold(strlen(seq));
		energy = energy_of_struct(seq, ss);
		// Prevent the outliers
		if (energy > 10) energy = 10;
		tot_energy += energy * msa->wgt[i];
	}
	free(seq);
	free(ss);
	return tot_energy / tot_weight;
}

double average_GC (ESL_MSA *msa, float tot_weight)
{
	double total_GC = 0;
	int    i,j;
	int    GC_count;
	int    count;
	double GC;
	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}
	for(i=0; i < msa->nseq; i++) {
		GC_count=0;
		count = 0;
		for (j=0; j < msa->alen; j++) {
			int nuc=msa->ax[i][j+1];
			if (esl_abc_XIsGap(abc,nuc)) {
				continue;
			}
			if (nuc==1 || nuc==2) { /* G or C */
				GC_count++;
			}
			count++;
		}
		GC =  (double)(GC_count) / (double)(count);
		total_GC += GC * msa->wgt[i];
	}
	return (total_GC / tot_weight);
}

int conserved_position (ESL_MSA *msa,float tot_weight,int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col)
{
	const int min_block_size = 4;
	const double gap_threshold = 0.6;
	const double cons_threshold = 0.7;
	int i,j;
	int block_size=0;
	int accumulated_block_size = 0;
	double freq[NUM_NUCS]={0,0,0,0};
	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}
	for (i=1; i <= msa->alen; i++) {
		int conserved = 0;
		int  gap_count = 0;
		double counted_tot_weight=0;
		double truncated_weight=0;
		double denom_weight;
		for (j=0; j < NUM_NUCS; j++) {
			freq[j]=0;
		}
		for (j=0; j < msa->nseq; j++) {
			int nuc=msa->ax[j][i];
			if (use_fragmentary) {
				/* ignore if either nuc is truncated by fragmentary */
				if (seq_fragment_start_col[j]!=-1 && i-1<seq_fragment_start_col[j]) {
					assert(esl_abc_XIsGap(abc,nuc));
					truncated_weight += msa->wgt[j];
					continue;
				}
				if (seq_fragment_end_col[j]!=-1 && i-1>seq_fragment_end_col[j]) {
					assert(esl_abc_XIsGap(abc,nuc));
					truncated_weight += msa->wgt[j];
					continue;
				}
			}
			if (esl_abc_XIsGap(abc,nuc)) {
				gap_count += msa->wgt[j];
			}
			else {
				if (esl_abc_XIsDegenerate(abc,nuc)) {
					// actually, let's skip degenerates
					continue;
				}
				/* handles both canonical and degenerate cases */
				if (esl_abc_DCount(abc,freq,nuc,msa->wgt[j])!=eslOK) { /* importantly, we already know it's not a gap */
					esl_fatal("esl_abc_DCount failed");
				}
			}
			counted_tot_weight += msa->wgt[j];
		}
		if (use_fragmentary) {
			denom_weight=counted_tot_weight;
			if (counted_tot_weight==0) {
				continue;
			}
			if ((gap_count+truncated_weight) / (denom_weight+truncated_weight) > gap_threshold) {
				/* if there's lots of truncated seqs, then that can artificially inflate conservation.  let's semi-treat those as gaps to guard against this. */
				continue;
			}
		}
		else {
			denom_weight=tot_weight;
		}
		if (gap_count / denom_weight > gap_threshold) {
			continue;
		}
		for (j=0; j < NUM_NUCS; j++) {
			if (freq[j] / denom_weight > cons_threshold && freq[j] >= 3) {
				block_size ++;
				conserved = 1;
			}
		}
		if (!conserved) {
			if (block_size >= min_block_size) {
				accumulated_block_size += block_size;
			}
			block_size = 0;
		}
	}
	if (block_size >= min_block_size) {
		accumulated_block_size += block_size;
	}
	return accumulated_block_size;
}

void bad_base_pair (ESL_MSA* msa, int *pt, double **bp_pr, double* ret_conflict_bp, double* ret_del_bp)
{
	int i,j;
	double conflict_bp = 0;
	double del_bp = 0;
	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}
	for (i=0; i < msa->alen; i++) {
		if (pt[i] > i) { /* only look at base-paired columsn (pt[i]>0) and only count each base pair once (pt[i]>i>=0) */
			for (j=0; j < msa->nseq; j++) {
				int l=msa->ax[j][i+1];
				int r=msa->ax[j][pt[i]+1];
				if (!IsCanonicalDigitalBpAndNotDegen(l,r)) {
					if (esl_abc_XIsGap(abc,l) && esl_abc_XIsGap(abc,r)) {
						del_bp += msa->wgt[j] * bp_pr[pt[i]+1][i+1];
					}
					else {
						conflict_bp += msa->wgt[j] * bp_pr[pt[i]+1][i+1];
					}
				}
			}
		}
	}
	*ret_conflict_bp=conflict_bp;
	*ret_del_bp=del_bp;
}

/* Adapted from esl_dst_XDiffMx in esl_distance.c, modified just to do truncation
 */
int
esl_dst_XPairId_fragmentary(const ESL_ALPHABET *abc, const ESL_DSQ *ax1, const ESL_DSQ *ax2, 
		double *opt_distance, int *opt_nid, int *opt_n, int use_fragmentary,int seq_fragment_start_col1,int seq_fragment_end_col1,int seq_fragment_start_col2,int seq_fragment_end_col2)
{
	int     status;
	int     idents;               /* total identical positions  */
	int     len1, len2;           /* lengths of seqs            */
	int     i;                    /* position in aligned seqs   */

	idents = len1 = len2 = 0;
	for (i = 1; ax1[i] != eslDSQ_SENTINEL && ax2[i] != eslDSQ_SENTINEL; i++) {
	
		if (use_fragmentary) {
			if (seq_fragment_start_col1!=-1 && i-1<seq_fragment_start_col1) {
				assert(esl_abc_XIsGap(abc,ax1[i])); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
				continue;
			}
			if (seq_fragment_end_col1!=-1 && i-1>seq_fragment_end_col1) {
				assert(esl_abc_XIsGap(abc,ax1[i]));
				continue;
			}
			if (seq_fragment_start_col2!=-1 && i-1<seq_fragment_start_col2) {
				assert(esl_abc_XIsGap(abc,ax2[i])); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
				continue;
			}
			if (seq_fragment_end_col2!=-1 && i-1>seq_fragment_end_col2) {
				assert(esl_abc_XIsGap(abc,ax2[i]));
				continue;
			}
		}
		
		if (esl_abc_XIsCanonical(abc, ax1[i])) len1++;
		if (esl_abc_XIsCanonical(abc, ax2[i])) len2++;

		if (esl_abc_XIsCanonical(abc, ax1[i]) && esl_abc_XIsCanonical(abc, ax2[i]) && ax1[i] == ax2[i])
			idents++;
	}
	if (len2 < len1) len1 = len2;

	if (ax1[i] != eslDSQ_SENTINEL || ax2[i] != eslDSQ_SENTINEL)
		ESL_XEXCEPTION(eslEINVAL, "strings not same length, not aligned");

	if (opt_distance != NULL)  *opt_distance = ( len1==0 ? 0. : (double) idents / (double) len1 );
	if (opt_nid      != NULL)  *opt_nid      = idents;
	if (opt_n        != NULL)  *opt_n        = len1;
	return eslOK;

	ERROR:
	if (opt_distance != NULL)  *opt_distance = 0.;
	if (opt_nid      != NULL)  *opt_nid      = 0;
	if (opt_n        != NULL)  *opt_n        = 0;
	return status;
}

/* Adapted from esl_dst_XDiffMx in esl_distance.c, modified just to do truncation
 */
int
esl_dst_XPairIdMx_fragmentary(const ESL_ALPHABET *abc,  ESL_DSQ **ax, int N, ESL_DMATRIX **ret_S, int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col)
{
	int status;
	ESL_DMATRIX *S = NULL;
	int i,j;

	if (( S = esl_dmatrix_Create(N,N) ) == NULL) goto ERROR;

	for (i = 0; i < N; i++) {
		S->mx[i][i] = 1.;
		for (j = i+1; j < N; j++) {
			status = esl_dst_XPairId_fragmentary(abc, ax[i], ax[j], &(S->mx[i][j]), NULL, NULL, use_fragmentary,seq_fragment_start_col[i],seq_fragment_end_col[i],seq_fragment_start_col[j],seq_fragment_end_col[j]);
			if (status != eslOK)
				ESL_XEXCEPTION(status, "Pairwise identity calculation failed at seqs %d,%d\n", i,j);
			S->mx[j][i] =  S->mx[i][j];
		}
	}
	if (ret_S != NULL) *ret_S = S; else esl_dmatrix_Destroy(S);
	return eslOK;

ERROR:
	if (S     != NULL)  esl_dmatrix_Destroy(S);
	if (ret_S != NULL) *ret_S = NULL;
	return status;
}

/* Adapted from esl_dst_XDiffMx in esl_distance.c, modified just to do truncation
 */
int
esl_dst_XDiffMx_fragmentary(const ESL_ALPHABET *abc, ESL_DSQ **ax, int N, ESL_DMATRIX **ret_D, int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col)
{
	int status;
	ESL_DMATRIX *D = NULL;
	int i,j;

	status = esl_dst_XPairIdMx_fragmentary(abc, ax, N, &D, use_fragmentary,seq_fragment_start_col,seq_fragment_end_col);
	if (status != eslOK) goto ERROR;

	for (i = 0; i < N; i++) {
		D->mx[i][i] = 0.;
		for (j = i+1; j < N; j++) {
			D->mx[i][j] = 1. - D->mx[i][j];
			D->mx[j][i] = D->mx[i][j];
		}
	}
	if (ret_D != NULL) *ret_D = D; else esl_dmatrix_Destroy(D);
	return eslOK;

ERROR:
	if (D     != NULL)  esl_dmatrix_Destroy(D);
	if (ret_D != NULL) *ret_D = NULL;
	return status;
}

/* Adapted from esl_msaweight_GSC in esl_msaweight.c, modified just to do truncation
 */
int
esl_msaweight_GSC_fragmentary(ESL_MSA *msa,int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col)
{
	ESL_DMATRIX *D = NULL;     /* distance matrix */
	ESL_TREE    *T = NULL;     /* UPGMA tree */
	double      *x = NULL;     /* storage per node, 0..N-2 */
	double       lw, rw;       /* total branchlen on left, right subtrees */
	double       lx, rx;	     /* distribution of weight to left, right side */
	int i;		     /* counter over nodes */
	int status;

	/* Contract checks
	*/
	ESL_DASSERT1( (msa       != NULL) );
	ESL_DASSERT1( (msa->nseq >= 1)    );
	ESL_DASSERT1( (msa->alen >= 1)    );
	ESL_DASSERT1( (msa->wgt  != NULL) );
	if (msa->nseq == 1) { msa->wgt[0] = 1.0; return eslOK; }

	/* GSC weights use a rooted tree with "branch lengths" calculated by
	* UPGMA on a fractional difference matrix - pretty crude.
	*/
	if (! (msa->flags & eslMSA_DIGITAL)) {
		assert(0); /* text case not implemented */
	} 
	else {
		if ((status = esl_dst_XDiffMx_fragmentary(msa->abc, msa->ax, msa->nseq, &D, use_fragmentary,seq_fragment_start_col,seq_fragment_end_col)) != eslOK) goto ERROR;
	}

	/* oi, look out here.  UPGMA is correct, but old squid library uses
	* single linkage, so for regression tests ONLY, we use single link. 
	*/
#ifdef  eslMSAWEIGHT_REGRESSION
	if ((status = esl_tree_SingleLinkage(D, &T)) != eslOK) goto ERROR; 
#else
	if ((status = esl_tree_UPGMA(D, &T)) != eslOK) goto ERROR; 
#endif
	esl_tree_SetCladesizes(T);	

	ESL_ALLOC(x, sizeof(double) * (T->N-1));

	/* Postorder traverse (leaves to root) to calculate the total branch
	* length under each internal node; store this in x[].  Remember the
	* total branch length (x[0]) for a future sanity check.
	*/
	for (i = T->N-2; i >= 0; i--) {
		x[i] = T->ld[i] + T->rd[i];
		if (T->left[i]  > 0) x[i] += x[T->left[i]];
		if (T->right[i] > 0) x[i] += x[T->right[i]];
	}

	/* Preorder traverse (root to leaves) to calculate the weights.  Now
	* we use x[] to mean, the total weight *above* this node that we will
	* apportion to the node's left and right children. The two
	* meanings of x[] never cross: every x[] beneath x[i] is still a
	* total branch length.
	*
	* Because the API guarantees that msa is returned unmodified in case
	* of an exception, and we're touching msa->wgt here, no exceptions
	* may be thrown from now on in this function.
	*/
	x[0] = 0;			/* initialize: no branch to the root. */
	for (i = 0; i <= T->N-2; i++) {
		lw = T->ld[i];   if (T->left[i]  > 0) lw += x[T->left[i]];
		rw = T->rd[i];   if (T->right[i] > 0) rw += x[T->right[i]];

		if (lw+rw == 0.) {
			/* A special case arises in GSC weights when all branch lengths in a subtree are 0.
			* In this case, all seqs in this clade should get equal weights, sharing x[i] equally.
			* So, split x[i] in proportion to cladesize, not to branch weight.
			*/
			if (T->left[i] > 0)  lx =  x[i] * ((double) T->cladesize[T->left[i]]  / (double) T->cladesize[i]);
			else                 lx =  x[i] / (double) T->cladesize[i];

			if (T->right[i] > 0) rx =  x[i] * ((double) T->cladesize[T->right[i]] / (double) T->cladesize[i]);
			else                 rx =  x[i] / (double) T->cladesize[i];
		} 
		else {/* normal case: x[i] split in proportion to branch weight. */
			lx = x[i] * lw/(lw+rw);
			rx = x[i] * rw/(lw+rw);
		}

		if (T->left[i]  <= 0) msa->wgt[-(T->left[i])] = lx + T->ld[i];
		else                  x[T->left[i]] = lx + T->ld[i];

		if (T->right[i] <= 0) msa->wgt[-(T->right[i])] = rx + T->rd[i];
		else                  x[T->right[i]] = rx + T->rd[i];
	} 

	/* Renormalize weights to sum to N.
	*/
	esl_vec_DNorm(msa->wgt, msa->nseq);
	esl_vec_DScale(msa->wgt, msa->nseq, (double) msa->nseq);
	msa->flags |= eslMSA_HASWGTS;

	free(x);
	esl_tree_Destroy(T);
	esl_dmatrix_Destroy(D);
	return eslOK;

ERROR:
	if (x != NULL) free(x);
	if (T != NULL) esl_tree_Destroy(T);
	if (D != NULL) esl_dmatrix_Destroy(D);
	return status;
}

typedef struct {
	int leftNucsToMove,rightNucsToMove;
} ShiftDistalMispair_OneSeqAndTermLoop;
typedef struct {
	int pairLeft,pairRight;
	ShiftDistalMispair_OneSeqAndTermLoop *perSeqInfo;
} ShiftDistalMispairs_TermLoop;
void ShiftDistalMispairsIntoTerminalLoops_DoMoving(ESL_MSA *new_msa,int seqnum,int src,int dest,int nucsToMove,int dir)
{
	int numMoved=0;
	while (numMoved<nucsToMove) {
		assert(src>=0 && src<new_msa->alen);
		if (!esl_abc_XIsGap(abc,new_msa->ax[seqnum][src+1])) {
			assert(dest>=0 && dest<new_msa->alen);
			new_msa->ax[seqnum][dest+1]=new_msa->ax[seqnum][src+1];
			new_msa->ax[seqnum][src+1]=esl_abc_XGetGap(abc);
			dest += dir;
			numMoved++;
		}
		src += dir;
	}
}
ESL_MSA * ShiftDistalMispairsIntoTerminalLoops (ESL_MSA *msa,CmfinderVars *cmfinderVars,int *seq_fragment_start_col,int *seq_fragment_end_col)
{
	int col;
	int prevOpenPair;
	int tl;

	if (!cmfinderVars->shiftDistalMispairsIntoTerminalLoops) {
		return msa;
	}

	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}
	
#if 0
	SaveMsa("/c/scratch/oak.sto",msa);
#endif

	/* remove any zero-column term loops, which are ridiculous, and will also screw up this code */
	for (col=1; col<msa->alen; col++) {
		if (pair_left(msa->ss_cons[col-1]) && pair_right(msa->ss_cons[col])) {
			msa->ss_cons[col-1]='.';
			msa->ss_cons[col]='.';
		}
	}
	
	int numTermLoops=0;
	ShiftDistalMispairs_TermLoop *termLoopInfoVector=(ShiftDistalMispairs_TermLoop *)MallocOrDie(sizeof(ShiftDistalMispairs_TermLoop)*1); /* just alloc something, so we can use realloc later */
	
	/* find terminal loops, and do analysis for each one*/
	prevOpenPair=-1;
	for (col=0; col<msa->alen; col++) {
		if (prevOpenPair!=-1) {
			if (pair_right(msa->ss_cons[col])) {
				
				ShiftDistalMispairs_TermLoop *termLoopInfo;
				
				/* here's a terminal loop */
				int pairLeft=prevOpenPair;
				int pairRight=col;
				int seqnum;
				
				numTermLoops++;
				termLoopInfoVector=(ShiftDistalMispairs_TermLoop *)realloc(termLoopInfoVector,sizeof(ShiftDistalMispairs_TermLoop)*numTermLoops);
				if (termLoopInfoVector==NULL) {
					esl_fatal("realloc failed");
				}
				termLoopInfo=&(termLoopInfoVector[numTermLoops-1]);
				
				termLoopInfo->pairLeft=pairLeft;
				termLoopInfo->pairRight=pairRight;
				termLoopInfo->perSeqInfo=(ShiftDistalMispair_OneSeqAndTermLoop *)MallocOrDie(sizeof(ShiftDistalMispair_OneSeqAndTermLoop)*msa->nseq);
				
				for (seqnum=0; seqnum<msa->nseq; seqnum++) {
					
					int i;
					int l=pairLeft;
					int r=pairRight;
					int leftNucsToMove=0;
					int rightNucsToMove=0;
					int noMorePairs=0;
					int numNucsInTermLoop=0;
					
					for (i=l+1; i<r; i++) {
						if (!esl_abc_XIsGap(abc,msa->ax[seqnum][i+1])) {
							numNucsInTermLoop++;
						}
					}
						
					while (!noMorePairs && l>=0 && r<msa->alen) {
						
						/* first check truncation.  if something is truncated, we don't proceed further, no matter what */
						if (cmfinderVars->use_fragmentary) {
							if (seq_fragment_start_col[seqnum]!=-1 && l<seq_fragment_start_col[seqnum]) {  /* i<j, so if j is left-truncated, then i will be also */
								assert(esl_abc_XIsGap(abc,msa->ax[seqnum][l+1])); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
								break;
							}
							if (seq_fragment_end_col[seqnum]!=-1 && r>seq_fragment_end_col[seqnum]) {
								assert(esl_abc_XIsGap(abc,msa->ax[seqnum][r+1])); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
								break;
							}
						}
						
						/* now check if we _have_ to move into the loop due to --min-term-loop-nucs */
						if (numNucsInTermLoop < cmfinderVars->minTermLoopNucsWithMovingIntoTerminalLoops) {
							/* the base pair doesn't matter, we must move this into the term loop */

							/* bookkeeping: we're adding nucs to the loop */
							if (!esl_abc_XIsGap(abc,msa->ax[seqnum][l+1])) {
								numNucsInTermLoop++;
							}
							if (!esl_abc_XIsGap(abc,msa->ax[seqnum][r+1])) {
								numNucsInTermLoop++;
							}
						}
						else {
					
							/* if --min-term-loop-nucs isn't an issue, we check the pair, and leave it if the bp is canonical */
							if (IsCanonicalDigitalBpAndNotDegen(msa->ax[seqnum][l+1],msa->ax[seqnum][r+1])) {
								/* we've got a good pair, so stop looking */
								break;
							}
						}
						
						if (!esl_abc_XIsGap(abc,msa->ax[seqnum][l+1])) {
							leftNucsToMove++;
						}
						if (!esl_abc_XIsGap(abc,msa->ax[seqnum][r+1])) {
							rightNucsToMove++;
						}
						/* look for the next pair, and also worry about non-paired columns.  we stop if we find any unpaired column with a nuc in it */

						while (l>=0 && !noMorePairs) {
							l--;
							if (pair_right(msa->ss_cons[l])) {
								/* went out of the terminal stem */
								noMorePairs=1;
								break;
							}
							else {
								if (pair_left(msa->ss_cons[l])) {
									/* another base pair to test */
									break;
								}
								else {
									if (!esl_abc_XIsGap(abc,msa->ax[seqnum][l+1])) {
										/* a nuc in an unpaired column aborts */
										noMorePairs=1;
									}
								}
							}
						}
						
						while (r<msa->alen && !noMorePairs) {
							r++;
							if (pair_left(msa->ss_cons[r])) {
								/* went out of the terminal stem */
								noMorePairs=1;
								break;
							}
							else {
								if (pair_right(msa->ss_cons[r])) {
									/* another base pair to test */
									break;
								}
								else {
									if (!esl_abc_XIsGap(abc,msa->ax[seqnum][r+1])) {
										/* a nuc in an unpaired column aborts */
										noMorePairs=1;
									}
								}
							}
						}
					}
					
					termLoopInfo->perSeqInfo[seqnum].leftNucsToMove=leftNucsToMove;
					termLoopInfo->perSeqInfo[seqnum].rightNucsToMove=rightNucsToMove;
				}
				
				/* but now we've consumed the open pair */
				prevOpenPair=-1;
			}
		}
		/* any time we see pair_left, we consider this as terminal loop, since it'd clobber any prev term loop */
		if (pair_left(msa->ss_cons[col])) {
			prevOpenPair=col;
		}
	}

	ColumnsToAdd *columnsToAdd=(ColumnsToAdd *)MallocOrDie(sizeof(ColumnsToAdd)*(numTermLoops*2));
	for (tl=0; tl<numTermLoops; tl++) {
		int max_leftNucsToMove=0;
		int max_rightNucsToMove=0;
		int seqnum;
		for (seqnum=0; seqnum<msa->nseq; seqnum++) {
			if (termLoopInfoVector[tl].perSeqInfo[seqnum].leftNucsToMove>max_leftNucsToMove) {
				max_leftNucsToMove=termLoopInfoVector[tl].perSeqInfo[seqnum].leftNucsToMove;
			}
			if (termLoopInfoVector[tl].perSeqInfo[seqnum].rightNucsToMove>max_rightNucsToMove) {
				max_rightNucsToMove=termLoopInfoVector[tl].perSeqInfo[seqnum].rightNucsToMove;
			}
		}
		
		columnsToAdd[tl*2+0].colAfterInserts=termLoopInfoVector[tl].pairLeft+1;
		columnsToAdd[tl*2+0].numCols=max_leftNucsToMove;
		columnsToAdd[tl*2+1].colAfterInserts=termLoopInfoVector[tl].pairRight;
		columnsToAdd[tl*2+1].numCols=max_rightNucsToMove;
	}

	ESL_MSA *new_msa;
	int *oldToNewColumnMap;
	int columnsAdded;
	int columnsToAddEntries=numTermLoops*2;
	AddColumnsToMsa (&new_msa,msa,columnsToAdd,columnsToAddEntries,&oldToNewColumnMap,&columnsAdded);

	if (columnsAdded==0) {
		printf("--bad-distal-pairs-to-loop : no change\n");
	}
	else {
		int seqnum;

		/* move any nucs that we're supposed to move */
		for (tl=0; tl<numTermLoops; tl++) {
			for (seqnum=0; seqnum<msa->nseq; seqnum++) {
				int nextLeftSrc,nextLeftDest;
				int nextRightSrc,nextRightDest;
				
				nextLeftSrc=oldToNewColumnMap[termLoopInfoVector[tl].pairLeft];
				nextLeftDest=oldToNewColumnMap[termLoopInfoVector[tl].pairLeft+1]-1; /* right-most column within the recently inserted columns */
				ShiftDistalMispairsIntoTerminalLoops_DoMoving(new_msa,seqnum,nextLeftSrc,nextLeftDest,termLoopInfoVector[tl].perSeqInfo[seqnum].leftNucsToMove,-1);
				
				nextRightSrc=oldToNewColumnMap[termLoopInfoVector[tl].pairRight];
				nextRightDest=oldToNewColumnMap[termLoopInfoVector[tl].pairRight-1]+1; /* left-most column within the just-inserted columns */
				ShiftDistalMispairsIntoTerminalLoops_DoMoving(new_msa,seqnum,nextRightSrc,nextRightDest,termLoopInfoVector[tl].perSeqInfo[seqnum].rightNucsToMove,+1);
			}
		}
		
		esl_msa_Destroy(msa);
		msa=new_msa;

		free(oldToNewColumnMap);
	}

	for (tl=0; tl<numTermLoops; tl++) {
		free(termLoopInfoVector[tl].perSeqInfo);
	}
	free(termLoopInfoVector);
	
	free(columnsToAdd);
	
	return msa;
}


int cmpstringp(const void *p1, const void *p2) /* copied from 'man qsort' */
{
  return strcmp(* (char * const *) p1, * (char * const *) p2);
}

int CalcQuasiNumSpecies (ESL_MSA *msa)
{
  int i,j;
  int numSeqId=0;
  char **hitIdList;

  /* make array of hitIds, converting NCBI-project-ID-type seqIds to just the project */
  hitIdList=(char **)MallocOrDie(msa->nseq * sizeof(char *));
  for (i=0; i<msa->nseq; i++) {
    hitIdList[i]=NULL;
    if (strlen(msa->sqname[i])<5) { /* definitely isn't NCBI project ID format */
    }
    else {
      for (j=0; j<4; j++) {
	if (!isalpha(msa->sqname[i][j])) {
	  break;
	}
      }
      if (j==4 && isdigit(msa->sqname[i][4])) {
	hitIdList[i]=(char *)MallocOrDie(5);
	strncpy(hitIdList[i],msa->sqname[i],4);
	hitIdList[i][4]=0;
      }
    }
    if (hitIdList[i]==NULL) {
      /* strip past slash to extract seqId, assuming it's in Rfam's hitId format (SEQID/START-END) */
      char *slash=strchr(msa->sqname[i],'/');
      if (slash==NULL) {
	/* weird, but okay */
	slash=msa->sqname[i]+strlen(msa->sqname[i]);
      }
      hitIdList[i]=(char *)MallocOrDie(slash-msa->sqname[i]+1);
      strncpy(hitIdList[i],msa->sqname[i],slash-msa->sqname[i]);
      hitIdList[i][slash-msa->sqname[i]]=0;
    }
  }

  /* count unique ids by first sorting, and then running through, comparing each entry with the previous one, since there isn't a convenient 'set' datastructure in standard C */
  qsort(hitIdList,msa->nseq,sizeof(char *),cmpstringp);
  for (i=0; i<msa->nseq; i++) {
    if (i==0) { /* no previous entry, so it's unique */
      numSeqId++;
    }
    else {
      if (strcmp(hitIdList[i-1],hitIdList[i])!=0) {
	numSeqId++;
      }
    }
  } 
  
  for (i=0; i<msa->nseq; i++) {
    free(hitIdList[i]);
  }
  free(hitIdList);

  return numSeqId;
}
