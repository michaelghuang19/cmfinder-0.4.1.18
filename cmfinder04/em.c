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

#include "esl_exponential.h"

enum M_step_traceback {
	M_step_traceback__impossible=0, /* the numbers don't matter, but it's convenient to assign them explicitly for debugging */
	M_step_traceback__min_hairpin=1,
	M_step_traceback__skip_gap_col_left=2,
	M_step_traceback__skip_gap_col_right=3,
	M_step_traceback__extend_pair=4,
	M_step_traceback__open_pair=5,
	M_step_traceback__emit_left=6,
	M_step_traceback__emit_right=7,
	M_step_traceback__bifurcate=8
};

void M_step_backtrace_sscons(ESL_MSA *msa,int *num_consensus_bp,int i,int j,double **pmat,int **pmat_traceback,double **zmat,int **zmat_traceback,int **zmat_traceback_mid)
{
	int is_pmat;
	int mid;
	int cont=1;
	
	is_pmat=0;
	while (cont) { /* terminates when traceback==M_step_traceback__min_hairpin, or after bifurcation */
		enum M_step_traceback traceback;
		if (is_pmat) {
			traceback=pmat_traceback[j][i];
		}
		else {
			traceback=zmat_traceback[j][i];
		}
		if (0) { printf("%d,%d,%d : %d\n",i,j,is_pmat,traceback); }
		switch (traceback) {
		case M_step_traceback__impossible:
			esl_fatal("got M_step_traceback__impossible");
			break;
		case M_step_traceback__min_hairpin:
			/* too small, so no more pairs */
			cont=0;
			break;
		case M_step_traceback__skip_gap_col_left:
		case M_step_traceback__emit_left:
			i++;
			break;
		case M_step_traceback__skip_gap_col_right:
		case M_step_traceback__emit_right:
			j--;
			break;
		case M_step_traceback__extend_pair:
		case M_step_traceback__open_pair:
			msa->ss_cons[i-1]='<'; /* everything here was 1-based, but the ss_cons line is 0-based */
			msa->ss_cons[j-1]='>';
			(*num_consensus_bp)++;
			i++;
			j--;
			if (traceback==M_step_traceback__open_pair) {
				is_pmat=0;
			}
			else {
				is_pmat=1;
			}
			break;
		case M_step_traceback__bifurcate:
			assert(is_pmat==0);
			mid=zmat_traceback_mid[j][i];
			M_step_backtrace_sscons(msa,num_consensus_bp,i,mid,pmat,pmat_traceback,zmat,zmat_traceback,zmat_traceback_mid);
			M_step_backtrace_sscons(msa,num_consensus_bp,mid+1,j,pmat,pmat_traceback,zmat,zmat_traceback,zmat_traceback_mid);
			cont=0;
			break;
		default:
			esl_fatal("internal error");
		}
	}
}

/* Dynamic programming to actually solve the sscons
* The basic recursion is:
*    Sij = max { Si+1,j  (emit left, no covariance)
*                Si,j-1  (emit right, no covariance)
*                Si+1,j-1 + xy[j][i].
*                max over mid: Si,mid + Smid+1,j (bifurcation)
*                }
* 
* But, we use matrices pmat for cases where the outer nucs are paired, 
* and zmat for non-paired outer nucs.
* These two matrices allow us to model the openPairCost
*/
void M_step_infer_sscons_dynamic_programming(ESL_MSA *msa,double **xy,int *is_consensus_col,float openPairCost,int *num_consensus_bp)
{
	double **pmat;
	double **zmat;
	int **zmat_traceback,**zmat_traceback_mid;
	int **pmat_traceback;
	int diff;
	int i,j;
	int mid;

	/* allocate arrays */
	zmat = DoubleAlloc2DArray(msa->alen+1);
	pmat = DoubleAlloc2DArray(msa->alen+1);
	
	zmat_traceback=IntAlloc2DArray(msa->alen+1);
	zmat_traceback_mid=IntAlloc2DArray(msa->alen+1);
	pmat_traceback=IntAlloc2DArray(msa->alen+1);

	/* Initialize to base case */
	for(diff = 0; diff <= MIN_HAIRPIN; diff++) {
		for (i = 1; (j = i + diff) <= msa->alen; i++) { /* changed the test to be '<=', but it was '<' in the original CMfinder code */
			zmat[j][i] = 0;
			zmat_traceback[j][i]=M_step_traceback__min_hairpin;
			pmat[j][i] = INT_NEGINFINITY;
			pmat_traceback[j][i]=M_step_traceback__impossible;
		}
	}

	/* here's the actual dynamic program */
	for (diff = MIN_HAIRPIN+1; diff < msa->alen; diff++) {
		for (i = 1; (j = i+diff) <= msa->alen; i++) {
			/* Ignore gap */
			if (!is_consensus_col[j-1]) {
				zmat[j][i] = zmat[j-1][i];
				pmat[j][i] = pmat[j-1][i];
				zmat_traceback[j][i]=M_step_traceback__skip_gap_col_right;
				pmat_traceback[j][i]=M_step_traceback__skip_gap_col_right;
			}
			else if (!is_consensus_col[i-1]) {
				zmat[j][i] = zmat[j][i+1];
				pmat[j][i] = pmat[j][i+1];
				zmat_traceback[j][i]=M_step_traceback__skip_gap_col_left;
				pmat_traceback[j][i]=M_step_traceback__skip_gap_col_left;
			}
			else {
				
				/* update pmat 
				 * since pmat means that the outer nucs are paired, we only have one case
				 */
				if (xy[j][i] > 0){
					pmat[j][i] =  xy[j][i]; /* any pairing has this cost */
					if (pmat[j-1][i+1] > zmat[j-1][i+1] - openPairCost) {
						pmat[j][i] += pmat[j-1][i+1]; /* end a previous pair */
						pmat_traceback[j][i]=M_step_traceback__extend_pair;
					}
					else {
						pmat[j][i] += zmat[j-1][i+1] - openPairCost; /* start a new helix */
						pmat_traceback[j][i]=M_step_traceback__open_pair;
					}
				}
				else {
					pmat[j][i] = INT_NEGINFINITY; /* evidence against the pairing, so say it's impossible */
					pmat_traceback[j][i]=M_step_traceback__impossible;
				}

				/* update zmat */
				zmat[j][i] = zmat[j-1][i]; /* assume right emit, to start */
				zmat_traceback[j][i]=M_step_traceback__emit_right;
				if (zmat[j][i+1] > zmat[j][i]) {
					zmat[j][i] = zmat[j][i+1]; /* left emit is better */
					zmat_traceback[j][i]=M_step_traceback__emit_left;
				}
				for (mid = i+1; mid < j-1; mid++) { /* consider bifurcations */
					if (zmat[mid][i] + zmat[j][mid+1] > zmat[j][i]) {
						zmat[j][i] = zmat[mid][i] + zmat[j][mid+1];
						zmat_traceback[j][i]=M_step_traceback__bifurcate;
						zmat_traceback_mid[j][i]=mid;
					}
				}
				/* If it is OK to base pair i,j */
				if (xy[j][i] > 0 && pmat[j-1][i+1] + xy[j][i] > zmat[j][i]){
					/* terminate a pairing run here.  
					 * note: these recursions allow us to stop a helix and then start it with an 
					 * immediately adjacent pair, but that's okay -- it's just 
					 * suboptimal because of openPairCost 
					 * 
					 * I think the advantage of doing the recursions this way is that 
					 * we don't have to add the possibility of coming from pmat in the recursions for zmat 
					 */
					zmat[j][i] = pmat[j-1][i+1] + xy[j][i];
					zmat_traceback[j][i]=M_step_traceback__extend_pair;
				}
			}
		}
	}
	
	/* and do traceback */
	*num_consensus_bp=0;
	M_step_backtrace_sscons(msa,num_consensus_bp,1,msa->alen,pmat,pmat_traceback,zmat,zmat_traceback,zmat_traceback_mid);

	/* and free memory */
	Free2DArray((void **)zmat,msa->alen+1);
	Free2DArray((void **)zmat_traceback,msa->alen+1);
	Free2DArray((void **)zmat_traceback_mid,msa->alen+1);
	Free2DArray((void **)pmat,msa->alen+1);
	Free2DArray((void **)pmat_traceback,msa->alen+1);
}

/* M_step_infer_sscons_and_rf is roughly like the original cmfinder function 'Automodelmaker'
 * code having to do with 'cluster' is removed, because it didn't seem to be used in the original cmfinder.
 */
void M_step_infer_sscons_and_rf (ESL_MSA *msa, ESL_MSA *text_msa,double **pxy,CmfinderVars *cmfinderVars,double openPairCost,int *num_consensus_bp)
{
	double **mxy;
	double **xy;
	int *is_consensus_col;	/* 0..alen-1 array; 0=insert col, 1=match col ; called 'matassign' in original cmfinder code */
	double   sum_weight=0;
	int i,j;
	int seqnum;
	int *seq_fragment_start_col,*seq_fragment_end_col;
	
	/* figure out start/end's due to truncation */
	seq_fragment_start_col=(int *)MallocOrDie(msa->nseq * sizeof(int));
	seq_fragment_end_col=(int *)MallocOrDie(msa->nseq * sizeof(int));
	CalcFragmentStartEnds (msa,cmfinderVars->use_fragmentary,seq_fragment_start_col,seq_fragment_end_col);

	/* remove existing structure, since we're going to infer the structure for the next iteration */
	AllocSSconsOrRfIfNotThere(&(msa->ss_cons),msa);
	for(i=0; i < msa->alen; i++) {
		msa->ss_cons[i] = '.';
	}

	/* calculate informative priors for all possible pairs, based on Boltzmann distribution */
	if (pxy==NULL) {
		if (cmfinderVars->avoidPartFunc) {
			pxy = DoubleAlloc2DArray(text_msa->alen+1);
			for(j = 0; j < text_msa->alen; j++) {
				for(i = 0; i < j; i++) {
					pxy[j+1][i+1]=0;
					pxy[j][i]=0;
				}
			}
		}
		else {
			if (cmfinderVars->columnOnlyBasePairProbs) {
				prxy_column_only(msa, NULL, cmfinderVars, seq_fragment_start_col, seq_fragment_end_col, &pxy);
			}
			else {
				prxy(text_msa, NULL, cmfinderVars, seq_fragment_start_col, seq_fragment_end_col, &pxy);
			}
		}
	}

	/* build Mxy matrix, which gives mutual information for all possible pairs */
	mixy(msa, cmfinderVars, &mxy);   /* mutual information matrix */

	merge(mxy, pxy, msa->alen, &xy);

	is_consensus_col=MallocOrDie(sizeof(int) * (msa->alen+1));

	sum_weight = 0;
	for (i = 0; i < msa->nseq; i++) {
		sum_weight += msa->wgt[i];
	}
	
	for (i = 0; i < msa->alen; i++) {
		double this_sum_weight=0;
		double gaps = 0;
		for (seqnum = 0;seqnum < msa->nseq; seqnum++) {
			if (cmfinderVars->use_fragmentary) {
				if (seq_fragment_start_col[seqnum]!=-1 && i<seq_fragment_start_col[seqnum]) {
					assert(esl_abc_CIsGap(abc,text_msa->aseq[seqnum][i])); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
					continue;
				}
				if (seq_fragment_end_col[seqnum]!=-1 && i>seq_fragment_end_col[seqnum]) {
					assert(esl_abc_CIsGap(abc,text_msa->aseq[seqnum][i])); /* just making sure I got the truncation right -- if it's truncated it should certainly be a gap */
					continue;
				}
			}
			this_sum_weight += msa->wgt[seqnum];
			if (esl_abc_XIsGap(abc,msa->ax[seqnum][i+1])) {
				gaps += msa->wgt[seqnum];
			}
		}
		is_consensus_col[i] = ( gaps / this_sum_weight > cmfinderVars->gapthreshold) ? 0 : 1;
	}
	
	/* set rf line from is_consensus_col */
	AllocSSconsOrRfIfNotThere(&(msa->rf),msa);
	for(i=0; i < msa->alen; i++) {
		if (is_consensus_col[i]) {
			msa->rf[i]='x';
		}
		else {
			msa->rf[i]='.';
		}
	}

	M_step_infer_sscons_dynamic_programming(msa,xy,is_consensus_col,openPairCost,num_consensus_bp);
	
	if (cmfinderVars->saveMsaAfterFirstMstep) {
		char *fileName=(char *)MallocOrDie(strlen(cmfinderVars->final_file)+100);
		sprintf(fileName,"%s.after-first-M-step.sto",cmfinderVars->final_file);
		SaveMsa(fileName,msa);
		printf("saved MSA after first M-step to file \"%s\"\n",fileName);
		free(fileName);
	}

	Free2DArray((void **)mxy,  msa->alen+1);
	Free2DArray((void **)pxy,  msa->alen+1);
	Free2DArray((void **)xy,  msa->alen+1);
	free(is_consensus_col);
	free(seq_fragment_start_col);
	free(seq_fragment_end_col);
}

CM_t * M_step (ESL_MSA *digital_msa, double **pxy, CmfinderVars *cmfinderVars,double openPairCost)
{
	char cmbuild_commandline_flags[COMMANDLINE_FLAGS_MAX];
	CmdlineBuilder cmdline;
	CM_t *cm; /* new covariance model */
	ESL_MSA *text_msa;
	ESL_STOPWATCH *timer_m_step=esl_stopwatch_Create();
	ESL_STOPWATCH *timer_cmbuild=esl_stopwatch_Create();
	esl_stopwatch_Start(timer_m_step);

	int num_consensus_bp;
	
	text_msa=MakeTextMsaClone(digital_msa);

        //printf("openPairCost=%g\n",openPairCost);
	//SaveMsa("just-before-M_step_infer_sscons_and_rf.sto",digital_msa);
	M_step_infer_sscons_and_rf(digital_msa,text_msa,pxy,cmfinderVars,openPairCost,&num_consensus_bp);

	if (num_consensus_bp < cmfinderVars->minPredictedBasePairs) {
		printf("only %d base pairs predicted, so stopping\n",num_consensus_bp);
		cm=NULL;
	}
	else {
	
		/* run cmbuild on the msa */
#if 0
		char more_compatible_with_older_cmfinder[]="--p56 --enone";
#endif
		esl_stopwatch_Start(timer_cmbuild);
		InitCmdlineBuilder(&cmdline,cmbuild_commandline_flags);
		PrintfCmdlineBuilder(&cmdline,"executable-dummy-token");
		if (esl_opt_GetBoolean(cmfinderVars->go, "--enone")) {
			PrintfCmdlineBuilder(&cmdline," --enone");
		}
		if (esl_opt_GetBoolean(cmfinderVars->go, "--p56")) {
			PrintfCmdlineBuilder(&cmdline," --p56");
		}
		PrintfCmdlineBuilder(&cmdline," --hand --wgiven -n CMfinder-internal");
		if (0) { printf("cmbuild flags: %s\n",cmbuild_commandline_flags); }
		cmbuild(cmbuild_commandline_flags,"CMfinder-internal.sto",digital_msa,&cm,cmfinderVars->hmmStatsToCalibrate);

		esl_stopwatch_Stop(timer_cmbuild);
		esl_stopwatch_Include(cmfinderVars->timer_cmbuild,timer_cmbuild);
	}
	
	esl_msa_Destroy(text_msa);
	
	esl_stopwatch_Destroy(timer_cmbuild);

	esl_stopwatch_Stop(timer_m_step);
	esl_stopwatch_Include(cmfinderVars->timer_m_step,timer_m_step);
	esl_stopwatch_Destroy(timer_m_step);
	
	return cm;
}

/* for --bg-score-evalue
 * similar to fit_histogram from cmcalibrate.c, but 
 * (1) doesn't take ESL_GETOPTS go or struct cfg, 
 * (2) takes CM_TOPHITS instead of 'scores' array, since that's more convenient for me
 * (3) doesn't dump to files based on that code, 
 * (4) takes 'tailp' as a parameter, since it can't read it from command-line flags, 
 * (5) doesn't need 'exp_mode', since that was just to get tailp from command-line flags,
 * (6) alloc's its own errbuf
 * (7) added ret_evdWorks flag (see comment under 'int evdWorks=1;' below)
 */
static int
fit_histogram_tophits (int *ret_evdWorks,float tailp,CM_TOPHITS *topHits, double *ret_mu, double *ret_lambda, int *ret_nrandhits, float *ret_tailp)
{
	int status;
	double mu;
	double lambda;
	double *xv;         /* raw data from histogram */
	int     n,z;  
	double  params[2];
	int     nrandhits; 
	int hitnum;
	int evdWorks=1; /* assume we can infer the EVD parameters, until we find an error.  In that case, back out with defensive code.  Most errors are unexpected and considered fatal, but some are anticipated cases where we're unable to infer an EVD. */

	ESL_HISTOGRAM *h = NULL;       /* histogram of scores */

	/* Initialize histogram; these numbers are guesses */
	if((h = esl_histogram_CreateFull(-100., 100., .1)) == NULL) return eslEMEM;    

	/* fill histogram */
	for (hitnum = 0; hitnum < topHits->N; hitnum++) {
		CM_HIT *hit=topHits->hit[hitnum];
		if ((hit->flags&CM_HIT_IS_REMOVED_DUPLICATE)!=0) {
			continue;
		}
		if((status = esl_histogram_Add(h, hit->score)) != eslOK) esl_fatal("fit_histogram(), esl_histogram_Add() call returned non-OK status: %d\n", status);
	}

	if (esl_histogram_GetTailByMass(h, tailp, &xv, &n, &z) != eslOK) { /* fit to right 'tailfit' fraction, 'tailfit' was determined in above block */
		esl_fatal("esl_histogram_GetTailByMass failed");
	}
	if(n <= 1) {
		evdWorks=0;
		printf("--bg-score-evalue: too few points in right tailfit.  #hits=%u, tailp=%g\n",(unsigned int)(topHits->N),tailp); /* avoid PRIu64 macro, since I'm not sure about portability */
	}
	
	if (evdWorks) {
		if (esl_exp_FitComplete(xv, n, &(params[0]), &(params[1])) != eslOK) {
			esl_fatal("esl_exp_FitComplete failed");
		}
		if (esl_histogram_SetExpectedTail(h, params[0], tailp, &esl_exp_generic_cdf, &params) != eslOK) {
			esl_fatal("esl_histogram_SetExpectedTail failed");
		}

		/* printf("# Exponential fit to %.7f%% tail: lambda = %f\n", tailp*100.0, params[1]); */
		mu = params[0];
		lambda = params[1];
		if(isnan(lambda)) esl_fatal("fit_histogram(), exp tail fit lambda is NaN, too few hits in histogram. Increase <x> with -L <x>.");
		if(isinf(lambda)) esl_fatal("fit_histogram(), exp tail fit lambda is inf, too few hits in histogram. Increase <x> with -L <x>.");
		nrandhits = h->n; /* total number of hits in the histogram */
	}

	/* print to output files if nec */
	esl_histogram_Destroy(h);

	*ret_mu     = mu;
	*ret_lambda = lambda;
	*ret_nrandhits = nrandhits;
	*ret_tailp = tailp;
	*ret_evdWorks=evdWorks;
	return eslOK;
}

double ScanCand_GetScoreThreshold (CmfinderVars *cmfinderVars,CM_t *input_cm,ScanMode parent_scanMode,ESL_SQCACHE *inputSeqCache)
{
	if (cmfinderVars->bgSeqCache->seq_count==0) {
		return cmfinderVars->cm_scoreThreshold;
	}
	else {
		int bgScoreNumNucsScanned=cmfinderVars->inputSeqSize*cmfinderVars->bgScoreSize;
		CM_t *cmsearch_cm;
		CM_TOPHITS *topHits;
		int hitnum;
		double *maxScorePerBgCopy;
		double cm_scoreThreshold;
		ScanMode scanMode=parent_scanMode;
		int i;
		ESL_HISTOGRAM *h=NULL; /* shut up compiler warnings */
		
		scanMode.need_alignments=0;
		
		if (cmfinderVars->bgScoreNonFrag) {
			scanMode.use_fragmentary=0;
		}
		
		cm_scoreThreshold=cmfinderVars->cm_scoreThreshold;
		maxScorePerBgCopy=(double *)MallocOrDie(sizeof(double)*cmfinderVars->bgScoreSize);
		for (i=0; i<cmfinderVars->bgScoreSize; i++) {
			maxScorePerBgCopy[i]=0;
		}
		assert(cmfinderVars->bgSeqCache->seq_count = inputSeqCache->seq_count * cmfinderVars->bgScoreSize); /* otherwise some of the bgCopyNum code doesn't work */
		
		cmsearch_cm=save_and_load_CM_via_memory_copying(input_cm,cmfinderVars); /* the hack to avoid having 'cmsearch' complain that it's already configured */
		cmsearch_wrapper (cmsearch_cm,cmfinderVars->bgSeqCache,&topHits,cmfinderVars,cmfinderVars->bgScoreScanBitThresh,scanMode);

		if (cmfinderVars->bgScoreEvalue>0) {
			if((h = esl_histogram_CreateFull(-100., 100., .1)) == NULL) esl_fatal("esl_histogram_CreateFull failed");
		}
		for (hitnum = 0; hitnum < topHits->N; hitnum++) {
			int bgCopyNum;
			CM_HIT *hit=topHits->hit[hitnum];
			if ((hit->flags&CM_HIT_IS_REMOVED_DUPLICATE)!=0) {
				continue;
			}

			if (cmfinderVars->bgScoreEvalue>0) {
				if(esl_histogram_Add(h, hit->score) != eslOK) esl_fatal("esl_histogram_Add failed");
			}
			
			bgCopyNum=hit->seq_idx % cmfinderVars->bgScoreSize;
			if (hit->score > maxScorePerBgCopy[bgCopyNum]) {
				maxScorePerBgCopy[bgCopyNum]=hit->score;
			}
			
			if (hit->score > cm_scoreThreshold) { /* can only increase the cm_scoreThreshold.  This is an issue if cmfinderVars->bgScoreScanBitThresh < cm_scoreThreshold, where some hits can be lower than cm_scoreThreshold */
				cm_scoreThreshold=hit->score;
			}
		}
		
		printf("--bg-score-size : max(score)=%lg\n",cm_scoreThreshold);

		if (0) {
			printf("per-copy bg-score stats:");
			for (i=0; i<cmfinderVars->bgScoreSize; i++) {
				printf(" %lg",maxScorePerBgCopy[i]);
			}
			printf("\n");
		}

		if (cmfinderVars->bgScoreEvalue>0) {
			
			double input_tailp=0.6;
			double tmp_mu=0, tmp_lambda=0; /* initialize to avoid compiler warnings, but too bad valgrind won't work any more -- it's up to us to make sure it's correct */
			int tmp_nrandhits=0;
			float tmp_tailp;
			ExpInfo_t *expInfo;
			double pvalue;
			int evdWorks=1; /* see comment for same-named variable in 'fit_histogram_tophits' */
			
			/* this from code from serial_master in cmcalibrate.c */
			if (fit_histogram_tophits (&evdWorks,input_tailp, topHits, &tmp_mu, &tmp_lambda, &tmp_nrandhits, &tmp_tailp) != eslOK) {
				cm_Fail("fit_histogram_tophits failed");
			}
			
			if (evdWorks) {
				double backCalcEvalue;
				double evdBasedScoreThreshold;
				
				expInfo=CreateExpInfo();
				SetExpInfo(expInfo, tmp_lambda, tmp_mu, (double)(bgScoreNumNucsScanned), tmp_nrandhits, tmp_tailp);
				
				/* from UpdateExpsForDBSize in stats.c, and its call from cm_pli_NewModel in cm_pipeline.c */
				expInfo->cur_eff_dbsize = ((double)(cmfinderVars->inputSeqSize) / expInfo->dbsize) * ((double) expInfo->nrandhits);
				
				/* (BTW, code in cm_pipeline.c) P = esl_exp_surv(sc, cm->expA[pli->fcyk_cm_exp_mode]->mu_extrap, cm->expA[pli->fcyk_cm_exp_mode]->lambda); */

				pvalue=cmfinderVars->bgScoreEvalue/expInfo->cur_eff_dbsize;  /* inverse of 'evalue=pvalue*expInfo->cur_eff_dbsize;' from cm_tophits_ComputeEvalues in cm_tophits.c */
				evdBasedScoreThreshold=esl_exp_invsurv(pvalue, expInfo->mu_extrap, expInfo->lambda);
				if (evdBasedScoreThreshold>cm_scoreThreshold) { /* never lower score threshold.  I think this is the conservative way to play it. */
					cm_scoreThreshold=evdBasedScoreThreshold;
				}

				backCalcEvalue=esl_exp_surv(cm_scoreThreshold,expInfo->mu_extrap,expInfo->lambda)*expInfo->cur_eff_dbsize;
				printf("bg info cm_scoreThreshold=%lg, pvalue=%lg, target E-value=%lg, back-calc'd E-value(there is 0<=p<=1 problems)=%lg\n",cm_scoreThreshold,pvalue,cmfinderVars->bgScoreEvalue,backCalcEvalue);
				
				free(expInfo);
			}
			else {
				/* just use the max(cm_scoreThreshold) from above */
			}
		}
		
		cm_tophits_Destroy(topHits);
		FreeCM(cmsearch_cm);
		free(maxScorePerBgCopy);
		
		return cm_scoreThreshold;
	}
}

void EliminateHitsByVariousCriteria(Cand **cand,int *ncand,ESL_SQCACHE *inputSeqCache,int max_cand,CmfinderVars *cmfinderVars)
{
	int d=0;
	if (cmfinderVars->eliminateIdenticalSeqs || cmfinderVars->eliminateIdenticalSubseqs || cmfinderVars->eliminateSeqsWithManyDegen) {
		int seqnum;
		int candnum,s2,c2;
		int nextcandnum;
		int *remove;
		ESL_STOPWATCH *timer_elim_iden_seqs=esl_stopwatch_Create();
		
		esl_stopwatch_Start(timer_elim_iden_seqs);

		remove=(int *)MallocOrDie(sizeof(int)*max_cand);
		for (seqnum=0; seqnum<inputSeqCache->seq_count; seqnum++) {

			/* initialize: don't remove anything in this seqnum */
			for (candnum=0; candnum < ncand[seqnum]; candnum++) {
				remove[candnum]=0;
			}
			
			/* eliminate seqs with too many N's */
			if (cmfinderVars->eliminateSeqsWithManyDegen) {
				ESL_SQ *sq=&(inputSeqCache->sq_list[seqnum]);
				for (candnum=0; candnum < ncand[seqnum]; candnum++) {
					Cand *c=&(cand[seqnum][candnum]);
					int numDegen;
					int start,stop,i;
					assert(c->seq_id==seqnum);
					
					start=c->start - cmfinderVars->flankingNucsForCountingDegen;
					stop=c->stop + cmfinderVars->flankingNucsForCountingDegen;
					if (start<1) {
						start=1;
					}
					if (stop>sq->L) {
						stop=sq->L;
					}
					
					numDegen=0;
					for (i=start; i<=stop; i++) {
						if (esl_abc_XIsDegenerate(abc,sq->dsq[i])) {
							numDegen++;
						}
					}
					
					if (numDegen > cmfinderVars->maxDegenPerSeq) {
						remove[candnum]=1;
					}
				}
			}
			
			/* eliminate identical seqs
			 * we could use a hash join to compare for equality more quickly, but that'd be a hassle to implement,
			 * and I doubt this will compare to the time required for the cmbuild/cmsearch/partition-function steps */
			if (cmfinderVars->eliminateIdenticalSeqs) {
				for (candnum=0; candnum < ncand[seqnum]; candnum++) {
					for (s2=seqnum; s2<inputSeqCache->seq_count; s2++) {
						for(c2=0; c2 < ncand[s2]; c2++) {
							if (seqnum<s2 || (seqnum==s2 && candnum<c2)) { /* impose ordering so we always leave a survivor */
								if (cand[seqnum][candnum].len==cand[s2][c2].len) {
									if (strcmp(cand[seqnum][candnum].seq,cand[s2][c2].seq)==0) {
										remove[candnum]=1;
										break;
									}
								}
							}
						}
						if (remove[candnum]) {
							break;
						}
					}
				}
			}
			
			/* eliminate seqs that are identical to subseqs of other seqs */
			if (cmfinderVars->eliminateIdenticalSubseqs) {
				for (candnum=0; candnum < ncand[seqnum]; candnum++) {
					for (s2=0; s2<inputSeqCache->seq_count; s2++) { /* in this case, we compare all against all, because a couple of lines down we'll require a strict order of 'len', which guarantees a survivor */
						for(c2=0; c2 < ncand[s2]; c2++) {
							if (cand[seqnum][candnum].len < cand[s2][c2].len) {  /* can only be subseq if its len is smaller, and exact match is taken care of above.  BTW, this also imposes an ordering; we can't delete both a subseq and the superseq. */
								if (strstr(cand[s2][c2].seq,cand[seqnum][candnum].seq)!=NULL) {  /* find within s2,c2 an occurrence of seqnum,candnum (params are opposite order to what I'd have expected) */
									remove[candnum]=1;
									if (d) {
										printf("substr: %s in %s\n",cand[seqnum][candnum].seq,cand[s2][c2].seq);
									}
									break;
								}
							}
						}
						if (remove[candnum]) {
							break;
						}
					}
				}
			}
			
			nextcandnum=0;
			for (candnum=0; candnum < ncand[seqnum]; candnum++) {
				if (remove[candnum]) {
					DeleteOneCand(&(cand[seqnum][candnum]));
				}
				else {
					if (nextcandnum!=candnum) {
						if (!cmfinderVars->dieOnUntestedCode) {
							assert(0); /* congrats: you've exercized this code */
						}
						cand[seqnum][nextcandnum]=cand[seqnum][candnum];
					}
					nextcandnum++;
				}
			}
			ncand[seqnum]=nextcandnum;
		}
		free(remove);

		esl_stopwatch_Stop(timer_elim_iden_seqs);
		esl_stopwatch_Include(cmfinderVars->timer_elim_iden_seqs,timer_elim_iden_seqs);
		esl_stopwatch_Destroy(timer_elim_iden_seqs);
	}
}

/* adapted from hit_sorter_for_overlap_removal (cm_tophits.c) */
int hit_sorter_by_position(const void *vh1, const void *vh2)
{
	CM_HIT *h1 = *((CM_HIT **) vh1);  /* don't ask. don't change. Don't Panic. */
	CM_HIT *h2 = *((CM_HIT **) vh2);
	
	/* key: cm_idx */
	if (h1->cm_idx != h2->cm_idx) {
		return h1->cm_idx - h2->cm_idx;
	}
	
	/* key: seq_idx */
	if (h1->seq_idx != h2->seq_idx) {
		return h1->seq_idx - h2->seq_idx;
	}
	
	/* key: in_rc , even tho CMfinder won't use rev-comp */
	if (h1->in_rc != h2->in_rc) {
		return h1->in_rc - h2->in_rc;
	}
	
	/* key: start (here's where it differs from infernal) */
	if (h1->start != h2->start) {
		return h1->start - h2->start;
	}
	
	/* they're equal, as far as we care */
	return 0;
}

void ScanCand_FilterNonFrag (ESL_SQCACHE **get_filterSeqCache,int **get_scannedSeqToActualSeq,int **get_scannedSeqToOriginalStart,CM_t *input_cm,ESL_SQCACHE *inputSeqCache,CmfinderVars *cmfinderVars,double cm_scoreThreshold,ScanMode scanMode)
{
	CM_t *filter_cmsearch_cm;
	CM_TOPHITS *topHits=NULL;
	int seqnum,start,stop;
	int prev_seqnum,prev_start,prev_stop;
	int hitnum,i;
	ESL_SQCACHE *filterSeqCache=NULL;
	int *scannedSeqToActualSeq=NULL;
	int *scannedSeqToOriginalStart=NULL;

	int debug_merge_overlap=0;
	
	filter_cmsearch_cm=save_and_load_CM_via_memory_copying(input_cm,cmfinderVars); /* the hack to avoid having 'cmsearch' complain that it's already configured */

	if (0) {
		SaveCM("test.cm",filter_cmsearch_cm); SaveSeqCacheAsFasta ("test.fasta",inputSeqCache);
	}
	scanMode.use_fragmentary=0; /* we're doing --filter-non-frag */
	cmsearch_wrapper (filter_cmsearch_cm,inputSeqCache,&topHits,cmfinderVars,cm_scoreThreshold,scanMode);

	/* adapted from cm_tophits_SortForOverlapRemoval (cm_tophits.c)*/
	for (hitnum = 0; hitnum < topHits->N; hitnum++) topHits->hit[hitnum] = topHits->unsrt + hitnum;
	if (topHits->N > 1) {
		qsort(topHits->hit, topHits->N, sizeof(CM_HIT *), hit_sorter_by_position);
	}
	topHits->is_sorted_by_evalue           = FALSE;
	topHits->is_sorted_by_position         = FALSE;
	topHits->is_sorted_for_overlap_removal = TRUE;

	filterSeqCache=(ESL_SQCACHE *)MallocOrDie(sizeof(ESL_SQCACHE));
	filterSeqCache->seq_count=0; /* we'll add the seqs as we go */
	
	if (topHits->N>0) {
		filterSeqCache->sq_list=(ESL_SQ *)MallocOrDie(sizeof(ESL_SQ)*topHits->N);  /* topHits->N is an upper bound; overlapping hits will make it smaller.  But I think the wasted memory will be negligable */
		scannedSeqToActualSeq=(int *)MallocOrDie(sizeof(int)*topHits->N);  /* ditto */
		scannedSeqToOriginalStart=(int *)MallocOrDie(sizeof(int)*topHits->N);  /* ditto */
		
		prev_seqnum=-1;
		for (hitnum = 0; hitnum < topHits->N; hitnum++) {
			ESL_SQ *src_sq;
			CM_HIT *hit=topHits->hit[hitnum];
			int add_prev;
			if ((hit->flags&CM_HIT_IS_REMOVED_DUPLICATE)!=0) {
				continue;
			}
			if (hit->score<cm_scoreThreshold) {
				assert(cmfinderVars->use_evalue); /* this is the only case in which a hit below the score threshold could be reported, since -E and -T in cmsearch are incompatible */
				continue;
			}
			seqnum=hit->seq_idx;
			assert(!hit->in_rc); /* CMfinder shouldn't be looking for rev-comp */
			assert(hit->start<=hit->stop);
			src_sq=&(inputSeqCache->sq_list[seqnum]);
			start=hit->start - cmfinderVars->filterNonFragPad;
			stop=hit->stop + cmfinderVars->filterNonFragPad;
			if (start<1) {
				start=1;
			}
			if (stop>src_sq->L) {
				stop=src_sq->L;
			}
			
			if (debug_merge_overlap) printf("--filter-non-frag : got hit %d/%d-%d --> %d-%d\n",(int)seqnum,(int)hit->start,(int)hit->stop,(int)start,(int)stop);
			
			/* validating the sort */
			if (hitnum>0) {
				assert(prev_seqnum<seqnum || (prev_seqnum==seqnum && prev_start<start) || (prev_seqnum==seqnum && prev_start==start && prev_stop<stop) || (prev_seqnum==seqnum && prev_start==start && prev_stop==stop)); /* 2nd-to-last case is because with cmfinderVars->filterNonFragPad , prev_start==start is quite possible if it bleeds into the start of the seq fragment.  last case is because it happened once, and maybe it's legit */
			}
			
			add_prev=0;
			if (prev_seqnum==-1) {
				// no previous hit to add, so nothing to do
			}
			else {
				if (seqnum==prev_seqnum) {
					if (start<=prev_stop) {
						if (debug_merge_overlap) printf("\tmerge with previous\n");
						
						/* merge this overlapping seq, so don't add the hit yet */
						start=prev_start; /* will get merged at the bottom of the loop, where we assign to prev_... */
					}
					else {
						// doesn't overlap, so add the previous
						add_prev=1;
					}
				}
				else {
					// different seq, so add the previous
					add_prev=1;
				}
			}
			if (add_prev) {
				ESL_SQ *sq=&(filterSeqCache->sq_list[filterSeqCache->seq_count]);
				assert(prev_seqnum>=0 && prev_seqnum<inputSeqCache->seq_count);
				scannedSeqToActualSeq[filterSeqCache->seq_count]=prev_seqnum;
				scannedSeqToOriginalStart[filterSeqCache->seq_count]=prev_start;
				/* we'll fill out the sq fields that come from the source seq later */
				sq->start=prev_start;
				sq->end=prev_stop;
				filterSeqCache->seq_count++;
			}
			
			prev_seqnum=seqnum;
			prev_start=start;
			prev_stop=stop;
		}
		assert(prev_seqnum>=0); /* since topHits->N, there should be something */
		
		/* add the last hit in */
		assert(prev_seqnum!=-1); // if there were no hits, we shouldn't be running this code
		ESL_SQ *sq=&(filterSeqCache->sq_list[filterSeqCache->seq_count]);
		scannedSeqToActualSeq[filterSeqCache->seq_count]=prev_seqnum;
		assert(prev_seqnum>=0 && prev_seqnum<inputSeqCache->seq_count);
		scannedSeqToOriginalStart[filterSeqCache->seq_count]=prev_start;
		/* we'll fill out the sq fields that come from the source seq later */
		sq->start=prev_start;
		sq->end=prev_stop;
		filterSeqCache->seq_count++;
		
		/* and now fill in the generic data that we get from the source seqs */
		for (i=0; i<filterSeqCache->seq_count; i++) {
			ESL_SQ *sq=&(filterSeqCache->sq_list[i]);
			int actualSeqnum=scannedSeqToActualSeq[i];
			ESL_SQ *src_sq=&(inputSeqCache->sq_list[actualSeqnum]);
			assert(actualSeqnum>=0 && actualSeqnum<inputSeqCache->seq_count);

			sq->L=src_sq->L;
			sq->abc=src_sq->abc;
			sq->acc=src_sq->acc;
			sq->desc=src_sq->desc;
			sq->idx=-1;
			sq->name=src_sq->name;
			sq->seq=NULL;
			sq->source=src_sq->source;
			sq->ss=NULL;
			
			sq->n=sq->end - sq->start + 1;
			
			sq->dsq=(ESL_DSQ *)MallocOrDie(sq->n+2);
			sq->dsq[0]=eslDSQ_SENTINEL;
			memcpy(sq->dsq+1,src_sq->dsq+sq->start,sq->n);
			sq->dsq[sq->n+1]=eslDSQ_SENTINEL;
		}
	}

	cm_tophits_Destroy(topHits);
	FreeCM(filter_cmsearch_cm);

	*get_filterSeqCache=filterSeqCache;
	*get_scannedSeqToActualSeq=scannedSeqToActualSeq;
	*get_scannedSeqToOriginalStart=scannedSeqToOriginalStart;
}
void ScanCand_FilterNonFrag_Free(ESL_SQCACHE *filterSeqCache,int *scannedSeqToActualSeq,int *scannedSeqToOriginalStart)
{
	if (filterSeqCache!=NULL) {
		int i;
		for (i=0; i<filterSeqCache->seq_count; i++) {
			free(filterSeqCache->sq_list[i].dsq);
		}
		if (filterSeqCache->seq_count>0) {
			free(filterSeqCache->sq_list);
		}
		free(filterSeqCache);
	}
	if (scannedSeqToActualSeq!=NULL) {
		free(scannedSeqToActualSeq);
	}
	if (scannedSeqToOriginalStart!=NULL) {
		free(scannedSeqToOriginalStart);
	}
}


/* ScanCand also incorporates code from CM_Search via 'cmsearch_wrapper'
 * because: unlike the original CMfinder code, I've decided to scan all sequences at once, 
 * since this is slightly more natural given the original cmsearch.c code.
 * Also, I've removed the 'range' functionality from the original CMfinder, so some
 * parts of code were removed, or simplified.
 */
void ScanCand (CM_t *input_cm,ESL_SQCACHE *inputSeqCache,int max_cand,Cand **cand, int *ncand,CmfinderVars *cmfinderVars)
{
	int hitnum,seqnum;
	int total_cand;
	CM_TOPHITS *topHits=NULL;
	double cm_scoreThreshold;
	CM_t *cmsearch_cm;
	ScanMode scanMode;
	int *scannedSeqToActualSeq;
	int *scannedSeqToOriginalStart;
	ESL_SQCACHE *filterSeqCache;
	ESL_SQCACHE *scanSeqCache;
	
	scanMode.force_inside_alg=0;
	scanMode.need_alignments=1;
	scanMode.disable_filters=0;
	scanMode.use_fragmentary=cmfinderVars->use_fragmentary;
	
	cm_scoreThreshold=ScanCand_GetScoreThreshold(cmfinderVars,input_cm,scanMode,inputSeqCache);

	cmsearch_cm=save_and_load_CM_via_memory_copying(input_cm,cmfinderVars); /* the hack to avoid having 'cmsearch' complain that it's already configured */
	if (cmfinderVars->use_evalue) {
		CM_t *cmcal_cm=save_and_load_CM_via_memory_copying(input_cm,cmfinderVars);
		cmcalibrate_viterbi(cmcal_cm,cmfinderVars,&(cmsearch_cm->expA)); /* cmcal_cm gets destroyed in this function */
		cmsearch_cm->flags |= CMH_EXPTAIL_STATS;
	}

	filterSeqCache=NULL; /* defaults are no filtering */
	scannedSeqToActualSeq=NULL;
	scannedSeqToOriginalStart=NULL;
	if (cmfinderVars->filterNonFrag) {
		ScanCand_FilterNonFrag(&filterSeqCache,&scannedSeqToActualSeq,&scannedSeqToOriginalStart,cmsearch_cm,inputSeqCache,cmfinderVars,cm_scoreThreshold,scanMode);
	}
	
	if (filterSeqCache!=NULL) {
		scanSeqCache=filterSeqCache;
	}
	else {
		scanSeqCache=inputSeqCache;
	}

	cmsearch_wrapper (cmsearch_cm,scanSeqCache,&topHits,cmfinderVars,cm_scoreThreshold,scanMode);

	/* initialize to zero cands per each seq */
	for (seqnum=0; seqnum<inputSeqCache->seq_count; seqnum++) {
		ncand[seqnum]=0;
	}
	
	if (cmfinderVars->bgScoreEvalue>0) {
		printf("--bg-score-evalue %lg : got %u hits\n",cmfinderVars->bgScoreEvalue,(unsigned int)(topHits->N));
	}

	/* for all the hits we got, add them as a cand of the appropriate seq (note: re-worked from CMfinder 0.3 code, since we now search all seqs at once) */
	for (hitnum = 0; hitnum < topHits->N; hitnum++) {
		int orig_start;
		CM_HIT *hit=topHits->hit[hitnum];
		if ((hit->flags&CM_HIT_IS_REMOVED_DUPLICATE)!=0) {
			continue;
		}
		if (hit->score<cm_scoreThreshold) {
			assert(cmfinderVars->use_evalue); /* this is the only case in which a hit below the score threshold could be reported, since -E and -T in cmsearch are incompatible */
			continue;
		}
		
		orig_start=1;
		seqnum=hit->seq_idx;
		if (scannedSeqToOriginalStart!=NULL) {
			orig_start=scannedSeqToOriginalStart[seqnum];
		}
		if (scannedSeqToActualSeq!=NULL) {
			seqnum=scannedSeqToActualSeq[seqnum];
		}
		
		if (ncand[seqnum]==maxNumCand) {
			/* reached maximum for seqnum, so ignore additional hits */
		}
		else {
			Cand *c=&(cand[seqnum][ncand[seqnum]]);
			ESL_DSQ *dsq=inputSeqCache->sq_list[seqnum].dsq;
			int i;
			assert(!hit->in_rc); /* CMfinder shouldn't be looking for rev-comp */
			assert(hit->start<=hit->stop);
			assert(hit->ad->sqfrom <= hit->ad->sqto); /* my understanding is that this should always be true (even if we had rev-comp hits) */
			
			c->weight=0;
			c->energy=0;

			if (cmfinderVars->use_evalue) {
				printf("evalue=%lg\n",hit->evalue);
				c->evalue=hit->evalue;
			}
			else {
				c->evalue=-1;
			}
			c->score=hit->score;
			c->cand_id=ncand[seqnum];
			c->seq_id=seqnum;
			c->start=hit->ad->sqfrom + orig_start-1;
			c->stop=hit->ad->sqto + orig_start-1;
			c->len=c->stop-c->start + 1;
			c->ad=hit->ad;
			hit->ad=NULL; /* detach the CM_ALIDISPLAY member from the hit */
			
			/* copied from Textize in esl_alphabet.c -- there isn't an easel function for a digital substring */
			for (i = 0; i < c->len; i++) {
				c->seq[i] = abc->sym[dsq[c->start+i]];  /* cand.start is already 1-based, so we don't need to add 1 */
			}
			c->seq[i] = '\0';
			c->ss[0]='\0'; /* no secondary structure yet */
			
			ncand[seqnum]++;
		}
	}
	
	/* For the hits in each seqnum,
	 * sort hits by score
	 * and take the best K cands from best-score-first, where K is limited by
	 * the number of hits we got (ncand[seqnum]), max_chosen_cand
	 * and the heuristics that we don't accepts scores < cmfinderVars->minCandScore bits, 
	 * and more than a factor of cmfinderVars->minCandScoreToBestScoreRatio worst than the best_score
	 */
	total_cand = 0;
	for (seqnum=0; seqnum<inputSeqCache->seq_count; seqnum++) {
		if (ncand[seqnum]>0) { /* code might get weird with allocations of 0 bytes, and there's nothing to do anyway if ncand[seqnum]==0 */
			Cand *tempCands = MallocOrDie(sizeof(Cand)* ncand[seqnum]);
			Cand **sort_cand = SortCand(cand[seqnum], ncand[seqnum], CompCandByScore);
			double best_score = sort_cand[0]->score;
			int j;
			for (j=0; j < ncand[seqnum] && j < max_cand; j++) {
				if (sort_cand[j]->score < best_score / cmfinderVars->minCandScoreToBestScoreRatio && sort_cand[j]->score < cmfinderVars->minCandScore) {
					break;
				}
				tempCands[j]=*(sort_cand[j]);
				tempCands[j].cand_id = j;
			}
			ncand[seqnum]=j;
			memcpy(cand[seqnum], tempCands, sizeof(Cand) * ncand[seqnum]);
			free(tempCands);
			free(sort_cand);
			
			total_cand += ncand[seqnum];
		}
	}

	EliminateHitsByVariousCriteria(cand,ncand,inputSeqCache,max_cand,cmfinderVars);
	
	/* if there's a big excess of candidates per sequence, then increase the cm_scoreThreshold, but only up to a limit of max_cm_scoreThreshold */
	if ((double)(total_cand) > (double)(inputSeqCache->seq_count) * cmfinderVars->overlyHighNumSeqToTotalCandRatio && cmfinderVars->cm_scoreThreshold < cmfinderVars->max_cm_scoreThreshold) {
		cmfinderVars->cm_scoreThreshold += cmfinderVars->cm_scoreThreshold_increment;
	}
	
	ScanCand_FilterNonFrag_Free(filterSeqCache,scannedSeqToActualSeq,scannedSeqToOriginalStart);
	cm_tophits_Destroy(topHits);
	FreeCM(cmsearch_cm);
}

double Compute_ZOOP_Weight(ESL_SQCACHE *inputSeqCache, int *ncand, Cand **cand, double *ret_totweight,CmfinderVars *cmfinderVars)
{
	static const double gamma0 = 0.3;  
	double  lambda;  
	double sum_prob_odd;
	double totscore =0 ;
	int    i,j;
	double totweight = 0;

	/* Calculate motif weight using ZOOP model*/
	for(i=0; i < inputSeqCache->seq_count; i++) {
		lambda  = cmfinderVars->zoop_gamma / inputSeqCache->sq_list[i].L;
		sum_prob_odd = 0;
		for (j = 0; j < ncand[i]; j++)
			sum_prob_odd += pow(2, cand[i][j].score);
		for(j = 0; j < ncand[i]; j++) {
			cand[i][j].weight =pow(2, cand[i][j].score)*lambda / (1 - cmfinderVars->zoop_gamma + sum_prob_odd * lambda);
			totweight += cand[i][j].weight;
			totscore += cand[i][j].score * cand[i][j].weight;
		}
	}
	cmfinderVars->zoop_gamma = 0;
	for(i=0; i < inputSeqCache->seq_count; i++) {
		for(j = 0; j < ncand[i]; j++) {
			cmfinderVars->zoop_gamma += cand[i][j].weight;
		}
	}
	cmfinderVars->zoop_gamma = (cmfinderVars->zoop_gamma + gamma0) / (inputSeqCache->seq_count + 1);
	*ret_totweight = totweight;
	return totscore;
}
typedef struct {
	int inputSeqNum,candNum;
} ScannedSeqToCand;
ESL_SQCACHE ESL_SQCACHE_FromCands__create (ESL_SQCACHE *inputSeqCache,int *ncand,Cand **cand,ScannedSeqToCand **ret_scannedSeqToCand,double minViterbiScoreToSkipInside)
{
	ESL_SQCACHE candSeqCache;
	int inputSeqNum,candNum;
	int outputSeqNum=0;
	int totalNumCand=GetTotalNumCand(inputSeqCache->seq_count,ncand,cand);
	candSeqCache.seq_count=totalNumCand; /* actually this is an upper bound on what we'll need; later we'll set the actual amount.  I expect the extra memory will be negligeable */
	candSeqCache.sq_list=(ESL_SQ *)MallocOrDie(sizeof(ESL_SQ)*totalNumCand);
	*ret_scannedSeqToCand=(ScannedSeqToCand *)MallocOrDie(sizeof(ScannedSeqToCand)*totalNumCand);
	for (inputSeqNum=0; inputSeqNum<inputSeqCache->seq_count; inputSeqNum++) {
		for (candNum=0; candNum<ncand[inputSeqNum]; candNum++) {
			Cand *c=&(cand[inputSeqNum][candNum]);
			ESL_SQ *sq=&(candSeqCache.sq_list[outputSeqNum]);
			ESL_SQ *src_sq=&(inputSeqCache->sq_list[inputSeqNum]);
			
			if (c->score >= minViterbiScoreToSkipInside) {
				continue; /* don't need to scan this with Inside alg, because it's Viterbi score is high enough. */
			}
			
			sq->L=src_sq->L;
			sq->abc=src_sq->abc;
			sq->acc=src_sq->acc;
			sq->desc=src_sq->desc;
			sq->idx=-1;
			sq->name=src_sq->name;
			sq->seq=NULL;
			sq->source=src_sq->source;
			sq->ss=NULL;
			
			sq->n=c->len;
			sq->start=c->start;
			sq->end=c->stop;
			
			sq->dsq=(ESL_DSQ *)MallocOrDie(c->len+2);
			sq->dsq[0]=eslDSQ_SENTINEL;
			memcpy(sq->dsq+1,src_sq->dsq+c->start,c->len);
			sq->dsq[c->len+1]=eslDSQ_SENTINEL;
			
			(*ret_scannedSeqToCand)[outputSeqNum].inputSeqNum=inputSeqNum;
			(*ret_scannedSeqToCand)[outputSeqNum].candNum=candNum;

			outputSeqNum++;
		}
	}
	candSeqCache.seq_count=outputSeqNum; /* we're not losing track of any mallocs this way, so it's okay */
	/* initialize other struct members just to avoid compiler warnings */
	candSeqCache.abc=abc;
	candSeqCache.hdr_size=0;
	candSeqCache.res_size=0;
	candSeqCache.filename=NULL;
	candSeqCache.format=0;
	candSeqCache.max_seq=candSeqCache.seq_count;
	candSeqCache.res_count=0;
	candSeqCache.header_mem=NULL;
	candSeqCache.residue_mem=NULL;
	return candSeqCache;
}
void ESL_SQCACHE_FromCands__destroy (ESL_SQCACHE *candSeqCache,ScannedSeqToCand *scannedSeqToCand)
{
	int i;
	for (i=0; i<candSeqCache->seq_count; i++) {
		free(candSeqCache->sq_list[i].dsq);
	}
	free(candSeqCache->sq_list);
	free(scannedSeqToCand);
}
double * FindBestScorePerCand(CM_TOPHITS *topHits,ScannedSeqToCand *scannedSeqToCand,int numCands)
{
	double *bestInsideScore=NULL;
	int hitnum,candSeqNum;

	/* we assume that the candSeqs might generate multiple hits in some cases (seems unlikely, but don't want to rule it out, 
	 * so we keep track of the best inside score */
	bestInsideScore=(double *)MallocOrDie(sizeof(double)*numCands);
	for (candSeqNum=0; candSeqNum<numCands; candSeqNum++) {
		bestInsideScore[candSeqNum]=DOUBLE_NEGINFINITY;
	}
	
	for (hitnum = 0; hitnum < topHits->N; hitnum++) {
		CM_HIT *hit=topHits->hit[hitnum];
		if ((hit->flags&CM_HIT_IS_REMOVED_DUPLICATE)!=0) {
			continue;
		}
		candSeqNum=hit->seq_idx;
		assert(!hit->in_rc); /* CMfinder shouldn't be looking for rev-comp */
		if (hit->score > bestInsideScore[candSeqNum]) {
			bestInsideScore[candSeqNum]=hit->score;
		}
	}
	
	return bestInsideScore;
}
double Compute_TCM_Weight(CM_t *cm, ESL_SQCACHE *inputSeqCache, int *ncand, Cand **cand, CmfinderVars *cmfinderVars,double *ret_totweight)
{
	/* variables from the original Compute_TCM_Weight function */
	CM_t *cmsearch_cm=save_and_load_CM_via_memory_copying(cm,cmfinderVars); /* the hack to avoid having 'cmsearch' complain that it's already configured */
	double lambda = 1.0 / cmfinderVars->DB_length;
	double totscore=0;
	double totweight = 0;
	double likelihood_ratio;
	double minViterbiScoreToSkipInside=sreLOG2(100.0/lambda); /* see bottom of page 40 of Zizhen's thesis */
	int i,j;
	
	/* we need the Inside scores of the candidates.  infernal's code has a nicely implemented pipeline, but it's not
	 * as obvious to me how to execute the internal functions, so I'd rather just use the pipeline, and just run it in inside mode.
	 * We'll use ESL_SQCACHE_FromCands__create to create a set of sequences to scan corresponding to the cand's, and
	 * then run the cmsearch pipeline on that.  If we originally use the Inside score, then we can simply use the existing score.
	 * (But I'd like to start off using Viterbi scores so I can compare to the old CMfinder results)
	 */
	int candSeqNum;
	CM_TOPHITS *topHits=NULL;
	double *bestInsideScore=NULL;
	ScannedSeqToCand *scannedSeqToCand=NULL;
	ESL_SQCACHE candSeqCache=ESL_SQCACHE_FromCands__create(inputSeqCache,ncand,cand,&scannedSeqToCand,minViterbiScoreToSkipInside);

	/* first, set all the cands that we don't have to do inside for */
	for (i=0; i<inputSeqCache->seq_count; i++) {
		for (j=0; j<ncand[i]; j++) {
			if (cand[i][j].score >= minViterbiScoreToSkipInside) {
				cand[i][j].weight = 1;
			}
		}
	}

	/* now set the ones that we did Inside for */
	if (candSeqCache.seq_count>0) {
		ScanMode scanMode;
		scanMode.use_fragmentary=cmfinderVars->use_fragmentary;
		scanMode.force_inside_alg=1; /* should always be 1, except when tested by using the Viterbi */
		scanMode.need_alignments=0;
		scanMode.disable_filters=1; /* if we use filters, they can sometimes reject a short cand seq, where the flanking sequence would increase the filter score allowing it to pass
		 * (see truncated_algs_notebook.doc entry for 2013-02-14).  Anyway, the filters and even envelope definition just restrict searches to subseqs, so once we already have
		 * a cand seq, there's no benefit to the filters.  we still want HMM/QDB-bands, but --max won't stop that */

#if 0
		/* compute Viterbi scores for cands, just for my debugging of the Inside scores */
		ScanMode vitScanMode=scanMode;
		vitScanMode.force_inside_alg=0;
		CM_TOPHITS *viterbiTopHits=NULL;
		double *bestViterbiRescore=NULL;
		cmsearch_wrapper (cmsearch_cm,&candSeqCache,&viterbiTopHits,cmfinderVars,cmfinderVars->cm_scoreThreshold,vitScanMode);
		bestViterbiRescore=FindBestScorePerCand(viterbiTopHits,scannedSeqToCand,candSeqCache.seq_count);
		cm_tophits_Destroy(viterbiTopHits);
		FreeCM(cmsearch_cm);
		cmsearch_cm=save_and_load_CM_via_memory_copying(cm,cmfinderVars); /* the hack to avoid having 'cmsearch' complain that it's already configured */
#endif

		cmsearch_wrapper (cmsearch_cm,&candSeqCache,&topHits,cmfinderVars,cmfinderVars->cm_scoreThreshold,scanMode); /* I assume Inside-score >= Viterbi-score, so we can use the same threshold */
		bestInsideScore=FindBestScorePerCand(topHits,scannedSeqToCand,candSeqCache.seq_count);

		/* phew, finally we get to the original TCM function */
		for (candSeqNum=0; candSeqNum<candSeqCache.seq_count; candSeqNum++) {
			i=scannedSeqToCand[candSeqNum].inputSeqNum;
			j=scannedSeqToCand[candSeqNum].candNum;
			assert(cand[i][j].score < minViterbiScoreToSkipInside); /* should have taken care of this above, in the non-Inside-alg case */
			if (bestInsideScore[candSeqNum] < cand[i][j].score) {
				/* assume that the hit got lost by the filters, so just warn and continue */
				if (!cmfinderVars->printedBestInsideScoreLessThanViterbi) {
					printf("encountered Inside score < Viterbi score.  Assuming something got lost in filters or envelope definition.  This is probably benign.\n");
					cmfinderVars->printedBestInsideScoreLessThanViterbi=1;
					/* leave cand[i][j].score as-is, just keep the Viterbi score */
					
					/*
					printf("bestInsideScore[candSeqNum] < cand[i][j].score (%lg<%lg)\n",bestInsideScore[candSeqNum],cand[i][j].score);
					fflush(stdout);
					assert(0);
					 */
				}
			}
			else {
				/* use Inside log likelihood. */
				cand[i][j].score = bestInsideScore[candSeqNum];
			}
			
			likelihood_ratio = sreEXP2(cand[i][j].score);
			cand[i][j].weight = likelihood_ratio * lambda / (1- lambda + likelihood_ratio * lambda);
		}

		cm_tophits_Destroy(topHits);
	}

	/* third, calculate totscore and totweight, regardless of how we set the weight (Viterbi or Inside) */
	for (i=0; i<inputSeqCache->seq_count; i++) {
		for (j=0; j<ncand[i]; j++) {
			totscore += cand[i][j].score * cand[i][j].weight;
			totweight += cand[i][j].weight;
		}
	}
	
	FreeCM(cmsearch_cm);
	free(bestInsideScore);
	ESL_SQCACHE_FromCands__destroy (&candSeqCache,scannedSeqToCand);

	*ret_totweight = totweight;
	return totscore;
}

void TrimDegenOnEndsOfMSA_AdjustPerSequenceSS (ESL_MSA *msa,int s,int p)
{
	if (msa->ss==NULL) {
		return;
	}
	if (msa->ss[s]==NULL) {
		return;
	}

	/* oops, it'll always be a gap, because the AddPerSequenceSsLines function won't put a pair with a degenerate nucleotide */
	if (esl_abc_CIsGap(abc,msa->ss[s][p])) {
		/* nothing to do */
	}
	else {
		/* remove the pair.  calling 'GetPairtable' is inefficient, but this will happen rarely, and should be trivial compared to other things, esp. since I'll generally use --max-degen-per-hit 2 */
		int *pt=GetPairtable(msa->ss[s]);
		assert(pt[p]>=0); /* otherwise it should have been a gap */
		msa->ss[s][p]='.';
		msa->ss[s][pt[p]]='.';
		free(pt);
	}
}

void TrimDegenOnEndsOfMSA (ESL_MSA *msa,Cand **best_msa_cand)
{
	int s,p;
	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}
	
	for (s=0; s<msa->nseq; s++) {
		
		/* from left to right */
		for (p=0; p<msa->alen; p++) {
			if (esl_abc_XIsDegenerate(abc,msa->ax[s][p+1])) {
				msa->ax[s][p+1]=esl_abc_XGetGap(abc);
				TrimDegenOnEndsOfMSA_AdjustPerSequenceSS (msa,s,p);
				if (best_msa_cand!=NULL) {
					best_msa_cand[s]->start++;
				}
			}
			else {
				if (esl_abc_XIsGap(abc,msa->ax[s][p+1])) {
					/* okay, continue looking for flanking */
				}
				else {
					/* we're done, we've gotten to a non-degen nuc on the end of the hit */
					break;
				}
			}
		}

		/* from right to left */
		for (p=msa->alen-1; p>=0; p--) {
			if (esl_abc_XIsDegenerate(abc,msa->ax[s][p+1])) {
				msa->ax[s][p+1]=esl_abc_XGetGap(abc);
				TrimDegenOnEndsOfMSA_AdjustPerSequenceSS (msa,s,p);
				if (best_msa_cand!=NULL) {
					best_msa_cand[s]->stop--;
				}
			}
			else {
				if (esl_abc_XIsGap(abc,msa->ax[s][p+1])) {
					/* okay, continue looking for flanking */
				}
				else {
					/* we're done, we've gotten to a non-degen nuc on the end of the hit */
					break;
				}
			}
		}
	}
}

void AddPerSequenceSsLines(ESL_MSA *best_msa,Parsetree_t **best_msa_tr,ESL_SQ **best_msa_sq,CM_t *cm)
{
	int i;
	int *dsqOffsetToMsaColumn=NULL;
	assert(best_msa->ss==NULL); /* no-one should have set it already, otherwise it'll get clobbered below, and we'll leak memory */
	best_msa->ss=(char **)MallocOrDie(sizeof(char*)*best_msa->nseq);
	dsqOffsetToMsaColumn=(int *)MallocOrDie(sizeof(int)*best_msa->alen+1); /* just alloc once with an upper bound of the memory we might need */
	for (i=0; i < best_msa->nseq; i++) {
		Parsetree_t *tr=best_msa_tr[i];
		int tidx,pos,aindex;
		char *ss;
		
		
		/* verify there aren't any gaps within the seq itself */
		for (pos=1; pos<=best_msa_sq[i]->L; pos++) {
			assert(!esl_abc_XIsGap(abc,best_msa_sq[i]->dsq[pos]));
		}
		
		/* now use gaps in MSA for this sequence to compute the dsqOffsetToMsaColumn map */
		if ((best_msa->flags & eslMSA_DIGITAL)==0) {
			esl_fatal("best_msa should be digital");
		}
		assert(best_msa->ax[i][0]==eslDSQ_SENTINEL); /* make sure it's in the expected digital format, with first nuc at position 1 */
		assert(best_msa->ax[i][best_msa->alen+1]==eslDSQ_SENTINEL);
		pos=1;
		for (aindex=1; aindex<=best_msa->alen; aindex++) {
			if (!esl_abc_XIsGap(abc,best_msa->ax[i][aindex])) {
				dsqOffsetToMsaColumn[pos]=aindex-1; /* -1 to get to the text string offset */
				pos++;
			}
		}
		
		ss=(char *)MallocOrDie(sizeof(char)*best_msa->alen+1);
		memset(ss,'.',best_msa->alen);
		for (tidx = 0; tidx < tr->n; tidx++) {  
			if(cm->sttype[tr->state[tidx]] == MP_st) {
				int left=tr->emitl[tidx];
				int right=tr->emitr[tidx];
				assert(left>=1 && left<=best_msa_sq[i]->L);
				assert(right>=1 && right<=best_msa_sq[i]->L);
				if (IsCanonicalDigitalBpAndNotDegen(best_msa_sq[i]->dsq[left],best_msa_sq[i]->dsq[right])) {
					ss[dsqOffsetToMsaColumn[left]]='<';
					ss[dsqOffsetToMsaColumn[right]]='>';
				}
			}
		}
		ss[best_msa->alen]=0;
		best_msa->ss[i]=ss;
	}
	free(dsqOffsetToMsaColumn);
}

void AddEvaluesToMsa(ESL_MSA *best_msa,Cand **best_msa_cand,CmfinderVars *cmfinderVars)
{
	int i;
	enum {MAX_EVALUE_STRING=64};
	char evalue_string[MAX_EVALUE_STRING];
	if (cmfinderVars->use_evalue) {
		for (i=0; i < best_msa->nseq; i++) {
			double evalue=best_msa_cand[i]->evalue;
			assert(evalue>=0); /* cmfinderVars->use_evalues should mean that we have valid E-values */
			sprintf(evalue_string,"%lg",evalue);
			if (esl_msa_AddGS(best_msa,"EVALUE", -1, i, evalue_string, -1) != eslOK) {
				esl_fatal("esl_msa_AddGS failed");
			}
		}
	}
	else {
		/* nothing to do, there's no evalues available */
	}
}

ESL_MSA * /*digital_msa */ E_step(CM_t *input_cm,ESL_SQCACHE *inputSeqCache,int *ncand,Cand **cand,double *ret_totscore,double ***ret_pxy,CmfinderVars *cmfinderVars)
{
	char             errbuf[eslERRBUFSIZE];
	int max_chosen_cand=cmfinderVars->max_cand_for_ScanCand;

	Cand        **chosen=NULL;
	ESL_SQ      **chosen_sq=NULL;       /* The chosen cands for each sequence        */
	char        **chosen_pp=NULL;
	double        *chosen_weight=NULL;
	Parsetree_t **chosen_tr=NULL;  /* tracebacks for each chosen cand              */

	Cand        **best_msa_cand = NULL;
	ESL_SQ      **best_msa_sq = NULL;       /* The chosen cands for each sequence        */
	char        **best_msa_pp=NULL; /* postProbAnnotationOfAlignment */
	char        **best_msa_pp_for_making_alignment=NULL;
	double       *best_msa_weight = NULL;
	Parsetree_t **best_msa_tr = NULL;  /* tracebacks for each chosen cand              */
	int           best_msa_nseq = 0;
	int nseq=inputSeqCache->seq_count;

	int      i, j;
	double   totscore = 0;
	double   totweight = 0;
	double   tot_weight = 0; /* used near the end of the function */
	int      chosen_idx=0;       /* The index of chosen candidate */
	ESL_MSA *digital_msa=NULL;

	ESL_STOPWATCH *timer_e_step=esl_stopwatch_Create();
	esl_stopwatch_Start(timer_e_step);

	ScanCand(input_cm,inputSeqCache,cmfinderVars->max_cand_for_ScanCand,cand,ncand,cmfinderVars);

	chosen = (Cand **) MallocOrDie( sizeof(Cand *) * nseq * max_chosen_cand);
	memset(chosen, 0, sizeof(Cand*) * nseq * max_chosen_cand);
	chosen_tr = (Parsetree_t **) MallocOrDie (sizeof(Parsetree_t *) *nseq * max_chosen_cand);
	chosen_pp = (char**) MallocOrDie (sizeof(char *) *nseq * max_chosen_cand);
	chosen_sq = (ESL_SQ **) MallocOrDie (sizeof(ESL_SQ *) *nseq * max_chosen_cand);

	totscore = 0;
	totweight = 0;
	chosen_idx = 0;

	if (cmfinderVars->do_zoop) {
		totscore=Compute_ZOOP_Weight(inputSeqCache, ncand, cand, &totweight,cmfinderVars);
	}
	else {
		totscore =  Compute_TCM_Weight(input_cm, inputSeqCache, ncand, cand, cmfinderVars, &totweight);
	}
	if (totscore > cmfinderVars->best_totscore) {
		cmfinderVars->best_totscore = totscore;
		best_msa_cand = (Cand **) MallocOrDie( sizeof(Cand *) * nseq);
		memset(best_msa_cand, 0, sizeof(Cand*) * nseq );
		best_msa_tr = (Parsetree_t **) MallocOrDie (sizeof(Parsetree_t *) * nseq);
		best_msa_pp = (char**) MallocOrDie (sizeof(char *) * nseq);
		best_msa_sq = (ESL_SQ **)MallocOrDie(sizeof(ESL_SQ *) * nseq);
	}

	for (i=0; i < nseq; i++) {    
		if (ncand[i] <= 0) continue;

		for(j=0; j <ncand[i] && j < max_chosen_cand; j++) {
			assert(cand[i][j].ad!=NULL); /* we should have taken it earlier in ScanCand */
			if (cm_alidisplay_Backconvert(input_cm,cand[i][j].ad,errbuf,&(chosen_sq[chosen_idx]),&(chosen_tr[chosen_idx]),&(chosen_pp[chosen_idx]))!=eslOK) {
				esl_fatal("cm_alidisplay_Backconvert failed: %s",errbuf);
			}
			
			/* don't need the .ad alidisplay member for this cand any more */
			cm_alidisplay_Destroy(cand[i][j].ad);
			cand[i][j].ad=NULL;

			chosen[chosen_idx] = &cand[i][j];

			if (j==0 && best_msa_cand) {
				int accept=1;
				
				if (cand[i][j].weight < cmfinderVars->min_acceptable_seq_weight) { /* NOTE: importantly, cand.weight does NOT reflect sequence weighting like GSC, which would make it much lower */
					/* reject this cand */
					/* printf("min_acceptable_seq_weight reject %s  , weight=%lg\n",cand[i][j].seq,cand[i][j].weight); */
					accept=0;
				}
				
				if (cand[i][j].score < cmfinderVars->min_seq_score_in_final_msa) {
					/* printf("min_seq_score_in_final_msa reject %s  , score=%lg\n",cand[i][j].seq,cand[i][j].score); */
					accept=0;
				}
				
				if (accept) {
					best_msa_cand[best_msa_nseq] =  &cand[i][j];
					best_msa_tr[best_msa_nseq] =  chosen_tr[chosen_idx];
					best_msa_sq[best_msa_nseq] = chosen_sq[chosen_idx];
					best_msa_pp[best_msa_nseq] = chosen_pp[chosen_idx];
					best_msa_nseq++;
				}
			}

			chosen_idx++;
		}
	}

	if (ret_totscore) {
		*ret_totscore = totscore;
	}

	tot_weight = 0;
	chosen_weight = (double*) malloc(sizeof(double) * chosen_idx);
	for(i=0; i < chosen_idx; i++) {
		chosen_weight[i] =  chosen[i]->weight;
		tot_weight += chosen[i]->weight;
	}
	if (tot_weight < cmfinderVars->min_acceptable_totweight) {
		printf("Too few sequences to predict any structure\n");
		cmfinderVars->cannotFindAcceptableStructure=1;
		digital_msa=NULL;
	}
	else {
		int *allmsa_seq_fragment_start_col,*allmsa_seq_fragment_end_col;

		if (Parsetrees2Alignment(input_cm,errbuf,abc,chosen_sq,chosen_weight,chosen_tr,chosen_pp,chosen_idx,NULL,NULL,TRUE,FALSE,&digital_msa)!=eslOK) {
			esl_fatal("Parsetrees2Alignment failed");
		}
		if ((digital_msa->flags & eslMSA_DIGITAL)==0) {
			esl_msa_Digitize(abc,digital_msa,errbuf);
		}
		MultiplyEMWeightsByWeightingStrategy (digital_msa,cmfinderVars->weightingStrategy,cmfinderVars->timer_weighting, cmfinderVars->use_fragmentary);
		TrimDegenOnEndsOfMSA (digital_msa,NULL);
	
		allmsa_seq_fragment_start_col=(int *)MallocOrDie(digital_msa->nseq * sizeof(int));
		allmsa_seq_fragment_end_col=(int *)MallocOrDie(digital_msa->nseq * sizeof(int));
		CalcFragmentStartEnds (digital_msa,cmfinderVars->use_fragmentary,allmsa_seq_fragment_start_col,allmsa_seq_fragment_end_col);
		digital_msa=ShiftDistalMispairsIntoTerminalLoops (digital_msa,cmfinderVars,allmsa_seq_fragment_start_col,allmsa_seq_fragment_end_col);

		if (best_msa_cand) {
			/* Save the best alignment so far */
			int *bestmsa_seq_fragment_start_col,*bestmsa_seq_fragment_end_col;
			char msaAuthor[MAX_MSA_AU];
			ESL_MSA *best_msa=NULL;
			
			/* calc tot_weight, set sequence descriptions and optionally set the hit names in the MSA to conform with CMfinder 0.3 */
			double tot_weight;
			best_msa_weight = (double*) malloc(sizeof(double) * best_msa_nseq);
			tot_weight =0;
			for (i=0; i < best_msa_nseq; i++) {
				if (cmfinderVars->putStartEndCoordsInHitIds) {
					/* this is the default for cm_alidisplay_Backconvert, so nothing to do */
					esl_fatal("oops, cmfinderVars->putStartEndCoordsInHitIds is problematic with TrimDegenOnEndsOfMSA, since the start/end coords go in, and _then_ the seqs get trimmed, and the coords have to be adjusted (although this is an issue only rarely, since more than a little bit of degen means that the seq gets killed).  We'd have to adjust TrimDegenOnEndsOfMSA so that it re-does the coords");
				}
				else {
					char *newname;
					char *oldname=best_msa_sq[i]->name;
					char *slash=strrchr(oldname,'/');
					assert(slash!=NULL); /* cm_alidisplay_Backconvert puts it in, so if it's not there, that's kinda weird and therefore ominous. */
					newname=(char *)MallocOrDie(slash-oldname+1);
					strncpy(newname,oldname,slash-oldname);
					newname[slash-oldname]=0;
					esl_sq_SetName(best_msa_sq[i],newname);
					free(newname);
				}
				best_msa_weight[i] =  best_msa_cand[i]->weight;
				tot_weight += best_msa_weight[i];
			}
			
			/* compute the alignment */
			if (1) {
				best_msa_pp_for_making_alignment=NULL;
			}
			else {
				best_msa_pp_for_making_alignment=best_msa_pp; /* this is how we could allow the user to have the #=GR ... PP lines (for posterior probability) */
			}
			if (Parsetrees2Alignment(input_cm,errbuf,abc,best_msa_sq,best_msa_weight,best_msa_tr,best_msa_pp_for_making_alignment,best_msa_nseq,NULL,NULL,TRUE,FALSE,&best_msa)!=eslOK) {
				esl_fatal("Parsetrees2Alignment failed");
			}
			if ((best_msa->flags & eslMSA_DIGITAL)==0) {
				esl_msa_Digitize(abc,best_msa,errbuf);
			}
			MultiplyEMWeightsByWeightingStrategy (best_msa,cmfinderVars->weightingStrategy,cmfinderVars->timer_weighting, cmfinderVars->use_fragmentary);
			bestmsa_seq_fragment_start_col=(int *)MallocOrDie(digital_msa->nseq * sizeof(int));
			bestmsa_seq_fragment_end_col=(int *)MallocOrDie(digital_msa->nseq * sizeof(int));
			CalcFragmentStartEnds (digital_msa,cmfinderVars->use_fragmentary,bestmsa_seq_fragment_start_col,bestmsa_seq_fragment_end_col);
			best_msa=ShiftDistalMispairsIntoTerminalLoops (best_msa,cmfinderVars,bestmsa_seq_fragment_start_col,bestmsa_seq_fragment_end_col);
			
			/* add per-sequence #=GR SS lines (note: Parsetrees2Alignment doesn't implement per-sequence SS annot, so we can't just add this to the individual sequence) */
			AddPerSequenceSsLines(best_msa,best_msa_tr,best_msa_sq,input_cm);
			
			TrimDegenOnEndsOfMSA (best_msa,best_msa_cand); /* note: it's important that we do this after AddPerSequenceSsLines, since that function counts non-gap characters to figure out where it is in the MSA */

			assert(best_msa_nseq==best_msa->nseq);
			for (i=0; i<best_msa->nseq; i++) {
				if (esl_msa_FormatSeqDescription(best_msa,i, "%d..%d\t%lg",best_msa_cand[i]->start,best_msa_cand[i]->stop,best_msa_cand[i]->score)!=eslOK) {
					esl_fatal("problem with esl_msa_FormatSeqDescription");
				}
			}
			
			AddEvaluesToMsa(best_msa,best_msa_cand,cmfinderVars);
			
			/* add other infor to msa */
			sprintf(msaAuthor,"CMfinder (infernal %s)",INFERNAL_VERSION);
			if (esl_msa_SetAuthor(best_msa,msaAuthor,-1)!=eslOK) {
				esl_fatal("esl_msa_SetAuthor failed");
			}
			
			/* keep this alignment as best-so-far */
			if (cmfinderVars->best_msa!=NULL) {
				esl_msa_Destroy(cmfinderVars->best_msa);
			}
			cmfinderVars->best_msa=best_msa;
		}
	}

	if (best_msa_cand) {
		/* and free memory */
		free(best_msa_cand);
		free(best_msa_sq);
		free(best_msa_weight);
		free(best_msa_tr);
		free(best_msa_pp);
	}

	if (ret_pxy) {
		*ret_pxy=NULL;
		/* investigate this stuff later.  It looks like the original CMfinder code is trying to calculate partition function probs
		* only once for overlapping candidates, and I think it's probably also caching partition function probs across iterations, since
		* the same thing gets calculated if the candidates are in the same or nearby locations.  The latter looks plausible given that the
		* variables are declared 'static'.  However, I'd rather get things basically working first. */
	}
	
	for (i = 0; i < chosen_idx; i++) {
		FreeParsetree(chosen_tr[i]);
		esl_sq_Destroy(chosen_sq[i]);
		free(chosen_pp[i]);
	} 
	free(chosen_pp);
	free(chosen_tr);
	free(chosen_weight);
	free(chosen_sq);
	free(chosen);
	
	esl_stopwatch_Stop(timer_e_step);
	esl_stopwatch_Include(cmfinderVars->timer_e_step,timer_e_step);
	esl_stopwatch_Destroy(timer_e_step);

	return digital_msa;
}
