/* cmfinder.c
* 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include "squid.h"		/* general sequence analysis library    */
#include "msa.h"                /* squid's multiple alignment i/o       */
#include "structs.h"		/* data structures, macros, #define's   */
#include "funcs.h"		/* external functions                   */
#include "version.h"            /* versioning info for Infernal         */
#include "prior.h"
#include "cand.h"
#include "global.h"
#include "histogram.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif



extern CM_t * Automodelmaker(MSA *msa, char** dsq,  int ** pxy,double gapthresh, int gapcost, int trim, int* constraint);
extern void prxy(char **aseq,int nseq,int alen,double **bp_pr,float *weights,int ***ret_prxy);
extern double* bppr_seq(char* seq);

extern int    do_local;
extern int    do_small;
extern int    do_binary;
extern int    do_banded;
extern int    do_hbanded;
extern int    do_hmm;
extern int    do_evd;
extern int    do_zoop;
extern int    do_tcm;
extern int    do_inside;

extern int    proximity;
extern int    DB_length;
extern double mu;
extern double lambda;
extern char  *cmfile;
extern char  *final_file; 

int    watson_crick=1;
int use_sums = 0;
float    hmm_threshold = 0;
float    cm_threshold = 0;
int debug_level = 0;
float      etarget = 1.46;


MSA * Parsetrees2AlignmentSS(CM_t *cm, char **dsq, SQINFO *sqinfo, float *wgt, Parsetree_t **tr, int nseq);
extern void actually_write_stockholm(FILE *fp, MSA *msa, int cpl);


void CM_save(CM_t *cm, char* cmfile)
{
	FILE      *cmfp;		  /* OUTPUT: fp to cvfile                      */  
	if ((cmfp = fopen(cmfile, "w")) == NULL)
		Die("Failed to open %s for writing", cmfile);
	CMFileWrite(cmfp, cm, do_binary);
	fclose(cmfp);  
}


void MSA_save(MSA *msa, char* msafile)
{
	FILE *fout;            
	if ((fout = fopen(msafile, "w")) == NULL)   Die("Can't write to output file ");      
#if 1
	actually_write_stockholm(fout, msa, 5000);
#else
	WriteStockholm(fout, msa);
#endif
	fclose(fout);
}


/*Zizhen: Change insert transition and emission scores at the root to 0*/
void
CMHackRootScores(CM_t *cm)
{
	int v, x,y;
	for (v = 0; v < cm->M; v++)
	{
		if (cm->stid[v]  == ROOT_S || cm->stid[v]  == ROOT_IL || cm->stid[v]  == ROOT_IR){
			for(x=0; x < cm->cnum[v]; x++){
				if (cm->sttype[cm->cfirst[v] + x] == IL_st || cm->sttype[cm->cfirst[v]+ x] == IR_st){
					cm->tsc[v][x]=0;
				}
				else{
					if (v != 0){
						for(y=0; y < cm->cnum[0]; y++){
							if (cm->cfirst[v] + x == cm->cfirst[0]+y){
								cm->tsc[v][x] = cm->tsc[0][y];
								break;
							}
						}
					}
				}
			}
		}
		else break;
	}
} 


HMM_Band* HMM_CM_Band(CM_t *cm, HMM* hmm, HMM_CM_Map* hmm_cm_map, char *dsq, int i, int j, int do_search)
{
	/* Step 1: Get HMM posteriors.*/
	float fwd_score, bck_score;
	HMM_Matrix *hmm_mx;        /* growable DP matrix for viterbi                       */
	HMM_Matrix *hmm_fwd;       /* growable DP matrix for forward                       */
	HMM_Matrix *hmm_bck;       /* growable DP matrix for backward                      */
	HMM_Matrix *hmm_posterior; /* growable DP matrix for posterior decode              */
	HMM_Band   *hmm_band;
	float hbandp = 0.0001;
	int  v;

	fwd_score = CP9Forward(dsq, i, j, hmm, &hmm_fwd);
	bck_score = CP9Backward(dsq, i, j, hmm, &hmm_bck);

	hmm_posterior = hmm_bck;
	CP9FullPosterior(dsq, i, j, hmm, hmm_fwd, hmm_bck, hmm_posterior);


	hmm_band = AllocCP9Bands(cm, hmm);

	if(use_sums)
		CP9_ifill_post_sums(hmm_posterior, i, j, hmm_band->hmm_M,
		hmm_band->isum_pn_m, hmm_band->isum_pn_i, hmm_band->isum_pn_d);


	/* match states */
	CP9_hmm_band_bounds(hmm_posterior->mmx, i, j, hmm_band->hmm_M,
		hmm_band->isum_pn_m, hmm_band->pn_min_m, hmm_band->pn_max_m,
		(1.-hbandp), HMMMATCH, use_sums, debug_level);
	/* insert states */
	CP9_hmm_band_bounds(hmm_posterior->imx, i, j, hmm_band->hmm_M,
		hmm_band->isum_pn_i, hmm_band->pn_min_i, hmm_band->pn_max_i,
		(1.-hbandp), HMMINSERT, use_sums, debug_level);
	/* delete states */
	CP9_hmm_band_bounds(hmm_posterior->dmx, i, j, hmm_band->hmm_M,
		hmm_band->isum_pn_d, hmm_band->pn_min_d, hmm_band->pn_max_d,
		(1.-hbandp), HMMDELETE, use_sums, debug_level);

	/* Step 3: HMM bands  ->  CM bands. */
	hmm2ij_bands(cm, hmm_cm_map, i, j, hmm_band->pn_min_m, hmm_band->pn_max_m, 
		hmm_band->pn_min_i, hmm_band->pn_max_i, hmm_band->pn_min_d, hmm_band->pn_max_d, 
		hmm_band->imin, hmm_band->imax, hmm_band->jmin, hmm_band->jmax, debug_level);

	if(do_search){
		relax_root_bands(hmm_band->imin, hmm_band->imax, hmm_band->jmin, hmm_band->jmax);		    
	}

	/* Use the CM bands on i and j to get bands on d, specific to j. */
	for(v = 0; v < cm->M; v++)
	{
		hmm_band->hdmin[v] = malloc(sizeof(int) * (hmm_band->jmax[v] - hmm_band->jmin[v] + 1));
		hmm_band->hdmax[v] = malloc(sizeof(int) * (hmm_band->jmax[v] - hmm_band->jmin[v] + 1));
	}
	ij2d_bands(cm, j - i + 1 , hmm_band->imin, hmm_band->imax, 
		hmm_band->jmin, hmm_band->jmax, hmm_band->hdmin, hmm_band->hdmax, -1);

	FreeCPlan9Matrix(hmm_fwd);
	FreeCPlan9Matrix(hmm_bck);
	return hmm_band;
}  

float CM_align(CM_t *cm, char *dsq, int L, int* dmin, int* dmax, 
			   HMM* hmm, HMM_CM_Map* hmm_cm_map, 
			   Parsetree_t **ret_tr)
{
	int   debug_level = 0;
	Parsetree_t    *tr;          /* parse trees for the sequences */
	float sc;
	if (hmm && do_hbanded){    
		HMM_Band* hmm_band = HMM_CM_Band(cm,  hmm, hmm_cm_map, dsq,1,L, 0);
		sc = CYKInside_b_jd(cm, dsq, L, 0, 1, L, &tr, 
			hmm_band->jmin, hmm_band->jmax, hmm_band->hdmin, hmm_band->hdmax, 
			hmm_band->safe_hdmin, hmm_band->safe_hdmax);
		*ret_tr = tr;
		FreeCP9Bands(hmm_band);
		free(hmm_band);
		return sc;
	}

	sc = CYKDivideAndConquer(cm, dsq, L, 0, 1, L, &tr, dmin, dmax);  
	*ret_tr= tr;
	return sc;
}



CM_t * M_step (MSA* msa, char **dsq, Prior_t *pri, float *null,int** pxy, float gapthreshold, int gapcost, HMM** ret_hmm, HMM_CM_Map** ret_map)
{ 
	CM_t      *cm;  
	int        i; 
	FILE      *cmfp;		  /* OUTPUT: fp to cvfile                      */
	float     *weight;
	float      tot_weight=0;

	cm = Automodelmaker( msa, dsq, pxy, gapthreshold, gapcost,1, NULL);
	if (cm == NULL){
		return NULL;
	}

	/* Convert to probabilities, and the global log-odds form
	* we save the model in.
	*/  

	for(i=0; i < msa->nseq; i++)
		tot_weight += msa->wgt[i];

	//float eff_nseq = CM_Eweight(cm, pri, (float) msa->nseq, etarget);
	//CMRescale(cm, eff_nseq / (float) msa->nseq);

	CMSetNullModel(cm, null);  
	PriorifyCM(cm, pri, NULL);  
	CMLogoddsify(cm);
	CMHackInsertScores(cm);	/* "TEMPORARY" fix for bad priors */  
	CM_save(cm, cmfile);  
	if (do_hmm){
		HMM  * hmm;
		HMM_CM_Map *map;
		if(!build_cp9_hmm(cm, &hmm, &map, FALSE, 0.0001, 0))
			Die("Couldn't build a CP9 HMM from the CM\n");           
		*ret_hmm = hmm;
		*ret_map = map;    
	}
	return cm;  
}

int TrimCand(CM_t *cm, Parsetree_t *tr, Cand* cand, char* cand_dsq)
{
	int tpos;
	int start = 1;
	int stop = cand->len;
	int length;
	int new_root=0;
	int i;
	int nxtl=-1;
	int nxtr=-1;
	char* temp;
	for(tpos =0 ; tpos < tr->n; tpos++){
		int v = tr ->state[tpos];
		if (cm->stid[v]  == ROOT_S ||cm->stid[v]  == ROOT_IR ||cm->stid[v]  == ROOT_IL ) {
			if (cm->stid[v]  == ROOT_IL){ 
				start = tr->emitl[tpos] +1;
			}    
			else if (cm->stid[v]  == ROOT_IR){
				stop = tr->emitr[tpos] -1;
			}          
			nxtl = tr->nxtl[tpos];
			nxtr = tr->nxtr[tpos];
		}
		else{
			new_root= tpos;      
			break;
		}
	}
	length = stop - start + 1;  
	if (length < cand->len){
		//printf("Before TrimCand start %d stop %d\n", 1, cand->len);
		//ParsetreeDump(stdout, tr, cm, cand_dsq);    
		//Trim parsetree
		for(tpos=new_root; tpos < tr->n; tpos++){
			int new_pos = tpos - new_root + 1;
			tr->emitl[new_pos] = tr->emitl[tpos] - start + 1;
			tr->emitr[new_pos] = tr->emitr[tpos] - start + 1;
			tr->state[new_pos] = tr->state[tpos];
			if ( tr->nxtl[tpos] == -1){
				tr->nxtl[new_pos] = -1;
			}
			else{
				tr->nxtl[new_pos] = tr->nxtl[tpos] - new_root;
			}
			if ( tr->nxtr[tpos] == -1){
				tr->nxtr[new_pos] = -1;
			}
			else{
				tr->nxtr[new_pos] = tr->nxtr[tpos] - new_root;
			}
			if ( tr->prv[tpos] == -1){
				tr->prv[new_pos] = -1;
			}
			else{
				tr->prv[new_pos] = tr->prv[tpos] - new_root;
			}
		}
		tr->n -= new_root;            
		tr->emitl[0] = 1;
		tr->emitr[0] = stop - start + 1;    
		if (nxtl != -1)    tr->nxtl[0] = nxtl- new_root;
		if (nxtr != -1)    tr->nxtr[0] = nxtr - new_root;

		if (cand_dsq){
			for(i=start; i <= stop; i++){
				cand_dsq[i - start + 1] = cand_dsq[i];
			}
			cand_dsq[stop+1] = cand_dsq[cand->len + 1];        
		}
		//printf("After TrimCand %d %d\n", start, stop);
		//ParsetreeDump(stdout, tr, cm, cand_dsq);

		//printf("Trim %d %d into %d %d\n", cand->start, cand->stop, cand->start+start-1, cand->start + stop-1);
		temp = MallocOrDie(sizeof(char) * (length + 1));
		strncpy(temp, cand->seq + start -1, length);
		temp[length]='\0';
		cand->start  += start - 1;
		cand->stop = cand->start + length - 1;  
		cand->len = length;
		strcpy(cand->seq, temp);
		return 1;
	}
	return 0;
}

/* get a sub matrix of a triangular matrix*/
double* get_submatrix(int start, int len, double* orig_tri_matrix)
{
	int i,j;
	double* new_tri_matrix = (double*) MallocOrDie(sizeof(double) * (TriIndex(len, len-1) - 1));
	for(j=1; j <  len ;j++)
		for(i=0; i < j; i++){
			int coor1= TriIndex(i,j);
			int coor2 = TriIndex(i + start, j+ start);
			new_tri_matrix[coor1] = orig_tri_matrix[coor2];  
		}
		return new_tri_matrix;  
}



double Compute_ZOOP_Weight(int nseq, SQINFO *sqinfo, int* ncand, Cand** cand, Range* range, double* ret_totweight)
{
	static const double gamma0 = 0.3;  
	static double gamma = 0.3;
	double  lambda;  
	double sum_prob_odd;
	double totscore =0 ;
	int    i,j;
	double totweight = 0;

	/* Calculate motif weight using ZOOP model*/
	for(i=0; i < nseq; i++){
		lambda  = gamma / sqinfo[i].len;    
		if (range) {
			lambda = gamma/ (range[i].stop - range[i].start + 1);
		}
		sum_prob_odd = 0;    
		for (j = 0; j < ncand[i]; j++) 
			sum_prob_odd += pow(2, cand[i][j].score); 	
		for(j = 0; j < ncand[i]; j++) {
			cand[i][j].weight =pow(2, cand[i][j].score)*lambda / (1 - gamma + sum_prob_odd * lambda);			
			totweight += cand[i][j].weight;
			totscore += cand[i][j].score * cand[i][j].weight;	
		}
	}
	gamma = 0;  
	for(i=0; i < nseq; i++) {    
		for(j = 0; j < ncand[i]; j++) {      
			gamma += cand[i][j].weight;      
		}    
	}  
	gamma  = (gamma + gamma0) / (nseq + 1);  
	*ret_totweight = totweight;
	return totscore;
}


double Compute_TCM_Weight(CM_t *cm, int nseq, SQINFO *sqinfo, int* ncand, Cand** cand, double* ret_totweight)
{
	double  lambda = 1.0 / (DB_length);
	double sum_prob_odd;
	double totscore=0;
	double totweight = 0;
	double likelihood_ratio;

	int    i,j;          

	/* Calculate motif weight using ZOOP model*/
	for(i=0; i < nseq; i++){
		for (j = 0; j < ncand[i]; j++) {
			if (cand[i][j].score > sreLOG2(100/lambda)){
				cand[i][j].weight = 1;
			}
			else{
				//Inside-outside log likelihood. 
				char* dsq = DigitizeSequence(cand[i][j].seq, cand[i][j].len);	
				cand[i][j].score = FInside(cm, dsq, cand[i][j].len, 1, cand[i][j].len, 0, NULL, NULL,NULL,NULL, 0);
				free(dsq);
				likelihood_ratio = sreEXP2(cand[i][j].score);
				cand[i][j].weight = likelihood_ratio * lambda / (1- lambda + likelihood_ratio * lambda);
			}
			totscore += cand[i][j].score * cand[i][j].weight;		  
			totweight += cand[i][j].weight;		  
		}
	}  
	*ret_totweight = totweight;
	return totscore;
}

void train_EVD(struct cm_s *cm)
{  
	int nsample= 20;
	int idx;
	char* seq;
	char* dsq;

	float score;
	float background[4]={0.25, 0.25, 0.25, 0.25};
	struct histogram_s *hist;
	int i;

	float sqlen =  1.5 * cm->nodes;
	hist = AllocHistogram(-200, 200, 100);  
	for (idx = 0; idx < nsample; idx++)
	{
		/* generate it */
		seq = RandomSequence(Alphabet, background, Alphabet_size, sqlen);
		dsq = DigitizeSequence(seq, sqlen);	
		score = CYKInsideScore(cm, dsq, 0, 1, sqlen, sqlen, NULL, NULL);      
		AddToHistogram(hist, score,1);
		free(dsq);
		free(seq);
	}  
	if (! ExtremeValueFitHistogram(hist, TRUE, 9999.)){
		if (! ExtremeValueFitHistogram(hist, FALSE, 9999.)){
			PrintASCIIHistogram(stdout, hist);
			Die("EVD fit failed");
		}
	}

	mu    = hist->param[EVD_MU];
	lambda  = hist->param[EVD_LAMBDA];  
	hmm_threshold =  (log(DB_length))/lambda + mu;
}

double Compute_Evalue_Weight(int nseq, SQINFO *sqinfo, int* ncand, Cand** cand, double* ret_totweight)
{
	int i,j;
	double totscore=0;
	double totweight = 0;
	static double min_weight = 3;

	printf("mu %f\t lambda %f\n", mu, lambda);
	for(i=0; i < nseq; i++){
		int length = sqinfo[i].len;
		if (length < DB_length){
			length = DB_length;
		}
		for (j = 0; j < ncand[i]; j++) {
			double pval= ExtremeValueP(cand[i][j].score, mu, lambda);      
			cand[i][j].weight = exp(-length * pval);      
			totweight += cand[i][j].weight;      
			totscore += cand[i][j].weight * cand[i][j].score;
		}
	}
	/* Inflate the weights if they are too small */
	if (totweight < min_weight){
		double scale = min_weight/ totweight;
		for(i=0; i < nseq; i++){
			for (j = 0; j < ncand[i]; j++) {
				cand[i][j].weight *= scale;
				if (cand[i][j].weight > 1) cand[i][j].weight = 1;	
			}
		}
		totscore *= scale;
	}
	*ret_totweight = totweight;
	return totscore;
}


MSA* E_step(struct cm_s *cm, 
			HMM* hmm,
			HMM_CM_Map* hmm_cm_map,
			int *dmin,
			int *dmax,
			int nseq, 
			SQINFO *sqinfo, 
			double * seq_weight, 
			int* ncand, Cand** cand, 
			int max_chosen_cand,
			Range* range,  
			double *ret_score, 
			int*** ret_pxy, 
			FILE* scan_fout,
			double*  best_totscore)
{  
	static const double gamma0 = 0.3;  
	static double gamma = 0.3;
	/* The computed partition functions for all segments in all segments.
	Introduced to save computational time*/
	static double*** all_pf_segment = NULL;
	static int*  num_pf_segment = NULL;
	static int** pf_segment_start = NULL;
	static int** pf_segment_len = NULL;

	Cand        **chosen;  
	char        **chosen_dsq;         /* two dim array of all cands for each sequence */
	SQINFO       *chosen_sqinfo;       /* The chosen cands for each sequence        */
	float        *chosen_weight;       
	Parsetree_t **chosen_tr;  /* tracebacks for each chosen cand              */

	Cand        **best_msa_cand = NULL;  
	char        **best_msa_dsq = NULL;         /* two dim array of all cands for each sequence */
	SQINFO       *best_msa_sqinfo = NULL;       /* The chosen cands for each sequence        */
	float        *best_msa_weight = NULL;       
	Parsetree_t **best_msa_tr = NULL;  /* tracebacks for each chosen cand              */
	int           best_msa_nseq = 0;

	int     i, j;  
	double  bestscore;  
	double  totscore = 0;    
	double  totweight = 0;    
	double  tot_weight = 0; /* used near the end of the function */
	int     chosen_idx=0;       /* The index of chosen candidate */
	Cand**  sort_cand;  
	int     free_dsq=0;
	MSA    *msa;
	double  lambda;  
	Fancyali_t      *ali;         /* alignment, formatted for display */   
	CMConsensus_t   *cons;	/* precalculated consensus info for display */

	chosen = (Cand **) malloc( sizeof(Cand *) * nseq * max_chosen_cand);
	memset(chosen, 0, sizeof(Cand*) * nseq * max_chosen_cand);  
	chosen_tr = (Parsetree_t **) malloc (sizeof(Parsetree_t *) *nseq * max_chosen_cand);
	chosen_dsq = (char**) malloc (sizeof(char *) *nseq * max_chosen_cand);

	totscore = 0; 
	totweight = 0;
	chosen_idx = 0;

	if (do_zoop){
		totscore=Compute_ZOOP_Weight(nseq, sqinfo, ncand, cand,range, &totweight);
	}
	else if (do_evd){
		totscore= Compute_Evalue_Weight(nseq, sqinfo, ncand, cand, &totweight);
	}
	else{
		totscore =  Compute_TCM_Weight(cm, nseq, sqinfo, ncand, cand, &totweight);
	}
	if (totscore > *best_totscore){
		*best_totscore = totscore;
		best_msa_cand = (Cand **) malloc( sizeof(Cand *) * nseq);
		memset(best_msa_cand, 0, sizeof(Cand*) * nseq );  
		best_msa_tr = (Parsetree_t **) malloc (sizeof(Parsetree_t *) *nseq);
		best_msa_dsq = (char**) malloc (sizeof(char *) *nseq);    
	}

	for(i=0; i < nseq; i++) {    
		if (ncand[i] <= 0) continue;
		bestscore = cand[i][0].score;

		if (scan_fout){
			fprintf(scan_fout, "sequence: %s\n", sqinfo[i].name);    
			cons = CreateCMConsensus(cm, 3.0, 1.0); 
		}

		for(j=0; j <ncand[i] && j < max_chosen_cand; j++){
			double old_cand_score=cand[i][j].score;
			double cand_score_change;
			char* cand_dsq = DigitizeSequence(cand[i][j].seq, cand[i][j].len);      
			cand[i][j].score= CM_align(cm, cand_dsq, cand[i][j].len, dmin, dmax, hmm, hmm_cm_map,  &chosen_tr[chosen_idx]);
			cand_score_change=old_cand_score-cand[i][j].score;
			if (fabs(cand_score_change)>1e-4) {
				int q=9;
			}

			chosen_dsq[chosen_idx] = cand_dsq;
			chosen[chosen_idx] = &cand[i][j];

			if (j==0 && best_msa_cand){
				best_msa_cand[best_msa_nseq] =  &cand[i][j];	
				best_msa_tr[best_msa_nseq] =  chosen_tr[chosen_idx];
				best_msa_dsq[best_msa_nseq] = chosen_dsq[chosen_idx];
				best_msa_nseq++;
			}

			/*
			printf("Cand %d-%d Score %.2f weight %f\n", 
			chosen[chosen_idx]->seq_id,j, chosen[chosen_idx]->score, chosen[chosen_idx]->weight);
			*/
			chosen_idx++;        
			//TrimCand(cm, chosen_tr[chosen_idx], &cand[i][j], cand_dsq);	            
			if (scan_fout){
				fprintf(scan_fout, "hit %-4d: %6d %6d %8.2f bits\n", j,
					cand[i][j].start, cand[i][j].stop, cand[i][j].score);	
				ali = CreateFancyAli(chosen_tr[chosen_idx], cm, 
					cons, cand_dsq);
				PrintFancyAli(scan_fout, ali);
				FreeFancyAli(ali);
			}      
		}
	}

	for(i=0; i < nseq; i++){
		/*Update the scan range */
		if (range && range[i].start > 0){
			for(j=0; j < ncand[i]; j++){
				if (cand[i][j].weight > 0.9){
					if (range[i].start != -1 && cand[i][j].start - proximity< range[i].start ){
						range[i].start = cand[i][j].start - proximity;
						if (range[i].start < 1) 
							range[i].start = 1;
					}
					if (range[i].stop != -1 && cand[i][j].stop + proximity > range[i].stop ){
						range[i].stop = cand[i][j].stop + proximity;
					}
					if ((range[i].start == -1 || range[j].stop == -1) && cand[i][j].weight > 0.9){
						range[i].start = cand[i][j].start - FLANK; 
						range[i].stop = cand[i][j].stop + FLANK; 
					}
				}
			}
		}
	}
	if (ret_score) *ret_score = totscore;  


	/* transform Cand data structure to SQINFO */
	chosen_sqinfo = Cand2Sqinfo(chosen, chosen_idx, sqinfo);  
	chosen_weight = (float*) malloc(sizeof(float) * chosen_idx);
	tot_weight = 0;
	for(i=0; i < chosen_idx; i++) {
		chosen_weight[i] =  chosen[i]->weight;
		tot_weight += chosen[i]->weight;
		if (seq_weight) chosen_weight[i] *= seq_weight[chosen[i]->seq_id];    
	}

	if (tot_weight < MIN_WEIGHT){
		Die("Too few sequences");
	}

	msa = Parsetrees2AlignmentSS(cm, chosen_dsq, chosen_sqinfo, chosen_weight, chosen_tr, chosen_idx);  

	/* Save the best alignment so far */
	if (best_msa_cand){
		MSA *best_msa;
		double tot_weight;
		best_msa_sqinfo = Cand2Sqinfo(best_msa_cand, best_msa_nseq, sqinfo);
		best_msa_weight = (float*) malloc(sizeof(float) * best_msa_nseq);
		tot_weight =0;
		for(i=0; i < best_msa_nseq; i++) {
			best_msa_weight[i] =  best_msa_cand[i]->weight;
			tot_weight += best_msa_weight[i];
			if (seq_weight) best_msa_weight[i] *= seq_weight[best_msa_cand[i]->seq_id];          
		}
		best_msa = Parsetrees2AlignmentSS(cm, best_msa_dsq, best_msa_sqinfo, best_msa_weight, best_msa_tr, best_msa_nseq);  
		if (final_file){
			MSA_save(best_msa, final_file);
		}
		MSAFree(best_msa);
		if (cmfile){
			CM_save(cm, cmfile);
		}
		free(best_msa_cand);
		free(best_msa_sqinfo);
		free(best_msa_weight);
		free(best_msa_tr);
		free(best_msa_dsq);
	}

	if (ret_pxy){
		/*update partition function segments*/
		double** bp_pr;
		if (all_pf_segment==NULL){
			all_pf_segment = (double***)MallocOrDie(sizeof(double**) * nseq);    
			num_pf_segment = (int*) MallocOrDie(sizeof(int) * nseq);    
			pf_segment_start = (int**) MallocOrDie(sizeof(int*) * nseq);    
			pf_segment_len = (int**) MallocOrDie(sizeof(int*) * nseq);    
			for(i= 0; i < nseq; i++){
				all_pf_segment[i] = (double**)MallocOrDie(sizeof(double*) * MAXCAND);    
				pf_segment_start[i] = (int*)MallocOrDie(sizeof(int) * MAXCAND);    
				pf_segment_len[i] = (int*)MallocOrDie(sizeof(int) * MAXCAND);    
				num_pf_segment[i] = 0;      
			}    
		}
		bp_pr= (double**) MallocOrDie(sizeof(double*) * chosen_idx);
		for(i=0; i < chosen_idx; i++){
			int seq_id = chosen[i]->seq_id;
			int start = chosen[i]->start;
			int stop = chosen[i]->stop;
			int match = -1;
			/* We would like to extract the partition function of this segment from existing partition function*/	
			for(j=0; j < num_pf_segment[seq_id]; j++){
				int loffset = start - pf_segment_start[seq_id][j];
				int roffset = pf_segment_len[seq_id][j] + pf_segment_start[seq_id][j] -1 - stop;      
				if (loffset >=0 && roffset >= 0 && loffset + roffset <= 10){	
					match = j;
					bp_pr[i] = get_submatrix(loffset, chosen[i]->len, all_pf_segment[seq_id][j]); 
					break;
				}
			} 
			if (match == -1){
				bp_pr[i] = bppr_seq(chosen[i]->seq);
				all_pf_segment[seq_id][num_pf_segment[seq_id]] = bp_pr[i];	
				pf_segment_start[seq_id][num_pf_segment[seq_id]] = chosen[i]->start;
				pf_segment_len[seq_id][num_pf_segment[seq_id]] = chosen[i]->len;		
				num_pf_segment[seq_id]++;	
				if (num_pf_segment[seq_id] >= MAXCAND)
					Die("Too  many segments");
			}

		}
		prxy(msa->aseq, msa->nseq, msa->alen, bp_pr, chosen_weight, ret_pxy);
		//Free2DArray((void**)bp_pr, chosen_idx);
	}

	for (i = 0; i < chosen_idx; i++) 
	{
		FreeParsetree(chosen_tr[i]);
		free(chosen_dsq[i]);
	} 
	free(chosen_dsq);
	free(chosen_tr);
	free(chosen_weight);
	free(chosen_sqinfo);
	free(chosen);  
	return (msa);
}


void CM_Scan_hit(CM_t *cm, HMM *hmm, HMM_CM_Map *hmm_cm_map, 
				 char *dsq, int len, int* dmin, int* dmax, int window, 
				 int *ret_nhits, int **ret_hitr, int **ret_hiti, int **ret_hitj, float **ret_hitsc)
{
	int    nhits;			/* number of hits in a seq */
	int   *hitr;			/* initial states for hits */
	int   *hiti;                  /* start positions of hits */
	int   *hitj;                  /* end positions of hits */
	float *hitsc;			/* scores of hits */
	int    i;
	int    x;
	if (do_hmm){      
		/* information on hits found with the CP9 HMM derived from the CM */
		int    hmm_nhits;		/* number of hits in a seq */
		int   *hmm_hitr;		/* initial states for hits */
		int   *hmm_hiti;              /* start positions of hits */
		int   *hmm_hitj;              /* end positions of hits */
		float *hmm_hitsc;		/* scores of hits */   
		HMM_Matrix *hmm_fwd;
		HMM_Matrix *hmm_bck;
		int    hmm_pad=0;
		float  fb_sc;      
		int    tmp_nhits;		 /* number of hits in a seq */
		int   *tmp_hitr;		 /* initial states for hits */
		int   *tmp_hiti;               /* start positions of hits */
		int   *tmp_hitj;               /* end positions of hits */
		float *tmp_hitsc;		 /* scores of hits */
		int alloc_nhits;

		fb_sc = CP9ForwardBackwardScan(dsq, 1, len, window, hmm, &hmm_fwd, &hmm_bck, 
			&hmm_nhits, &hmm_hitr, &hmm_hiti, &hmm_hitj, &hmm_hitsc, hmm_threshold,
			hmm_pad);    
		FreeCPlan9Matrix(hmm_fwd);
		FreeCPlan9Matrix(hmm_bck);    

		alloc_nhits = 10;
		nhits = 0; /* number of CM hits is set to 0, but we'll check each HMM hit with the CM,
				   * and this number may grow */
		hitr  = MallocOrDie(sizeof(int)   * alloc_nhits);
		hitj  = MallocOrDie(sizeof(int)   * alloc_nhits);
		hiti  = MallocOrDie(sizeof(int)   * alloc_nhits);
		hitsc = MallocOrDie(sizeof(float) * alloc_nhits);
		for (i = 0; i < hmm_nhits; i++){
			//printf("HMM hit %d: %d-%d sc %f\n", i, hmm_hiti[i], hmm_hitj[i], hmm_hitsc[i]);
			if (do_hbanded){
				HMM_Band* hmm_band = HMM_CM_Band(cm,  hmm, hmm_cm_map, dsq, hmm_hiti[i], hmm_hitj[i],1);
				CYKBandedScan_jd(cm, dsq, hmm_band->jmin, hmm_band->jmax, hmm_band->hdmin, hmm_band->hdmax, 
					hmm_hiti[i], hmm_hitj[i], window, 
					&tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, cm_threshold);
				FreeCP9Bands(hmm_band);
				free(hmm_band);
			}
			else{
				if (do_banded) 
					CYKBandedScan(cm, dsq, dmin, dmax, hmm_hiti[i], hmm_hitj[i], window, 			
					&tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, cm_threshold);
				else{
					CYKScan(cm, dsq, hmm_hiti[i], hmm_hitj[i], window, 
						&tmp_nhits, &tmp_hitr, &tmp_hiti, &tmp_hitj, &tmp_hitsc, cm_threshold);
				}
			}
			for (x = 0; x < tmp_nhits; x++)
			{
				hitr[nhits] = tmp_hitr[x];
				hiti[nhits] = tmp_hiti[x];
				hitj[nhits] = tmp_hitj[x];
				hitsc[nhits] = tmp_hitsc[x];
				nhits++;
				if (nhits == alloc_nhits) {
					hitr  = ReallocOrDie(hitr,  sizeof(int)   * (alloc_nhits + 10));
					hitj  = ReallocOrDie(hitj,  sizeof(int)   * (alloc_nhits + 10));
					hiti  = ReallocOrDie(hiti,  sizeof(int)   * (alloc_nhits + 10));
					hitsc = ReallocOrDie(hitsc, sizeof(float) * (alloc_nhits + 10));
					alloc_nhits += 10;
				}
			}
			free(tmp_hitr);
			free(tmp_hiti);
			free(tmp_hitj);
			free(tmp_hitsc);
		}    
		free(hmm_hitr);
		free(hmm_hiti);
		free(hmm_hitj);
		free(hmm_hitsc);
	}
	else{
		if (do_banded){
			CYKBandedScan(cm, dsq,  dmin, dmax, 1, len, window, 
				&nhits, &hitr, &hiti, &hitj, &hitsc, cm_threshold);
		}
		else{
			CYKScan(cm, dsq, 1, len, window, 
				&nhits, &hitr, &hiti, &hitj, &hitsc, cm_threshold);    
		}
	}  

	*ret_nhits = nhits;
	*ret_hitr = hitr;
	*ret_hiti = hiti;
	*ret_hitj = hitj;
	*ret_hitsc = hitsc;  
}


int  CM_Search(CM_t *cm, HMM *hmm, HMM_CM_Map *hmm_cm_map, 
			   char *dsq, char* rseq, int seq_id, int max_chosen_cand, 
			   int* dmin, int* dmax, int window, Cand *cand, Range* range)
{
	int    j;
	int    nhits;			/* number of hits in a seq */
	int   *hitr;			/* initial states for hits */
	int   *hiti;                  /* start positions of hits */
	int   *hitj;                  /* end positions of hits */
	float *hitsc;			/* scores of hits */
	double max_score=0;
	int    seq_len = strlen(rseq);    
	int    left_bound  = 1;
	int    right_bound = seq_len;
	int    x,i;
	int sub_len;
	char* sub_dsq;
	float max_sc;

	if (range){
		left_bound = range->start;
		right_bound = range->stop;
		if (right_bound > seq_len) right_bound = seq_len;
	}

	sub_len = right_bound - left_bound + 1;
	sub_dsq = DigitizeSequence(rseq + left_bound - 1, sub_len);
	max_sc = 0;
	if (sub_len < seq_len / 4 ){    
		CM_Scan_hit(cm, hmm, hmm_cm_map, sub_dsq, sub_len, dmin, dmax, window, 
			&nhits, &hitr, &hiti, &hitj, &hitsc);        
		for(i=0; i < nhits; i++){
			if (hitsc[i] > max_sc){
				max_sc = hitsc[i];
			}
		}
		if ( max_sc < 30){      
			free(hitr);
			free(hiti);
			free(hitj);
			free(hitsc);
		}
		free(sub_dsq);
	}  
	if ( max_sc < 30){
		char* dsq =  DigitizeSequence(rseq, seq_len);
		CM_Scan_hit(cm, hmm, hmm_cm_map, dsq, seq_len, dmin, dmax, window, 
			&nhits, &hitr, &hiti, &hitj, &hitsc);
		left_bound = 1;
		right_bound = seq_len;
		free(dsq);
	}

	if (nhits > 0) {
		Cand** sort_cand;
		float best_score;
		Cand* temp;
		if (nhits > MAXCAND) nhits = MAXCAND;
		memset(cand, 0, sizeof(Cand) * nhits);
		for(j=0; j < nhits; j++) {	
			cand[j].score =  hitsc[j];
			cand[j].cand_id = j;
			cand[j].seq_id = seq_id;
			cand[j].start =  hiti[j] + left_bound - 1;
			cand[j].stop =   hitj[j] + left_bound - 1;
			cand[j].len =  cand[j].stop - cand[j].start + 1;
			strncpy(cand[j].seq, rseq + cand[j].start - 1, cand[j].len);
			cand[j].seq[cand[j].len] = '\0';	  	
			cand[j].ss[0]='\0';          //no secondary structure yet
		}
		sort_cand = SortCand(cand, nhits, CompCandByScore);
		best_score = sort_cand[0]->score;      
		temp = MallocOrDie(sizeof(Cand)* nhits);
		for(j=0; j < nhits && j < max_chosen_cand; j++){
			if (sort_cand[j]->score < best_score / 2 && sort_cand[j]->score < 10) break;	
			memcpy(&temp[j], sort_cand[j], sizeof(Cand));	
			temp[j].cand_id = j;
		}
		nhits = j;
		memcpy(cand, temp, sizeof(Cand) * nhits);
		free(temp);      
		free(sort_cand);
	}  
	//Free bad candidate here
	free(hitr);
	free(hiti);
	free(hitj);
	free(hitsc);
	return nhits;
}


void ScanCand(CM_t *cm, HMM* hmm, HMM_CM_Map* hmm_cm_map, 	      
			  int nseq, char** dsq, char** rseqs, 
			  int window, int max_cand, int* dmin, int* dmax, 
			  Cand **cand, int *ncand, Range* range)
{
	int i,j; 
	double  **gamma;		/* cumulative distribution p(len <= n) for state v */
	double   bandp = 0.0001;	/* tail loss probability for banding */
	int total_cand;

	if (do_local) {
		ConfigLocal(cm, 0.5, 0.5);  
		CMLogoddsify(cm);
		if (do_hmm){
			CPlan9SWConfig(hmm, 0.5, 0.5);
			CP9Logoddsify(hmm);        
		}
	}
	CMHackInsertScores(cm);	/* "TEMPORARY" fix for bad priors */  
	//CMHackRootScores(cm);	/* ignore the flanking sequences of the real hit */  

	total_cand = 0;
	for(i=0; i < nseq; i++){
		Range *this_range = NULL;
		if (range && range[i].start != -1){
			this_range = &range[i];
		}
		ncand[i] =0;
		ncand[i]= CM_Search(cm, hmm, hmm_cm_map, 
			dsq[i], rseqs[i], i, max_cand, 
			dmin, dmax, window, cand[i], this_range);
		total_cand += ncand[i];
	}
	if(total_cand > nseq * 5 && cm_threshold < 30){
		cm_threshold += 10;
		hmm_threshold += 10;
	}  
}


void ExtendCand(char* rseq, Cand* cand, int flank_left, int flank_right){

	int seq_len = strlen(rseq);
	if (cand->len + flank_left + flank_right >= MAXLENGTH){
		Die("Candidate too long!");
	}
	cand->start = cand->start > flank_left ? cand->start - flank_left: 1;
	cand->stop  = cand->stop + flank_right <= seq_len? cand->stop + flank_right : seq_len ;	  
	cand->len   = cand->stop - cand->start + 1;
	strncpy(cand->seq, rseq + cand->start - 1, cand->len);
	cand->seq[cand->len] = '\0';	  	
}


int CompareCand(int nseq, Cand** old_cand, int* old_ncand, Cand** cand, int* ncand, int** seq_delta, double** score_delta)
{
	int i,j;
	//Compare Old candidate and new candidate, compute seq_delta and score_delta

	for(i=0; i < nseq; i++){
		for(j=0; j < ncand[i]; j++){	
			int k;
			int overlap=0;
			int min_overlap, max_overlap, min, max;
			int match = -1;
			seq_delta[i][j] = cand[i][j].len;
			for(k=0; k < old_ncand[i]; k++){
				int olap = Overlap(&old_cand[i][k], &cand[i][j], &min_overlap, &max_overlap, &min, &max);	
				if (olap > overlap){
					overlap = olap;
					match = k;
				}
			}
			if (match >= 0){
				Overlap(&old_cand[i][match], &cand[i][j], &min_overlap, &max_overlap, &min, &max);	
				seq_delta[i][j] = abs(min - min_overlap) >  abs(max - max_overlap)? abs(min - min_overlap) :  abs(max - max_overlap);
				score_delta[i][j] = (abs(cand[i][j].score - old_cand[i][match].score) + 1) / (cand[i][j].score + 1);		  
			}
			else{
				score_delta[i][j] =1;
				seq_delta[i][j] = FLT_MAX;
			}
		}
	}
}


/* Parsetrees2Alignment 
* Copy From Infernal package 
*/
MSA *
Parsetrees2AlignmentSS(CM_t *cm, char **dsq, SQINFO *sqinfo, float *wgt, 
					   Parsetree_t **tr, int nseq)
{
	MSA         *msa;          /* multiple sequence alignment */
	CMEmitMap_t *emap;         /* consensus emit map for the CM */
	int          i;            /* counter over traces */
	int          v, nd;        /* state, node indices */
	int          cpos;         /* counter over consensus positions (0)1..clen */
	int         *matuse;       /* TRUE if we need a cpos in mult alignment */
	int         *iluse;        /* # of IL insertions after a cpos for 1 trace */
	int         *eluse;        /* # of EL insertions after a cpos for 1 trace */
	int         *iruse;        /* # of IR insertions after a cpos for 1 trace */
	int         *maxil;        /* max # of IL insertions after a cpos */
	int         *maxel;        /* max # of EL insertions after a cpos */
	int         *maxir;        /* max # of IR insertions after a cpos */
	int	      *matmap;       /* apos corresponding to a cpos */
	int         *ilmap;        /* first apos for an IL following a cpos */
	int         *elmap;        /* first apos for an EL following a cpos */
	int         *irmap;        /* first apos for an IR following a cpos */
	int          alen;	     /* length of msa in columns */
	int          apos;	     /* position in an aligned sequence in MSA */
	int          rpos;	     /* position in an unaligned sequence in dsq */
	int          tpos;         /* position in a parsetree */
	int          el_len;	     /* length of an EL insertion in residues */
	CMConsensus_t *con;        /* consensus information for the CM */
	int          prvnd;	     /* keeps track of previous node for EL */

	emap = CreateEmitMap(cm);

	matuse = malloc(sizeof(int)*(emap->clen+1));   
	iluse  = malloc(sizeof(int)*(emap->clen+1));   
	eluse  = malloc(sizeof(int)*(emap->clen+1));   
	iruse  = malloc(sizeof(int)*(emap->clen+1));   
	maxil  = malloc(sizeof(int)*(emap->clen+1));   
	maxel  = malloc(sizeof(int)*(emap->clen+1));   
	maxir  = malloc(sizeof(int)*(emap->clen+1));   
	matmap = malloc(sizeof(int)*(emap->clen+1));   
	ilmap  = malloc(sizeof(int)*(emap->clen+1));   
	elmap  = malloc(sizeof(int)*(emap->clen+1));   
	irmap  = malloc(sizeof(int)*(emap->clen+1));   

	for (cpos = 0; cpos <= emap->clen; cpos++) 
	{
		matuse[cpos] = 0;
		maxil[cpos] = maxel[cpos] = maxir[cpos] = 0;
		ilmap[cpos] = elmap[cpos] = irmap[cpos] = 0;
	}

	/* Look at all the traces; find maximum length of
	* insert needed at each of the clen+1 possible gap
	* points. (There are three types of insert, IL/EL/IR.)
	* Also find whether we don't need some of the match
	* (consensus) columns.
	*/
	for (i = 0; i < nseq; i++) 
	{
		for (cpos = 0; cpos <= emap->clen; cpos++) 
			iluse[cpos] = eluse[cpos] = iruse[cpos] = 0;

		for (tpos = 0; tpos < tr[i]->n; tpos++)
		{
			v  = tr[i]->state[tpos];
			if (cm->sttype[v] == EL_st) nd = prvnd;
			else                        nd = cm->ndidx[v];

			switch (cm->sttype[v]) {
	  case MP_st: 
		  matuse[emap->lpos[nd]] = 1;
		  matuse[emap->rpos[nd]] = 1;
		  break;
	  case ML_st:
		  matuse[emap->lpos[nd]] = 1;
		  break;
	  case MR_st:
		  matuse[emap->rpos[nd]] = 1;
		  break;
	  case IL_st:
		  iluse[emap->lpos[nd]]++;
		  break;
	  case IR_st:		
		  /* remember, convention on rpos is that IR precedes this
		  * cpos. Make it after the previous cpos, hence the -1. 
		  */
		  iruse[emap->rpos[nd]-1]++;
		  break;
	  case EL_st:
		  el_len = tr[i]->emitr[tpos] - tr[i]->emitl[tpos] + 1;
		  eluse[emap->epos[nd]] = el_len;
		  /* not possible to have >1 EL in same place; could assert this */
		  break;
			}

			prvnd = nd;
		} /* end looking at trace i */

		for (cpos = 0; cpos <= emap->clen; cpos++) 
		{
			if (iluse[cpos] > maxil[cpos]) maxil[cpos] = iluse[cpos];
			if (eluse[cpos] > maxel[cpos]) maxel[cpos] = eluse[cpos];
			if (iruse[cpos] > maxir[cpos]) maxir[cpos] = iruse[cpos];
		}
	} /* end calculating lengths used by all traces */


	/* Now we can calculate the total length of the multiple alignment, alen;
	* and the maps ilmap, elmap, and irmap that turn a cpos into an apos
	* in the multiple alignment: e.g. for an IL that follows consensus position
	* cpos, put it at or after apos = ilmap[cpos] in aseq[][].
	* IR's are filled in backwards (3'->5') and rightflushed.
	*/
	alen = 0;
	for (cpos = 0; cpos <= emap->clen; cpos++)
	{
		if (matuse[cpos]) {
			matmap[cpos] = alen; 
			alen++;
		} else 
			matmap[cpos] = -1;

		ilmap[cpos] = alen; alen += maxil[cpos];
		elmap[cpos] = alen; alen += maxel[cpos];
		alen += maxir[cpos]; irmap[cpos] = alen-1; 
	}

	/* We're getting closer.
	* Now we can allocate for the MSA.
	*/
	msa = MSAAlloc(nseq, alen);
	msa->nseq = nseq;
	msa->alen = alen;
	msa->ss   = (char **)malloc(sizeof(char *) * nseq);    
	for (i = 0; i < nseq; i++)
	{	
		msa->ss[i] = (char *)malloc(sizeof(char) * (alen + 1));      
	}

	for (i = 0; i < nseq; i++){    
		/* Initialize the aseq with all pads '.' (in insert cols) 
		* and deletes '-' (in match cols).
		*/



		for (apos = 0; apos < alen; apos++){	
			msa->aseq[i][apos] = '.';
			msa->ss[i][apos] = '.';
		}

		for (cpos = 0; cpos <= emap->clen; cpos++)
			if (matmap[cpos] != -1) msa->aseq[i][matmap[cpos]] = '-';


		msa->aseq[i][alen] = '\0';
		msa->ss[i][alen] = '\0';

		/* Traverse this guy's trace, and place all his
		* emitted residues.
		*/
		for (cpos = 0; cpos <= emap->clen; cpos++)
			iluse[cpos] = iruse[cpos] = 0;

		for (tpos = 0; tpos < tr[i]->n; tpos++) 
		{
			v  = tr[i]->state[tpos];	 

			if (cm->sttype[v] == EL_st) nd = prvnd;
			else                        nd = cm->ndidx[v];

			switch (cm->sttype[v]) {
	  case MP_st:
		  {
			  char lc,rc;
			  int  l,r;	      
			  cpos = emap->lpos[nd];
			  apos = matmap[cpos];
			  rpos = tr[i]->emitl[tpos];
			  msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
			  msa->ss[i][apos] = '<';	    
			  l = apos;
			  lc= msa->aseq[i][l];

			  cpos = emap->rpos[nd];
			  apos = matmap[cpos];
			  rpos = tr[i]->emitr[tpos];
			  msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
			  msa->ss[i][apos] = '>';	    
			  r  = apos;	      
			  rc= msa->aseq[i][r];

			  if(watson_crick && !IsBasePair(lc,rc)) {
				  msa->ss[i][l]='-';
				  msa->ss[i][r]='-';
			  }	      
		  }	    
		  break;	    
	  case ML_st:
		  cpos = emap->lpos[nd];
		  apos = matmap[cpos];
		  rpos = tr[i]->emitl[tpos];
		  msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
		  break;

	  case MR_st:
		  cpos = emap->rpos[nd];
		  apos = matmap[cpos];
		  rpos = tr[i]->emitr[tpos];
		  msa->aseq[i][apos] = Alphabet[(int) dsq[i][rpos]];
		  break;

	  case IL_st:
		  cpos = emap->lpos[nd];
		  apos = ilmap[cpos] + iluse[cpos];
		  rpos = tr[i]->emitl[tpos];
		  msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
		  iluse[cpos]++;
		  break;

	  case EL_st: 
		  /* we can assert eluse[cpos] always == 0 when we enter,
		  * because we can only have one EL insertion event per 
		  * cpos. If we ever decide to regularize (split) insertions,
		  * though, we'll want to calculate eluse in the rpos loop.
		  */
		  cpos = emap->epos[nd]; 
		  apos = elmap[cpos]; 
		  for (rpos = tr[i]->emitl[tpos]; rpos <= tr[i]->emitr[tpos]; rpos++)
		  {
			  msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
			  apos++;
		  }
		  break;

	  case IR_st: 
		  cpos = emap->rpos[nd]-1;  /* -1 converts to "following this one" */
		  apos = irmap[cpos] - iruse[cpos];  /* writing backwards, 3'->5' */
		  rpos = tr[i]->emitr[tpos];
		  msa->aseq[i][apos] = tolower((int) Alphabet[(int) dsq[i][rpos]]);
		  iruse[cpos]++;
		  break;

	  case D_st:
		  if (cm->stid[v] == MATP_D || cm->stid[v] == MATL_D) 
		  {
			  cpos = emap->lpos[nd];
			  if (matuse[cpos]) msa->aseq[i][matmap[cpos]] = '-';
		  }
		  if (cm->stid[v] == MATP_D || cm->stid[v] == MATR_D) 
		  {
			  cpos = emap->rpos[nd];
			  if (matuse[cpos]) msa->aseq[i][matmap[cpos]] = '-';
		  }
		  break;

			} /* end of the switch statement */
			prvnd = nd;
		} /* end traversal over trace i. */

		/* Here is where we could put some insert-regularization code
		* a la HMMER: reach into each insert, find a random split point,
		* and shove part of it flush-right. But, for now, don't bother.
		*/

	} /* end loop over all parsetrees */


	/* Gee, wasn't that easy?
	* Add the rest of the ("optional") information to the MSA.
	*/
	con = CreateCMConsensus(cm, 3.0, 1.0);

	/* "author" info */
	msa->au   = malloc(sizeof(char) * (strlen(RELEASE)+10));
	sprintf(msa->au, "CMfinder %s", RELEASE);
	if (wgt != NULL) msa->flags |= MSA_SET_WGT;  
	for (i = 0; i < nseq; i++)
	{
		msa->sqname[i] = sre_strdup(sqinfo[i].name, -1);

		msa->sqlen[i]  = sqinfo[i].len;
		if (sqinfo[i].flags & SQINFO_ACC)
			MSASetSeqAccession(msa, i, sqinfo[i].acc);
		if (sqinfo[i].flags & SQINFO_DESC)
			MSASetSeqDescription(msa, i, sqinfo[i].desc);
		if (wgt == NULL) msa->wgt[i] = 1.0;
		else    {
			msa->wgt[i] = wgt[i];
		}
	}

	/* Construct the secondary structure consensus line, msa->ss_cons:
	*       IL, IR are annotated as .
	*       EL is annotated as ~
	*       and match columns use the structure code.
	* Also the primary sequence consensus/reference coordinate system line,
	* msa->rf.
	*/
	msa->ss_cons = malloc(sizeof(char) * (alen+1));
	msa->rf = malloc(sizeof(char) * (alen+1));
	for (cpos = 0; cpos <= emap->clen; cpos++) 
	{
		if (matuse[cpos]) 
		{ /* CMConsensus is off-by-one right now, 0..clen-1 relative to cpos's 1..clen */

			/* bug i1, xref STL7 p.12. Before annotating something as a base pair,
			* make sure the paired column is also present.
			*/
			if (con->ct[cpos-1] != -1 && matuse[con->ct[cpos-1]+1] == 0) {
				msa->ss_cons[matmap[cpos]] = '.';
				msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
			} else {
				msa->ss_cons[matmap[cpos]] = con->cstr[cpos-1];	
				msa->rf[matmap[cpos]]      = con->cseq[cpos-1];
			}
		}
		if (maxil[cpos] > 0) 
			for (apos = ilmap[cpos]; apos < ilmap[cpos] + maxil[cpos]; apos++)
			{
				msa->ss_cons[apos] = '.';
				msa->rf[apos] = '.';
			}
			if (maxel[cpos] > 0)
				for (apos = elmap[cpos]; apos < elmap[cpos] + maxel[cpos]; apos++)
				{
					msa->ss_cons[apos] = '~';
					msa->rf[apos] = '~';
				}
				if (maxir[cpos] > 0)	/* remember to write backwards */
					for (apos = irmap[cpos]; apos > irmap[cpos] - maxir[cpos]; apos--)
					{
						msa->ss_cons[apos] = '.';
						msa->rf[apos] = '.';
					}
	}
	msa->ss_cons[alen] = '\0';
	msa->rf[alen] = '\0';

	FreeCMConsensus(con);
	FreeEmitMap(emap);
	free(matuse);
	free(iluse);
	free(eluse);
	free(iruse);
	free(maxil);
	free(maxel);
	free(maxir);
	free(matmap);
	free(ilmap);
	free(elmap);
	free(irmap);
	return msa;
}
