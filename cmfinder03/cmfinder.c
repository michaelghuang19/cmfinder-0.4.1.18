/* cmfinder.c
* 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

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

static char usage[]  = "\
					   Usage: cmfinder [-options] <seqfile in> <cmfile output> \n\
					   where options are:\n\
					   -c <candidate file>: the candidate file \n\
					   -a <align file>    : the initial motif alignment \n\
					   -i <cm file>       : the initial covariance model\n\
					   -o <align file>    : the output motif structural alignment in stockholm format \n\
					   -v verbose         : print intermediate alignments \n	\
					   -h                 : print short help and version info\n\
					   ";

static char experts[] = "\
						Expert, in development, or infrequently used options are:\n		\
						--g <gap threshold> : the gap threshold to determine the conserved column\n\
						--hmm               : apply HMM filter \n\
						--cmzasha           : apply cmzasha filter \n				\
						--update            : Update instead of scanning for new candidates at each iteration \n\
						--informat <s>: specify that input alignment is in format <s>\n\
						--fragmentary : account for fragmentary input sequences\n\
						";


static struct opt_s OPTIONS[] = {
	{ "-c", TRUE, sqdARG_STRING},
	{ "-a", TRUE, sqdARG_STRING},
	{ "-i", TRUE, sqdARG_STRING},
	{ "-o", TRUE, sqdARG_STRING},
	{ "-r", TRUE, sqdARG_STRING},
	{ "-w", TRUE, sqdARG_STRING},
	{ "-l", TRUE, sqdARG_INT},
	{ "-v", TRUE, sqdARG_NONE },
	{ "-h", TRUE, sqdARG_NONE },
	{ "--t",        FALSE, sqdARG_STRING},
	{ "--g",        FALSE, sqdARG_FLOAT},
	{ "--hmm",      FALSE, sqdARG_NONE},
	{ "--hbanded",  FALSE, sqdARG_NONE},
	{ "--local",    FALSE, sqdARG_NONE},
	{ "--scan",     FALSE, sqdARG_STRING},
	{ "--update",   FALSE, sqdARG_NONE },
	{ "--informat", FALSE, sqdARG_STRING},
	{ "--fragmentary", FALSE, sqdARG_NONE}
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static char banner[] = "cmfinder: learning motif covariance model for unaligned sequences\n";

extern void train_EVD(struct cm_s *cm);
extern void CM_save(CM_t *cm, char* cmfile);
extern int save_hit(int left, int right, double score);
extern int init_hits(char* s);
extern int return_hits(Cand** ret_cand, int * ret_ncand);
extern CM_t * Automodelmaker(MSA *msa, char** dsq,  int ** pxy,double gapthresh, int gapcost, int trim, int** constraint);
extern double* bppr_seq(char* seq);
extern CM_t * M_step (MSA* msa, char **dsq, Prior_t *pri, float *null,int** pxy, float gapthreshold, int gapcost, HMM** ret_hmm, HMM_CM_Map** ret_map);
extern void ScanCand(CM_t *cm, HMM* hmm, HMM_CM_Map* hmm_cm_map, 	      
					 int nseq, char** dsq, char** rseqs, 
					 int window, int max_cand, int* dmin, int* dmax, 
					 Cand **cand, int *ncand, Range* range);
extern MSA* E_step(struct cm_s *cm, 
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
				   double *bestcmscore
				   );



int    do_local = 0;
int    do_small = 1;
int    do_binary = 0;
int    do_banded = 1;
int    do_hmm = 0;
int    do_hbanded = 0;
int    do_update = 0;
int    do_evd = 0; 
int    do_zoop = 0;
int    do_tcm = 1;
int    do_inside = 0;
int    use_fragmentary=0;

int    proximity=30;
int    DB_length = 100000;       
double mu;
double lambda;
char  *cmfile=NULL;
char  *final_file = NULL;   /* files contains the final alignment of motifs */

extern int    do_evd;
extern int    do_zoop;
extern int    do_tcm;





int
main(int argc, char **argv)
{
	int        format;              /* alifile format                            */
	char      *seqfile=NULL;        /* training sequence file                    */
	char      *alifile=NULL;        /* file contain the initial alignment of selected cand */
	char      *candfile = NULL;     /* file contain candidate motifs             */
	CMFILE    *in_cmfp;	          /* open CM file */

	char 	    **rseqs;	          /* training sequences                        */
	char      **dsq;                /* Digitized training sequences              */
	SQINFO    *sqinfo;		  /* array of sqinfo structures for rseqs      */
	int  	    nseq;		  /* number of seqs */                           

	Cand      **cand;        /* all the candidates for all the sequences */
	int       *ncand;	   /* the number of cands of each sequence */     

	CM_t     *cm      = NULL;       /* current model                        */

	double    bestcmscore;     /* The score of the best cm */  
	double    totscore;	     /* summed scores over training seqs          */
	double    oldscore;	     /* previous totscore for old model           */
	double    delta;	     /* fractional change in scores for iteration */
	int       iteration;	     /* iteration number we're on                 */
	double    gapthreshold=0.6;/* gap threshold to determine the conserved columns */  
	double    gapcost     =100;/* gap penalty in M step */  
	int       max_iterations;
	double    threshold;		/* fractional tolerance, test for convergence   */   

	char      *in_cmfile = NULL;		/* file containing input model       */
	Prior_t   *pri       = NULL;           /* mixture Dirichlet prior structure */

	double    *seq_weight= NULL;    
	MSA       *msa       = NULL;  
	char     **msa_dsq   = NULL;  
	float      null[4];  
	double     nt_count = 0;  
	int      **pxy = NULL;  /* partition function log ratio */  

	char  *optname;                /* name of option found by Getopt()        */
	char  *optarg;                 /* argument found by Getopt()              */
	int    optind;                 /* index in argv[]                         */	
	int    verbose=0;
	int    i, j;
	int    window        = 200;    
	int    safe_windowlen        = 200;    /* for Banded window computational*/
	int   *dmin=NULL;
	int   *dmax=NULL;
	double **gamma=NULL;
	float  beta              = (float)0.0001;

	Range     *scan_range = NULL;
	char      *range_file   = NULL;
	HMM       *hmm=NULL;  
	HMM_CM_Map *hmm_cm_map = NULL;
	char      *scan_file = NULL;
	FILE*     scan_fout=NULL;

	/*Parse command line */	

#ifdef MEMDEBUG
	unsigned long histid1, histid2, orig_size, current_size;
#endif

	/*********************************************** 
	* Parse commnd line
	***********************************************/

	threshold            = 0.02;		/* default: 1% */
	max_iterations       = 100;
	format               = MSAFILE_STOCKHOLM;  

	while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))  {
			if      (strcmp(optname, "-c") == 0)        candfile  = optarg;
			else if (strcmp(optname, "-a") == 0)        alifile   = optarg;
			else if (strcmp(optname, "-i") == 0)        in_cmfile = optarg;
			else if (strcmp(optname, "-o") == 0)        final_file= optarg;
			else if (strcmp(optname, "-r") == 0)        range_file= optarg;
			else if (strcmp(optname, "-l") == 0)        DB_length= atoi(optarg);
			else if (strcmp(optname, "-v") == 0)        verbose   = 1;
			else if (strcmp(optname, "--hmm")==0)       do_hmm    = 1;
			else if (strcmp(optname, "--hbanded")==0)   do_hbanded = 1;
			else if (strcmp(optname, "--local")==0)     do_local = 1;
			else if (strcmp(optname, "--scan")==0)      scan_file = optarg;
			else if (strcmp(optname, "--update")==0)    do_update = 1;
			else if (strcmp(optname, "--g") == 0)       gapthreshold   = atof(optarg);
			else if (strcmp(optname, "--fragmentary") == 0) use_fragmentary=1;
			else if (strcmp(optname, "-w") == 0)        {      
				if ( strcmp(optarg, "EVD") == 0) do_evd = 1;
				else if ( strcmp(optarg, "ZOOP") == 0) do_zoop = 1;
				else if ( strcmp(optarg, "TCM") == 0) do_tcm = 1;
			}
			else if (strcmp(optname, "--informat") == 0){
				format = String2SeqfileFormat(optarg);
				if (format == MSAFILE_UNKNOWN) 
					Die("unrecognized sequence file format \"%s\"", optarg);
				if (! IsAlignmentFormat(format))
					Die("%s is an unaligned format, can't read as an alignment", optarg);
			}
			else if (strcmp(optname, "-h") == 0) {
				puts(banner);
				puts(usage);
				puts(experts);
				exit(EXIT_SUCCESS);
			}    
	}  
	if (argc - optind < 2) Die("%s\n", usage);

	seqfile = argv[optind++];  
	cmfile  = argv[optind++];  

	/*********************************************** 
	* Get sequence data
	***********************************************/
	/* read the training seqs from file */
	if (! ReadMultipleRseqs(seqfile, SQFILE_FASTA, &rseqs, &sqinfo, &nseq))
		Die("Failed to read any sequences from file %s", seqfile);

	/* Preprocess */
	for (i = 0; i < nseq; i++){
		PrepareSequence(rseqs[i],use_fragmentary); /* randomly sample from degenerate nucs */
	}   
	dsq = (char **) malloc(sizeof(char *) * nseq);
	for(i=0; i < Alphabet_size; i++) null[i]=0;  
	for (i = 0; i < nseq; i++) {    
		dsq[i] = DigitizeSequence(rseqs[i], sqinfo[i].len);
		for(j=1; j <= sqinfo[i].len; j++) {      
			if (dsq[i][j] < 0 || dsq[i][j] >= Alphabet_size) 
				Die("Seq %d position %d: Invalid letter", i, j, dsq[i][j]);      
		}    
	}  
	CMDefaultNullModel(null);

	if (range_file){
		scan_range = ReadCluster(range_file, sqinfo, nseq) ;
	}

	/* Get prior */
	pri = Prior_Default();
	/*Construct Initial CM */
	if (in_cmfile){
		/* Read CM File*/
		if ((in_cmfp = CMFileOpen(in_cmfile, NULL)) == NULL)
			Die("Failed to open covariance model save file %s\n%s\n", in_cmfile, usage);
		if (! CMFileRead(in_cmfp, &cm))
			Die("Failed to read a CM from %s -- file corrupt?\n", in_cmfile);
		if (cm == NULL) 
			Die("%s empty?\n", in_cmfile);
		CMFileClose(in_cmfp);    
		/* Construct CM */
		CMLogoddsify(cm);    
		CM_save(cm, cmfile);
		if (do_hmm){
			if(!build_cp9_hmm(cm, &hmm, &hmm_cm_map, FALSE, (float)0.0001, 0))
				Die("Couldn't build a CP9 HMM from the CM\n");
		}
	}  
	else if (alifile){
		/* Read initial alignment file*/
		MSAFILE     *afp = NULL;        /* file handle of initial alignment          */    
		if ((afp = MSAFileOpen(alifile, format, NULL)) == NULL)
			Die("Alignment file %s could not be opened for reading", alifile);
		if ((msa = MSAFileRead(afp)) != NULL)
		{
			for (i = 0; i < msa->nseq; i++){
				PrepareSequence(msa->aseq[i],use_fragmentary);
			}	
			/* Estimate CM */
			MSAFileClose(afp);	
			msa_dsq = DigitizeAlignment(msa->aseq, msa->nseq, msa->alen);
			cm = M_step(msa, msa_dsq, pri, null, NULL, (float)gapthreshold, 50, &hmm, &hmm_cm_map);	
			if (cm== NULL){
				Die("No secondary structure detected");
			}
			CM_save(cm, cmfile);

			Free2DArray((void **)msa_dsq, msa->nseq);      
			MSAFree(msa);
		}
		else Die("Fail to read alignment file %s", alifile);
	}
	else {
		Die("No initial Covariance Model or initial alignment is available");
	}

	if (do_evd) train_EVD(cm);

	cand = (Cand **) malloc(sizeof(Cand*) * nseq);
	ncand = (int *) malloc(sizeof(int) * nseq);    
	for(i =0 ; i < nseq; i++)
		cand[i] = (Cand *)malloc(sizeof(Cand) *MAXCAND);        

	for(i=0; i< nseq; i++)  ncand[i] = 0;  

	bestcmscore = NEGINFINITY;  
	oldscore = NEGINFINITY;  
	iteration = 0;

	/* Now Scan the whole sequence looking to candidates. */    
	if (do_banded){
		safe_windowlen = window * 2;
		while(!(BandCalculationEngine(cm, safe_windowlen, beta, 0, &dmin, &dmax, &gamma, do_local)))
		{
			FreeBandDensities(cm, gamma);
			free(dmin);
			free(dmax);
			safe_windowlen *= 2;
		}
		FreeBandDensities(cm, gamma);
		window = dmax[0];      
	}
	ScanCand(cm, hmm, hmm_cm_map, nseq, dsq, rseqs, window, 10, dmin, dmax, cand, ncand, scan_range);        

	while( iteration < max_iterations ) {    
		iteration++;
		if (verbose)
			printf("Iteration %4d  : model of %d nodes\n ", iteration, cm->nodes);    

		totscore = 0.0;           
		msa = E_step(cm, hmm, hmm_cm_map, dmin, dmax, 
			nseq, sqinfo, seq_weight, 
			ncand, cand, 10, scan_range, &totscore, &pxy,NULL, &bestcmscore);              


		FreeCM(cm);
		if (hmm) FreeCPlan9(hmm);
		if (hmm_cm_map) FreeCP9Map(hmm_cm_map);

		if (totscore < -10) {
			Die("Bad alignment %f!", totscore);	
		} 

		delta = (totscore - oldscore) / fabs(totscore);   

		if(verbose){
			printf("score %.3f, delta %.3f\n", totscore / (double) nseq, delta);    
			/*WriteStockholm(stdout, msa);    */
		}   

		/* If we've converged, stop.
		* Else, make a new model from the alignment. 
		*/	
		if (delta < threshold )
		{		
			/* we've converged. Free traces and break out of iteration loop.*/
			break;	  
		}
		if (do_evd && delta > 0.3) train_EVD(cm);
		oldscore = totscore;    	
		msa_dsq = DigitizeAlignment(msa->aseq, msa->nseq, msa->alen);

		cm = M_step(msa, msa_dsq, pri, null, pxy, (float)gapthreshold, gapcost, &hmm, &hmm_cm_map);	
		if (cm == NULL){
			Die("No secondary structure detected");
		}

		if (do_banded){
			if (dmax) free(dmax);
			if (dmin) free(dmin);
			safe_windowlen = window * 2;
			while(!(BandCalculationEngine(cm, safe_windowlen * 2, beta, 0, &dmin, &dmax, &gamma, do_local)))
			{
				FreeBandDensities(cm, gamma);
				free(dmin);
				free(dmax);
				safe_windowlen *= 2;
			}
			FreeBandDensities(cm, gamma);
			window = dmax[0];
		}
		ScanCand(cm, hmm, hmm_cm_map, nseq, dsq, rseqs, window, 10, dmin, dmax, cand, ncand, scan_range);                    

		Free2DArray((void **)msa_dsq, msa->nseq);      
		MSAFree(msa); 
	}            

	for(i=0; i < nseq;i++)  {
		free(cand[i]);  
	}
	free(cand);
	free(ncand);

	if (verbose)
		printf("Final score %.3f\n", bestcmscore / (double) nseq);  

	if (scan_range){
		free(scan_range);
	}

	if (dmax) free(dmax);
	if (dmin) free(dmin);    
	Prior_Destroy(pri);
	MSAFree(msa);
	for (i = 0; i < nseq; i++)
		FreeSequence(rseqs[i], &(sqinfo[i]));
	Free2DArray((void**)dsq, nseq);
	free(sqinfo);

	return 0;

}
