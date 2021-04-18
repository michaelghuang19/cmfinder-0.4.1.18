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
#include "esl_random.h"
#include "esl_randomseq.h"

#ifdef HMMER_THREADS
#include "esl_threads.h"
#endif

#include "hmmer.h"

#include "infernal.h"

#include "cand.h"
#include "cmfinder.h"

#define CMFINDER_PACKAGE_VERSION "0.4.1.18"

#ifndef _POSIX_VERSION
// supply a body for this function so that it builds under MinGW
char *ctime_r(const time_t *timep, char *buf)
{
	strcpy(buf,"Wed Jun 30 21:49:08 2012\n"); // we need to give it the proper form, otherwise the in-memory CM files can't be parsed properly
	return buf;
}
#endif

#define EFFOPTS "--eent,--enone"                        /* Exclusive options for effective sequence number calculation */
#define DEGENOPTS "--degen-rand,--degen-keep"           /* options for dealing with degenerate nucleotides */
static ESL_OPTIONS cmfinder_options[] = {
	/* name           type            default  env   range      toggles  reqs   incomp     help                                     docgroup*/
	{ "-h",           eslARG_NONE,    FALSE,   NULL, NULL,      NULL,    NULL,  NULL,      "show brief help on version and usage",  1 },
	{ "-a",           eslARG_INFILE,  NULL,    NULL, NULL,      NULL,    NULL,  NULL,      "input alignment file (.sto)",           1},
	{ "-o",           eslARG_OUTFILE, NULL,    NULL, NULL,      NULL,    NULL,  NULL,      "output alignment file (.sto)",          1},
    { "--version",    eslARG_NONE,    FALSE,   NULL, NULL,      NULL,    NULL,  NULL,      "show version", 1},
	
	{ "--degen-rand", eslARG_NONE,    FALSE,   NULL, NULL,      DEGENOPTS,NULL, NULL,      "randomize degenerate nucs like CMfinder 0.3", 2},
	{ "--degen-keep", eslARG_NONE,    "default",NULL,NULL,      DEGENOPTS,NULL, NULL,      "keep degenerate nucs and marginalize", 2},
	{ "--fragmentary",eslARG_NONE,    FALSE,   NULL, NULL,      NULL,    NULL,  NULL,      "allow truncated hits (independent of --degen-X, unlike pscore)", 2},
	{ "--non-frag-avg-bppr",eslARG_NONE,NULL,  NULL, NULL,      NULL,    NULL,  NULL,      "ignore --fragmentary for calculating average base pair probs", 2},
	
	{ "--wgsc",        eslARG_NONE,    NULL,    NULL,NULL,      NULL,    NULL,  NULL,      "use GSC alg to weight sequences for redundancy",2},
	{ "--wpb",        eslARG_NONE,    NULL,    NULL, NULL,      NULL,    NULL,  NULL,      "use PB alg to weight sequences for redundancy",2},
	
	{ "--ints-like-03",eslARG_NONE,   FALSE,   NULL, NULL,      NULL,    NULL,  NULL,      "use ints for mutual info and partition func, like CMfinder 0.3 did", 2},
    { "--avoid-part-func",eslARG_NONE,FALSE,   NULL, NULL,      NULL,    NULL,  NULL,      "don't use partition function to predict structure, only use mutual-info"},
    { "--column-only-base-pair-probs",eslARG_NONE,FALSE,NULL,NULL,NULL,  NULL,  NULL,      "don't use partition function to predict structure, but use per-column-only method just based on freqs"},
	{ "--min-seq-weight",eslARG_REAL, "0.01",  NULL, NULL,      NULL,    NULL,  NULL,      "eliminate seqs from MSA whose TCM weight is below this value", 2},
	{ "--no-elim-iden-seq",eslARG_NONE,FALSE,  NULL, NULL,      NULL,    NULL,  NULL,      "don't eliminate identical sequences as candidate members",2},
	{ "--no-elim-iden-subseq",eslARG_NONE,FALSE,NULL,NULL,      NULL,    NULL,  NULL,      "don't eliminate identical sub-sequences of other sequences as candidate members",2},
	{ "--allow-untested-code",eslARG_NONE,FALSE,NULL,NULL,      NULL,    NULL,  NULL,      "run code that was never exercized in tests; don't abort to allow testing",2},
	{ "--min-cand-score-in-final",eslARG_REAL,"0",NULL,NULL,    NULL,    NULL,  NULL,      "min cmsearch score to put a seq into the saved MSA.",2},
    { "--min-base-pairs",eslARG_INT,  "3",     NULL, NULL,      NULL,    NULL,  NULL,      "if fewer than this number of base pairs are predicted in the M-step, then abort"},
	{ "--bg-score-size",eslARG_INT,   "0",     NULL, NULL,      NULL,    NULL,  NULL,      "create this many randomized seqs for each input seq to get background score, below which candidates are rejected",2},
	{ "--bg-score-evalue",eslARG_REAL, "-1",   NULL, NULL,      NULL,    NULL,  NULL,      "try to get an EVD from --bg-score-size, and apply this E-value",2},
	{ "--bg-score-scan-thresh",eslARG_REAL,"0",NULL, NULL,      NULL,    NULL,  NULL,      "bitscore threshold (-T in cmsearch) to use for scanning --bg-score-size data.",2},
	{ "--bg-score-non-frag",eslARG_NONE,FALSE, NULL, NULL,      NULL,    NULL,  NULL,      "prevent --bg-score-size scans from using fragmentary modes -- force --nontrunc",2},
	{ "--filter-non-frag",eslARG_NONE, FALSE,  NULL, NULL,      NULL,    NULL,  NULL,      "first run cmsearch with --notrunc, then run normal cmsearch only on the hits from --notrunc",2},
	{ "--filter-non-frag-pad",eslARG_INT,"20", NULL, NULL,      NULL,    NULL,  NULL,      "with --filter-non-frag, add this many nucs on 5' and 3' side of the non-frag hits",2},
	{ "--max-degen-per-hit",eslARG_INT,NULL,   NULL, NULL,      NULL,    NULL,  NULL,      "eliminate hits with this many degen nucs or more",2},
	{ "--max-degen-flanking-nucs",eslARG_INT,NULL,NULL,NULL,    NULL,    NULL,  NULL,      "consider this many nucs beyond the 5' and 3' ends of each hit in counting degen nucs for --max-degen-per-hit",2},
	{ "--bad-distal-pairs-to-loop",eslARG_NONE,FALSE,NULL,NULL, NULL,    NULL,  NULL,      "shift non-canon pairs in distal part of stems to the terminal loop",2},
	{ "--bad-distal-pairs-to-loop-only",eslARG_NONE,FALSE,NULL,NULL,NULL,NULL,  NULL,      "just run the input msa (-a flag) thru --bad-distal-pairs-to-loop-test and save to output msa (-o flag)",2},
	{ "--min-term-loop-nucs",eslARG_INT,NULL,  NULL, NULL,      NULL,    NULL,  NULL,      "only with --bad-distal-pairs-to-loop, move even good base pairs into loop if there are fewer than this many nucs in term loop.  But leave it alone if bp is truncated (with --fragmentary)",2},
    { "--open-pair-cost",eslARG_REAL,NULL,NULL,NULL,NULL,NULL,NULL,"set openPairCost in M-step"},
    { "--open-pair-cost-initial",eslARG_REAL,NULL,NULL,NULL,NULL,NULL,NULL,"set openPairCost in first M-step"},

	{ "--seed",         eslARG_INT,   "0",     NULL,"n>=0",     NULL,    NULL,  NULL,      "set random number generator's seed to <n>",2},
	
	{ "--evalue",     eslARG_REAL,    NULL,    NULL, NULL,      NULL,    NULL,  NULL,      "use this E-value in ScanCand, in addition to a threshold (note: implies running internal cmcalibrate, which will be very slow)", 2},
	
	{ "--create-file-on-success",eslARG_OUTFILE,NULL,NULL,NULL, NULL,    NULL,  NULL,      "create this file, empty, upon successful completion of the program, for convenience elsewhere",2},
	{ "--save-after-first-M-step",eslARG_NONE,NULL,NULL,NULL,   NULL,    NULL,  NULL,      "for debugging.  program exits after saving the file",2},
	{ "--save-in-progress",eslARG_NONE,NULL,   NULL, NULL,      NULL,    NULL,  NULL,      "for debugging, save MSA's as we processed",2},
	{ "--timer-append",eslARG_OUTFILE,NULL,    NULL, NULL,      NULL,    NULL,  NULL,      "append timing stats to tab-delimited file",2},
	
	/* parameters that are essentially just for cmcalibrate */
	{ "--tailn",      eslARG_INT,     NULL,    NULL, NULL,      NULL,    NULL,  NULL,      "pass --gtailn or --ltailn as appropriate to cmcalibrate (default: accept cmcalibrate's default)", 3},
	
	/* parameters that are essentially just for cmsearch */
	{ "--local",      eslARG_NONE,    FALSE,   NULL, NULL,      NULL,    NULL,  NULL,      "local mode, i.e. don't pass -g to internal cmsearch",3},
	{ "--noF4F5",     eslARG_NONE,    FALSE,   NULL, NULL,      NULL,    NULL,  NULL,      "set --noF4 and --noF5 (env def) to avoid glocal hmm", 3},
	{ "--max",        eslARG_NONE,    FALSE,   NULL, NULL,      NULL,    NULL,  "--noF4F5","pass --max to cmsearch (and skip calibrations in cmbuild)", 3},
	{ "--amaa",       eslARG_NONE,    FALSE,   NULL, NULL,      NULL,    NULL,  NULL,      "use maximal-alignment accuracy in cmsearch, i.e. don't pass --acyk", 3},
#ifdef HMMER_THREADS
	{ "--cpu",        eslARG_INT,     "0",     NULL,"n>=-1",     NULL,    NULL,  NULL,      "pass to internal cmsearch/cmcalibrate.  --cpu -1 means use all CPUs (default is --cpu 0, which is single-threaded)",       3 },
#endif
	
	/* parameters that just go to cmbuild */
	{ "--p56",       eslARG_NONE,    NULL,    NULL, NULL,         NULL,      NULL,    "--prior", "use the default prior from Infernal v0.56 through v1.0.2",     4 },
	{ "--prior",     eslARG_INFILE,  NULL,    NULL, NULL,         NULL,      NULL,  "--rsearch", "read priors from file <f>",                                      4 }, /* not implemented yet */
	{ "--eent",    eslARG_NONE, "default",    NULL,   NULL, EFFOPTS,     NULL,   NULL, "adjust eff seq # to achieve relative entropy target",           4 },
	{ "--enone",   eslARG_NONE,     FALSE,    NULL,   NULL, EFFOPTS,     NULL,   NULL, "no effective seq # weighting: just use nseq",                   4 },
	{ "--EmN",     eslARG_INT,    "200", NULL, "n>0",   NULL,  NULL, NULL,        "number of sampled seqs to use for p7 local MSV calibration",    4 }, /* not implemented yet */
	{ "--EvN",     eslARG_INT,    "200", NULL, "n>0",   NULL,  NULL, NULL,        "number of sampled seqs to use for p7 local Vit calibration",    4 }, /* not implemented yet */
	{ "--ElfN",    eslARG_INT,    "200", NULL, "n>0",   NULL,  NULL, NULL,        "number of sampled seqs to use for p7 local Fwd calibration",    4 }, /* not implemented yet */
	{ "--EgfN",    eslARG_INT,    "200", NULL, "n>0",   NULL,  NULL, NULL,        "number of sampled seqs to use for p7 glocal Fwd calibration",   4 }, /* not implemented yet */
	
	/* flags for summarize functionality */
	{ "--summarize", eslARG_NONE,    NULL,    NULL, NULL,         NULL,      NULL,    NULL,  "perform functionality like 'summarize' program.  commandline has the .sto file", 5},
	{ "--summarize-gsc", eslARG_NONE,NULL,    NULL, NULL,         NULL,      NULL,    NULL,  "use GSC alg weighting for --summarize", 5},
	{ "--summarize-save-msa",eslARG_OUTFILE,NULL,NULL,NULL,       NULL,      NULL,    NULL,  "save MSA used by --summarize, for debugging modifications on loading",5},
	{ "--summarize-no-de-tag", eslARG_NONE,NULL,    NULL, NULL,         NULL,      NULL,    NULL,  "skip calculations that require #=GR ... DE tags", 5},
	

  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
/*
 * EXIT CODE: 0=success, 1=failure, 2=couldn't produce an acceptable alignment based on algorithmic criteria (i.e. everything's working properly, the input just doesn't have a good motif)
 */
 
static char usage[]="[options] <input-sto-file>\nOR --summarize [options] <input-sto-file>";

ESL_ALPHABET *abc=NULL;

/* this function is analogous to global.c/PrepareSequence in older CMfinder
 */
void RandomizeDegen_digital (ESL_DSQ *dsq,int len)
{
	int64_t pos;
	for (pos=1; pos<=len; pos++) {
		ESL_DSQ ch=dsq[pos];
		if (esl_abc_XIsDegenerate(abc,ch)) {
			esl_fatal("not allowing degenerate nucs for now, while I get the basics working");
		}
	}
}
void RandomizeDegen (char *seq)
{
	char *sym;
	for (sym = seq; *sym != '\0'; sym++) {
		*sym = toupper((int)*sym);
		if (*sym == 'T') *sym= 'U';
		if (esl_abc_CIsGap(abc,*sym)) {
			continue;
		}
		
		if (esl_abc_CIsDegenerate(abc,*sym)) {
			esl_fatal("not allowing degenerate nucs for now while I get the basics working");
		}
	}
}
void RandomizeDegenInAllSeqs (ESL_SQCACHE *inputSeqCache)
{
	int64_t seqnum;
	ESL_SQ *sq;
	for (seqnum=0; seqnum<inputSeqCache->seq_count; seqnum++) {
		sq=&(inputSeqCache->sq_list[seqnum]);
		assert(esl_sq_IsDigital(sq));
		RandomizeDegen_digital(sq->dsq,sq->L);
	}
}

void RandomizeDegenInMsaDigital(ESL_MSA *msa)
{
	int seqnum;
	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}
	for (seqnum=0; seqnum<msa->nseq; seqnum++) {
		RandomizeDegen_digital(msa->ax[seqnum],msa->alen);
	}
}
void RandomizeDegenInMsaText(ESL_MSA *msa)
{
	int seqnum;
	if ((msa->flags & eslMSA_DIGITAL)!=0) {
		esl_fatal("msa should not be digital");
	}
	for (seqnum=0; seqnum<msa->nseq; seqnum++) {
		RandomizeDegen(msa->aseq[seqnum]);
	}
}

void MultiplyEMWeightsByWeightingStrategy (ESL_MSA *msa,WeightingStrategy ws,ESL_STOPWATCH *timer_weighting,int use_fragmentary)
{
	int *seq_fragment_start_col, *seq_fragment_end_col;
	ESL_STOPWATCH *timer_this=esl_stopwatch_Create();
	esl_stopwatch_Start(timer_this);
	
	seq_fragment_start_col=(int *)MallocOrDie(msa->nseq * sizeof(int));
	seq_fragment_end_col=(int *)MallocOrDie(msa->nseq * sizeof(int));
	CalcFragmentStartEnds (msa,use_fragmentary,seq_fragment_start_col,seq_fragment_end_col);
	
	switch (ws) {
	case WS_EM_weights_only:
		/* assume they're okay.  they should be set from eslx_msafile_read, and if the alignment doesn't have #=GS ... WT lines, then it'll be set to 1.0 */
		break;
	case WS_uniform_weights:
		esl_vec_DSet(msa->wgt, msa->nseq, 1.);
		break;
	case WS_GSC:
		esl_msaweight_GSC_fragmentary(msa,use_fragmentary,seq_fragment_start_col,seq_fragment_end_col); /* clobbers whatever weights were there before */
		break;
	case WS_GSC_and_EM:
	case WS_PB_and_EM:
		{
			int i;
			double *wgt=(double *)MallocOrDie(sizeof(double)*msa->nseq);
			memcpy(wgt,msa->wgt,sizeof(double)*msa->nseq);
			if (ws==WS_GSC_and_EM) {
				esl_msaweight_GSC_fragmentary(msa,use_fragmentary,seq_fragment_start_col,seq_fragment_end_col);
			}
			else {
				if (ws==WS_PB_and_EM) {
					esl_msaweight_PB(msa);
				}
				else {
					assert(0);
				}
			}
			for (i=0; i<msa->nseq; i++) {
				msa->wgt[i] *= wgt[i];
			}
			free(wgt);
		}
		break;
	default:
		assert(0);
	}

	free(seq_fragment_start_col);
	free(seq_fragment_end_col);

	esl_stopwatch_Stop(timer_this);
	if (timer_weighting!=NULL) {
		esl_stopwatch_Include(timer_weighting,timer_this);
	}
	esl_stopwatch_Destroy(timer_this);
}

ESL_MSA *LoadMSAFromFile(char *alifile,WeightingStrategy ws,ESL_STOPWATCH *timer_weighting,int use_fragmentary)
{
	ESL_MSA *msa=NULL;
	ESL_MSA *temp_msa=NULL;
	int nali=0;
	int status;
	int fmt=eslMSAFILE_UNKNOWN;
	ESLX_MSAFILE *afp=NULL;
	if((status = eslx_msafile_Open(&abc, alifile, NULL, fmt, NULL, &afp)) != eslOK) { 
		eslx_msafile_OpenFailure(afp, status);
	}
	while ((status = eslx_msafile_Read(afp, &temp_msa)) != eslEOF) {
		if (status != eslOK) {
			eslx_msafile_ReadFailure(afp, status);
		}
		msa=temp_msa;
		nali++;
		if (nali>1) {
			cm_Fail("alignment file %s has more than one alignment in it, which CMfinder does not allow",alifile);
		}
	}
	if (msa==NULL) {
		cm_Fail("didn't get any alignment from %s",alifile);
	}

	TrimDegenOnEndsOfMSA (msa,NULL);
	MultiplyEMWeightsByWeightingStrategy (msa,ws,timer_weighting,use_fragmentary);
	
	eslx_msafile_Close(afp);
	
	return msa;
}

char *PrintInProgressFileName (const char *final_file,int iteration,const char *stepDescription,const char *fileExt)
{
	char *filename=(char *)MallocOrDie(strlen(final_file)+100);
	sprintf(filename,"%s.inprogress.%d.%s.%s",final_file,iteration,stepDescription,fileExt);
	printf("writing in-progress MSA to file %s\n",filename);
	return filename;
}
void SaveInProgressMsa(const char *final_file,int iteration,const char *stepDescription,ESL_MSA *msa)
{
	char *filename=PrintInProgressFileName (final_file,iteration,stepDescription,"sto");
	SaveMsa(filename,msa);
	free(filename);
}
void SaveInProgressCm(const char *final_file,int iteration,const char *stepDescription,CM_t *cm)
{
	char *filename=PrintInProgressFileName (final_file,iteration,stepDescription,"cm");
        if (cm!=NULL) {
          SaveCM(filename,cm);
        }
        else {
          printf("(not saving %s because cm is NULL)\n",filename);
        }
        free(filename);
}

void InitInfernal (void)
{
	/* union of initializations from cmbuild.c and cmsearch.c */
	init_ilogsum();
	FLogsumInit();
	p7_FLogsumInit();
	impl_Init();
}

DegenStrategy GetDegenStrategy(ESL_GETOPTS *go)
{
	DegenStrategy degenStrategy=DegenKeep;
	if (esl_opt_GetBoolean(go, "--degen-rand")) {
		degenStrategy=DegenRand;
	}
	if (esl_opt_GetBoolean(go, "--degen-keep")) {
		degenStrategy=DegenKeep;
	}
	return degenStrategy;
}

/*
 * code copied from the L=0, -d case of function 'seq_shuffling' in esl_shuffle.c
 */
void BgSeqCache_Create(ESL_SQCACHE **ret_bgSeqCache,ESL_SQCACHE *dbSeqs,CmfinderVars *cmfinderVars)
{
	ESL_SQCACHE *bgSeqCache;
	int inputSeqnum,copynum;
	int targetSeqnum;
	ESL_SQ *shuff=esl_sq_Create();
	
	bgSeqCache=(ESL_SQCACHE *)MallocOrDie(sizeof(ESL_SQCACHE));
	
	bgSeqCache->seq_count=dbSeqs->seq_count * cmfinderVars->bgScoreSize; /* actually this is an upper bound on what we'll need; later we'll set the actual amount.  I expect the extra memory will be negligeable */
	if (bgSeqCache->seq_count==0) {
		/* no need to do any more */
	}
	else {
		bgSeqCache->sq_list=(ESL_SQ *)MallocOrDie(sizeof(ESL_SQ)*bgSeqCache->seq_count);
		
		targetSeqnum=0;
		for (inputSeqnum=0; inputSeqnum<dbSeqs->seq_count; inputSeqnum++) {

			ESL_SQ *src_sq=&(dbSeqs->sq_list[inputSeqnum]);
			ESL_SQ *text_sq;
			
			/* make text-version of input sequence */
			text_sq=esl_sq_Create();
			esl_sq_Copy(src_sq,text_sq);
			esl_sq_Textize(text_sq);
			
			/* alloc memory for shuffled version */
			esl_sq_GrowTo(shuff, src_sq->n); /* make sure shuff can hold sq */	  
			shuff->n = src_sq->n;

			for (copynum=0; copynum<cmfinderVars->bgScoreSize; copynum++) {
				ESL_SQ *sq=&(bgSeqCache->sq_list[targetSeqnum]);
				
				if (esl_rsq_CShuffleDP(cmfinderVars->randomness, text_sq->seq, shuff->seq) != eslOK) {
					esl_fatal("esl_sq_FormatName failed");
				}
				if (esl_sq_FormatName(shuff, "%s-shuffled-%d", src_sq->name, copynum) != eslOK) {
					esl_fatal("esl_sq_FormatName failed");
				}
				
				sq->L=src_sq->L;
				sq->abc=src_sq->abc;
				sq->acc="dummy-acc";
				sq->desc="dummy-desc";
				sq->idx=-1;
				sq->source=NULL;
				sq->ss=NULL;
				sq->n=sq->L;
				sq->start=1;
				sq->end=sq->L;

				sq->name=(char *)MallocOrDie(strlen(shuff->name)+1);
				strcpy(sq->name,shuff->name);
				sq->seq=NULL;  /* we need a digitial seq for scanning */
				sq->dsq=(ESL_DSQ *)MallocOrDie(sizeof(ESL_DSQ)*(sq->L+2));
				if (esl_abc_Digitize(abc, shuff->seq, sq->dsq) != eslOK) {
					esl_fatal("esl_abc_Digitize failed");
				}
				
				targetSeqnum++;
			}
			
			esl_sq_Destroy(text_sq);
		}
	}
	
	esl_sq_Destroy(shuff);
	
	*ret_bgSeqCache=bgSeqCache;
}
void BgSeqCache_Destroy(ESL_SQCACHE *bgSeqCache)
{
	int i;
	for (i=0; i<bgSeqCache->seq_count; i++) {
		free(bgSeqCache->sq_list[i].name);
		free(bgSeqCache->sq_list[i].dsq);
	}
	if (bgSeqCache->seq_count>0) {
		free(bgSeqCache->sq_list);
	}
	free(bgSeqCache);
}

int compare_sq_len (const void *voidx,const void *voidy)
{
	const ESL_SQ *x=(const ESL_SQ*)voidx;
	const ESL_SQ *y=(const ESL_SQ*)voidy;
	if (x->L==y->L) {
		return 0;
	}
	if (x->L < y->L) {
		return +1;
	}
	else {
		return -1;
	}
}
void SortSqcacheLongestToShortest(ESL_SQCACHE *dbSeqs) /* doesn't reduce parallel cmsearch time, at least in the one case I tried */
{
	qsort(dbSeqs->sq_list,dbSeqs->seq_count,sizeof(ESL_SQ),compare_sq_len);
}

int main_cmfinder(ESL_GETOPTS *go)
{
	int exitCode=1;
	int needSeqFile,numExpectedCommandlineParams;
	char *alifile=NULL;
	char *seqfile=NULL;
	ESL_MSA *digital_msa=NULL;
	ESL_SQCACHE *inputSeqCache=NULL;
	ESL_SQCACHE *bgSeqCache=NULL;
	Cand **cand;
	int i;
	int *ncand;
	double **pxy = NULL;  /* partition function log ratio */  
	CmfinderVars alloc_cmfinderVars;
	CmfinderVars *cmfinderVars=&alloc_cmfinderVars;
	const char *createFileNameOnSuccess=NULL;
	FILE *timerAppendFile=NULL;
	typedef struct TimerInfo_ {
		ESL_STOPWATCH **timer;
		const char *desc;
	} TimerInfo;
	TimerInfo timerInfoList[]={
		{&(cmfinderVars->timer_overall),"overall"},
		{&(cmfinderVars->timer_partition_function),"partition_function"},
		{&(cmfinderVars->timer_m_step),"m_step"},
		{&(cmfinderVars->timer_e_step),"e_step"},
		{&(cmfinderVars->timer_cmsearch_viterbi),"cmsearch_viterbi"},
		{&(cmfinderVars->timer_cmsearch_inside),"cmsearch_inside"},
		{&(cmfinderVars->timer_cmbuild),"cmbuild"},
		{&(cmfinderVars->timer_elim_iden_seqs),"elim_iden_seqs"},
		{&(cmfinderVars->timer_weighting),"weighting"},
		{&(cmfinderVars->timer_cmcalibrate),"cmcalibrate"}
	};

	cmfinderVars->gapthreshold=0.6;
	/* the original cmfinder code calls this 'gapcost', but I think this is a misnomer.
	 * Also, in the original cmfinder.c code, it first uses cost=50 (which is hard-coded), then cost=100.  
	 * I think this is an oversight, and won't change the results much anyway 
	 * but I want to have to option to emulate the old behavior, hence the two separate vars. */
	cmfinderVars->openPairCost=100;
	cmfinderVars->openPairCost_initial=50;
	cmfinderVars->max_cand_for_ScanCand=10; /* I'm not sure why this is hard coded, and why it's not equal to maxNumCand */
	cmfinderVars->cm_scoreThreshold=0;
	cmfinderVars->max_cm_scoreThreshold=30;
	cmfinderVars->overlyHighNumSeqToTotalCandRatio=5;
        cmfinderVars->minPredictedBasePairs=3;
	cmfinderVars->cm_scoreThreshold_increment=10;
	cmfinderVars->minCandScore=10;
	cmfinderVars->minCandScoreToBestScoreRatio=2;
	cmfinderVars->best_totscore = DOUBLE_NEGINFINITY;  
	cmfinderVars->oldscore = DOUBLE_NEGINFINITY;  
	cmfinderVars->iteration = 0;
	cmfinderVars->max_iteration=100;
	cmfinderVars->min_acceptable_totscore=-10;
	cmfinderVars->min_acceptable_totweight=2.5;
	cmfinderVars->weightingStrategy=WS_EM_weights_only;
	cmfinderVars->cmsearch_mxsize=4000;
	cmfinderVars->cmsearch_smxsize=4000;
	cmfinderVars->DB_length = 100000;	
	cmfinderVars->do_zoop=0;
	cmfinderVars->zoop_gamma=0.3;
	cmfinderVars->printedBestInsideScoreLessThanViterbi=0;
	cmfinderVars->convergence_threshold=0.02;
	cmfinderVars->best_msa=NULL;
	cmfinderVars->putStartEndCoordsInHitIds=0; /* compatibility with CMfinder 0.3 */
	cmfinderVars->use_fragmentary=0;
	cmfinderVars->filter_noF4F5=0;
	cmfinderVars->filter_max=0;
	cmfinderVars->use_maxAlignAccuracy=0;
	cmfinderVars->do_local=0;
	cmfinderVars->useIntsForMiAndPrLikeOldCmfinder=0;
        cmfinderVars->avoidPartFunc=0;
        cmfinderVars->columnOnlyBasePairProbs=0;
	cmfinderVars->min_acceptable_seq_weight=0.01;
	cmfinderVars->eliminateIdenticalSeqs=1;
	cmfinderVars->eliminateIdenticalSubseqs=1;
	cmfinderVars->dieOnUntestedCode=1;
	cmfinderVars->min_seq_score_in_final_msa=0;
	cmfinderVars->numCpu=0;
	cmfinderVars->bgScoreSize=0;
	cmfinderVars->bgScoreEvalue=-1;
	cmfinderVars->use_evalue=0;
	cmfinderVars->cmcalibrate_tailn=-1;
	cmfinderVars->bgScoreScanBitThresh=0;
	cmfinderVars->bgScoreNonFrag=0;
	cmfinderVars->filterNonFrag=0;
	cmfinderVars->filterNonFragPad=10;
	cmfinderVars->saveMsaAfterFirstMstep=0;
	cmfinderVars->saveInProgress=0;
	cmfinderVars->nonFragmentaryAvgBppr=0;
	cmfinderVars->eliminateSeqsWithManyDegen=0;
	cmfinderVars->maxDegenPerSeq=3;
	cmfinderVars->flankingNucsForCountingDegen=0;
	cmfinderVars->shiftDistalMispairsIntoTerminalLoops=0;
	cmfinderVars->minTermLoopNucsWithMovingIntoTerminalLoops=0;

	for (i=0; i<sizeof(timerInfoList)/sizeof(TimerInfo); i++) {
		*(timerInfoList[i].timer)=esl_stopwatch_Create();
	}
	esl_stopwatch_Start(cmfinderVars->timer_overall);

	InitInfernal();

	needSeqFile=1;
	if (esl_opt_GetBoolean(go, "--bad-distal-pairs-to-loop-only")) {
		needSeqFile=0;
	}
	
	numExpectedCommandlineParams=(needSeqFile?1:0);
	if (esl_opt_ArgNumber(go) != numExpectedCommandlineParams) {
		int i;
		for (i=1; i<esl_opt_ArgNumber(go); i++) {
			printf("arg #%d: \"%s\"\n",i,esl_opt_GetArg(go,i));
		}
		cm_Fail("Incorrect number of command line arguments.");
	}
	if (needSeqFile) {
		seqfile=esl_opt_GetArg(go,1);
		if (seqfile==NULL) { cm_Fail("didn't get seqfile argument on commandline");	}
	}

	alifile=esl_opt_GetString(go, "-a");
	cmfinderVars->final_file=esl_opt_GetString(go, "-o");

	if (alifile==NULL) {
		esl_fatal("you must give an input alignment file with the -a option");
	}
	if (cmfinderVars->final_file==NULL) {
		esl_fatal("you must give an output file with the -o option");
	}

	cmfinderVars->degenStrategy=GetDegenStrategy(go);
	if (esl_opt_GetBoolean(go, "--fragmentary")) {
		cmfinderVars->use_fragmentary=1;
	}
	if (esl_opt_GetBoolean(go, "--noF4F5")) {
		cmfinderVars->filter_noF4F5=1;
	}
	if (esl_opt_GetBoolean(go, "--max")) {
		cmfinderVars->filter_max=1;
	}
	if (esl_opt_GetBoolean(go, "--amaa")) {
		cmfinderVars->use_maxAlignAccuracy=1;
	}
	if (esl_opt_GetBoolean(go, "--local")) {
		cmfinderVars->do_local=1;
	}
	if (esl_opt_GetBoolean(go, "--avoid-part-func")) {
                cmfinderVars->avoidPartFunc=1;
        }
	if (esl_opt_GetBoolean(go, "--column-only-base-pair-probs")) {
                cmfinderVars->columnOnlyBasePairProbs=1;
        }
	if (esl_opt_GetBoolean(go, "--ints-like-03")) {
		cmfinderVars->useIntsForMiAndPrLikeOldCmfinder=1;
	}
	if (esl_opt_IsOn(go, "--min-seq-weight")) {
		cmfinderVars->min_acceptable_seq_weight=esl_opt_GetReal(go,"--min-seq-weight");
	}
	if (esl_opt_IsOn(go, "--min-base-pairs")) {
		cmfinderVars->minPredictedBasePairs=esl_opt_GetInteger(go,"--min-base-pairs");
	}
	if (esl_opt_GetBoolean(go, "--no-elim-iden-seq")) {
		cmfinderVars->eliminateIdenticalSeqs=0;
	}
	if (esl_opt_GetBoolean(go, "--no-elim-iden-subseq")) {
		cmfinderVars->eliminateIdenticalSubseqs=0;
	}
	if (esl_opt_GetBoolean(go, "--allow-untested-code")) {
		cmfinderVars->dieOnUntestedCode=0;
	}
	if (esl_opt_IsOn(go, "--min-cand-score-in-final")) {
		cmfinderVars->min_seq_score_in_final_msa=esl_opt_GetReal(go,"--min-cand-score-in-final");
	}
	if (esl_opt_IsOn(go, "--open-pair-cost")) {
		cmfinderVars->openPairCost=esl_opt_GetReal(go,"--open-pair-cost");
	}
	if (esl_opt_IsOn(go, "--open-pair-cost-initial")) {
		cmfinderVars->openPairCost_initial=esl_opt_GetReal(go,"--open-pair-cost-initial");
	}
#ifdef HMMER_THREADS
	if (esl_opt_IsOn(go ,"--cpu")) {
		cmfinderVars->numCpu=esl_opt_GetInteger(go,"--cpu");
		if (cmfinderVars->numCpu<0) {
			if (esl_threads_CPUCount(&(cmfinderVars->numCpu)) != eslOK) {
				esl_fatal("esl_threads_CPUCount failed");
			}
		}
	}
#endif
	if (esl_opt_IsOn(go, "--bg-score-size")) {
		cmfinderVars->bgScoreSize=esl_opt_GetInteger(go,"--bg-score-size");
	}
	if (esl_opt_IsOn(go, "--evalue")) {
		cmfinderVars->use_evalue=1;
		cmfinderVars->evalue=esl_opt_GetReal(go,"--evalue");
	}
	if (esl_opt_IsOn(go, "--tailn")) {
		cmfinderVars->cmcalibrate_tailn=esl_opt_GetInteger(go,"--tailn");
	}
	if (esl_opt_IsOn(go,"--create-file-on-success")) {
		createFileNameOnSuccess=esl_opt_GetString(go,"--create-file-on-success");
	}
	if (esl_opt_IsOn(go, "--bg-score-evalue")) {
		cmfinderVars->bgScoreEvalue=esl_opt_GetReal(go,"--bg-score-evalue");
	}
	if (esl_opt_IsOn(go, "--bg-score-scan-thresh")) {
		cmfinderVars->bgScoreScanBitThresh=esl_opt_GetReal(go, "--bg-score-scan-thresh");
	}
	if (esl_opt_GetBoolean(go, "--bg-score-non-frag")) {
		cmfinderVars->bgScoreNonFrag=1;
	}
	if (esl_opt_GetBoolean(go, "--filter-non-frag")) {
		if (!cmfinderVars->use_fragmentary) {
			esl_fatal("--filter-non-frag used without --fragmentary : that doesn't make sense");
		}
		cmfinderVars->filterNonFrag=1;
	}
	if (esl_opt_IsOn(go, "--filter-non-frag-pad")) {
		cmfinderVars->filterNonFragPad=esl_opt_GetInteger(go,"--filter-non-frag-pad");
	}
	if (esl_opt_GetBoolean(go, "--save-after-first-M-step")) {
		cmfinderVars->saveMsaAfterFirstMstep=1;
	}
	if (esl_opt_GetBoolean(go, "--save-in-progress")) {
		cmfinderVars->saveInProgress=1;
	}
	if (esl_opt_GetBoolean(go, "--non-frag-avg-bppr")) {
		cmfinderVars->nonFragmentaryAvgBppr=1;
	}
	if (esl_opt_GetBoolean(go, "--wgsc")) {
		cmfinderVars->weightingStrategy=WS_GSC_and_EM;
	}
	if (esl_opt_GetBoolean(go, "--wpb")) {
		cmfinderVars->weightingStrategy=WS_PB_and_EM;
	}
	if (esl_opt_IsOn(go,"--timer-append")) {
		const char *fn=esl_opt_GetString(go,"--timer-append");
		timerAppendFile=fopen(fn,"at");
		if (timerAppendFile==NULL) {
			esl_fatal("cannot open file %s",fn);
		}
	}
	if (esl_opt_IsOn(go, "--max-degen-per-hit")) {
		cmfinderVars->eliminateSeqsWithManyDegen=1;
		cmfinderVars->maxDegenPerSeq=esl_opt_GetInteger(go,"--max-degen-per-hit");
	}
	if (esl_opt_IsOn(go, "--max-degen-flanking-nucs")) {
		if (!esl_opt_IsOn(go, "--max-degen-per-hit")) {
			esl_fatal("--max-degen-flanking-nucs doesn't make sense without --max-degen-per-hit");
		}
		cmfinderVars->flankingNucsForCountingDegen=esl_opt_GetInteger(go,"--max-degen-flanking-nucs");
	}
	if (esl_opt_GetBoolean(go, "--bad-distal-pairs-to-loop")) {
		cmfinderVars->shiftDistalMispairsIntoTerminalLoops=1;
	}
	if (esl_opt_IsOn(go, "--min-term-loop-nucs")) {
		cmfinderVars->minTermLoopNucsWithMovingIntoTerminalLoops=esl_opt_GetInteger(go,"--min-term-loop-nucs");
	}
			
	for (i=0; i<CM_p7_NEVPARAM; i++) {
		cmfinderVars->hmmStatsToCalibrate[i]=1;
	}
	if (cmfinderVars->filter_max) {
		for (i=0; i<CM_p7_NEVPARAM; i++) {
			cmfinderVars->hmmStatsToCalibrate[i]=0;
		}
	}
	if (cmfinderVars->filter_noF4F5) {
		cmfinderVars->hmmStatsToCalibrate[CM_p7_GFMU]=0;
		cmfinderVars->hmmStatsToCalibrate[CM_p7_GFLAMBDA]=0;
	}
	
	cmfinderVars->go=go;
	
	cmfinderVars->randomness = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));

	if (needSeqFile) {
		/* load input seqs into cache, and create array of their lengths, for use with the cm_pipeline code */
		if (esl_sqfile_Cache(abc,seqfile,eslSQFILE_FASTA,NULL,&inputSeqCache) != eslOK) {
			esl_fatal("Problem with opening input sequence file seqfile=%s",seqfile);
		}
		cmfinderVars->inputSeqSize=0;
		for (i=0; i<inputSeqCache->seq_count; i++) {
			assert(inputSeqCache->sq_list[i].n==inputSeqCache->sq_list[i].L); /* no subseq weirdness */
			cmfinderVars->inputSeqSize += inputSeqCache->sq_list[i].n;
		}
		if (cmfinderVars->degenStrategy==DegenRand) {
			RandomizeDegenInAllSeqs (inputSeqCache);
		}
		if (cmfinderVars->numCpu>1) {
			SortSqcacheLongestToShortest(inputSeqCache);
		}

		BgSeqCache_Create(&bgSeqCache,inputSeqCache,cmfinderVars);
		cmfinderVars->bgSeqCache=bgSeqCache;
	}

	/* older CMfinder has the ability to get the starting point from a .cm file, but we're hardcoded to starting from a .sto file */
	digital_msa=LoadMSAFromFile(alifile,cmfinderVars->weightingStrategy,cmfinderVars->timer_weighting, cmfinderVars->use_fragmentary); /* seems to be loaded in digital */
	if ((digital_msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("loaded MSA in file %s which is not in digital mode (in internal results returned from the easel library).  try converting to .sto format");
	}
	if (cmfinderVars->degenStrategy==DegenRand) {
		RandomizeDegenInMsaDigital(digital_msa);
	}
	
	if (esl_opt_GetBoolean(go, "--bad-distal-pairs-to-loop-only")) {
		ESL_MSA *new_msa;
		int *seq_fragment_start_col,*seq_fragment_end_col;
		
		seq_fragment_start_col=(int *)MallocOrDie(digital_msa->nseq * sizeof(int));
		seq_fragment_end_col=(int *)MallocOrDie(digital_msa->nseq * sizeof(int));
		CalcFragmentStartEnds (digital_msa,cmfinderVars->use_fragmentary,seq_fragment_start_col,seq_fragment_end_col);
		cmfinderVars->shiftDistalMispairsIntoTerminalLoops=1;
		new_msa=ShiftDistalMispairsIntoTerminalLoops (digital_msa,cmfinderVars,seq_fragment_start_col,seq_fragment_end_col);
		SaveMsa(cmfinderVars->final_file,new_msa);
		printf("completed --bad-distal-pairs-to-loop-only (saved to %s)\n",cmfinderVars->final_file);
		exit(0);
	}

	CM_t *cm = M_step(digital_msa, NULL /* pxy */,cmfinderVars,cmfinderVars->openPairCost_initial);
	if (cmfinderVars->saveMsaAfterFirstMstep) {
		printf("stopping after M-step\n");
		exit(0);
	}
	if (cmfinderVars->saveInProgress) {
		SaveInProgressMsa(cmfinderVars->final_file,cmfinderVars->iteration,"after-M-step",digital_msa);
                SaveInProgressCm(cmfinderVars->final_file,cmfinderVars->iteration,"after-M-step",cm);
	}
	
	if (cm==NULL) {
          printf("fyi: cm returned from M_step is NULL, so no acceptable structures\n");
		cmfinderVars->cannotFindAcceptableStructure=1;
	}
	else {
		
		InitCand(&cand,&ncand,inputSeqCache->seq_count);
		
		cmfinderVars->iteration=0;
		cmfinderVars->cannotFindAcceptableStructure=0;
		while (cmfinderVars->iteration < cmfinderVars->max_iteration) {
			double totscore,delta;
			
			cmfinderVars->iteration++;
			
			totscore=0;
			esl_msa_Destroy(digital_msa); /* about to make a new one */
			digital_msa=NULL;
			digital_msa = E_step(cm, inputSeqCache, ncand, cand, &totscore, &pxy, cmfinderVars);
			FreeCM(cm); /* don't need this any more */

			if (cmfinderVars->saveInProgress) {
				SaveInProgressMsa(cmfinderVars->final_file,cmfinderVars->iteration,"after-E-step",digital_msa);
			}

			if (totscore < cmfinderVars->min_acceptable_totscore) {
				printf("Bad alignment (totscore=%f)\n", totscore);
				cmfinderVars->cannotFindAcceptableStructure=1;
			}
			
			if (cmfinderVars->cannotFindAcceptableStructure) {
				break;
			}

			/* If we've converged, stop. */
			delta = (totscore - cmfinderVars->oldscore) / fabs(totscore);
			
			if (1) {
				printf("score %g, delta %g\n", totscore / (double)(inputSeqCache->seq_count), delta);
			}

			if (delta < cmfinderVars->convergence_threshold) {
				printf("converged in iteration %d (delta=%g)\n",cmfinderVars->iteration,delta);
				break;
			}
			cmfinderVars->oldscore=totscore;
			
			/* Else, make a new model from the alignment. */

			cm = M_step(digital_msa, pxy,cmfinderVars,cmfinderVars->openPairCost);
			if (cm==NULL) {
                          printf("fyi: cm returned from M_step is NULL, so no acceptable structures\n");
				cmfinderVars->cannotFindAcceptableStructure=1;
				break;
			}

			if (cmfinderVars->saveInProgress) {
				SaveInProgressMsa(cmfinderVars->final_file,cmfinderVars->iteration,"after-M-step",digital_msa);
				SaveInProgressCm(cmfinderVars->final_file,cmfinderVars->iteration,"after-M-step",cm);
			}
		}

		DeleteCand(cand,ncand,inputSeqCache->seq_count);
	}
	
	if (cmfinderVars->cannotFindAcceptableStructure) {
		printf("failed to find acceptable structure prediction after %d iterations\n",cmfinderVars->iteration);
		exitCode=2;
	}
	else {
		if (cmfinderVars->best_msa!=NULL && cmfinderVars->final_file!=NULL) {
			SaveMsa(cmfinderVars->final_file,cmfinderVars->best_msa);
		}
		exitCode=0;
	}
	if (cmfinderVars->best_msa!=NULL) {
		esl_msa_Destroy(cmfinderVars->best_msa);
		cmfinderVars->best_msa=NULL;
	}
	
	if (createFileNameOnSuccess!=NULL) {
		FILE *f=fopen(createFileNameOnSuccess,"w");
		if (f==NULL) {
			esl_fatal("can't open %s",createFileNameOnSuccess);
		}
		fclose(f);
	}

	esl_stopwatch_Stop(cmfinderVars->timer_overall);

	if (timerAppendFile!=NULL) {
		for (i=0; i<sizeof(timerInfoList)/sizeof(TimerInfo); i++) {
			ESL_STOPWATCH *timer=*(timerInfoList[i].timer);
			fprintf(timerAppendFile,"cmfinder04-%s\t%lg\t%lg\n",timerInfoList[i].desc,timer->user+timer->sys,timer->elapsed);
		}
		fclose(timerAppendFile);
	}
	
	printf("timers:\n");
	for (i=0; i<sizeof(timerInfoList)/sizeof(TimerInfo); i++) {
		char buf[256];
		sprintf(buf,"timer %s ",timerInfoList[i].desc);
		esl_stopwatch_Display(stdout,*(timerInfoList[i].timer),buf);
	}
	for (i=0; i<sizeof(timerInfoList)/sizeof(TimerInfo); i++) {
		esl_stopwatch_Destroy(*(timerInfoList[i].timer));
	}
	
	esl_getopts_Destroy(go);
	if (digital_msa!=NULL) esl_msa_Destroy(digital_msa);
	if (inputSeqCache!=NULL) esl_sqfile_Free(inputSeqCache);
	if (bgSeqCache!=NULL) BgSeqCache_Destroy(bgSeqCache);
	esl_randomness_Destroy(cmfinderVars->randomness);
	
	return exitCode;
}

int main_summarize(ESL_GETOPTS *go)
{
	ESL_MSA *digital_msa,*text_msa;
	ESL_STOPWATCH *timer_weighting_dummy=NULL;
	char *alifile;
	WeightingStrategy weightingStrategy;
	DegenStrategy degenStrategy;
	double nullModel[NUM_NUCS] = {0.25, 0.25, 0.25, 0.25};
	double tot_weight,e,m,avg_bppr,avg_bp,avg_score,avg_seq_id,avg_seq_len,avg_energy,avg_GC;
	double conf_bp=0,del_bp=0;
	int cons_pos;
	int i;
	int *pt;
	double **bp_pr=NULL;
	double *any_bp_pr=NULL;
	double yaoFormulaScore;
	double numSpecies;
	int *seq_fragment_start_col,*seq_fragment_end_col;
	int use_fragmentary=0;
	FILE *timerAppendFile=NULL;
	ESL_STOPWATCH *timer_summarize=esl_stopwatch_Create();
	esl_stopwatch_Start(timer_summarize);
	int skipDeTag=0;

	weightingStrategy=WS_EM_weights_only;
	degenStrategy=DegenKeep;
	
	if (esl_opt_ArgNumber(go)                 != 1)     { cm_Fail("Incorrect number of command line arguments."); }
	alifile=esl_opt_GetArg(go,1);
	if (alifile==NULL) {
		esl_fatal("can't get alifile");
	}

	
	if (esl_opt_GetBoolean(go, "--summarize-no-de-tag")) {
	  skipDeTag=1;
	}
	if (esl_opt_GetBoolean(go, "--summarize-gsc")) {
		weightingStrategy=WS_GSC;
	}
	degenStrategy=GetDegenStrategy(go);
	
	if (esl_opt_GetBoolean(go, "--fragmentary")) {
		use_fragmentary=1;
	}
	
	if (esl_opt_IsOn(go,"--timer-append")) {
		const char *fn=esl_opt_GetString(go,"--timer-append");
		timerAppendFile=fopen(fn,"at");
		if (timerAppendFile==NULL) {
			esl_fatal("cannot open file %s",fn);
		}
	}

	digital_msa=LoadMSAFromFile(alifile,weightingStrategy,timer_weighting_dummy,use_fragmentary); /* seems to be loaded in digital */
	if ((digital_msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("loaded MSA in file %s which is not in digital mode (in internal results returned from the easel library).  try converting to .sto format");
	}
	if (degenStrategy==DegenRand) {
		RandomizeDegenInMsaDigital(digital_msa);
	}
	
	if (esl_opt_IsOn(go, "--summarize-save-msa")) {
		const char *fn=esl_opt_GetString(go,"--summarize-save-msa");
		SaveMsa(fn,digital_msa);
	}
	
	if (digital_msa->ss_cons==NULL) {
		esl_fatal("input MSA does not have an SS_cons line");
	}
	if (!skipDeTag) {
	  if (digital_msa->ss==NULL) {
	    esl_fatal("input MSA does not have per-sequence #=GR ... SS lines (this program could implement the cons_ss2ss function from summarize.c in CMfinder 0.3, but this isn't done)");
	  }
	}
	
	/* figure out start/end's due to truncation */
	seq_fragment_start_col=(int *)MallocOrDie(digital_msa->nseq * sizeof(int));
	seq_fragment_end_col=(int *)MallocOrDie(digital_msa->nseq * sizeof(int));
	CalcFragmentStartEnds (digital_msa,use_fragmentary,seq_fragment_start_col,seq_fragment_end_col);
	
	text_msa=MakeTextMsaClone(digital_msa);
	Avg_bppr(text_msa, NULL, use_fragmentary, seq_fragment_start_col, seq_fragment_end_col, &bp_pr, &any_bp_pr);
	pt = GetPairtable(digital_msa->ss_cons);
	
	tot_weight=0;
	for (i = 0; i < digital_msa->nseq; i++) {
		tot_weight += digital_msa->wgt[i];
	}
	
	for (i=0; i<NUM_NUCS; i++) {
		if (nullModel[i]<=0 || nullModel[i]>=1) {
			esl_fatal("nullModel probabilities must be greater than 0 and less than 1.  Equality to zero causes degenerate cases in some code.");
		}
	}
	
	e=entropy( digital_msa,nullModel,tot_weight);
	m= summarize_mxy(digital_msa, nullModel, tot_weight,pt);
	avg_bppr= weighted_base_pair(digital_msa,pt,bp_pr);
	avg_bp= average_base_pair(digital_msa,pt);
	if (skipDeTag) {
	  avg_score=0;
	  avg_energy=0;
	}
	else {
	  avg_score= average_score(digital_msa,tot_weight);
	  avg_energy = average_energy(text_msa, tot_weight);  
	}
	avg_seq_id = average_seq_id(digital_msa,tot_weight,use_fragmentary,seq_fragment_start_col,seq_fragment_end_col);
	avg_seq_len = average_seq_len(digital_msa,tot_weight);
	avg_GC = average_GC(digital_msa, tot_weight);  
	cons_pos = conserved_position(digital_msa, tot_weight,use_fragmentary,seq_fragment_start_col,seq_fragment_end_col);
	bad_base_pair(digital_msa, pt, bp_pr, &conf_bp, &del_bp);
	numSpecies=(double)(CalcQuasiNumSpecies(digital_msa));

	/* formula from Zizhen Yao's dissertation */
	yaoFormulaScore=sqrt((cons_pos+0.2)*avg_bppr/avg_seq_id)*numSpecies*(1+log((double)(digital_msa->nseq)/numSpecies));


	if (skipDeTag) {
	  fprintf(stderr,"WARNING: skipped calculations of \"Score\" and \"Energy\" due to --summarize-no-de-tag\n");
	}
	printf("Num=%d\t Weight=%.2f\t Len=%.1f\t Score=%.2f\t Entropy=%.2f\t MI=%.2f\t BP=%.2f\t BP.org=%.2f\t Seq_id=%.2f\t Energy=%.2f\t GC=%.2f\t Conserved_pos=%d\t Conf_bp=%.2f\t Del_bp=%.2f\t yaoFormulaScore=%f\t numSpecies=%f\n",
	       digital_msa->nseq, tot_weight, avg_seq_len, avg_score, e, m, avg_bppr, avg_bp, avg_seq_id, avg_energy, avg_GC, cons_pos, conf_bp, del_bp,yaoFormulaScore,numSpecies);

	esl_stopwatch_Stop(timer_summarize);
	if (timerAppendFile!=NULL) {
		fprintf(timerAppendFile,"summarize\t%lg\t%lg\n",timer_summarize->user+timer_summarize->sys,timer_summarize->elapsed);
		fclose(timerAppendFile);
	}
	esl_stopwatch_Destroy(timer_summarize);

	free(seq_fragment_start_col);
	free(seq_fragment_end_col);
    Free2DArray((void **)bp_pr, text_msa->alen+1);
	free(any_bp_pr);
	free(pt);
	esl_msa_Destroy(digital_msa);
	esl_msa_Destroy(text_msa);

	return 0;
}

int main (int argc,char *argv[])
{
	ESL_GETOPTS *go     = NULL;
	int status;

	abc=esl_alphabet_Create(eslRNA); /* not sure why the distinction between RNA and DNA */

	if ((go = esl_getopts_Create(cmfinder_options))   == NULL) { cm_Fail("Internal failure creating options object");}
	if (esl_opt_ProcessEnvironment(go)         != eslOK)  { esl_fatal("Failed to process environment: %s\n", go->errbuf); }
	if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { esl_fatal("Failed to parse command line: %s\n", go->errbuf); }
	if (esl_opt_VerifyConfig(go)               != eslOK)  { esl_fatal("Failed to parse command line: %s\n", go->errbuf); }

	if (esl_opt_GetBoolean(go, "-h")) {
		const int textwidth=16000; // yup, that's a bit big, but I don't want to shorten things
		const int indentation=2;
		esl_usage(stdout, argv[0], usage);
		
		puts("\nBasic options:");
		esl_opt_DisplayHelp(stdout, go, 1, indentation, textwidth);
		puts("\nGeneral cmfinder options:");
		esl_opt_DisplayHelp(stdout, go, 2, indentation, textwidth);
		puts("\noptions related to internal cmcalibrate/cmsearch:");
		esl_opt_DisplayHelp(stdout, go, 3, indentation, textwidth);
		puts("\noptions related to internal cmbuild:");
		esl_opt_DisplayHelp(stdout, go, 4, indentation, textwidth);
		puts("\noptions for --summarizea:");
		esl_opt_DisplayHelp(stdout, go, 5, indentation, textwidth);
		
		exit(0);
	}

	if (esl_opt_GetBoolean(go,"--version")) {
		printf("CMFINDER_PACKAGE_VERSION=%s.\n",CMFINDER_PACKAGE_VERSION);
		exit(0);
	}

	if (esl_opt_GetBoolean(go, "--summarize")) {
		status=main_summarize(go);
	}
	else {
		status=main_cmfinder(go);		
	}
	
	esl_alphabet_Destroy(abc);
	
	return status;
}
