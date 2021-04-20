extern ESL_ALPHABET *abc;

#define INT_NEGINFINITY  -9999999 /* these are for 'mixy' function, which uses ints.  hopefully I can convert it to floats/doubles at some point. */
#define INT_POSINFINITY  9999999

#define DOUBLE_NEGINFINITY INT_NEGINFINITY

#define NUM_NUCS 4

#define MIN_HAIRPIN 3
#define INT_MINSCORE -400

#define pair_left(c) (c == '(' || c == '<' || c == '{' || c == '[')
#define pair_right(c) (c == ')' || c == '>' || c == '}' || c == ']')

/* functions copied from infernal-0.72 */
#define sreLOG2(x)  ((x) > 0 ? log(x) * 1.44269504 : IMPOSSIBLE)
#define sreEXP2(x)  (exp((x) * 0.69314718 ))

#define MAX_SQ_DESC 128 /* should be plenty */
#define MAX_MSA_AU 128

typedef enum {
	WS_EM_weights_only,
	WS_uniform_weights,
	WS_GSC_and_EM,
	WS_PB_and_EM,
	WS_GSC
} WeightingStrategy;
typedef enum {
	DegenRand,DegenKeep
} DegenStrategy;

typedef struct {
	ESL_GETOPTS *go;
	WeightingStrategy weightingStrategy;
	DegenStrategy degenStrategy;
	int do_local;
	int use_maxAlignAccuracy;
	int useIntsForMiAndPrLikeOldCmfinder;
        int avoidPartFunc;
        int columnOnlyBasePairProbs;
	int nonFragmentaryAvgBppr;
	int filter_noF4F5,filter_max;
	int use_fragmentary;
	int do_zoop;
	int iteration,max_iteration,cannotFindAcceptableStructure;
	int max_cand_for_ScanCand;
	double gapthreshold;
	char *final_file;
	int inputSeqSize;
	double min_acceptable_totscore; /* was hardcoded as -10 in cmfinder.c */
	double min_acceptable_totweight; /* was #defined as MIN_WEIGHT in global.h */
	double convergence_threshold; /* was 'threshold' variable in 'main' function in cmfinder.c */
	int printedBestInsideScoreLessThanViterbi;
	double min_acceptable_seq_weight; /* new thing, especially because --fragmentary can lead to poorly scoring seqs */
	double min_seq_score_in_final_msa;
	int eliminateIdenticalSeqs;
	int eliminateIdenticalSubseqs;
	int eliminateSeqsWithManyDegen,maxDegenPerSeq,flankingNucsForCountingDegen;
	int numCpu;
	int bgScoreSize;
	ESL_SQCACHE *bgSeqCache;
	double bgScoreEvalue; /* <0 means disabled */
	double bgScoreScanBitThresh;
	int bgScoreNonFrag; /* prevent bgscore from allowing fragmentary cmsearches */
	int shiftDistalMispairsIntoTerminalLoops,minTermLoopNucsWithMovingIntoTerminalLoops;
	ESL_RANDOMNESS *randomness;
	
	int filterNonFrag,filterNonFragPad;

	double openPairCost,openPairCost_initial; /* was gapcost in old cmfinder */
	double zoop_gamma; /* was a static variable defined within the Compute_ZOOP_Weight function, initialized to 0.3 */
	double best_totscore,oldscore;
	double overlyHighNumSeqToTotalCandRatio; /* hardcoded to 5 in ScanCand (em.c) in old cmfinder */
	float cm_scoreThreshold; /* was cm_threshold in old cmfinder */
	float cm_scoreThreshold_increment; /* hardcoded to 10 in ScanCand (em.c) in old cmfinder */
	float max_cm_scoreThreshold; /* hardcoded to 30 in ScanCand (em.c) in old cmfinder */
	double minCandScore; /* hardcoded to 10 in ScanCand (em.c) in old cmfinder */
	double minCandScoreToBestScoreRatio; /* hardcoded to 2 in ScanCand (em.c) in old cmfinder */
        int minPredictedBasePairs;
	int DB_length; /* set to 100000 in old cmfinder.c, but there was a command-line flag to use a different value */
	char hmmStatsToCalibrate[CM_p7_NEVPARAM];
	ESL_MSA *best_msa;
	int putStartEndCoordsInHitIds; /* makes hits in output MSA like INPUTSEQID/START-END, which is default in infernal 1.1 and similar to Rfam, but is different from what CMfinder 0.3 did */

	double cmsearch_mxsize,cmsearch_smxsize; /* options for cmsearch */
	
	int use_evalue;
	double evalue;
	int cmcalibrate_tailn;

	int saveMsaAfterFirstMstep,saveInProgress;
	int dieOnUntestedCode;
	/* timers */
	ESL_STOPWATCH *timer_overall,*timer_partition_function,*timer_m_step,*timer_e_step,*timer_cmsearch_viterbi,*timer_cmsearch_inside,*timer_cmbuild,*timer_elim_iden_seqs,*timer_weighting,*timer_cmcalibrate;
} CmfinderVars;

void *MallocOrDie (size_t n);

/* the function 'save_and_load_CM_via_memory' is a trick to avoid problems that cmbuild calls configure_cm and cmsearch does as well, and
 * the function can't be called twice.  Rather than figure out how to make the cmbuild and cmsearch code compatible, I'd rather just
 * simulate saving then loading the CM (i.e. simulating running the cmbuild and then cmsearch UNIX commands), which'll definitely work.
 * By using RAM, I'd guess it won't take that much extra time to do this hack.
 */
CM_t * /*new CM copy*/ save_and_load_CM_via_memory_copying (CM_t *in_cm,CmfinderVars *cmfinderVars);
void save_and_load_CM_via_memory_freeing_old (CM_t **ref_cm,CmfinderVars *cmfinderVars);

int** IntAlloc2DArray(int acol);
double** DoubleAlloc2DArray(int acol);
void Free2DArray(void **p, int dim1);
int TriIndex(int i,int j);

int IsCanonicalDigitalBpAndNotDegen (int nuc1,int nuc2);
char *remove_gap(char* seq,int** ret_idx_map);
void Avg_bppr(ESL_MSA *text_msa,double **bp_pr,int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col,double ***ret_avg_bppr,double **ret_avg_any_bppr);  /* Pr(not pair)=1-avg_any_bppr */
void prxy_column_only(ESL_MSA *msa,double **bp_pr,CmfinderVars *cmfinderVars,int *seq_fragment_start_col,int *seq_fragment_end_col,double ***ret_prxy);
void prxy(ESL_MSA *text_msa,double **bp_pr,CmfinderVars *cmfinderVars,int *seq_fragment_start_col,int *seq_fragment_end_col,double ***ret_prxy);
void mixy(ESL_MSA *digital_msa,CmfinderVars *cmfinderVars,double ***ret_mxy);
void merge(double **mxy, double **pxy, int alen, double ***ret_xy);
double entropy(ESL_MSA *msa, double *nullModel, double tot_weight);
double summarize_mxy(ESL_MSA *msa, double *nullModel, double tot_weight,int *pt);
double weighted_base_pair (ESL_MSA* msa,int *pt,double **bp_pr);
double average_base_pair(ESL_MSA* msa,int *pt);
double average_score(ESL_MSA *msa, float tot_weight);
double average_seq_id(ESL_MSA *msa, float tot_weight,int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col);
double average_seq_len(ESL_MSA *msa, float tot_weight);
double average_energy(ESL_MSA  *text_msa, float tot_weight);
double average_GC (ESL_MSA *msa, float tot_weight);
int conserved_position (ESL_MSA *msa, float tot_weight,int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col);
void bad_base_pair (ESL_MSA* msa, int *pt, double **bp_pr, double* ret_conflict_bp, double* ret_del_bp);
void CalcFragmentStartEnds (ESL_MSA *digital_msa,int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col);
void SaveInProgressMsa(const char *final_file,int iteration,const char *stepDescription,ESL_MSA *msa);
void TrimDegenOnEndsOfMSA (ESL_MSA *msa,Cand **best_msa_cand);
int esl_msaweight_GSC_fragmentary(ESL_MSA *msa,int use_fragmentary,int *seq_fragment_start_col,int *seq_fragment_end_col);
ESL_MSA * ShiftDistalMispairsIntoTerminalLoops (ESL_MSA *msa,CmfinderVars *cmfinderVars,int *seq_fragment_start_col,int *seq_fragment_end_col);
int CalcQuasiNumSpecies (ESL_MSA *msa);

void MultiplyEMWeightsByWeightingStrategy (ESL_MSA *msa,WeightingStrategy ws,ESL_STOPWATCH *timer_weighting,int use_fragmentary);
CM_t * M_step (ESL_MSA *digital_msa, double **pxy,CmfinderVars *cmfinderVars,double openPairCost);
ESL_MSA * /*digital_msa */ E_step(CM_t *cm,ESL_SQCACHE *inputSeqCache,int *ncand,Cand **cand,double *ret_totscore,double ***ret_pxy,CmfinderVars *cmfinderVars);
void ScanCand(CM_t *cm,ESL_SQCACHE *inputSeqCache,int max_cand,Cand **cand, int *ncand,CmfinderVars *cmfinderVars);

/* cmbuild and functions in from_cmbuild.c are derived from src/cmbuild.c in infernal-1.1rc2
 * changes:
 * - input is an ESL_MSA struct or CM_t* and/or ESL_SQCACHE.
 */
/* passing commang line flags:
 * it's kind of a hassle to print up the string with optional args, but it's also a hassle to make a fake argc/argv with converting numbers to strings, 
 * and then remembering to free the memory later.  esl_optargs internally stores optional arguments as strings (not as ints or doubles), so
 * using my own functions to set everything still requires conversion to strings (although I can ask the esl_optargs code to store the memory).
 * I could also modify esl_optargs or remove infernal's dependence on it, but it seems easier to use everything as a black box. 
 */
#define COMMANDLINE_FLAGS_MAX 1024 /* should be plenty, and the true max is certainly finite (although I haven't calculated it) */
typedef struct _CmdlineBuilder {
	char *buf;
	int pos,maxpos;
} CmdlineBuilder;
void InitCmdlineBuilder (CmdlineBuilder *cmdline,char *buf);
void PrintfCmdlineBuilder (CmdlineBuilder *cmdline,const char *format, ...);
void cmbuild (const char *cmbuild_commandline_flags,char *alifile,ESL_MSA *msa,CM_t **ret_cm,char *hmmStatsToCalibrate);
void cmsearch(const char *cmsearch_commandline_flags,CM_t *cm,ESL_SQCACHE *dbSeqs,ESL_ALPHABET *abc,CM_TOPHITS **ret_th,int calc_evalues);
void cmcalibrate (const char *cmcalibrate_commandline_flags,CM_t *cm,int *expModesToCalibrate,ExpInfo_t ***ret_expA);

typedef struct {
	int force_inside_alg,need_alignments,disable_filters,use_fragmentary;
} ScanMode;
void cmsearch_wrapper (CM_t *cm,ESL_SQCACHE *dbSeqs,CM_TOPHITS **ret_th,CmfinderVars *cmfinderVars,double cm_scoreThreshold,ScanMode scanMode);
void cmcalibrate_viterbi (CM_t *cm,CmfinderVars *cmfinderVars,ExpInfo_t ***ret_expA);

void SaveMsa (const char *filename,ESL_MSA *msa);
ESL_MSA *MakeTextMsaClone(ESL_MSA *digital_msa);
void SaveSeqCacheAsFasta (const char *fileName,ESL_SQCACHE *seqs);
void SaveCM(const char *fileName,CM_t *cm);
void AllocSSconsOrRfIfNotThere(char **p,ESL_MSA *msa);
void AllocSSconsOrRf(char **p,ESL_MSA *msa);
typedef struct {
	int colAfterInserts; /* insert between colAfterInserts-1 and colAfterInserts */
	int numCols;
} ColumnsToAdd;
void AddColumnsToMsa (ESL_MSA **get_new_msa,ESL_MSA *msa,ColumnsToAdd *columnsToAdd,int numColumnsToAddEntries,int **get_oldToNewColumnMap,int *get_columnsAdded);
