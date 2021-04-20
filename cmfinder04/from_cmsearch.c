/* derived from cmsearch.c */
#include "esl_config.h"
#include "p7_config.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_getopts.h"
#include "esl_scorematrix.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_ssi.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

/* note: HAVE_MPI code has been removed */

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"

#include "infernal.h"
#include "cand.h"
#include "cmfinder.h"

#define CMSEARCH_MAX_RESIDUE_COUNT 100000 /* put we'll need to remove this from thread_loop, if I decide to use it */

typedef struct {
#ifdef HMMER_THREADS
	pthread_mutex_t  *takeSeqMutex;	/* coordinating who gets the next sequence */
	ESL_SQCACHE      *dbSeqs;
	volatile int *global_seqnum;
#endif /*HMMER_THREADS*/
	CM_PIPELINE      *pli;         /* work pipeline                           */
	CM_TOPHITS       *th;          /* top hit results                         */
	CM_t             *cm;          /* a covariance model                      */
	P7_BG            *bg;          /* null models                             */
	P7_OPROFILE      *om;          /* optimized query profile HMM             */
	P7_PROFILE       *gm;          /* generic   query profile HMM                       */
	P7_PROFILE       *Rgm;         /* generic   query profile HMM for 5' truncated hits */
	P7_PROFILE       *Lgm;         /* generic   query profile HMM for 3' truncated hits */
	P7_PROFILE       *Tgm;         /* generic   query profile HMM for 5' and 3' truncated hits */
	P7_MSVDATA       *msvdata;     /* MSV/SSV specific data structure */
	float            *p7_evparam;  /* [0..CM_p7_NEVPARAM] E-value parameters */
	float             smxsize;     /* max size (Mb) of allowable scan mx (only relevant if --nohmm or --max) */
	/* added for simpler CMfinder thread parallelization */
} WORKER_INFO;

#define REPOPTS     "-E,-T,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--incE,--incT,--cut_ga,--cut_nc,--cut_tc"
#define FMODEOPTS   "--FZ,--hmmonly,--rfam,--mid,--nohmm,--max"
#define TIMINGOPTS  "--timeF1,--timeF2,--timeF3,--timeF4,--timeF5,--timeF6"

/* ** Large sets of options are InCompatible With (ICW) --max, --nohmm,
* --mid, --rfam, --FZ, Previously (before these were commented out) I
* used this defines in the 'incompatible with' field of the
* esl_getopts definition, but they're too long and cause a error
* message buffer overflow, so I now check and enforce each
* incompatibility within process_commandline() below, and (perhaps
* confusingly) the 'incompatible with' field is empty for these
* options which are actually incompatible with a lot of other
* options. 
*
* #define ICWMAX   "--nohmm,--mid,--default,--rfam,--FZ,--noF1,--noF2,--noF3,--noF4,--noF6,--doF1b,--noF2b,--noF3b,--noF4b,--doF5b,--F1,--F1b,--F2,--F2b,--F3,--F3b,--F4,--F4b,--F5,--F6,--ftau,--fsums,--fqdb,--fbeta,--fnonbanded,--nocykenv,--cykenvx,--tau,--sums,--nonbanded,--rt1,--rt2,--rt3,--ns,--maxtau,--anytrunc"
* #define ICWNOHMM "--max,--mid,--default,--rfam,--FZ,--noF1,--noF2,--noF3,--noF4,--doF1b,--noF2b,--noF3b,--noF4b,--doF5b,--F1,--F1b,--F2,--F2b,--F3,--F3b,--F4,--F4b,--F5,--ftau,--fsums,--tau,--sums,--rt1,--rt2,--rt3,--ns,--maxtau,--anytrunc"
* #define ICWMID   "--max,--nohmm,--default,--rfam,--FZ,--noF1,--noF2,--noF3,--doF1b,--noF2b,--F1,--F1b,--F2,--F2b"
* #define ICWDF    "--max,--nohmm,--mid,--rfam,--FZ"
* #define ICWRFAM  "--max,--nohmm,--mid,--default,--FZ"
* #define ICWFZ    "--max,--nohmm,--mid,--default,--rfam"
*/

/* previously tested HAVE_MPI */
#define CPUOPTS     NULL
#define MPIOPTS     NULL

static ESL_OPTIONS options[] = {
	/* name           type      default  env  range     toggles   reqs   incomp            help                                                          docgroup*/
	{ "-h",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "show brief help on version and usage",                          1 },
	{ "-g",           eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  "--hmmonly",     "configure CM for glocal alignment [default: local]",            1 },
	{ "-Z",           eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL,  NULL,            "set search space size in *Mb* to <x> for E-value calculations", 1 },
	{ "--devhelp",    eslARG_NONE,   NULL,  NULL, NULL,    NULL,  NULL,  NULL,            "show list of otherwise hidden developer/expert options",        1 },
	/* Control of output */
	/* name           type         default   env  range toggles   reqs   incomp           help                                                            docgroup*/
	{ "-o",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "direct output to file <f>, not stdout",                        2 },
	{ "-A",           eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save multiple alignment of all significant hits to file <s>",  2 },
	{ "--tblout",     eslARG_OUTFILE, NULL, NULL, NULL,    NULL,  NULL,  NULL,            "save parseable table of hits to file <s>",                     2 },
	{ "--acc",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "prefer accessions over names in output",                       2 },
	{ "--noali",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "don't output alignments, so output is smaller",                2 },
	{ "--notextw",    eslARG_NONE,    NULL, NULL, NULL,    NULL,  NULL, "--textw",        "unlimit ASCII text output line width",                         2 },
	{ "--textw",      eslARG_INT,    "120", NULL, "n>=120",NULL,  NULL, "--notextw",      "set max width of ASCII text output lines",                     2 },
	{ "--verbose",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "report extra information; mainly useful for debugging",        2 },
	/* Control of reporting thresholds */
	/* name           type         default   env  range toggles   reqs   incomp           help                                                            docgroup*/
	{ "-E",           eslARG_REAL,  "10.0", NULL, "x>0",   NULL,  NULL,  REPOPTS,         "report sequences <= this E-value threshold in output",         3 },
	{ "-T",           eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  REPOPTS,         "report sequences >= this score threshold in output",           3 },
	/* Control of inclusion (significance) thresholds */
	/* name           type         default   env  range toggles   reqs   incomp           help                                                            docgroup*/
	{ "--incE",       eslARG_REAL,  "0.01", NULL, "x>0",   NULL,  NULL,  INCOPTS,         "consider sequences <= this E-value threshold as significant",  4 },
	{ "--incT",       eslARG_REAL,   FALSE, NULL, NULL,    NULL,  NULL,  INCOPTS,         "consider sequences >= this score threshold as significant",    4 },
	/* Model-specific thresholding for both reporting and inclusion */
	/* name           type         default   env  range toggles   reqs   incomp           help                                                            docgroup*/
	{ "--cut_ga",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use CM's GA gathering cutoffs as reporting thresholds",        5 },
	{ "--cut_nc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use CM's NC noise cutoffs as reporting thresholds",            5 },
	{ "--cut_tc",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  THRESHOPTS,      "use CM's TC trusted cutoffs as reporting thresholds",          5 },
	/* Control of filtering mode/acceleration level */
	/* name           type         default   env  range toggles   reqs   incomp                   help                                                              docgroup*/
	{ "--max",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL, /* see ** above */ "turn all heuristic filters off (slow)",                          6 },
	{ "--nohmm",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL, /* see ** above */ "skip all HMM filter stages, use only CM (slow)",                 6 },
	{ "--mid",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL, /* see ** above */ "skip first two HMM filter stages (SSV & Vit)",                   6 },
	{ "--default",    eslARG_NONE,"default",NULL, NULL,    NULL,  NULL,  NULL, /* see ** above */ "default: run search space size-dependent pipeline",              6 },
	{ "--rfam",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL, /* see ** above */ "set heuristic filters at Rfam-level (fast)",                     6 },
	{ "--hmmonly",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL, /* see ** above */ "use HMM only, don't use a CM at all",                            6 },
	{ "--FZ",         eslARG_REAL,    NULL, NULL, NULL,    NULL,  NULL,  NULL, /* see ** above */ "set filters to defaults used for a search space of size <x> Mb", 6 },
	{ "--Fmid",       eslARG_REAL,  "0.02", NULL, NULL,    NULL,"--mid", NULL,                    "with --mid, set P-value threshold for HMM stages to <x>",        6 },
	/* Other options */
	/* name           type         default   env  range toggles   reqs   incomp                help                                                            docgroup*/
	{ "--notrunc",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,                 "do not allow truncated hits at sequence termini",              7 },
	{ "--anytrunc",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,"-g,--notrunc",         "allow truncated hits anywhere within sequences",               7 },
	{ "--nonull3",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,                 "turn off the NULL3 post hoc additional null model",            7 },
	{ "--mxsize",     eslARG_REAL,  "128.", NULL, "x>0.1", NULL,  NULL,  NULL,                 "set max allowed size of alignment DP matrices to <x> Mb",      7 },
	{ "--smxsize",    eslARG_REAL,  "128.", NULL, "x>0.1", NULL,  NULL,  NULL,                 "set max allowed size of search DP matrices to <x> Mb",         7 },
	{ "--cyk",        eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,                 "use scanning CM CYK algorithm, not Inside in final stage",     7 },
	{ "--acyk",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,                 "align hits with CYK, not optimal accuracy",                    7 },
	{ "--wcx",        eslARG_REAL,   FALSE, NULL, "x>=1.25",NULL, NULL,"--nohmm,--qdb,--fqdb", "set W (expected max hit len) as <x> * cm->clen (model len)",   7 },
	{ "--toponly",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,                 "only search the top strand",                                   7 },
	{ "--bottomonly", eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,                 "only search the bottom strand",                                7 },
	{ "--tformat",    eslARG_STRING,  NULL, NULL, NULL,    NULL,  NULL,  NULL,                 "assert target <seqdb> is in format <s>: no autodetection",     7 },
	{ "--glist",      eslARG_INFILE,  NULL, NULL, NULL,    NULL,  NULL,  NULL,                 "BOGUS OPTION, NEVER ALLOWED",    999 },
#ifdef HMMER_THREADS 
	{ "--cpu",        eslARG_INT, NULL,"INFERNAL_NCPU","n>=0",NULL,  NULL,  CPUOPTS,      "number of parallel CPU workers to use for multithreads",       7 },
#endif

	/* All options below are developer options, only shown if --devhelp invoked */
	/* Options for precise control of each stage of the CM filter pipeline */
	/* name           type         default  env   range  toggles  reqs  incomp            help                                                      docgroup*/
	{ "--noF1",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--doF1b",        "skip the HMM SSV filter stage",                              101 },
	{ "--noF2",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--noF2b", NULL,          "skip the HMM Viterbi filter stage",                          101 },
	{ "--noF3",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--noF3b", NULL,          "skip the HMM Forward filter stage",                          101 },
	{ "--noF4",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,"--noF4b", NULL,          "skip the HMM glocal Forward filter stage",                   101 },
	{ "--noF6",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,             "skip the CM CYK filter stage",                               101 },
	{ "--doF1b",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,             "turn on  the HMM SSV composition bias filter",               101 },
	{ "--noF2b",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,             "turn off the HMM Vit composition bias filter",               101 },
	{ "--noF3b",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,             "turn off the HMM Fwd composition bias filter",               101 },
	{ "--noF4b",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,             "turn off the HMM glocal Fwd composition bias filter",        101 },
	{ "--doF5b",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,             "turn on  the HMM per-envelope composition bias filter",      101 },
	{ "--F1",         eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL, "--noF1",         "Stage 1 (SSV) threshold:         promote hits w/ P <= <x>",  101 },
	{ "--F1b",        eslARG_REAL,   FALSE, NULL, "x>0",   NULL,"--doF1b", NULL,          "Stage 1 (MSV) bias threshold:    promote hits w/ P <= <x>",  101 },
	{ "--F2",         eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL, "--noF2",         "Stage 2 (Vit) threshold:         promote hits w/ P <= <x>",  101 },
	{ "--F2b",        eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL, "--noF2b",        "Stage 2 (Vit) bias threshold:    promote hits w/ P <= <x>",  101 },
	{ "--F3",         eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL, "--noF3",         "Stage 3 (Fwd) threshold:         promote hits w/ P <= <x>",  101 },
	{ "--F3b",        eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL, "--noF3b",        "Stage 3 (Fwd) bias threshold:    promote hits w/ P <= <x>",  101 },
	{ "--F4",         eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL, "--noF4",         "Stage 4 (gFwd) glocal threshold: promote hits w/ P <= <x>",  101 },
	{ "--F4b",        eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL, "--noF4b",        "Stage 4 (gFwd) glocal bias thr:  promote hits w/ P <= <x>",  101 },
	{ "--F5",         eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL, NULL,             "Stage 5 (env defn) threshold:    promote hits w/ P <= <x>",  101 },
	{ "--F5b",        eslARG_REAL,   FALSE, NULL, "x>0",   NULL,"--doF5b", NULL,          "Stage 5 (env defn) bias thr:     promote hits w/ P <= <x>",  101 },
	{ "--F6",         eslARG_REAL,   FALSE, NULL, "x>0",   NULL,  NULL, "--noF6",         "Stage 6 (CYK) threshold:         promote hits w/ P <= <x>",  101 },
	/* Options for precise control of each stage of the HMM-only filter pipeline */
	/* name          type         default  env  range  toggles   reqs  incomp            help                                                         docgroup*/
	{ "--hmmmax",     eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--hmmF1,--hmmF2,--hmmF3,--hmmnobias", "in HMM-only mode, turn off all filters",  102 },
	{ "--hmmF1",      eslARG_REAL,  "0.02", NULL, "x>0",   NULL,  NULL, "--nohmmonly",    "in HMM-only mode, set stage 1 (SSV) P value threshold to <x>", 102 },
	{ "--hmmF2",      eslARG_REAL,  "1e-3", NULL, "x>0",   NULL,  NULL, "--nohmmonly",    "in HMM-only mode, set stage 2 (Vit) P value threshold to <x>", 102 },
	{ "--hmmF3",      eslARG_REAL,  "1e-5", NULL, "x>0",   NULL,  NULL, "--nohmmonly",    "in HMM-only mode, set stage 3 (Fwd) P value threshold to <x>", 102 },
	{ "--hmmnobias",  eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--nohmmonly",    "in HMM-only mode, turn off the bias composition filter",       102 },
	{ "--hmmnonull2", eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--nohmmonly",    "in HMM-only mode, turn off the null2 score correction",        102 },
	{ "--nohmmonly",  eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, "--hmmmax",       "never run HMM-only mode, not even for models with 0 basepairs",102 },
	/* Options for precise control of HMM envelope definition */
	/* name           type          default  env range toggles    reqs  incomp            help                                                      docgroup*/
	{ "--rt1",        eslARG_REAL,  "0.25", NULL, NULL,    NULL,  NULL, "--nohmm,--max",  "set domain/envelope definition rt1 parameter as <x>",        103 },
	{ "--rt2",        eslARG_REAL,  "0.10", NULL, NULL,    NULL,  NULL, "--nohmm,--max",  "set domain/envelope definition rt2 parameter as <x>",        103 },
	{ "--rt3",        eslARG_REAL,  "0.20", NULL, NULL,    NULL,  NULL, "--nohmm,--max",  "set domain/envelope definition rt3 parameter as <x>",        103 },
	{ "--ns",         eslARG_INT,   "200",  NULL, NULL,    NULL,  NULL, "--nohmm,--max",  "set number of domain/envelope tracebacks to <n>",            103 },
	/* Options for precise control of the CYK filter round of searching */
	/* name           type          default  env range      toggles     reqs  incomp            help                                                      docgroup*/
	{ "--ftau",       eslARG_REAL, "1e-4",  NULL, "1E-18<x<1", NULL,    NULL, "--fqdb",   "set HMM band tail loss prob for CYK filter to <x>",             104 },
	{ "--fsums",      eslARG_NONE,  FALSE,  NULL, NULL,        NULL,    NULL, "--fqdb",   "w/--fhbanded use posterior sums (widens bands)",                104 },
	{ "--fqdb",       eslARG_NONE,  FALSE,  NULL, NULL,        NULL,    NULL,   NULL,     "use QDBs in CYK filter round, not HMM bands",                   104 },
	{ "--fbeta",      eslARG_REAL, "1e-7",  NULL, "1E-18<x<1", NULL,    NULL,   NULL,     "set tail loss prob for CYK filter QDB calculation to <x>",      104 },
	{ "--fnonbanded", eslARG_NONE,  FALSE,  NULL, NULL,        NULL,    NULL,"--ftau,--fsums,--fqdb,--fbeta","do not use any bands for CYK filter round",  104 },
	{ "--nocykenv",   eslARG_NONE,  FALSE,  NULL, NULL,        NULL,    NULL, "--max",     "do not redefine envelopes after stage 6 based on CYK hits",    104 },
	{ "--cykenvx",    eslARG_INT,    "10",  NULL, "n>=1",      NULL,    NULL, "--max",     "CYK envelope redefinition threshold multiplier, <n> * F6",     104 },
	/* Options for precise control of the final round of searching */
	/* name           type          default  env range      toggles     reqs  incomp   help                                                      docgroup*/
	{ "--tau",        eslARG_REAL, "5e-6",  NULL, "1E-18<x<1", NULL,    NULL,"--qdb",  "set HMM band tail loss prob for final round to <x>",               105 },
	{ "--sums",       eslARG_NONE,  FALSE,  NULL, NULL,        NULL,    NULL,"--qdb",  "w/--hbanded use posterior sums (widens bands)",                    105 },
	{ "--qdb",        eslARG_NONE,  FALSE,  NULL, NULL,        NULL,    NULL,   NULL,  "use QDBs (instead of HMM bands) in final Inside round",            105 },
	{ "--beta",       eslARG_REAL,"1e-15",  NULL, "1E-18<x<1", NULL,    NULL,   NULL,  "set tail loss prob for final Inside QDB calculation to <x>",       105 },
	{ "--nonbanded",  eslARG_NONE,  FALSE,  NULL, NULL,        NULL,    NULL,"--tau,--sums,--qdb,--beta", "do not use QDBs or HMM bands in final Inside round of CM search", 105 },
	/* Options for timing individual pipeline stages */
	/* name          type         default  env  range  toggles   reqs  incomp            help                                                  docgroup*/
	{ "--timeF1",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,       "abort after Stage 1 SSV; for timing expts",          106 },
	{ "--timeF2",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,       "abort after Stage 2 Vit; for timing expts",          106 },
	{ "--timeF3",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,       "abort after Stage 3 Fwd; for timing expts",          106 },
	{ "--timeF4",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,       "abort after Stage 4 glocal Fwd; for timing expts",   106 },
	{ "--timeF5",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,       "abort after Stage 5 envelope def; for timing expts", 106 },
	{ "--timeF6",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, TIMINGOPTS,       "abort after Stage 6 CYK; for timing expts",          106 },
	/* Other expert options */
	/* name          type          default   env  range toggles   reqs  incomp            help                                                             docgroup*/
	{ "--nogreedy",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "do not resolve hits with greedy algorithm, use optimal one",    107 },
	{ "--cp9noel",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  "-g",            "turn off local ends in cp9 HMMs",                               107 },
	{ "--cp9gloc",    eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  "-g,--cp9noel",  "configure cp9 HMM in glocal mode",                              107 },
	{ "--null2",      eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,            "turn on null 2 biased composition HMM score corrections",       107 },
	{ "--maxtau",     eslARG_REAL,  "0.05", NULL,"0<x<0.5",NULL,  NULL,  NULL,            "set max tau <x> when tightening HMM bands",                     107 },
	{ "--seed",       eslARG_INT,    "181", NULL, "n>=0",  NULL,  NULL,  NULL,            "set RNG seed to <n> (if 0: one-time arbitrary seed)",           107 },

	/* flags added for CMfinder */
	{ "--noF5",       eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL,  NULL,          "skip the HMM glocal envelope definition",                   101 },
	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* struct cfg_s : "Global" application configuration shared by all threads/processes
* 
* This structure is passed to routines within main.c, as a means of semi-encapsulation
* of shared data amongst different parallel processes (threads or MPI processes).
*/
struct cfg_s {
	char            *dbfile;           /* target sequence database file                   */
	char            *cmfile;           /* query HMM file                                  */
	int64_t          Z;                /* database size, in nucleotides                   */
	enum cm_zsetby_e Z_setby;          /* how Z was set: CM_ZSETBY_SSIINFO, CM_ZSETBY_OPTION, CM_ZSETBY_FILEINFO */
	int              do_mpi;           /* TRUE if we're doing MPI parallelization         */
	int              nproc;            /* how many MPI processes, total                   */
	int              my_rank;          /* who am I, in 0..nproc-1                         */
};

static char banner[] = "search CM(s) against a sequence database";

static int serial_master(ESL_GETOPTS *go, struct cfg_s *cfg,CM_t *cm,ESL_SQCACHE *dbSeqs,ESL_ALPHABET *abc,CM_TOPHITS **ret_th,int calc_evalues);
static int serial_loop(WORKER_INFO *info, ESL_SQCACHE *dbSeqs);

#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000
int  thread_loop(WORKER_INFO *info, ESL_THREADS *obj, ESL_SQCACHE *dbSeqs);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/

static void process_commandline(const char *cmsearch_commandline_flags, ESL_GETOPTS **ret_go);
static int  output_header(FILE *ofp, const ESL_GETOPTS *go, char *cmfile, char *seqfile, int ncpus);

/* Functions to avoid code duplication for common tasks */
int          open_dbfile(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_SQFILE **ret_dbfp);
int          dbsize_and_seq_lengths(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_SQFILE **dbfp_ptr, char *errbuf, int64_t **ret_srcL, int64_t *ret_nseqs);
static WORKER_INFO *create_info(const ESL_GETOPTS *go);
static int          clone_info(ESL_GETOPTS *go, WORKER_INFO *src_info, WORKER_INFO *dest_infoA, int dest_infocnt, char *errbuf);
static void         free_info(WORKER_INFO *info);
static int          configure_cm(WORKER_INFO *info);
static int          setup_hmm_filter(ESL_GETOPTS *go, WORKER_INFO *info);

void 
cmsearch(const char *cmsearch_commandline_flags,CM_t *cm,ESL_SQCACHE *dbSeqs,ESL_ALPHABET *abc,CM_TOPHITS **ret_th,int calc_evalues)
{
	int              status   = eslOK;
	int added_expA_to_cm=0;

	ESL_GETOPTS     *go  = NULL;    /* command line processing                 */
	struct cfg_s     cfg;         /* configuration data                      */
	
	
	/* MAJOR HACK.  Pretend that the E-value statistics have been calculated,
	 * otherwise code in cm_pipeline.c complains.  It's CMfinder's responsibility
	 * to make sure that the exponents aren't actually needed.
	 * However, I commented out the line to initialize the expA data, so that
	 * valgrind will detect errors.
	 */
	if ((cm->flags & CMH_EXPTAIL_STATS)==0) {
		int i;
		cm->flags |= CMH_EXPTAIL_STATS;
		added_expA_to_cm=1;
		cm->expA=(ExpInfo_t **)(MallocOrDie(sizeof(ExpInfo_t *) * EXP_NMODES));
#if 1 /* don't initialize, so valgrind detects errors.  Note: we can't just set it to NULL, since then the cm_Clone command (called within cmsearch) will deference NULL */
		for(i = 0; i < EXP_NMODES; i++) cm->expA[i] = (ExpInfo_t *)MallocOrDie(sizeof(ExpInfo_t));
#else
		for(i = 0; i < EXP_NMODES; i++) cm->expA[i] = CreateExpInfo();
#endif
	}

	/* Initialize what we can in the config structure (without knowing the alphabet yet)
	*/
	cfg.cmfile     = NULL;
	cfg.dbfile     = NULL;

	cfg.do_mpi     = FALSE;               /* this gets reset below, if we init MPI */
	cfg.nproc      = 0;                   /* this gets reset below, if we init MPI */
	cfg.my_rank    = 0;                   /* this gets reset below, if we init MPI */

	process_commandline(cmsearch_commandline_flags, &go);

	status = serial_master(go, &cfg,cm,dbSeqs,abc,ret_th,calc_evalues);

	esl_getopts_Destroy(go);
	
	if (status!=eslOK) {
		esl_fatal("cmsearch function had problem");
	}
	
	/* put CM back into original state with respect to CMH_EXPTAIL_STATS */
	if (added_expA_to_cm) {
		int i;
		cm->flags &= (~CMH_EXPTAIL_STATS);
		for(i = 0; i < EXP_NMODES; i++) free(cm->expA[i]);
		free(cm->expA);
		cm->expA=NULL;
	}
}

/* serial_master()
* The serial version of cmsearch.
* For each query CM in <cmfile> search the database for hits.
* 
* A master can only return if it's successful. All errors are handled 
* immediately and fatally with cm_Fail().
*/
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg,CM_t *cm,ESL_SQCACHE *dbSeqs,ESL_ALPHABET *abc,CM_TOPHITS **ret_th,int calc_evalues)
{
	FILE            *ofp      = stdout;            /* results output file (-o)                        */
	FILE            *afp      = NULL;              /* alignment output file (-A)                      */
	FILE            *tblfp    = NULL;              /* output stream for tabular hits (--tblout)       */
	ESL_STOPWATCH   *w        = NULL;              /* timing one query model                          */
	ESL_STOPWATCH   *mw       = NULL;              /* timing all query models                         */
	int              textw    = 0;
	int              status   = eslOK;
	int              sstatus  = eslOK;
	int              i;
	int              cm_idx;
	double           eZ;                           /* effective database size */

	int              ncpus    = 0;

	WORKER_INFO     *tinfo;                        /* the template info, cloned to make the worked info */
	WORKER_INFO     *info          = NULL;         /* the worker info */
	int              infocnt       = 0;            /* number of worker infos */

	int              nbps;                         /* number of basepairs in current CM */

#ifdef HMMER_THREADS
	ESL_SQ_BLOCK    *block    = NULL;
	ESL_THREADS     *threadObj= NULL;
	ESL_WORK_QUEUE  *queue    = NULL;
#endif
	char             errbuf[eslERRBUFSIZE];

	w  = esl_stopwatch_Create();
	mw = esl_stopwatch_Create();
	esl_stopwatch_Start(mw);

	if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
	else                                     textw = esl_opt_GetInteger(go, "--textw");

	/* Open the results output files */
	if (esl_opt_IsOn(go, "-o"))           { if ((ofp       = fopen(esl_opt_GetString(go, "-o"),          "w")) == NULL) cm_Fail("Failed to open output file %s for writing\n",         esl_opt_GetString(go, "-o")); }
	if (esl_opt_IsOn(go, "-A"))           { if ((afp       = fopen(esl_opt_GetString(go, "-A"),          "w")) == NULL) cm_Fail("Failed to open alignment file %s for writing\n", esl_opt_GetString(go, "-A")); }
	if (esl_opt_IsOn(go, "--tblout"))     { if ((tblfp     = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL) cm_Fail("Failed to open tabular output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }

#ifdef HMMER_THREADS
	/* initialize thread data */
	if (esl_opt_IsOn(go, "--cpu")) ncpus = esl_opt_GetInteger(go, "--cpu");
	else                                   esl_threads_CPUCount(&ncpus);
	if (ncpus > 0) {
		threadObj = esl_threads_Create(&pipeline_thread);
		queue = esl_workqueue_Create(ncpus * 2);
	}
#endif

	infocnt = (ncpus == 0) ? 1 : ncpus;
	ESL_ALLOC(info, sizeof(WORKER_INFO) * infocnt);

	/* One-time initializations after alphabet <abc> becomes known */
	output_header(ofp, go, cfg->cmfile, cfg->dbfile, ncpus);

	if(esl_opt_IsUsed(go, "-Z")) { /* -Z enabled, use that size */
		cfg->Z       = (int64_t) (esl_opt_GetReal(go, "-Z") * 1000000.); 
		cfg->Z_setby = CM_ZSETBY_OPTION; 
	}

	for (i = 0; i < infocnt; ++i)    {
		info[i].pli          = NULL;
		info[i].th           = NULL;
		info[i].cm           = NULL;
		info[i].om           = NULL;
		info[i].bg           = NULL;
		ESL_ALLOC(info[i].p7_evparam, sizeof(float) * CM_p7_NEVPARAM);
	}

#ifdef HMMER_THREADS    
	for (i = 0; i < ncpus * 2; ++i) {
		block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abc);
		if (block == NULL)           esl_fatal("Failed to allocate sequence block");

		status = esl_workqueue_Init(queue, block);
		if (status != eslOK)         esl_fatal("Failed to add block to work queue");
	}
#endif

	/* In original cmsearch code: 
	 * "Outer loop: over each query CM in <cmfile>.", 
	 * but now there's only 1 CM */
	cm_idx=1;
	esl_stopwatch_Start(w);

	/* create a new template info, and point it to the cm we just read */
	tinfo = create_info(go);
	tinfo->cm = cm;

	/* sanity check: we better have a filter HMM */
	if(! (tinfo->cm->flags & CMH_FP7)) cm_Fail("no filter HMM was read for CM: %s\n", tinfo->cm->name);

	/* NOTE: disabled check in original cmsearch.c to verify that the CM has E-values, 
	 * since in CMfinder, we can't afford cmcalibrate on every iteration
	 */
	 
	nbps = CMCountNodetype(tinfo->cm, MATP_nd);

	fprintf(ofp, "Query:       %s  [CLEN=%d]\n", tinfo->cm->name, tinfo->cm->clen);
	if (tinfo->cm->acc)  fprintf(ofp, "Accession:   %s\n", tinfo->cm->acc);
	if (tinfo->cm->desc) fprintf(ofp, "Description: %s\n", tinfo->cm->desc);

	/* configure the CM (this builds QDBs if nec) and setup HMM filters 
	* (we need to do this before clone_info()). We need a pipeline to 
	* do this only b/c we need pli->cm_config_opts.
	*/
	tinfo->pli = cm_pipeline_Create(go, abc, tinfo->cm->clen, 100, cfg->Z, cfg->Z_setby, CM_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
	if((status = configure_cm(tinfo))         != eslOK) cm_Fail(tinfo->pli->errbuf);
	if((status = setup_hmm_filter(go, tinfo)) != eslOK) cm_Fail(tinfo->pli->errbuf);
	/* Clone all data structures in tinfo into the WORKER_INFO's in info */
	if((status = clone_info(go, tinfo, info, infocnt, errbuf)) != eslOK) cm_Fail(errbuf);

	/* Create processing pipeline and hit list */
#ifdef HMMER_THREADS
	volatile int global_seqnum=0;
	pthread_mutex_t  takeSeqMutex;
	if (pthread_mutex_init(&takeSeqMutex, NULL) != 0) esl_fatal("mutex init failed");
#endif
	for (i = 0; i < infocnt; ++i) {
		info[i].th   = cm_tophits_Create();
		info[i].pli  = cm_pipeline_Create(go, abc, tinfo->cm->clen, 100, cfg->Z, cfg->Z_setby, CM_SEARCH_SEQS); /* L_hint = 100 is just a dummy for now */
		if((status = cm_pli_NewModel(info[i].pli, CM_NEWMODEL_CM, info[i].cm, info[i].cm->clen, info[i].cm->W, nbps,
			info[i].om, info[i].bg, info[i].p7_evparam, info[i].om->max_length, cm_idx-1, NULL)) != eslOK) { 
				cm_Fail(info[i].pli->errbuf);
		}

#ifdef HMMER_THREADS
		info[i].global_seqnum=&global_seqnum;
		info[i].takeSeqMutex=&takeSeqMutex;
		info[i].dbSeqs=dbSeqs;
		if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
	}

#ifdef HMMER_THREADS
	if (ncpus > 0)  {
		sstatus = thread_loop(info, threadObj, dbSeqs);
	}
	else {
		sstatus = serial_loop(info, dbSeqs);
	}
	pthread_mutex_destroy(&takeSeqMutex);
#else
	sstatus = serial_loop(info, dbSeqs);
#endif
	switch(sstatus) {
	case eslEFORMAT:
		esl_fatal("Parse failed in serial_loop");
		break;
	case eslEOF:
		/* do nothing */
		break;
	default:
		esl_fatal("Unexpected error in serial_loop");
	}

	if (calc_evalues) {
		/* we need to re-compute e-values before merging (when list will be sorted) */
		for (i = 0; i < infocnt; ++i) { 
			if(info[i].pli->do_hmmonly_cur) eZ = info[i].pli->Z / (float) info[i].om->max_length;
			else                 	      eZ = info[i].cm->expA[info[i].pli->final_cm_exp_mode]->cur_eff_dbsize;
			cm_tophits_ComputeEvalues(info[i].th, eZ, 0);
		}
	}

	/* merge the search results */
	for (i = 1; i < infocnt; ++i) {
		cm_tophits_Merge(info[0].th,   info[i].th);
		cm_pipeline_Merge(info[0].pli, info[i].pli);
		free_info(&info[i]);
	}

	/* Sort by sequence index/position and remove duplicates */
	cm_tophits_SortForOverlapRemoval(info[0].th);
	if((status = cm_tophits_RemoveOverlaps(info[0].th, errbuf)) != eslOK) cm_Fail(errbuf);

	/* Resort by score and enforce threshold */
	cm_tophits_Threshold(info[0].th, info[0].pli);

	/* tally up total number of hits and target coverage */
	for (i = 0; i < info->th->N; i++) {
		if ((info[0].th->hit[i]->flags & CM_HIT_IS_REPORTED) || (info[0].th->hit[i]->flags & CM_HIT_IS_INCLUDED)) { 
			info[0].pli->acct[info[0].th->hit[i]->pass_idx].n_output++;
			info[0].pli->acct[info[0].th->hit[i]->pass_idx].pos_output += abs(info[0].th->hit[i]->stop - info[0].th->hit[i]->start) + 1;
		}
	}

	/* Output */
	cm_tophits_Targets(ofp, info[0].th, info[0].pli, textw);
	fprintf(ofp, "\n\n");

	if(info[0].pli->show_alignments) {
		if((status = cm_tophits_HitAlignments(ofp, info[0].th, info[0].pli, textw)) != eslOK) cm_Fail("Out of memory");
		fprintf(ofp, "\n\n");
		if(info[0].pli->be_verbose) { 
			cm_tophits_HitAlignmentStatistics(ofp, info[0].th, 
				(info[0].pli->cm_align_opts & CM_ALIGN_HBANDED), 
				(info[0].pli->cm_align_opts & CM_ALIGN_CYK),
				info[0].pli->final_tau);
			fprintf(ofp, "\n\n");
		}
	}

	if (tblfp != NULL) { 
		cm_tophits_TabularTargets(tblfp, info[0].cm->name, info[0].cm->acc, info[0].th, info[0].pli, (cm_idx == 1)); 
		fflush(tblfp);
	}
	esl_stopwatch_Stop(w);
	cm_pli_Statistics(ofp, info[0].pli, w);

	/* Output the results in an MSA (-A option) */
	if (afp) {
		ESL_MSA *msa = NULL;
		if((status = cm_tophits_Alignment(info[0].cm, info[0].th, errbuf, &msa)) == eslOK) { 
			if(msa != NULL) { 
				if (textw > 0) eslx_msafile_Write(afp, msa, eslMSAFILE_STOCKHOLM);
				else           eslx_msafile_Write(afp, msa, eslMSAFILE_PFAM);
				fprintf(ofp, "# Alignment of %d hits satisfying inclusion thresholds saved to: %s\n", msa->nseq, esl_opt_GetString(go, "-A"));
			}
			else { 
				fprintf(ofp, "# No hits satisfy inclusion thresholds; no alignment saved\n");
			}
			esl_msa_Destroy(msa);
		}
		else { /* status != eslOK */
			cm_Fail(errbuf);
		}
	}
	fprintf(ofp, "//\n");
	
	/* steal the tophits data structure for the caller to look at */
	*ret_th=info[0].th;
	info[0].th=NULL;

	tinfo->cm=NULL; /* leave this input CM_t*cm to the caller */
	free_info(tinfo);
	free(tinfo);
	free_info(&(info[0]));

#ifdef HMMER_THREADS
	if (ncpus > 0) {
		esl_workqueue_Reset(queue);
		while (esl_workqueue_Remove(queue, (void **) &block) == eslOK) {
			esl_sq_DestroyBlock(block);
		}
		esl_workqueue_Destroy(queue);
		esl_threads_Destroy(threadObj);
	}
#endif

	/* Terminate outputs... any last words?
	*/
	if (tblfp)         cm_tophits_TabularTail(tblfp,    "cmsearch", CM_SEARCH_SEQS, cfg->cmfile, cfg->dbfile, go);
	fprintf(ofp, "[ok]\n");

	esl_stopwatch_Stop(mw);
	if(esl_opt_GetBoolean(go, "--verbose")) esl_stopwatch_Display(stdout, mw, "Total runtime:");

	free(info);

	if(w  != NULL) esl_stopwatch_Destroy(w);
	if(mw != NULL) esl_stopwatch_Destroy(mw);

	if (ofp != stdout) fclose(ofp);
	if (afp)           fclose(afp);
	if (tblfp)         fclose(tblfp);

	return eslOK;

ERROR:
	return eslFAIL;
}


/* serial_loop(): 
* 
* Note: unlike the original cmsearch.c code,
* - we read from dbSeqs instead of from a file
* - we don't bother with reading CMSEARCH_MAX_RESIDUE_COUNT at a time; this seems to be a cmsearch.c-specific thing, so if you can read a seq into RAM, then you can scan it
* - I deleted the rev-comp code, since we'll never use it
* - turns out we don't need *srcL (it seems the code only needs to know the lengths of seqs in advance for something to do with reading blocks of long seqs in windows)
*/
static int
serial_loop(WORKER_INFO *info, ESL_SQCACHE *dbSeqs)
{
	int       status;
	int seqnum=0;
	int64_t   seq_idx = 0;
	ESL_SQ   *dbsq;

	for (seqnum=0; seqnum<dbSeqs->seq_count; seqnum++) {
		seq_idx=seqnum+1;
		
		dbsq=&(dbSeqs->sq_list[seqnum]);
		if (dbsq->L==0) {
			esl_fatal("sorry, altered cmsearch code does not allow zero-length sequences in input.  please remove it in input, or add code to cmfinder.c to remove zero-len sequences in inputSeqCache upon loading");
		}

		cm_pli_NewSeq(info->pli, dbsq, seq_idx-1);

		if (info->pli->do_top) { 
			if((status = cm_Pipeline(info->pli, info->cm->offset, info->om, info->bg, info->p7_evparam, info->msvdata, dbsq, info->th, FALSE, /* FALSE: not in reverse complement */
				NULL, &(info->gm), &(info->Rgm), &(info->Lgm), &(info->Tgm), &(info->cm))) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n%s\n", status, info->pli->errbuf);
			cm_pipeline_Reuse(info->pli); /* prepare for next search */
		}
	}

	return eslEOF;
}

#ifdef HMMER_THREADS

/*  The original cmsearch code won't work here without modification, because we're not reading from a file.
 * For simplicity, I've just re-done the serial_loop, but accessing sequences with a condition variable.
 */
int
thread_loop(WORKER_INFO *info, ESL_THREADS *obj, ESL_SQCACHE *dbSeqs)
{
	esl_threads_WaitForStart(obj);

	esl_threads_WaitForFinish(obj);

	return eslEOF;
}

/* pipeline_thread()
* 
* threads use mutexes and choose the next job -- no need for a master dishing out work
* code adapted from serial_loop
*/

static void 
pipeline_thread(void *arg)
{
	int       status;
	int64_t   seq_idx = 0;
	ESL_SQ   *dbsq;
	int workeridx;
	WORKER_INFO   *info;
	ESL_THREADS   *obj;
	ESL_SQCACHE *dbSeqs;

#ifdef HAVE_FLUSH_ZERO_MODE
	/* In order to avoid the performance penalty dealing with sub-normal
	* values in the floating point calculations, set the processor flag
	* so sub-normals are "flushed" immediately to zero.
	* On OS X, need to reset this flag for each thread
	* (see TW notes 05/08/10 for details)
	*/
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

	obj = (ESL_THREADS *) arg;
	esl_threads_Started(obj, &workeridx);
	info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);
	dbSeqs=info->dbSeqs;

	while (1) {
		int done;
		int seqnum;

		if (pthread_mutex_lock (info->takeSeqMutex) != 0) esl_fatal("mutex lock failed");
		if (*(info->global_seqnum)==dbSeqs->seq_count) {
			seqnum=-1;
			done=1;
		}
		else {
			done=0;
			seqnum=*(info->global_seqnum);
			(*(info->global_seqnum))++;
		}
		if (pthread_mutex_unlock  (info->takeSeqMutex) != 0) esl_fatal("mutex unlock failed");
		if (done) {
			break;
		}
		
		seq_idx=seqnum+1;
		
		dbsq=&(dbSeqs->sq_list[seqnum]);
		if (dbsq->L==0) {
			esl_fatal("sorry, altered cmsearch code does not allow zero-length sequences in input.  please remove it in input, or add code to cmfinder.c to remove zero-len sequences in inputSeqCache upon loading");
		}

		cm_pli_NewSeq(info->pli, dbsq, seq_idx-1);

		if (info->pli->do_top) { 
			if((status = cm_Pipeline(info->pli, info->cm->offset, info->om, info->bg, info->p7_evparam, info->msvdata, dbsq, info->th, FALSE, /* FALSE: not in reverse complement */
				NULL, &(info->gm), &(info->Rgm), &(info->Lgm), &(info->Tgm), &(info->cm))) != eslOK) cm_Fail("cm_pipeline() failed unexpected with status code %d\n%s\n", status, info->pli->errbuf);
			cm_pipeline_Reuse(info->pli); /* prepare for next search */
		}
	}

	esl_threads_Finished(obj, workeridx);
	
	return;
}
#endif   /* HMMER_THREADS */

/* process_commandline()
* 
* Processes the commandline, filling in fields in <cfg> and creating and returning
* an <ESL_GETOPTS> options structure. The help page (cmsearch -h) is formatted
* here.
*/
static void
process_commandline(const char *cmsearch_commandline_flags, ESL_GETOPTS **ret_go)
{
	ESL_GETOPTS *go     = NULL;
	int          do_dev = FALSE; /* set to TRUE if --devhelp used */

	if ((go = esl_getopts_Create(options))     == NULL)     cm_Fail("Internal failure creating options object");
	/*if (esl_opt_ProcessEnvironment(go)         != eslOK)  { printf("Failed to process environment: %s\n", go->errbuf); goto ERROR; }*/
	/*if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf);  goto ERROR; }*/
	if (esl_opt_ProcessSpoof(go,cmsearch_commandline_flags) != eslOK) { esl_fatal("Failed to parse command line: %s\n", go->errbuf); }
	if (esl_opt_VerifyConfig(go)               != eslOK)  { printf("Failed to parse command line: %s\n", go->errbuf);  goto ERROR; }

	/* help format: */
	do_dev = esl_opt_GetBoolean(go, "--devhelp") ? TRUE : FALSE;
	if (esl_opt_GetBoolean(go, "-h") || do_dev) {
		esl_fatal("umm, since this code is embedded within CMfinder, I don't support help");
	}

	if (esl_opt_ArgNumber(go)                  != 0)     { printf("Incorrect number of command line arguments -- you should give CM as a CM_t* and dbfile as ESL_SEQCACHE*.  optargs=%d",esl_opt_ArgNumber(go));     goto ERROR; }

	/* Check for incompatible option combinations I don't know how to disallow with esl_getopts */

	/* --glist never works, it only exists because cmscan has it, and consequently cm_pipeline_Create() expects it to */
	if (esl_opt_IsOn(go, "--glist")) { 
		puts("Failed to parse command line: --glist is an invalid cmsearch option (it only works with cmscan)");
		goto ERROR;
	}    

	/* --beta only makes sense with --qdb, --nohmm or --max */
	if (esl_opt_IsUsed(go, "--beta") && (! esl_opt_GetBoolean(go, "--qdb")) && 
		(! esl_opt_GetBoolean(go, "--nohmm")) && (! esl_opt_GetBoolean(go, "--max"))) { 
			puts("Failed to parse command line: --beta only makes sense in combination with --qdb, --nohmm or --max");
			goto ERROR;
	}    

	/* --fbeta only makes sense with --fqdb or --nohmm */
	if (esl_opt_IsUsed(go, "--fbeta") && (! esl_opt_GetBoolean(go, "--fqdb")) && (! esl_opt_GetBoolean(go, "--nohmm"))) { 
		puts("Failed to parse command line: --fbeta only makes sense in combination with --fqdb or --nohmm");
		goto ERROR;
	}    

	/* If we'll be using QDBs for both the CYK filter and final round,
	* make sure that beta (--beta <x>) is <= filter round beta (--fbeta <x>) 
	* There's two ways we'll need both sets of QDBs: 
	* 1. --nohmm and neither of --fnonbanded --nonbanded (1st half of ugly if below)
	* 2. --qdb and --fqdb (2nd half of ugly if below)
	*/
	if ((esl_opt_GetBoolean(go, "--nohmm") && (! esl_opt_GetBoolean(go, "--fnonbanded")) && (! esl_opt_GetBoolean(go, "--nonbanded"))) || 
		(esl_opt_GetBoolean(go, "--qdb") && esl_opt_GetBoolean(go, "--fqdb"))) {     
			if(esl_opt_IsUsed(go, "--beta") && esl_opt_IsUsed(go, "--fbeta")) { 
				if((esl_opt_GetReal(go, "--beta") - esl_opt_GetReal(go, "--fbeta")) > 1E-20) { 
					puts("Failed to parse command line: with --nohmm --fbeta <x1> --beta <x2>, <x1> must be >= <x2>\n");
					goto ERROR;
				}
			}
			else if(esl_opt_IsUsed(go, "--beta")) { 
				if((esl_opt_GetReal(go, "--beta") - esl_opt_GetReal(go, "--fbeta")) > 1E-20) { 
					printf("Failed to parse command line: with --nohmm --beta <x> (not in combination with --fbeta), <x> must be <= %g\n", esl_opt_GetReal(go, "--fbeta"));
					goto ERROR;
				}
			}
			else if(esl_opt_IsUsed(go, "--fbeta")) { 
				if((esl_opt_GetReal(go, "--beta") - esl_opt_GetReal(go, "--fbeta")) > 1E-20) { 
					printf("Failed to parse command line: with --nohmm --fbeta <x> (not in combination with --beta), <x> must be >= %g\n", esl_opt_GetReal(go, "--beta"));
					goto ERROR;
				}
			}
	}

	/* Finally, check for incompatible option combinations I *do* know
	* how to disallow with esl_getopts, but that would require an error
	* message like: "Option 'x' is incompatible with options
	* y1,y2,y3,y4....yn", where there's so many y's that the message is
	* truncated because errbuf runs out of space. As a workaround we
	* laboriously check for all incompatible options of that type here.
	*/
	if(esl_opt_IsUsed(go, "--max")) { 
		if(esl_opt_IsUsed(go, "--nohmm"))      { puts("Failed to parse command line: Option --max is incompatible with option --nohmm");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--mid"))        { puts("Failed to parse command line: Option --max is incompatible with option --mid");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--rfam"))       { puts("Failed to parse command line: Option --max is incompatible with option --rfam");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--FZ"))         { puts("Failed to parse command line: Option --max is incompatible with option --FZ");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF1"))       { puts("Failed to parse command line: Option --max is incompatible with option --noF1");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF2"))       { puts("Failed to parse command line: Option --max is incompatible with option --noF2");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF3"))       { puts("Failed to parse command line: Option --max is incompatible with option --noF3");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF4"))       { puts("Failed to parse command line: Option --max is incompatible with option --noF4");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF6"))       { puts("Failed to parse command line: Option --max is incompatible with option --noF6");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--doF1b"))      { puts("Failed to parse command line: Option --max is incompatible with option --doF1b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF2b"))      { puts("Failed to parse command line: Option --max is incompatible with option --noF2b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF3b"))      { puts("Failed to parse command line: Option --max is incompatible with option --noF3b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF4b"))      { puts("Failed to parse command line: Option --max is incompatible with option --noF4b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--doF5b"))      { puts("Failed to parse command line: Option --max is incompatible with option --doF5b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--F1"))         { puts("Failed to parse command line: Option --max is incompatible with option --F1");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F1b"))        { puts("Failed to parse command line: Option --max is incompatible with option --F1b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F2"))         { puts("Failed to parse command line: Option --max is incompatible with option --F2");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F2b"))        { puts("Failed to parse command line: Option --max is incompatible with option --F2b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F3"))         { puts("Failed to parse command line: Option --max is incompatible with option --F3");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F3b"))        { puts("Failed to parse command line: Option --max is incompatible with option --F3b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F4"))         { puts("Failed to parse command line: Option --max is incompatible with option --F4");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F4b"))        { puts("Failed to parse command line: Option --max is incompatible with option --F4b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F5"))         { puts("Failed to parse command line: Option --max is incompatible with option --F5");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F6"))         { puts("Failed to parse command line: Option --max is incompatible with option --F6");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--ftau"))       { puts("Failed to parse command line: Option --max is incompatible with option --ftau");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--fsums"))      { puts("Failed to parse command line: Option --max is incompatible with option --fsums");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--fqdb"))       { puts("Failed to parse command line: Option --max is incompatible with option --fqdb");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--fbeta"))      { puts("Failed to parse command line: Option --max is incompatible with option --fbeta");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--fnonbanded")) { puts("Failed to parse command line: Option --max is incompatible with option --fnonbanded"); goto ERROR; }
		if(esl_opt_IsUsed(go, "--nocykenv"))   { puts("Failed to parse command line: Option --max is incompatible with option --nocykenv");   goto ERROR; }
		if(esl_opt_IsUsed(go, "--cykenvx"))    { puts("Failed to parse command line: Option --max is incompatible with option --cykenvx");    goto ERROR; }
		if(esl_opt_IsUsed(go, "--tau"))        { puts("Failed to parse command line: Option --max is incompatible with option --tau");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--sums"))       { puts("Failed to parse command line: Option --max is incompatible with option --sums");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--nonbanded"))  { puts("Failed to parse command line: Option --max is incompatible with option --nonbanded");  goto ERROR; }
		if(esl_opt_IsUsed(go, "--rt1"))        { puts("Failed to parse command line: Option --max is incompatible with option --rt1");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--rt2"))        { puts("Failed to parse command line: Option --max is incompatible with option --rt2");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--rt3"))        { puts("Failed to parse command line: Option --max is incompatible with option --rt3");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--ns"))         { puts("Failed to parse command line: Option --max is incompatible with option --ns");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--maxtau"))     { puts("Failed to parse command line: Option --max is incompatible with option --maxtau");     goto ERROR; }
		if(esl_opt_IsUsed(go, "--anytrunc"))   { puts("Failed to parse command line: Option --max is incompatible with option --anytrunc");   goto ERROR; }
	}
	if(esl_opt_IsUsed(go, "--nohmm")) { 
		if(esl_opt_IsUsed(go, "--max"))        { puts("Failed to parse command line: Option --nohmm is incompatible with option --max");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--mid"))        { puts("Failed to parse command line: Option --nohmm is incompatible with option --mid");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--rfam"))       { puts("Failed to parse command line: Option --nohmm is incompatible with option --rfam");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--FZ"))         { puts("Failed to parse command line: Option --nohmm is incompatible with option --FZ");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF1"))       { puts("Failed to parse command line: Option --nohmm is incompatible with option --noF1");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF2"))       { puts("Failed to parse command line: Option --nohmm is incompatible with option --noF2");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF3"))       { puts("Failed to parse command line: Option --nohmm is incompatible with option --noF3");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF4"))       { puts("Failed to parse command line: Option --nohmm is incompatible with option --noF4");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--doF1b"))      { puts("Failed to parse command line: Option --nohmm is incompatible with option --doF1b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF2b"))      { puts("Failed to parse command line: Option --nohmm is incompatible with option --noF2b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF3b"))      { puts("Failed to parse command line: Option --nohmm is incompatible with option --noF3b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF4b"))      { puts("Failed to parse command line: Option --nohmm is incompatible with option --noF4b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--doF5b"))      { puts("Failed to parse command line: Option --nohmm is incompatible with option --doF5b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--F1"))         { puts("Failed to parse command line: Option --nohmm is incompatible with option --F1");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F1b"))        { puts("Failed to parse command line: Option --nohmm is incompatible with option --F1b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F2"))         { puts("Failed to parse command line: Option --nohmm is incompatible with option --F2");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F2b"))        { puts("Failed to parse command line: Option --nohmm is incompatible with option --F2b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F3"))         { puts("Failed to parse command line: Option --nohmm is incompatible with option --F3");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F3b"))        { puts("Failed to parse command line: Option --nohmm is incompatible with option --F3b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F4"))         { puts("Failed to parse command line: Option --nohmm is incompatible with option --F4");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F4b"))        { puts("Failed to parse command line: Option --nohmm is incompatible with option --F4b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F5"))         { puts("Failed to parse command line: Option --nohmm is incompatible with option --F5");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--ftau"))       { puts("Failed to parse command line: Option --nohmm is incompatible with option --ftau");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--fsums"))      { puts("Failed to parse command line: Option --nohmm is incompatible with option --fsums");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--tau"))        { puts("Failed to parse command line: Option --nohmm is incompatible with option --tau");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--sums"))       { puts("Failed to parse command line: Option --nohmm is incompatible with option --sums");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--rt1"))        { puts("Failed to parse command line: Option --nohmm is incompatible with option --rt1");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--rt2"))        { puts("Failed to parse command line: Option --nohmm is incompatible with option --rt2");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--rt3"))        { puts("Failed to parse command line: Option --nohmm is incompatible with option --rt3");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--ns"))         { puts("Failed to parse command line: Option --nohmm is incompatible with option --ns");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--maxtau"))     { puts("Failed to parse command line: Option --nohmm is incompatible with option --maxtau");     goto ERROR; }
		if(esl_opt_IsUsed(go, "--anytrunc"))   { puts("Failed to parse command line: Option --nohmm is incompatible with option --anytrunc");   goto ERROR; }
	}
	if(esl_opt_IsUsed(go, "--mid")) { 
		if(esl_opt_IsUsed(go, "--max"))      { puts("Failed to parse command line: Option --mid is incompatible with option --max");   goto ERROR; }
		if(esl_opt_IsUsed(go, "--nohmm"))    { puts("Failed to parse command line: Option --mid is incompatible with option --nohmm"); goto ERROR; }
		if(esl_opt_IsUsed(go, "--rfam"))     { puts("Failed to parse command line: Option --mid is incompatible with option --rfam");  goto ERROR; }
		if(esl_opt_IsUsed(go, "--FZ"))       { puts("Failed to parse command line: Option --mid is incompatible with option --FZ");    goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF1"))     { puts("Failed to parse command line: Option --mid is incompatible with option --noF1");  goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF2"))     { puts("Failed to parse command line: Option --mid is incompatible with option --noF2");  goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF3"))     { puts("Failed to parse command line: Option --mid is incompatible with option --noF3");  goto ERROR; }
		if(esl_opt_IsUsed(go, "--doF1b"))    { puts("Failed to parse command line: Option --mid is incompatible with option --doF1b"); goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF2b"))    { puts("Failed to parse command line: Option --mid is incompatible with option --noF2b"); goto ERROR; }
		if(esl_opt_IsUsed(go, "--F1"))       { puts("Failed to parse command line: Option --mid is incompatible with option --F1");    goto ERROR; }
		if(esl_opt_IsUsed(go, "--F1b"))      { puts("Failed to parse command line: Option --mid is incompatible with option --F1b");   goto ERROR; }
		if(esl_opt_IsUsed(go, "--F2"))       { puts("Failed to parse command line: Option --mid is incompatible with option --F2");    goto ERROR; }
		if(esl_opt_IsUsed(go, "--F2b"))      { puts("Failed to parse command line: Option --mid is incompatible with option --F2b");   goto ERROR; }
	}
	if(esl_opt_IsUsed(go, "--default")) { 
		if(esl_opt_IsUsed(go, "--max"))   { puts("Failed to parse command line: Option --default is incompatible with option --max");   goto ERROR; }
		if(esl_opt_IsUsed(go, "--nohmm")) { puts("Failed to parse command line: Option --default is incompatible with option --nohmm"); goto ERROR; }
		if(esl_opt_IsUsed(go, "--rfam"))  { puts("Failed to parse command line: Option --default is incompatible with option --rfam");  goto ERROR; }
		if(esl_opt_IsUsed(go, "--FZ"))    { puts("Failed to parse command line: Option --default is incompatible with option --FZ");    goto ERROR; }
	}
	if(esl_opt_IsUsed(go, "--rfam")) { 
		if(esl_opt_IsUsed(go, "--max"))     { puts("Failed to parse command line: Option --rfam is incompatible with option --max");     goto ERROR; }
		if(esl_opt_IsUsed(go, "--nohmm"))   { puts("Failed to parse command line: Option --rfam is incompatible with option --nohmm");   goto ERROR; }
		if(esl_opt_IsUsed(go, "--default")) { puts("Failed to parse command line: Option --rfam is incompatible with option --default"); goto ERROR; }
		if(esl_opt_IsUsed(go, "--FZ"))      { puts("Failed to parse command line: Option --rfam is incompatible with option --FZ");      goto ERROR; }
	}
	if(esl_opt_IsUsed(go, "--FZ")) { 
		if(esl_opt_IsUsed(go, "--max"))     { puts("Failed to parse command line: Option --FZ is incompatible with option --max");     goto ERROR; }
		if(esl_opt_IsUsed(go, "--nohmm"))   { puts("Failed to parse command line: Option --FZ is incompatible with option --nohmm");   goto ERROR; }
		if(esl_opt_IsUsed(go, "--default")) { puts("Failed to parse command line: Option --FZ is incompatible with option --default"); goto ERROR; }
		if(esl_opt_IsUsed(go, "--rfam"))    { puts("Failed to parse command line: Option --FZ is incompatible with option --rfam");    goto ERROR; }
	}
	if(esl_opt_IsUsed(go, "--hmmonly")) { 
		if(esl_opt_IsUsed(go, "--max"))        { puts("Failed to parse command line: Option --hmmonly is incompatible with option --max");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--nohmm"))      { puts("Failed to parse command line: Option --hmmonly is incompatible with option --nohmm");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--mid"))        { puts("Failed to parse command line: Option --hmmonly is incompatible with option --mid");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--rfam"))       { puts("Failed to parse command line: Option --hmmonly is incompatible with option --rfam");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--FZ"))         { puts("Failed to parse command line: Option --hmmonly is incompatible with option --FZ");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF1"))       { puts("Failed to parse command line: Option --hmmonly is incompatible with option --noF1");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF2"))       { puts("Failed to parse command line: Option --hmmonly is incompatible with option --noF2");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF3"))       { puts("Failed to parse command line: Option --hmmonly is incompatible with option --noF3");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF4"))       { puts("Failed to parse command line: Option --hmmonly is incompatible with option --noF4");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF6"))       { puts("Failed to parse command line: Option --hmmonly is incompatible with option --noF6");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--doF1b"))      { puts("Failed to parse command line: Option --hmmonly is incompatible with option --doF1b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF2b"))      { puts("Failed to parse command line: Option --hmmonly is incompatible with option --noF2b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF3b"))      { puts("Failed to parse command line: Option --hmmonly is incompatible with option --noF3b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--noF4b"))      { puts("Failed to parse command line: Option --hmmonly is incompatible with option --noF4b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--doF5b"))      { puts("Failed to parse command line: Option --hmmonly is incompatible with option --doF5b");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--F1"))         { puts("Failed to parse command line: Option --hmmonly is incompatible with option --F1");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F1b"))        { puts("Failed to parse command line: Option --hmmonly is incompatible with option --F1b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F2"))         { puts("Failed to parse command line: Option --hmmonly is incompatible with option --F2");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F2b"))        { puts("Failed to parse command line: Option --hmmonly is incompatible with option --F2b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F3"))         { puts("Failed to parse command line: Option --hmmonly is incompatible with option --F3");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F3b"))        { puts("Failed to parse command line: Option --hmmonly is incompatible with option --F3b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F4"))         { puts("Failed to parse command line: Option --hmmonly is incompatible with option --F4");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F4b"))        { puts("Failed to parse command line: Option --hmmonly is incompatible with option --F4b");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--F5"))         { puts("Failed to parse command line: Option --hmmonly is incompatible with option --F5");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--F6"))         { puts("Failed to parse command line: Option --hmmonly is incompatible with option --F6");         goto ERROR; }
		if(esl_opt_IsUsed(go, "--ftau"))       { puts("Failed to parse command line: Option --hmmonly is incompatible with option --ftau");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--fsums"))      { puts("Failed to parse command line: Option --hmmonly is incompatible with option --fsums");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--fqdb"))       { puts("Failed to parse command line: Option --hmmonly is incompatible with option --fqdb");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--fbeta"))      { puts("Failed to parse command line: Option --hmmonly is incompatible with option --fbeta");      goto ERROR; }
		if(esl_opt_IsUsed(go, "--fnonbanded")) { puts("Failed to parse command line: Option --hmmonly is incompatible with option --fnonbanded"); goto ERROR; }
		if(esl_opt_IsUsed(go, "--nocykenv"))   { puts("Failed to parse command line: Option --hmmonly is incompatible with option --nocykenv");   goto ERROR; }
		if(esl_opt_IsUsed(go, "--cykenvx"))    { puts("Failed to parse command line: Option --hmmonly is incompatible with option --cykenvx");    goto ERROR; }
		if(esl_opt_IsUsed(go, "--tau"))        { puts("Failed to parse command line: Option --hmmonly is incompatible with option --tau");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--sums"))       { puts("Failed to parse command line: Option --hmmonly is incompatible with option --sums");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--qdb"))        { puts("Failed to parse command line: Option --hmmonly is incompatible with option --qdb");        goto ERROR; }
		if(esl_opt_IsUsed(go, "--beta"))       { puts("Failed to parse command line: Option --hmmonly is incompatible with option --beta");       goto ERROR; }
		if(esl_opt_IsUsed(go, "--nonbanded"))  { puts("Failed to parse command line: Option --hmmonly is incompatible with option --nonbanded");  goto ERROR; }
		if(esl_opt_IsUsed(go, "--maxtau"))     { puts("Failed to parse command line: Option --hmmonly is incompatible with option --maxtau");     goto ERROR; }
		if(esl_opt_IsUsed(go, "--anytrunc"))   { puts("Failed to parse command line: Option --hmmonly is incompatible with option --anytrunc");   goto ERROR; }
		if(esl_opt_IsUsed(go, "--mxsize"))     { puts("Failed to parse command line: Option --hmmonly is incompatible with option --mxsize");     goto ERROR; }
		if(esl_opt_IsUsed(go, "--smxsize"))    { puts("Failed to parse command line: Option --hmmonly is incompatible with option --smxsize");    goto ERROR; }
		if(esl_opt_IsUsed(go, "--nonull3"))    { puts("Failed to parse command line: Option --hmmonly is incompatible with option --nonull3");    goto ERROR; }
		if(esl_opt_IsUsed(go, "--nohmmonly"))  { puts("Failed to parse command line: Option --hmmonly is incompatible with option --nohmmonly");  goto ERROR; }
		if(esl_opt_IsUsed(go, "--timeF4"))     { puts("Failed to parse command line: Option --hmmonly is incompatible with option --timeF4");     goto ERROR; }
		if(esl_opt_IsUsed(go, "--timeF5"))     { puts("Failed to parse command line: Option --hmmonly is incompatible with option --timeF5");     goto ERROR; }
		if(esl_opt_IsUsed(go, "--timeF6"))     { puts("Failed to parse command line: Option --hmmonly is incompatible with option --timeF6");     goto ERROR; }
		if(esl_opt_IsUsed(go, "--nogreedy"))   { puts("Failed to parse command line: Option --hmmonly is incompatible with option --nogreedy");   goto ERROR; }
		if(esl_opt_IsUsed(go, "--cp9noel"))    { puts("Failed to parse command line: Option --hmmonly is incompatible with option --cp9noel");    goto ERROR; }
		if(esl_opt_IsUsed(go, "--cp9gloc"))    { puts("Failed to parse command line: Option --hmmonly is incompatible with option --cp9gloc");    goto ERROR; }
		if(esl_opt_IsUsed(go, "--null2"))      { puts("Failed to parse command line: Option --hmmonly is incompatible with option --null2");      goto ERROR; }
	}

	*ret_go = go;
	return;

ERROR:  /* all errors handled here are user errors, so be polite.  */
	esl_fatal("there was a problem within parse_commandline in from_cmsearch.c (the embedding of cmsearch.c within CMfinder)");
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *cmfile, char *seqfile, int ncpus)
{
	cm_banner(ofp, go->argv[0], banner);

	fprintf(ofp, "# query CM file:                         %s\n", cmfile);
	fprintf(ofp, "# target sequence database:              %s\n", seqfile);
	if (esl_opt_IsUsed(go, "-g"))           fprintf(ofp, "# CM configuration:                      glocal\n");
	if (esl_opt_IsUsed(go, "-Z"))           fprintf(ofp, "# database size is set to:               %.1f Mb\n",        esl_opt_GetReal(go, "-Z"));
	if (esl_opt_IsUsed(go, "-o"))           fprintf(ofp, "# output directed to file:               %s\n",             esl_opt_GetString(go, "-o"));
	if (esl_opt_IsUsed(go, "-A"))           fprintf(ofp, "# MSA of significant hits saved to file: %s\n",             esl_opt_GetString(go, "-A"));
	if (esl_opt_IsUsed(go, "--tblout"))     fprintf(ofp, "# tabular output of hits:                %s\n",             esl_opt_GetString(go, "--tblout"));
	if (esl_opt_IsUsed(go, "--acc"))        fprintf(ofp, "# prefer accessions over names:          yes\n");
	if (esl_opt_IsUsed(go, "--noali"))      fprintf(ofp, "# show alignments in output:             no\n");
	if (esl_opt_IsUsed(go, "--notextw"))    fprintf(ofp, "# max ASCII text line length:            unlimited\n");
	if (esl_opt_IsUsed(go, "--textw"))      fprintf(ofp, "# max ASCII text line length:            %d\n",             esl_opt_GetInteger(go, "--textw"));
	if (esl_opt_IsUsed(go, "--verbose"))    fprintf(ofp, "# verbose output mode:                   on\n");
	if (esl_opt_IsUsed(go, "-E"))           fprintf(ofp, "# sequence reporting threshold:          E-value <= %g\n",  esl_opt_GetReal(go, "-E"));
	if (esl_opt_IsUsed(go, "-T"))           fprintf(ofp, "# sequence reporting threshold:          score >= %g\n",    esl_opt_GetReal(go, "-T"));
	if (esl_opt_IsUsed(go, "--incE"))       fprintf(ofp, "# sequence inclusion threshold:          E-value <= %g\n",  esl_opt_GetReal(go, "--incE"));
	if (esl_opt_IsUsed(go, "--incT"))       fprintf(ofp, "# sequence inclusion threshold:          score >= %g\n",    esl_opt_GetReal(go, "--incT"));
	if (esl_opt_IsUsed(go, "--cut_ga"))     fprintf(ofp, "# model-specific thresholding:           GA cutoffs\n");
	if (esl_opt_IsUsed(go, "--cut_nc"))     fprintf(ofp, "# model-specific thresholding:           NC cutoffs\n");
	if (esl_opt_IsUsed(go, "--cut_tc"))     fprintf(ofp, "# model-specific thresholding:           TC cutoffs\n");
	if (esl_opt_IsUsed(go, "--max"))        fprintf(ofp, "# Max sensitivity mode:                  on [all heuristic filters off]\n");
	if (esl_opt_IsUsed(go, "--nohmm"))      fprintf(ofp, "# CM-only mode:                          on [HMM filters off]\n");
	if (esl_opt_IsUsed(go, "--mid"))        fprintf(ofp, "# HMM MSV and Viterbi filters:           off\n");
	if (esl_opt_IsUsed(go, "--rfam"))       fprintf(ofp, "# Rfam pipeline mode:                    on [strict filtering]\n");
	if (esl_opt_IsUsed(go, "--FZ"))         fprintf(ofp, "# Filters set as if DB size in Mb is:    %f\n", esl_opt_GetReal(go, "--FZ"));
	if (esl_opt_IsUsed(go, "--Fmid"))       fprintf(ofp, "# HMM Forward filter thresholds set to:  %g\n", esl_opt_GetReal(go, "--Fmid"));
	if (esl_opt_IsUsed(go, "--hmmonly"))    fprintf(ofp, "# HMM-only mode (for all models):        on [CM will not be used]\n");
	if (esl_opt_IsUsed(go, "--notrunc"))    fprintf(ofp, "# truncated sequence detection:          off\n");
	if (esl_opt_IsUsed(go, "--anytrunc"))   fprintf(ofp, "# allowing truncated sequences anywhere: on\n");
	if (esl_opt_IsUsed(go, "--nonull3"))    fprintf(ofp, "# null3 bias corrections:                off\n");
	if (esl_opt_IsUsed(go, "--mxsize"))     fprintf(ofp, "# maximum DP alignment matrix size:      %.1f Mb\n", esl_opt_GetReal(go, "--mxsize"));
	if (esl_opt_IsUsed(go, "--smxsize"))    fprintf(ofp, "# maximum DP search matrix size:         %.1f Mb\n", esl_opt_GetReal(go, "--smxsize"));
	if (esl_opt_IsUsed(go, "--cyk"))        fprintf(ofp, "# use CYK for final search stage         on\n");
	if (esl_opt_IsUsed(go, "--acyk"))       fprintf(ofp, "# use CYK to align hits:                 on\n");
	if (esl_opt_IsUsed(go, "--wcx"))        fprintf(ofp, "# W set as <x> * cm->clen:               <x>=%g\n", esl_opt_GetReal(go, "--wcx"));
	if (esl_opt_IsUsed(go, "--toponly"))    fprintf(ofp, "# search top-strand only:                on\n");
	if (esl_opt_IsUsed(go, "--bottomonly")) fprintf(ofp, "# search bottom-strand only:             on\n");
	if (esl_opt_IsUsed(go, "--tformat"))    fprintf(ofp, "# targ <seqdb> format asserted:          %s\n", esl_opt_GetString(go, "--tformat"));
	/* Developer options, only shown to user if --devhelp used */
	if (esl_opt_IsUsed(go, "--noF1"))       fprintf(ofp, "# HMM MSV filter:                        off\n");
	if (esl_opt_IsUsed(go, "--noF2"))       fprintf(ofp, "# HMM Vit filter:                        off\n");
	if (esl_opt_IsUsed(go, "--noF3"))       fprintf(ofp, "# HMM Fwd filter:                        off\n");
	if (esl_opt_IsUsed(go, "--noF4"))       fprintf(ofp, "# HMM glocal Fwd filter:                 off\n");
	if (esl_opt_IsUsed(go, "--noF6"))       fprintf(ofp, "# CM CYK filter:                         off\n");
	if (esl_opt_IsUsed(go, "--doF1b"))      fprintf(ofp, "# HMM MSV biased comp filter:            on\n");
	if (esl_opt_IsUsed(go, "--noF2b"))      fprintf(ofp, "# HMM Vit biased comp filter:            off\n");
	if (esl_opt_IsUsed(go, "--noF3b"))      fprintf(ofp, "# HMM Fwd biased comp filter:            off\n");
	if (esl_opt_IsUsed(go, "--noF4b"))      fprintf(ofp, "# HMM gFwd biased comp filter:           off\n");
	if (esl_opt_IsUsed(go, "--doF5b"))      fprintf(ofp, "# HMM per-envelope biased comp filter:   on\n");
	if (esl_opt_IsUsed(go, "--F1"))         fprintf(ofp, "# HMM MSV filter P threshold:            <= %g\n", esl_opt_GetReal(go, "--F1"));
	if (esl_opt_IsUsed(go, "--F1b"))        fprintf(ofp, "# HMM MSV bias P threshold:              <= %g\n", esl_opt_GetReal(go, "--F1b"));
	if (esl_opt_IsUsed(go, "--F2"))         fprintf(ofp, "# HMM Vit filter P threshold:            <= %g\n", esl_opt_GetReal(go, "--F2"));
	if (esl_opt_IsUsed(go, "--F2b"))        fprintf(ofp, "# HMM Vit bias P threshold:              <= %g\n", esl_opt_GetReal(go, "--F2b"));
	if (esl_opt_IsUsed(go, "--F3"))         fprintf(ofp, "# HMM Fwd filter P threshold:            <= %g\n", esl_opt_GetReal(go, "--F3"));
	if (esl_opt_IsUsed(go, "--F3b"))        fprintf(ofp, "# HMM Fwd bias P threshold:              <= %g\n", esl_opt_GetReal(go, "--F3b"));
	if (esl_opt_IsUsed(go, "--F4"))         fprintf(ofp, "# HMM glocal Fwd filter P threshold:     <= %g\n", esl_opt_GetReal(go, "--F4"));
	if (esl_opt_IsUsed(go, "--F4b"))        fprintf(ofp, "# HMM glocal Fwd bias P threshold:       <= %g\n", esl_opt_GetReal(go, "--F4b"));
	if (esl_opt_IsUsed(go, "--F5"))         fprintf(ofp, "# HMM env defn filter P threshold:       <= %g\n", esl_opt_GetReal(go, "--F5"));
	if (esl_opt_IsUsed(go, "--F5b"))        fprintf(ofp, "# HMM env defn bias   P threshold:       <= %g\n", esl_opt_GetReal(go, "--F5b"));
	if (esl_opt_IsUsed(go, "--F6"))         fprintf(ofp, "# CM CYK filter P threshold:             <= %g\n", esl_opt_GetReal(go, "--F6"));

	if (esl_opt_IsUsed(go, "--hmmmax"))     fprintf(ofp, "# max sensitivity mode   (HMM-only):     on [all heuristic filters off]\n");
	if (esl_opt_IsUsed(go, "--hmmF1"))      fprintf(ofp, "# HMM MSV filter P threshold (HMM-only)  <= %g\n", esl_opt_GetReal(go, "--hmmF1"));
	if (esl_opt_IsUsed(go, "--hmmF2"))      fprintf(ofp, "# HMM Vit filter P threshold (HMM-only)  <= %g\n", esl_opt_GetReal(go, "--hmmF2"));
	if (esl_opt_IsUsed(go, "--hmmF3"))      fprintf(ofp, "# HMM Fwd filter P threshold (HMM-only)  <= %g\n", esl_opt_GetReal(go, "--hmmF3"));
	if (esl_opt_IsUsed(go, "--hmmnobias"))  fprintf(ofp, "# HMM MSV biased comp filter (HMM-only)  off\n");
	if (esl_opt_IsUsed(go, "--hmmnonull2")) fprintf(ofp, "# null2 bias corrections (HMM-only):     off\n");
	if (esl_opt_IsUsed(go, "--nohmmonly"))  fprintf(ofp, "# HMM-only mode for 0 basepair models:   no\n");

	if (esl_opt_IsUsed(go, "--rt1"))        fprintf(ofp, "# domain definition rt1 parameter        %g\n", esl_opt_GetReal(go, "--rt1"));
	if (esl_opt_IsUsed(go, "--rt2"))        fprintf(ofp, "# domain definition rt2 parameter        %g\n", esl_opt_GetReal(go, "--rt2"));
	if (esl_opt_IsUsed(go, "--rt3"))        fprintf(ofp, "# domain definition rt3 parameter        %g\n", esl_opt_GetReal(go, "--rt3"));
	if (esl_opt_IsUsed(go, "--ns"))         fprintf(ofp, "# number of envelope tracebacks sampled  %d\n", esl_opt_GetInteger(go, "--ns"));

	if (esl_opt_IsUsed(go, "--ftau"))       fprintf(ofp, "# tau parameter for CYK filter stage:    %g\n", esl_opt_GetReal(go, "--ftau"));
	if (esl_opt_IsUsed(go, "--fsums"))      fprintf(ofp, "# posterior sums (CYK filter stage):     on\n");
	if (esl_opt_IsUsed(go, "--fqdb"))       fprintf(ofp, "# QDBs (CYK filter stage)                on\n");
	if (esl_opt_IsUsed(go, "--fbeta"))      fprintf(ofp, "# beta parameter for CYK filter stage:   %g\n", esl_opt_GetReal(go, "--fbeta"));
	if (esl_opt_IsUsed(go, "--fnonbanded")) fprintf(ofp, "# no bands (CYK filter stage)            on\n");
	if (esl_opt_IsUsed(go, "--nocykenv"))   fprintf(ofp, "# CYK envelope redefinition:             off\n");
	if (esl_opt_IsUsed(go, "--cykenvx"))    fprintf(ofp, "# CYK envelope redefn P-val multiplier:  %d\n", esl_opt_GetInteger(go, "--cykenvx"));   

	if (esl_opt_IsUsed(go, "--tau"))        fprintf(ofp, "# tau parameter for final stage:         %g\n", esl_opt_GetReal(go, "--tau"));
	if (esl_opt_IsUsed(go, "--sums"))       fprintf(ofp, "# posterior sums (final stage):          on\n");
	if (esl_opt_IsUsed(go, "--qdb"))        fprintf(ofp, "# QDBs (final stage)                     on\n");
	if (esl_opt_IsUsed(go, "--beta"))       fprintf(ofp, "# beta parameter for final stage:        %g\n", esl_opt_GetReal(go, "--beta"));
	if (esl_opt_IsUsed(go, "--nonbanded"))  fprintf(ofp, "# no bands (final stage)                 on\n");

	if (esl_opt_IsUsed(go, "--timeF1"))     fprintf(ofp, "# abort after Stage 1 MSV (for timing)   on\n");
	if (esl_opt_IsUsed(go, "--timeF2"))     fprintf(ofp, "# abort after Stage 2 Vit (for timing)   on\n");
	if (esl_opt_IsUsed(go, "--timeF3"))     fprintf(ofp, "# abort after Stage 3 Fwd (for timing)   on\n");
	if (esl_opt_IsUsed(go, "--timeF4"))     fprintf(ofp, "# abort after Stage 4 gFwd (for timing)  on\n");
	if (esl_opt_IsUsed(go, "--timeF5"))     fprintf(ofp, "# abort after Stage 5 env defn (for timing) on\n");
	if (esl_opt_IsUsed(go, "--timeF6"))     fprintf(ofp, "# abort after Stage 6 CYK (for timing)   on\n");

	if (esl_opt_IsUsed(go, "--nogreedy"))   fprintf(ofp, "# greedy CM hit resolution:              off\n");
	if (esl_opt_IsUsed(go, "--cp9noel"))    fprintf(ofp, "# CP9 HMM local ends:                    off\n");
	if (esl_opt_IsUsed(go, "--cp9gloc"))    fprintf(ofp, "# CP9 HMM configuration:                 glocal\n");
	if (esl_opt_IsUsed(go, "--null2"))      fprintf(ofp, "# null2 bias corrections:                on\n");
	if (esl_opt_IsUsed(go, "--maxtau"))     fprintf(ofp, "# max tau during band tightening:        %g\n", esl_opt_GetReal(go, "--maxtau"));
	if (esl_opt_IsUsed(go, "--seed"))  {
		if (esl_opt_GetInteger(go, "--seed") == 0) fprintf(ofp, "# random number seed:                    one-time arbitrary\n");
		else                                       fprintf(ofp, "# random number seed set to:             %d\n", esl_opt_GetInteger(go, "--seed"));
	}
	/* if --max or --nohmm, truncated ends will be turned off 
	* in cm_pipeline_Create() (hopefully this will go away soon when D&C
	* truncated alignment is fixed).
	*/
	if(! (esl_opt_IsUsed(go, "--notrunc"))) { 
		if(esl_opt_IsUsed(go, "--max"))        {  fprintf(ofp, "# truncated hit detection:               off [due to --max]\n"); }
		if(esl_opt_IsUsed(go, "--nohmm"))      {  fprintf(ofp, "# truncated hit detection:               off [due to --nohmm]\n"); }
	}
	/* output number of processors being used, always (this differs from H3 which only does this if --cpu) */
	int output_ncpu = FALSE;
#ifdef HMMER_THREADS
	if (! output_ncpu)                       {  fprintf(ofp, "# number of worker threads:              %d%s\n", ncpus, (esl_opt_IsUsed(go, "--cpu") ? " [--cpu]" : "")); output_ncpu = TRUE; }
#endif 
	if (! output_ncpu)                       {  fprintf(ofp, "# number of worker threads:              0 [serial mode; threading unavailable]\n"); }
	fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n");
	return eslOK;
}

/* Function:  open_dbfile
* Synopsis:  Open the database file.
* Incept:    EPN, Mon Jun  6 09:13:08 2011
*
* Returns: eslOK on success. ESL_SQFILE *dbfp in *ret_dbfp.
*          Upon error, fills errbuf with error message and
*          returns error status code.
*/
int
open_dbfile(ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_SQFILE **ret_dbfp)
{
	int status;
	int dbfmt = eslSQFILE_UNKNOWN; /* format of dbfile                                */

	if (esl_opt_IsOn(go, "--tformat")) {
		dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
		if (dbfmt == eslSQFILE_UNKNOWN) ESL_FAIL(eslEINVAL, errbuf, "%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
	}
	status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, ret_dbfp);
	if      (status == eslENOTFOUND) ESL_FAIL(status, errbuf, "Failed to open sequence file %s for reading\n",          cfg->dbfile);
	else if (status == eslEFORMAT)   ESL_FAIL(status, errbuf, "Sequence file %s is empty or misformatted\n",            cfg->dbfile);
	else if (status == eslEINVAL)    ESL_FAIL(status, errbuf, "Can't autodetect format of a stdin or .gz seqfile");
	else if (status != eslOK)        ESL_FAIL(status, errbuf, "Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  

	if (esl_sqio_IsAlignment((*ret_dbfp)->format)) { 
		/* file is an alignment format, we can't deal with that since it may be interleaved 
		* and esl_sq_ReadBlock() can't handle that since it uses ReadWindow() to read 
		* possibly non-full length subsequences and isn't implemented to handle alignments.
		*/
		cm_Fail("%s autodetected as an alignment file; unaligned sequence file, like FASTA, is required\n", cfg->dbfile);
	}

	return eslOK;
}

/* Function:  dbsize_and_seq_lengths()
* Synopsis:  Determine size of the database by reading it
*            (but not storing it), as well as sequence lengths.
* Incept:    EPN, Mon Jun  6 09:17:32 2011
*
* 
* Returns:   eslOK on success: 
*               cfg->Z is set
*               array of seq lengths returned in *ret_srcL
*               # seqs returned in *ret_nseqs
*            eslEFORMAT if database file is screwy.
*/
int
dbsize_and_seq_lengths(ESL_GETOPTS *go, struct cfg_s *cfg, ESL_SQFILE **dbfp_ptr, char *errbuf, int64_t **ret_srcL, int64_t *ret_nseqs)
{
	int       status;
	ESL_SQ   *sq = NULL;
	int64_t   nres = 0;     /* total number of residues */
	/* variables used to store lengths of all target sequences */
	int64_t  *srcL        = NULL;    /* [0..nseqs-1] full length of each target sequence read */
	int64_t   nseqs       = 0;       /* total number of sequences */
	int64_t   nalloc_srcL = 0;       /* current allocation size of srcL */
	int       alloc_srcL  = 10000;   /* chunk size to increase allocation by for srcL */
	int64_t   i;                     /* counter */       
	char     *tmp_filename = NULL;   /* name of sqfile, used only if gzipped */
	int       tmp_fmt;               /* fmt of sqfile, used only if gzipped */

	/* we'll only use this if seqfile is gzipped */
	ESL_ALLOC(tmp_filename, sizeof(char) * (strlen((*dbfp_ptr)->filename) + 1));

	sq = esl_sq_Create();
	while ((status = esl_sqio_ReadInfo(*dbfp_ptr, sq)) == eslOK) { 
		nres += sq->L;
		if(nseqs == nalloc_srcL) { /* reallocate */
			nalloc_srcL += alloc_srcL;
			ESL_REALLOC(srcL, sizeof(int64_t) * nalloc_srcL);
			for(i = nalloc_srcL-alloc_srcL; i < nalloc_srcL; i++) srcL[i] = -1; /* initialize */
		}      
		srcL[nseqs++] = sq->L;
		esl_sq_Reuse(sq);
	}
	if(status != eslEOF) goto ERROR; 

	/* if we get here we've successfully read entire file */
	cfg->Z = nres;
	if((! esl_opt_GetBoolean(go, "--toponly")) && 
		(! esl_opt_GetBoolean(go, "--bottomonly"))) { 
			cfg->Z *= 2; /* we're searching both strands */
	}
	cfg->Z_setby = CM_ZSETBY_FILEREAD;

	if((*dbfp_ptr)->data.ascii.do_gzip == TRUE) { 
		/* file is gzipped, close it and reopen it.  
		we know we successfully opened it the first time, so a
		failure to reopen is an exception, not a user-reportable
		normal error. ENOTFOUND is the only normal error;
		EFORMAT error can't occur because we know the format and
		don't use autodetection.
		*/
		strcpy(tmp_filename, (*dbfp_ptr)->filename);
		tmp_fmt = (*dbfp_ptr)->format;
		esl_sqfile_Close(*dbfp_ptr);
		status = esl_sqfile_Open(tmp_filename, tmp_fmt, NULL, dbfp_ptr);
		if      (status == eslENOTFOUND) ESL_EXCEPTION(eslENOTFOUND, "failed to reopen alignment file");
		else if (status != eslOK)        ESL_FAIL(status, errbuf, "unexpected error when reopening sequence file");
	}
	else { 
		/* file is not gzipped, easier case */
		esl_sqfile_Position((*dbfp_ptr), 0);
	}

	if(tmp_filename != NULL) free(tmp_filename);
	if(sq           != NULL) esl_sq_Destroy(sq);

	*ret_srcL  = srcL;
	*ret_nseqs = nseqs;

	return eslOK;

ERROR: 
	if((*dbfp_ptr)->data.ascii.do_gzip == TRUE) { /* same as above, close and reopen */
		strcpy(tmp_filename, (*dbfp_ptr)->filename);
		tmp_fmt = (*dbfp_ptr)->format;
		esl_sqfile_Close((*dbfp_ptr));
		status = esl_sqfile_Open(tmp_filename, tmp_fmt, NULL, dbfp_ptr);
		if      (status == eslENOTFOUND) ESL_EXCEPTION(eslENOTFOUND, "failed to reopen alignment file");
		else if (status != eslOK)        ESL_XFAIL(status, errbuf, "unexpected error when reopening sequence file");
	}
	else { 
		/* file is not gzipped, easier case */
		esl_sqfile_Position((*dbfp_ptr), 0);
	}
	if(tmp_filename != NULL) free(tmp_filename);
	if(sq           != NULL) esl_sq_Destroy(sq);

	*ret_srcL  = NULL;
	*ret_nseqs = 0;
	return status;
}

/* Function:  create_info
* Incept:    EPN, Mon Jun  6 14:52:29 2011
*
* Purpose:  Create (allocate and initalize) a WORKER_INFO object.
*
* Returns:  The new WORKER_INFO is returned. NULL is returned upon an error.
*/
WORKER_INFO *
create_info(const ESL_GETOPTS *go)
{ 
	int status;
	WORKER_INFO *info = NULL;

	ESL_ALLOC(info, sizeof(WORKER_INFO));
	info->pli     = NULL;
	info->th      = NULL;
	info->cm      = NULL;
	info->gm      = NULL;
	info->Rgm     = NULL;
	info->Lgm     = NULL;
	info->Tgm     = NULL;
	info->om      = NULL;
	info->bg      = NULL;
	info->msvdata = NULL;
	ESL_ALLOC(info->p7_evparam, sizeof(float) * CM_p7_NEVPARAM);
	info->smxsize = esl_opt_GetReal(go, "--smxsize");
	return info;

ERROR: 
	if(info != NULL) free(info);
	return NULL;
}

/* Function:  clone_info
* Incept:    EPN, Mon Jun  6 10:24:21 2011
*
* Purpose: Given a template WORKER_INFO <tinfo>, clone it into the
*          <dest_infocnt> WORKER_INFOs in <dest_infoA>. After cloning the CM,
*          configure it.
*
* Returns: <eslOK> on success.
*          <eslEMEM> if out of memory, errbuf filled.
*/
int
clone_info(ESL_GETOPTS *go, WORKER_INFO *src_info, WORKER_INFO *dest_infoA, int dest_infocnt, char *errbuf)
{ 
	int status;
	int i;

	for (i = 0; i < dest_infocnt; ++i) {
		if((status = cm_Clone(src_info->cm, errbuf, &(dest_infoA[i].cm))) != eslOK) return status;
		if((dest_infoA[i].gm  = p7_profile_Clone(src_info->gm)) == NULL) goto ERROR;
		dest_infoA[i].Rgm = dest_infoA[i].Lgm = dest_infoA[i].Tgm = NULL; /* changed below if nec */
		if(src_info->Rgm != NULL) { if((dest_infoA[i].Rgm = p7_profile_Clone(src_info->Rgm)) == NULL) goto ERROR; }
		if(src_info->Lgm != NULL) { if((dest_infoA[i].Lgm = p7_profile_Clone(src_info->Lgm)) == NULL) goto ERROR; }
		if(src_info->Tgm != NULL) { if((dest_infoA[i].Tgm = p7_profile_Clone(src_info->Tgm)) == NULL) goto ERROR; }
		if((dest_infoA[i].om  = p7_oprofile_Clone(src_info->om)) == NULL) goto ERROR;
		if((dest_infoA[i].bg  = p7_bg_Create(src_info->bg->abc)) == NULL) goto ERROR;
		if(dest_infoA[i].p7_evparam == NULL) ESL_ALLOC(dest_infoA[i].p7_evparam, sizeof(float) * CM_p7_NEVPARAM);
		esl_vec_FCopy(src_info->cm->fp7_evparam, CM_p7_NEVPARAM, dest_infoA[i].p7_evparam);
		dest_infoA[i].msvdata = p7_hmm_MSVDataClone(src_info->msvdata, src_info->om->abc->Kp);
	}
	return eslOK;

ERROR:
	ESL_FAIL(status, errbuf, "clone_info(): out of memory");
}

/* Function:  free_info
* Incept:    EPN, Mon Jun  6 10:46:27 2011
*
* Purpose:  Free a WORKER_INFO object.
*
* Returns: void. Dies immediately upon an error.
*/
void
free_info(WORKER_INFO *info)
{ 
	if(info->pli        != NULL) cm_pipeline_Destroy(info->pli, info->cm); info->pli        = NULL;
	if(info->th         != NULL) cm_tophits_Destroy(info->th);             info->th         = NULL;
	if(info->cm         != NULL) FreeCM(info->cm);                         info->cm         = NULL;
	if(info->om         != NULL) p7_oprofile_Destroy(info->om);            info->om         = NULL;
	if(info->gm         != NULL) p7_profile_Destroy(info->gm);             info->gm         = NULL;
	if(info->Rgm        != NULL) p7_profile_Destroy(info->Rgm);            info->Rgm        = NULL;
	if(info->Lgm        != NULL) p7_profile_Destroy(info->Lgm);            info->Lgm        = NULL;
	if(info->Tgm        != NULL) p7_profile_Destroy(info->Tgm);            info->Tgm        = NULL;
	if(info->bg         != NULL) p7_bg_Destroy(info->bg);                  info->bg         = NULL;
	if(info->p7_evparam != NULL) free(info->p7_evparam);                   info->p7_evparam = NULL;
	if(info->msvdata    != NULL) p7_hmm_MSVDataDestroy(info->msvdata);     info->msvdata    = NULL;

	return;
}

/* Function:  configure_cm()
* Incept:    EPN, Mon Jun  6 12:06:38 2011
*
* Purpose: Given a WORKER_INFO <info> with a valid CM just read from
*          a file and config_opts from <info->pli>, configure the CM.
*          We use pli->cm_config_opts created in cm_pipeline_Create()
*          so that cmscan and cmsearch create config_opts identically
*          based on their common esl_getopts object.
*
*          Also create all the structures related to the p7 HMM
*          filter.
*
* Returns: eslOK on success. Upon an error, fills errbuf with
*          message and returns appropriate error status code.
*/
int
configure_cm(WORKER_INFO *info)
{ 
	int status;
	float reqMb = 0.;
	int   check_fcyk_beta;    /* TRUE if we need to check if cm's beta1 == info->pli->fcyk_beta */
	int   check_final_beta;   /* TRUE if we need to check if cm's beta2 == info->pli->final_beta */
	int   W_from_cmdline;     /* -1 (W not set on cmdline) unless pli->use_wcx is TRUE, which means --wcx was enabled */

	if(info->pli->cm_config_opts & CM_CONFIG_SCANMX)   reqMb += cm_scan_mx_SizeNeeded(info->cm, TRUE, TRUE);
	if(info->pli->cm_config_opts & CM_CONFIG_TRSCANMX) reqMb += cm_tr_scan_mx_SizeNeeded(info->cm, TRUE, TRUE);
	if(reqMb > info->smxsize) { 
		ESL_FAIL(eslERANGE, info->pli->errbuf, "search will require %.2f Mb > %.2f Mb limit.\nIncrease limit with --smxsize, or don't use --max,--nohmm,--qdb,--fqdb.", reqMb, info->smxsize);
	}

	/* cm_pipeline_Create() sets configure/align options in pli->cm_config_opts, pli->cm_align_opts */
	info->cm->config_opts = info->pli->cm_config_opts;
	info->cm->align_opts  = info->pli->cm_align_opts;

	/* check if we need to recalculate QDBs prior to building the scan matrix in cm_Configure() */
	check_fcyk_beta  = (info->pli->fcyk_cm_search_opts  & CM_SEARCH_QDB) ? TRUE : FALSE;
	check_final_beta = (info->pli->final_cm_search_opts & CM_SEARCH_QDB) ? TRUE : FALSE;
	if((status = CheckCMQDBInfo(info->cm->qdbinfo, info->pli->fcyk_beta, check_fcyk_beta, info->pli->final_beta, check_final_beta)) != eslOK) { 
		info->cm->config_opts   |= CM_CONFIG_QDB;
		info->cm->qdbinfo->beta1 = info->pli->fcyk_beta;
		info->cm->qdbinfo->beta2 = info->pli->final_beta;
	}
	/* else we don't have to change (*opt_cm)->qdbinfo->beta1/beta2 */

	W_from_cmdline = info->pli->do_wcx ? (int) (info->cm->clen * info->pli->wcx) : -1; /* -1: use W from CM file */
	if((status = cm_Configure(info->cm, info->pli->errbuf, W_from_cmdline)) != eslOK) return status;

	return eslOK;
}

/* Function:  setup_hmm_filter()
* Incept:    EPN, Mon Jun  6 11:31:42 2011
*
* Purpose:  Given a WORKER_INFO <info> with a valid non-configured
*           CM (just read from a file), set up the HMM 
*           filter related data in <info>.
*
*           Note: this is separate from cm_Configure() only 
*           because we need to call this before we call
*           clone_info(), but we have to call cm_Configure()
*           after we call clone_info() (only because we 
*           can't clone a configured CM, just a non-configured
*           one).
*
* Returns: <eslOK> on success. Dies immediately upon an error.
*/
int
setup_hmm_filter(ESL_GETOPTS *go, WORKER_INFO *info)
{ 
	int do_trunc_ends = (esl_opt_GetBoolean(go, "--notrunc") || esl_opt_GetBoolean(go, "--anytrunc")) ? FALSE : TRUE;

	/* set up the HMM filter-related structures */
	info->gm = p7_profile_Create (info->cm->fp7->M, info->cm->abc);
	info->om = p7_oprofile_Create(info->cm->fp7->M, info->cm->abc);
	info->bg = p7_bg_Create(info->cm->abc);
	p7_ProfileConfig(info->cm->fp7, info->bg, info->gm, 100, p7_LOCAL);  /* 100 is a dummy length for now; and MSVFilter requires local mode */
	p7_oprofile_Convert(info->gm, info->om);                             /* <om> is now p7_LOCAL, multihit */
	/* clone gm into Tgm before putting it into glocal mode */
	if(do_trunc_ends) { 
		info->Tgm = p7_profile_Clone(info->gm);
	}
	/* after om has been created, convert gm to glocal, to define envelopes in cm_pipeline() */
	p7_ProfileConfig(info->cm->fp7, info->bg, info->gm, 100, p7_GLOCAL);

	if(do_trunc_ends) { 
		/* create Rgm, Lgm, and Tgm specially-configured profiles for defining envelopes around 
		* hits that may be truncated 5' (Rgm), 3' (Lgm) or both (Tgm). */
		info->Rgm = p7_profile_Clone(info->gm);
		info->Lgm = p7_profile_Clone(info->gm);
		/* info->Tgm was created when gm was still in local mode above */
		/* we cloned Tgm from the while profile was still locally configured, above */
		p7_ProfileConfig5PrimeTrunc         (info->Rgm, 100);
		p7_ProfileConfig3PrimeTrunc         (info->cm->fp7, info->Lgm, 100);
		p7_ProfileConfig5PrimeAnd3PrimeTrunc(info->Tgm, 100);
	}
	else { 
		info->Rgm = NULL;
		info->Lgm = NULL;
		info->Tgm = NULL;
	}

	/* copy E-value parameters */
	esl_vec_FCopy(info->cm->fp7_evparam, CM_p7_NEVPARAM, info->p7_evparam); 

	/* compute msvdata */
	info->msvdata = p7_hmm_MSVDataCreate(info->om, FALSE);

	return eslOK;
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

