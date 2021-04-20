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

#define CONOPTS "--fast,--hand,--rsearch"                      /* Exclusive options for model construction                    */
#define WGTOPTS "--wpb,--wgsc,--wblosum,--wnone,--wgiven"      /* Exclusive options for relative weighting                    */
#define EFFOPTS "--eent,--enone,--eset"                        /* Exclusive options for effective sequence number calculation */
#define CLUSTOPTS "--ctarget,--cmaxid,--call,--corig,--cdump"  /* options for clustering the input aln and building a CM from each cluster */
ESL_OPTIONS cmbuild_options[]={
	/* name           type      default  env  range     toggles      reqs       incomp  help  docgroup*/
	{ "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "show brief help on version and usage",                     1 },
	{ "-n",        eslARG_STRING, NULL,  NULL, NULL,      NULL,      NULL,        NULL, "name the CM(s) <s>, (only if single aln in file)",         1 },
	{ "-F",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,        NULL, "force; allow overwriting of <cmfile_out>",                 1 },
	{ "-o",        eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,        NULL, "direct summary output to file <f>, not stdout",            1 },
	{ "-O",        eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,        NULL, "resave consensus/insert column annotated MSA to file <f>", 1 },
	{ "--devhelp", eslARG_NONE,   NULL,  NULL, NULL,      NULL,      NULL,        NULL, "show list of otherwise hidden developer/expert options",   1 },

	/* Expert model construction options */
	/* name          type            default  env  range       toggles       reqs        incomp  help  docgroup*/
	{ "--fast",      eslARG_NONE,"default",   NULL, NULL,      CONOPTS,      NULL,         NULL, "assign cols w/ >= symfrac residues as consensus",                2 },
	{ "--hand",      eslARG_NONE,    FALSE,   NULL, NULL,      CONOPTS,      NULL,         NULL, "use reference coordinate annotation to specify consensus",       2 },
	{ "--symfrac",   eslARG_REAL,    "0.5",   NULL, "0<=x<=1",    NULL,      NULL,         NULL, "fraction of non-gaps to require in a consensus column [0..1]",   2 },
	{ "--noss",      eslARG_NONE,    FALSE,   NULL, NULL,         NULL,      NULL,         NULL, "ignore secondary structure annotation in input alignment",       2 },
	{ "--rsearch",   eslARG_INFILE,  NULL,    NULL, NULL,      CONOPTS,      NULL,      "--p56", "use RSEARCH parameterization with RIBOSUM matrix file <f>",      2 }, 

	/* Other model construction options */
	/* name          type            default  env  range       toggles       reqs        incomp  help  docgroup*/
	{ "--null",      eslARG_INFILE,  NULL,    NULL, NULL,         NULL,      NULL,  "--rsearch", "read null (random sequence) model from file <f>",                3 },
	{ "--prior",     eslARG_INFILE,  NULL,    NULL, NULL,         NULL,      NULL,  "--rsearch", "read priors from file <f>",                                      3 },
	/* below are only shown with --devhelp */
	{ "--betaW",     eslARG_REAL,    "1E-7",  NULL, "x>1E-18",    NULL,      NULL,         NULL, "set tail loss prob for calc'ing W (max size of a hit) to <x>", 103 },
	{ "--beta1",     eslARG_REAL,    "1E-7",  NULL, "x>1E-18",    NULL,      NULL,         NULL, "set tail loss prob for calc'ing tighter set of QDBs to <x>",   103 },
	{ "--beta2",     eslARG_REAL,    "1E-15", NULL, "x>1E-18",    NULL,      NULL,         NULL, "set tail loss prob for calc'ing looser  set of QDBs to <x>",   103 },
	{ "--informat",  eslARG_STRING,  NULL,    NULL, NULL,         NULL,      NULL,         NULL, "specify input alignment is in format <s> (Stockholm or Pfam)", 103 },
	{ "--v1p0",      eslARG_NONE,    FALSE,   NULL, NULL,         NULL,      NULL,         NULL,  "parameterize CM using methods from Infernal v1.0.2",          103 },
	{ "--p56",       eslARG_NONE,    NULL,    NULL, NULL,         NULL,      NULL,    "--prior", "use the default prior from Infernal v0.56 through v1.0.2",     103 },
	{ "--noh3pri",   eslARG_NONE,    NULL,    NULL, NULL,         NULL,      NULL,"--v1p0,--p56","do not use the hmmer3 DNA prior for zero basepair models",     103 },
	{ "--iins",      eslARG_NONE,    FALSE,   NULL, NULL,         NULL,      NULL,         NULL, "allow informative insert emissions, do not zero them",         103 },
	{ "--iflank",    eslARG_NONE,    FALSE,   NULL, NULL,         NULL,      NULL,         NULL, "learn ROOT_IL/ROOT_IR transitions for 5'/3' flanking residues",103 },
	{ "--nobalance", eslARG_NONE,    FALSE,   NULL, NULL,         NULL,      NULL,         NULL, "don't rebalance the CM; number in strict preorder",            103 },
	{ "--nodetach",  eslARG_NONE,    FALSE,   NULL, NULL,         NULL,      NULL,         NULL, "do not 'detach' one of two inserts that model same column",    103 },
	{ "--elself",    eslARG_REAL,    "0.94",  NULL, "0<=x<=1",    NULL,      NULL,         NULL, "set EL self transition prob to <x>",                           103 },
	{ "--n2omega",   eslARG_REAL,    "0.000015258791",NULL,"x>0", NULL,      NULL,         NULL, "set prior probability of null2 model as <x>",                  103 }, 
	{ "--n3omega",   eslARG_REAL,    "0.000015258791",NULL,"x>0", NULL,      NULL,         NULL, "set prior probability of null3 model as <x>",                  103 }, 

	/* Alternate relative sequence weighting strategies */
	/* name        type         default   env  range     toggles         reqs  incomp  help  docgroup*/
	{ "--wpb",     eslARG_NONE,"default", NULL, NULL,    WGTOPTS,        NULL, NULL, "Henikoff position-based weights",                   4 },
	{ "--wgsc",    eslARG_NONE,     NULL, NULL, NULL,    WGTOPTS,        NULL, NULL, "Gerstein/Sonnhammer/Chothia tree weights",          4 },
	{ "--wnone",   eslARG_NONE,     NULL, NULL, NULL,    WGTOPTS,        NULL, NULL, "don't do any relative weighting; set all to 1",     4 },
	{ "--wgiven",  eslARG_NONE,     NULL, NULL, NULL,    WGTOPTS,        NULL, NULL, "use weights as given in MSA file",                  4 },
	{ "--wblosum", eslARG_NONE,     NULL, NULL, NULL,    WGTOPTS,        NULL, NULL, "Henikoff simple filter weights",                    4 },
	{ "--wid",     eslARG_REAL,   "0.62", NULL,"0<=x<=1",   NULL, "--wblosum", NULL, "for --wblosum: set identity cutoff",                4 },

	/* Alternate effective sequence weighting strategies */
	/* name        type            default    env     range toggles      reqs   incomp  help  docgroup*/
	{ "--eent",    eslARG_NONE, "default",    NULL,   NULL, EFFOPTS,     NULL,   NULL, "adjust eff seq # to achieve relative entropy target",           5 },
	{ "--enone",   eslARG_NONE,     FALSE,    NULL,   NULL, EFFOPTS,     NULL,   NULL, "no effective seq # weighting: just use nseq",                   5 },
	{ "--ere",     eslARG_REAL,      NULL,    NULL,  "x>0",    NULL, "--eent",   NULL, "for --eent: set CM target relative entropy to <x>",             5 },
	{ "--eset",    eslARG_REAL,      NULL,    NULL, "x>=0", EFFOPTS,     NULL,   NULL, "set eff seq # for all models to <x>",                           5 },
	{ "--eminseq", eslARG_REAL,     "0.1",    NULL, "x>=0",    NULL, "--eent",   NULL, "for --eent: set minimum effective sequence number to <x>",      5 },
	{ "--ehmmre",  eslARG_REAL,      NULL,    NULL,  "x>0",    NULL, "--eent",   NULL, "for --eent: set minimum HMM relative entropy to <x>",           5 }, 
	{ "--esigma",  eslARG_REAL,    "45.0",    NULL,  "x>0",    NULL, "--eent",   NULL, "for --eent: set sigma param to <x>",                            5 },

	/* Options controlling filter p7 HMM construction */
	/* name         type           default  env  range toggles  reqs  incomp    help  docgroup*/
	{ "--p7ere",    eslARG_REAL,     NULL, NULL, NULL, NULL,    NULL, "--p7ml", "for the filter p7 HMM, set minimum rel entropy/posn to <x>",   6 },
	{ "--p7ml",     eslARG_NONE,    FALSE, NULL, NULL, NULL,    NULL,     NULL, "define the filter p7 HMM as the ML p7 HMM",                    6 },
	/* below are only shown with --devhelp */
	{ "--p7prior",  eslARG_INFILE,   NULL, NULL, NULL, NULL,    NULL, "--p7ml", "read p7 prior for the filter HMM from file <f>",             106 },
	{ "--p7hprior", eslARG_NONE,     NULL, NULL, NULL, NULL,    NULL, "--p7ml", "use HMMER's default p7 prior, not Infernal's p7 prior",      106 },
	{ "--p7hemit",  eslARG_NONE,    FALSE, NULL, NULL, NULL,    NULL, "--p7ml", "use HMMER emission priors for filter HMM",                   106 }, 

	/* Options controlling filter p7 HMM calibration */
	/* name        type         default  env   range toggles   reqs  incomp       help  docgroup*/
	{ "--EmN",     eslARG_INT,    "200", NULL, "n>0",   NULL,  NULL, NULL,        "number of sampled seqs to use for p7 local MSV calibration",    7 },
	{ "--EvN",     eslARG_INT,    "200", NULL, "n>0",   NULL,  NULL, NULL,        "number of sampled seqs to use for p7 local Vit calibration",    7 },
	{ "--ElfN",    eslARG_INT,    "200", NULL, "n>0",   NULL,  NULL, NULL,        "number of sampled seqs to use for p7 local Fwd calibration",    7 },
	{ "--EgfN",    eslARG_INT,    "200", NULL, "n>0",   NULL,  NULL, NULL,        "number of sampled seqs to use for p7 glocal Fwd calibration",   7 },
	/* below are only shown with --devhelp */
	{ "--Elftp",   eslARG_REAL, "0.055", NULL, "x>0.",  NULL,  NULL, NULL,        "fit p7 local fwd exp tail to <f> fraction of scoring dist",   107 },
	{ "--Egftp",   eslARG_REAL, "0.065", NULL, "x>0.",  NULL,  NULL, NULL,        "fit p7 glocal fwd exp tail to <f> fraction of scoring dist",  107 },
	{ "--Ereal",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,        "sample realistic, not iid genomic seqs, for p7 calibration",  107 },
	{ "--Enull3",  eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,        "use null3 correction in p7 calibrations",                     107 },
	{ "--Ebias",   eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,        "use bias correction in p7 calibrations",                      107 },
	{ "--Efitlam", eslARG_NONE,   FALSE, NULL, NULL,    NULL,  NULL, NULL,        "fit lambda, don't use a fixed one near 0.693",                107 },
	{ "--Elcmult", eslARG_REAL,   "2.0", NULL, "x>0.",  NULL,  NULL, NULL,        "length of seqs to search for local stats is <x> * cm->clen",  107 },
	{ "--Egcmult", eslARG_REAL,   "2.0", NULL, "x>0.",  NULL,  NULL, NULL,        "length of seqs to search for glocal stats is <x> * cm->clen", 107 },
	{ "--ElL",     eslARG_INT,     NULL, NULL, "n>0",   NULL,  NULL, "--Elcmult", "length of seqs to search for local stats is <n>",             107 },
	{ "--EgL",     eslARG_INT,     NULL, NULL, "n>0",   NULL,  NULL, "--Egcmult", "length of seqs to search for glocal stats is <n>",            107 },

	/* Refining the input alignment */
	/* name          type            default  env  range    toggles      reqs         incomp  help  docgroup*/
	{ "--refine",    eslARG_OUTFILE,   NULL, NULL, NULL,    NULL,       NULL,           NULL, "refine input aln w/Expectation-Maximization, save to <f>",          8 },
	{ "-l",          eslARG_NONE,     FALSE, NULL, NULL,    NULL, "--refine",           NULL, "w/--refine, configure model for local alignment [default: global]", 8 },
	{ "--gibbs",     eslARG_NONE,     FALSE, NULL, NULL,    NULL, "--refine",           NULL, "w/--refine, use Gibbs sampling instead of EM",                      8 },
	{ "--seed",      eslARG_INT,        "0", NULL, "n>=0",  NULL,  "--gibbs",           NULL, "w/--gibbs, set RNG seed to <n> (if 0: one-time arbitrary seed)",    8 },
	{ "--cyk",       eslARG_NONE,     FALSE, NULL, NULL,    NULL, "--refine",           NULL, "w/--refine, use CYK instead of optimal accuracy",                   8 },
	{ "--notrunc",   eslARG_NONE,     FALSE, NULL, NULL,    NULL, "--refine",           NULL, "w/--refine, do not use truncated alignment algorithm",              8 },
	/* below are only shown with --devhelp */
	{ "--sub",       eslARG_NONE,     FALSE, NULL, NULL,    NULL, "--refine", "--notrunc,-l", "w/--refine, use sub CM for columns b/t HMM start/end points",     108 },
	{ "--nonbanded", eslARG_NONE,     FALSE, NULL, NULL,    NULL, "--refine",           NULL, "do not use bands to accelerate alignment with --refine",          108 },
	{ "--indi",      eslARG_NONE,     FALSE, NULL, NULL,    NULL, "--refine",           NULL, "print individual sequence scores during MSA refinement",          108 },
	{ "--fins",      eslARG_NONE,     FALSE, NULL, NULL,    NULL, "--refine",           NULL, "w/--refine, flush inserts left/right in alignments",              108 },
	{ "--tau",       eslARG_REAL,    "1E-7", NULL, "0<x<1", NULL, "--refine",  "--nonbanded", "set tail loss prob for HMM bands to <x>",                         108 },
	{ "--mxsize",    eslARG_REAL,  "2048.0", NULL, "x>0.",  NULL, "--refine",           NULL, "set maximum allowable DP matrix size to <x> Mb",                  108 },
	{ "--rdump",     eslARG_OUTFILE  , NULL, NULL, NULL,    NULL, "--refine",           NULL, "w/--refine, print all intermediate alignments to <f>",            108 },

	/* All options below are developer options, only shown if --devhelp invoked */
	/* Developer verbose output options */
	/* name        type          default  env   range toggles reqs  incomp help  docgroup*/
	{ "--verbose", eslARG_NONE,    FALSE, NULL, NULL,   NULL, NULL, NULL,  "be verbose with output",                              109 },
	{ "--cfile",   eslARG_OUTFILE,  NULL, NULL, NULL,   NULL, NULL, NULL,  "save count vectors to file <f>",                      109 },
	{ "--efile",   eslARG_OUTFILE,  NULL, NULL, NULL,   NULL, NULL, NULL,  "save emission score information to file <f>",         109 },
	{ "--tfile",   eslARG_OUTFILE,  NULL, NULL, NULL,   NULL, NULL, NULL,  "dump individual sequence parsetrees to file <f>",     109 },
	{ "--cmtbl",   eslARG_OUTFILE,  NULL, NULL, NULL,   NULL, NULL, NULL,  "save tabular description of CM topology to file <f>", 109 },
	{ "--emap",    eslARG_OUTFILE,  NULL, NULL, NULL,   NULL, NULL, NULL,  "save consensus emit map to file <f>",                 109 },
	{ "--gtree",   eslARG_OUTFILE,  NULL, NULL, NULL,   NULL, NULL, NULL,  "save tree description of master tree to file <f>",    109 },
	{ "--gtbl",    eslARG_OUTFILE,  NULL, NULL, NULL,   NULL, NULL, NULL,  "save tabular description of master tree to file <f>", 109 },

	/* Building multiple CMs after clustering input MSA */
	/* name        type            default env   range      toggles reqs  incomp    help  docgroup*/
	{ "--ctarget", eslARG_INT,     NULL,   NULL, "n>0" ,    NULL,   NULL, "--call", "build (at most) <n> CMs by partitioning MSA into <n> clusters", 110 },
	{ "--cmaxid",  eslARG_REAL,    NULL,   NULL, "0.<x<1.", NULL,   NULL, "--call", "max fractional id b/t 2 clusters is <x>, each cluster -> CM",   110 }, 
	{ "--call",    eslARG_NONE,    FALSE,  NULL, NULL,      NULL,   NULL,     NULL, "build a separate CM from every seq in MSA",                     110 },
	{ "--corig",   eslARG_NONE,    FALSE,  NULL, NULL,      NULL,   NULL,     NULL, "build an additional CM from the original, full MSA",            110 }, 
	{ "--cdump",   eslARG_OUTFILE, NULL,   NULL, NULL,      NULL,   NULL,     NULL, "dump the MSA for each cluster (CM) to file <f>",                110 },

	/* Developer options related to experimental local begin/end modes */
	/* name        type          default env   range      toggles reqs  incomp       help  docgroup*/
	{ "--pbegin",  eslARG_REAL,  "0.05", NULL, "0<x<1",   NULL,   NULL,  NULL,       "set aggregate local begin prob to <x>", 111 },
	{ "--pend",    eslARG_REAL,  "0.05", NULL, "0<x<1",   NULL,   NULL,  NULL,       "set aggregate local end prob to <x>",   111 },
	{ "--pebegin", eslARG_NONE,   FALSE, NULL, NULL,      NULL,   NULL,  "--pbegin", "set all local begins as equiprobable",  111 },
	{ "--pfend",   eslARG_REAL,   NULL,  NULL, "0<x<1",   NULL,   NULL,  "--pend",   "set all local end probs to <x>",        111 },

	/* options motivated by CMfinder */
	/* none, at the moment */

	{  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/* entire struct from cmbuild.c */
struct cfg_s {
	FILE         *ofp;		/* output file (default is stdout) */

	char         *alifile;	/* name of the alignment file we're building CMs from  */
	int           fmt;		/* format code for alifile */
	ESLX_MSAFILE *afp;            /* open alifile  */
	ESL_ALPHABET *abc;		/* digital alphabet */

	char         *cmfile;         /* file to write CM to                    */
	FILE         *cmoutfp;        /* CM output file handle                  */

	char         *postmsafile;	/* optional file to resave annotated MSAs to  */
	FILE         *postmsafp;	/* open <postmsafile>, or NULL */

	float        *null;		/* null model                              */
	Prior_t      *pri;		/* mixture Dirichlet prior for the CM     */
	Prior_t      *pri_zerobp;	/* mixture Dirichlet prior for any CMs with 0 basepairs */

	fullmat_t    *fullmat;        /* if --rsearch, the full RIBOSUM matrix */

	int           be_verbose;	/* standard verbose output, as opposed to one-line-per-CM summary */
	int           nali;		/* which # alignment this is in file */
	int           nnamed;		/* number of alignments that had their own names */
	int           ncm_total;      /* which # CM this is that we're constructing (we may build > 1 per file) */
	ESL_RANDOMNESS *r;            /* source of randomness, only created if --gibbs enabled */

	/* optional files used for building additional filter p7 HMMs */
	P7_BG        *fp7_bg;         /* background model for additional P7s */
	P7_BUILDER   *fp7_bld;        /* the P7_BUILDER */

	/* optional output files */
	FILE         *cfp;            /* for --cfile */
	FILE         *escfp;          /* for --efile */
	FILE         *tblfp;          /* for --cmtbl */
	FILE         *efp;            /* for --emap */
	FILE         *gfp;            /* for --gtree */
	FILE         *gtblfp;         /* for --gtbl */
	FILE         *tfp;            /* for --tfile */
	FILE         *cdfp;           /* if --cdump, output file handle for dumping clustered MSAs */
	FILE         *refinefp;       /* if --refine, output file handle for dumping refined MSAs */
	FILE         *rdfp;           /* if --rfile, output file handle for dumping intermediate MSAs during iterative refinement */
};

extern ESL_OPTIONS cmbuild_options[];

extern P7_PRIOR * cm_p7_prior_CreateNucleic(void);
extern int set_msa_name(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
extern int process_build_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t **ret_mtr, Parsetree_t ***ret_msa_tr,char *hmmStatsToCalibrate);
extern int refine_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *orig_cm, ESL_MSA *input_msa, Parsetree_t **input_msa_tr, CM_t **ret_cm, ESL_MSA **ret_msa, Parsetree_t **ret_mtr, Parsetree_t ***ret_trA, int *ret_niter);
extern int output_result(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int msaidx, int cmidx, ESL_MSA *msa, CM_t *cm, Parsetree_t *mtr, Parsetree_t **tr);
extern int check_and_clean_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
extern int    set_relative_weights(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
extern int build_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int do_print, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t **ret_mtr, Parsetree_t ***ret_msa_tr);
extern int annotate(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm);
extern int set_model_cutoffs(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm);
extern int set_effective_seqnumber(const ESL_GETOPTS *go, const struct cfg_s *cfg,char *errbuf, ESL_MSA *msa, CM_t *cm, const Prior_t *pri);
extern int parameterize(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int do_print, CM_t *cm, const Prior_t *prior, float msa_nseq);
extern int configure_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int iter);
extern int set_consensus(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm);
extern int build_and_calibrate_p7_filter(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm, int use_mlp7_as_filter,char *hmmStatsToCalibrate);
extern int print_countvectors(const struct cfg_s *cfg, char *errbuf, CM_t *cm);
extern double set_target_relent(const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int clen, int nbps);
extern double version_1p0_default_target_relent(const ESL_ALPHABET *abc, int clen, double eX);
extern int flatten_insert_emissions(CM_t *cm);
extern P7_PRIOR * p7_prior_Read(FILE *fp);

P7_PRIOR *cm_p7_prior_CreateNucleic(void)
{
	P7_PRIOR *pri = NULL;
	int q;
	int status;

	int num_comp = 4;

	static double defmq[5] = { 0.079226, 0.259549, 0.241578, 0.419647 };
	static double defm[4][4] = {
		{ 1.294511, 0.400028, 6.579555, 0.509916}, 
		{ 0.090031, 0.028634, 0.086396, 0.041186},
		{ 0.158085, 0.448297, 0.114815, 0.394151},
		{ 1.740028, 1.487773, 1.565443, 1.947555}
	};

	ESL_ALLOC(pri, sizeof(P7_PRIOR));
	pri->tm = pri->ti = pri->td = pri->em = pri->ei = NULL;

	pri->tm = esl_mixdchlet_Create(1, 3);  // match transitions; single component; 3 params
	pri->ti = esl_mixdchlet_Create(1, 2);  // insert transitions; single component; 2 params
	pri->td = esl_mixdchlet_Create(1, 2);  // delete transitions; single component; 2 params
	pri->em = esl_mixdchlet_Create(num_comp, 4); // match emissions; X component; 4 params
	pri->ei = esl_mixdchlet_Create(1, 4); // insert emissions; single component; 4 params

	if (pri->tm == NULL || pri->ti == NULL || pri->td == NULL || pri->em == NULL || pri->ei == NULL) goto ERROR;

	/* Transition priors: taken from hmmer's p7_prior.c::p7_prior_CreateNucleic() */
	/* Roughly, learned from rmark benchmark - hand-beautified (trimming overspecified significant digits)
	*/
	pri->tm->pq[0]       = 1.0;
	pri->tm->alpha[0][0] = 2.0; // TMM
	pri->tm->alpha[0][1] = 0.1; // TMI
	pri->tm->alpha[0][2] = 0.1; // TMD

	pri->ti->pq[0]       = 1.0;
	pri->ti->alpha[0][0] = 0.06; // TIM
	pri->ti->alpha[0][1] = 0.2; // TII

	pri->td->pq[0]       = 1.0;
	pri->td->alpha[0][0] = 0.1; // TDM
	pri->td->alpha[0][1] = 0.2; // TDD

	/* Match emission priors  */
	for (q = 0; q < num_comp; q++)
	{
		pri->em->pq[q] = defmq[q];
		esl_vec_DCopy(defm[q], 4, pri->em->alpha[q]);
	}

	/* Insert emission priors. */
	pri->ei->pq[0] = 1.0;
	esl_vec_DSet(pri->ei->alpha[0], 4, 1.0);

	return pri;

ERROR:
	if (pri != NULL) p7_prior_Destroy(pri);
	return NULL;
}

int
set_msa_name(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf, ESL_MSA *msa)
{
	char *name = NULL;
	int   status;

	if (cfg->nali == 1)  /* first (only?) MSA in file: */
	{
		if(esl_opt_IsUsed(go, "-n")) 
		{  
			if((status = esl_msa_SetName(msa, esl_opt_GetString(go, "-n"), -1)) != eslOK) return status;
		}
		else if (msa->name != NULL) 
		{ 
			cfg->nnamed++;
		}
		else if (cfg->afp->bf->filename) 
		{ 
			if ((status = esl_FileTail(cfg->afp->bf->filename, TRUE, &name)) != eslOK) return status; /* TRUE=nosuffix */	  
			if ((status = esl_msa_SetName(msa, name, -1))                    != eslOK) return status;
			free(name);
		}
		else ESL_FAIL(eslEINVAL, errbuf, "Failed to set model name: msa has no name, no msa filename, and no -n");
	}
	else 
	{
		if (esl_opt_IsUsed(go, "-n")) ESL_FAIL(eslEINVAL, errbuf, "Oops. Wait. You can't use -n with an alignment database.");
		else if (msa->name != NULL)   cfg->nnamed++;
		else                          ESL_FAIL(eslEINVAL, errbuf, "Oops. Wait. I need name annotation on each alignment in a multi MSA file; failed on #%d", cfg->nali);

		/* special kind of failure: the *first* alignment didn't have a name, and we used the filename to
		* construct one; now that we see a second alignment, we realize this was a boo-boo*/
		if (cfg->nnamed != cfg->nali) ESL_FAIL(eslEINVAL, errbuf, "Oops. Wait. I need name annotation on each alignment in a multi MSA file; first MSA didn't have one");
	}
	return eslOK;
}

int
process_build_workunit(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t **ret_mtr, Parsetree_t ***ret_msa_tr,char *hmmStatsToCalibrate)
{
	int      status;
	uint32_t checksum = 0;      /* checksum calculated for the input MSA. cmalign --mapali verifies against this. */
	CM_t    *cm = NULL;         /* the CM */
	int      pretend_cm_is_hmm; /* TRUE if we will use special HMM-like parameterization because this CM has 0 basepairs */
	Prior_t *pri2use = NULL;    /* cfg->pri or cfg->pri_zerobp (the latter if CM has no basepairs) */

	if ((status =  check_and_clean_msa          (go, cfg, errbuf, msa))                                 != eslOK) goto ERROR;
	if ((status =  esl_msa_Checksum             (msa, &checksum))                                       != eslOK) ESL_FAIL(status, errbuf, "Failed to calculate checksum"); 
	if ((status =  set_relative_weights         (go, cfg, errbuf, msa))                                 != eslOK) goto ERROR;
	if ((status =  build_model                  (go, cfg, errbuf, TRUE, msa, &cm, ret_mtr, ret_msa_tr)) != eslOK) goto ERROR;

	if((CMCountNodetype(cm, MATP_nd) == 0)     && 
		(! esl_opt_GetBoolean(go, "--v1p0"))    && 
		(! esl_opt_GetBoolean(go, "--p56"))     && 
		(! esl_opt_GetBoolean(go, "--noh3pri")) && 
		(! esl_opt_IsUsed    (go, "--prior"))) { 
			pretend_cm_is_hmm = TRUE;
	}
	else { 
		pretend_cm_is_hmm = FALSE;
	}
	pri2use = (pretend_cm_is_hmm) ? cfg->pri_zerobp : cfg->pri;

	if ((status =  annotate                     (go, cfg, errbuf, msa, cm))                             != eslOK) goto ERROR;
	if ((status =  set_model_cutoffs            (go, cfg, errbuf, msa, cm))                             != eslOK) goto ERROR;
	if ((status =  set_effective_seqnumber      (go, cfg, errbuf, msa, cm, pri2use))                    != eslOK) goto ERROR;
	if ((status =  parameterize                 (go, cfg, errbuf, TRUE, cm, pri2use, msa->nseq))        != eslOK) goto ERROR;
	if ((status =  configure_model              (go, cfg, errbuf, cm, 1))                               != eslOK) goto ERROR;
	if ((status =  set_consensus                (go, cfg, errbuf, cm))                                  != eslOK) goto ERROR;
	/* if <pretend_cm_is_hmm> we'll set the CM's filter p7 HMM as its maximum likelihood HMM */
	if ((status =  build_and_calibrate_p7_filter(go, cfg, errbuf, msa, cm, pretend_cm_is_hmm,hmmStatsToCalibrate))          != eslOK) goto ERROR;

	cm->checksum = checksum;
	cm->flags   |= CMH_CHKSUM;

	*ret_cm = cm;
	return eslOK;

ERROR:
	if(cm != NULL) FreeCM(cm);
	*ret_cm = NULL;
	return status;
}

int
check_and_clean_msa(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa)
{
	int status;
	ESL_STOPWATCH *w = NULL;

	if (cfg->be_verbose) {
		w = esl_stopwatch_Create();
		esl_stopwatch_Start(w);
		fprintf(cfg->ofp, "%-40s ... ", "Checking MSA");  
		fflush(cfg->ofp); 
	}

	if (esl_opt_GetBoolean(go, "--hand") && msa->rf == NULL)      ESL_FAIL(eslFAIL, errbuf, "--hand used, but alignment #%d has no reference coord annotation", cfg->nali);
	if (esl_opt_GetBoolean(go, "--noss")) { /* --noss: if SS_cons exists, strip all BPs from it; if it doesn't create it with zero bps */
		if(msa->ss_cons == NULL) { ESL_ALLOC(msa->ss_cons, sizeof(char) * (msa->alen+1)); msa->ss_cons[msa->alen] = '\0'; }
		memset(msa->ss_cons,  '.', msa->alen);
	}
	if (msa->ss_cons == NULL)                                     ESL_FAIL(eslFAIL, errbuf, "Alignment #%d has no consensus structure annotation, and --noss not used.", cfg->nali);
	if (! clean_cs(msa->ss_cons, msa->alen, (! cfg->be_verbose))) ESL_FAIL(eslFAIL, errbuf, "Failed to parse consensus structure annotation in alignment #%d", cfg->nali);

	if ( esl_opt_IsOn(go, "--rsearch")) { 
		if(msa->nseq != 1) ESL_FAIL(eslEINCOMPAT, errbuf,"with --rsearch option, all of the input alignments must have exactly 1 sequence");
		/* We can't have ambiguous bases in the MSA, only A,C,G,U will do. The reason is that rsearch_CMProbifyEmissions() expects each
		* cm->e prob vector to have exactly 1.0 count for exactly 1 singlet or base pair. If we have ambiguous residues we'll have a 
		* fraction of a count for more than one residue/base pair for some v. 
		* ribosum_MSA_resolve_degeneracies() replaces ambiguous bases with most likely compatible base */
		ribosum_MSA_resolve_degeneracies(cfg->fullmat, msa); /* cm_Fails() if some error is encountered */
	}

	/* MSA better have a name, we named it before */
	if(msa->name == NULL) ESL_FAIL(eslEINCONCEIVABLE, errbuf, "MSA is nameless, (we thought we named it...) shouldn't happen");

	if (cfg->be_verbose) { 
		fprintf(cfg->ofp, "done.  ");
		esl_stopwatch_Stop(w);
		esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
	}
	if(w != NULL) esl_stopwatch_Destroy(w);
	return eslOK;

ERROR: 
	ESL_FAIL(status, errbuf, "out of memory");
	return status; /* NOT REACHED */
}

int
set_relative_weights(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa)
{
	ESL_STOPWATCH *w = NULL;

	if (cfg->be_verbose) {
		w = esl_stopwatch_Create();
		esl_stopwatch_Start(w);
		fprintf(cfg->ofp, "%-40s ... ", "Relative sequence weighting");  
		fflush(cfg->ofp); 
	}

	if      (esl_opt_GetBoolean(go, "--wnone"))                  esl_vec_DSet(msa->wgt, msa->nseq, 1.);
	else if (esl_opt_GetBoolean(go, "--wgiven"))                 ;
	else if (esl_opt_GetBoolean(go, "--wpb"))                    esl_msaweight_PB(msa);
	else if (esl_opt_GetBoolean(go, "--wgsc"))                   esl_msaweight_GSC(msa);
	else if (esl_opt_GetBoolean(go, "--wblosum"))                esl_msaweight_BLOSUM(msa, esl_opt_GetReal(go, "--wid"));

	if (cfg->be_verbose) { 
		fprintf(cfg->ofp, "done.  ");
		esl_stopwatch_Stop(w);
		esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
	}

	if(w != NULL) esl_stopwatch_Destroy(w);
	return eslOK;
}

int
build_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int do_print, ESL_MSA *msa, CM_t **ret_cm, Parsetree_t **ret_mtr, Parsetree_t ***ret_msa_tr)
{
	int status;
	Parsetree_t     **tr;
	Parsetree_t     *mtr;
	int idx;
	CM_t *cm;
	ESL_STOPWATCH *w = NULL;
	int use_rf;
	int use_wts;
	int* used_el = NULL;
	int pretend_cm_is_hmm; /* TRUE if we will use special HMM-like parameterization because this CM has 0 basepairs */

	if (cfg->be_verbose && do_print) {
		w = esl_stopwatch_Create();
		esl_stopwatch_Start(w);
		fprintf(cfg->ofp, "%-40s ... ", "Constructing model "); 
		fflush(cfg->ofp);
	}

	use_rf  = (esl_opt_GetBoolean(go, "--hand")) ? TRUE : FALSE;
	use_wts = (use_rf || esl_opt_GetBoolean(go, "--v1p0")) ? FALSE : TRUE;
	if((status = HandModelmaker(msa, errbuf, use_rf, 
		FALSE, /* use_el: never when building a model */
		use_wts, (1. - esl_opt_GetReal(go, "--symfrac")), &cm, &mtr)) != eslOK) return status;

	/* set the CM's null model, if rsearch mode, use the bg probs used to calc RIBOSUM */
	if( esl_opt_IsOn(go, "--rsearch")) {
		CMSetNullModel(cm, cfg->fullmat->g);
	}
	else {
		CMSetNullModel(cm, cfg->null);
	}

	/* if we're using RSEARCH emissions (--rsearch) set the flag */
	if(esl_opt_GetString(go, "--rsearch") != NULL) cm->flags |= CM_RSEARCHEMIT;

	/* rebalance CM */
	if(! esl_opt_GetBoolean(go, "--nobalance"))
	{
		CM_t *new = NULL;
		if((status = CMRebalance(cm, errbuf, &new)) != eslOK) return status;
		FreeCM(cm);
		cm = new;
	}
	pretend_cm_is_hmm = ((CMCountNodetype(cm, MATP_nd) > 0)   || 
		esl_opt_GetBoolean(go, "--noh3pri")  || 
		esl_opt_GetBoolean(go, "--v1p0")) ? FALSE : TRUE;

	/* get counts */
	ESL_ALLOC(tr, sizeof(Parsetree_t *) * (msa->nseq));
	/* define the used_el array as all FALSE values, this means
	* Transmogrify will never create a parsetree with an EL, even if
	* the alignment seems to imply it. This is not ideal, but
	* currently necessary. If we wanted to allow ELs that existed in
	* the MSA (e.g. if the alignment were generated by cmalign) we'd
	* probably want to ensure that no EL columns in msa->rf get
	* defined as match columns, but that requires a significant change
	* to how match/insert columns are defined that I'm not ready to
	* make right now (very close to 1.1 release). We'll be forced to
	* revisit this when/if a jackhmmer analog in infernal is
	* ever implemented - in that case we will want to handle EL
	* emissions (more) correctly.
	*/
	ESL_ALLOC(used_el, sizeof(int) * (msa->alen+1));
	esl_vec_ISet(used_el, msa->alen+1, FALSE);
	for (idx = 0; idx < msa->nseq; idx++) {
		if((status = Transmogrify(cm, errbuf, mtr, msa->ax[idx], used_el, msa->alen, &(tr[idx]))) != eslOK) return status;
		if(pretend_cm_is_hmm) { 
			if((status = cm_parsetree_Doctor(cm, errbuf, tr[idx], NULL, NULL)) != eslOK) return status;
		}
		ParsetreeCount(cm, tr[idx], msa->ax[idx], msa->wgt[idx]);
		/*ParsetreeDump(cfg->ofp, tr[idx], cm, msa->ax[idx]);*/
	}
	free(used_el);
	cm->nseq     = msa->nseq;
	cm->eff_nseq = msa->nseq;

	if(cfg->be_verbose && do_print) { 
		fprintf(cfg->ofp, "done.  ");
		esl_stopwatch_Stop(w);
		esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
	}

	/* Set transition counts into ROOT_IL and ROOT_IR to 0, we don't
	* learn those counts from the alignment, unless --v1p0 (b/c we used
	* to in versions up to v1.0.2) or --iflank (which turns this
	* specific behavior off). The emission scores for these states will
	* be zeroed later so we don't touch them.
	*/
	if((! esl_opt_GetBoolean(go, "--v1p0")) && (! esl_opt_GetBoolean(go, "--iflank"))) { 
		if((status = cm_zero_flanking_insert_counts(cm, errbuf)) != eslOK) return status;
	}

	/* ensure the dual insert states we will detach were populated with 0 counts */
	if(!(esl_opt_GetBoolean(go, "--nodetach")))
	{
		if(cfg->be_verbose && do_print) { 
			esl_stopwatch_Start(w);
			fprintf(cfg->ofp, "%-40s ... ", "Finding and checking dual inserts");
		}
		cm_find_and_detach_dual_inserts(cm, 
			TRUE,   /* Do check (END_E-1) insert states have 0 counts */
			FALSE); /* Don't detach the states yet, wait til CM is priorified */
		if (cfg->be_verbose && do_print) {
			fprintf(cfg->ofp, "done.  ");
			esl_stopwatch_Stop(w);
			esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
		}
	}

	/* create the emitmap */
	if(cm->emap == NULL) cm->emap = CreateEmitMap(cm);

	/* set the EL self transition probability */
	cm->el_selfsc = sreLOG2(esl_opt_GetReal(go, "--elself"));

	/* set the beta parameters, these will be used to calculate W and QDBs that get stored in the CM file */
	cm->beta_W         = esl_opt_GetReal(go, "--betaW");
	cm->qdbinfo->beta1 = esl_opt_GetReal(go, "--beta1");
	cm->qdbinfo->beta2 = esl_opt_GetReal(go, "--beta2");

	/* set the cm->null2_omega and cm->null3_omega parameters */
	if(esl_opt_IsUsed(go, "--n2omega")) { /* user set --n2omega, use that */
		cm->null2_omega = esl_opt_GetReal(go, "--n2omega");
	}
	else { /* user didn't set --n2omega, definition of cm->null2_omega depends on whether --p56 was set or not */
		cm->null2_omega = ((esl_opt_GetBoolean(go, "--p56") == TRUE) ? V1P0_NULL2_OMEGA : esl_opt_GetReal(go, "--n2omega"));
	}
	if(esl_opt_IsUsed(go, "--n3omega")) { /* user set --n3omega, use that */
		cm->null3_omega = esl_opt_GetReal(go, "--n3omega");
	}
	else { /* user didn't set --n3omega, definition of cm->null3_omega depends on whether --p56 was set or not */
		cm->null3_omega = ((esl_opt_GetBoolean(go, "--p56") == TRUE) ? V1P0_NULL3_OMEGA : esl_opt_GetReal(go, "--n3omega"));
	}

	/* Before converting to probabilities, save a count vector file, if asked.
	* Used primarily for making data files for training priors.
	*/
	if (cfg->cfp != NULL) { 
		if ((status = print_countvectors(cfg, errbuf, cm)) != eslOK) goto ERROR;
	}

	*ret_cm  = cm;
	if(ret_mtr == NULL) FreeParsetree(mtr);
	else *ret_mtr = mtr;
	if(ret_msa_tr == NULL) {
		for(idx = 0; idx < msa->nseq; idx++)
			FreeParsetree(tr[idx]);
		free(tr);
		tr = NULL;
	}
	else *ret_msa_tr = tr;

	if(w != NULL) esl_stopwatch_Destroy(w);

	return eslOK;

ERROR:
	return status;
}

int
annotate(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm)
{
	int status = eslOK;
	ESL_STOPWATCH *w = NULL;

	if (cfg->be_verbose) {
		w = esl_stopwatch_Create();
		esl_stopwatch_Start(w);
		fprintf(cfg->ofp, "%-40s ... ", "Transferring MSA annotation");
		fflush(cfg->ofp);
	}

	if ((status = cm_SetName           (cm, msa->name))                    != eslOK)  ESL_XFAIL(status, errbuf, "Unable to set name for CM");
	if ((status = cm_SetAccession      (cm, msa->acc))                     != eslOK)  ESL_XFAIL(status, errbuf, "Failed to record MSA accession");
	if ((status = cm_SetDescription    (cm, msa->desc))                    != eslOK)  ESL_XFAIL(status, errbuf, "Failed to record MSA description");
	if ((status = cm_AppendComlog      (cm, go->argc, go->argv, FALSE, 0)) != eslOK)  ESL_XFAIL(status, errbuf, "Failed to record command log");
	if ((status = cm_SetCtime          (cm))                               != eslOK)  ESL_XFAIL(status, errbuf, "Failed to record timestamp");

	if (cfg->be_verbose) { 
		fprintf(cfg->ofp, "done.  ");
		esl_stopwatch_Stop(w);
		esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
	}

	if(w != NULL) esl_stopwatch_Destroy(w);
	return eslOK;

ERROR:
	if (cfg->be_verbose) { 
		fprintf(cfg->ofp, "FAILED.  ");
		esl_stopwatch_Stop(w);
		esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
	}
	if(w != NULL) esl_stopwatch_Destroy(w);
	return status;
}

int
set_model_cutoffs(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm)
{
	int cutoff_was_set = FALSE;
	if(msa->cutset[eslMSA_TC1]) { 
		cm->tc = msa->cutoff[eslMSA_TC1];
		cm->flags |= CMH_TC;
		cutoff_was_set = TRUE;
	}
	if(msa->cutset[eslMSA_GA1]) { 
		cm->ga = msa->cutoff[eslMSA_GA1];
		cm->flags |= CMH_GA;
		cutoff_was_set = TRUE;
	}
	if(msa->cutset[eslMSA_NC1]) { 
		cm->nc = msa->cutoff[eslMSA_NC1];
		cm->flags |= CMH_NC;
		cutoff_was_set = TRUE;
	}
	return eslOK;
}

int
set_effective_seqnumber(const ESL_GETOPTS *go, const struct cfg_s *cfg,
						char *errbuf, ESL_MSA *msa, CM_t *cm, const Prior_t *pri)
{
	int status;
	double neff;
	int used_hmm_etarget = FALSE;
	ESL_STOPWATCH *w = NULL;

	if(cfg->be_verbose) { 
		w = esl_stopwatch_Create();
		esl_stopwatch_Start(w);
		fprintf(cfg->ofp, "%-40s ... ", "Set effective sequence number");
		fflush(cfg->ofp);
	}

	if((esl_opt_GetBoolean(go, "--enone")) || ( esl_opt_IsOn(go, "--rsearch")))
	{
		neff = msa->nseq;
		if(cfg->be_verbose) fprintf(cfg->ofp, "done.  ");
	}
	else if(esl_opt_IsOn(go, "--eset")) 
	{
		neff = esl_opt_GetReal(go, "--eset");
		if(cfg->be_verbose) fprintf(cfg->ofp, "done.  ");
		cm->eff_nseq = neff;
		cm_Rescale(cm, neff / (float) msa->nseq);
	}
	else if (esl_opt_GetBoolean(go, "--eent") == TRUE)
	{
		double etarget; 
		double hmm_etarget; 
		double hmm_re;
		int clen = 0;
		int nd;
		for(nd = 0; nd < cm->nodes; nd++) { 
			if(cm->ndtype[nd] == MATP_nd) clen += 2;
			else if(cm->ndtype[nd] == MATL_nd) clen += 1;
			else if(cm->ndtype[nd] == MATR_nd) clen += 1;
		}
		if(esl_opt_GetBoolean(go, "--v1p0")) { 
			/* determine etarget with default method used by Infernal version 1.0-->1.0.2 */
			etarget = version_1p0_default_target_relent(cm->abc, clen, 6.0);
		}
		else { 
			etarget = set_target_relent(go, cm->abc, clen, CMCountNodetype(cm, MATP_nd));
		}

		status = cm_EntropyWeight(cm, pri, etarget, esl_opt_GetReal(go, "--eminseq"), FALSE, &hmm_re, &neff);
		/* if --ehmmre <x> enabled, ensure HMM relative entropy per match column is at least <x>, if not,
		* recalculate neff so HMM relative entropy of <x> is achieved.
		*/
		if( esl_opt_IsOn(go, "--ehmmre")) { 
			hmm_etarget = esl_opt_GetReal(go, "--ehmmre"); 
			if(hmm_re < hmm_etarget) { 
				status = cm_EntropyWeight(cm, pri, hmm_etarget, esl_opt_GetReal(go, "--eminseq"), TRUE, &hmm_re, &neff); /* TRUE says: pretend model is an HMM for entropy weighting */
				if      (status == eslEMEM) ESL_FAIL(status, errbuf, "memory allocation failed");
				else if (status != eslOK)   ESL_FAIL(status, errbuf, "internal failure in entropy weighting algorithm");
				used_hmm_etarget = TRUE;
			}
		}
		if      (status == eslEMEM) ESL_FAIL(status, errbuf, "memory allocation failed");
		else if (status != eslOK)   ESL_FAIL(status, errbuf, "internal failure in entropy weighting algorithm");
		cm->eff_nseq = neff;
		cm_Rescale(cm, neff / (float) msa->nseq);

		if(cfg->be_verbose) { 
			if(used_hmm_etarget) fprintf(cfg->ofp, "done.  ");
			else                 fprintf(cfg->ofp, "done.  ");
			esl_stopwatch_Stop(w);
			esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
		}
	}
	if(w != NULL) esl_stopwatch_Destroy(w);

	return eslOK;
}

int
parameterize(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, int do_print, CM_t *cm, const Prior_t *prior, float msa_nseq)
{
	int status; 
	ESL_STOPWATCH *w = NULL;

	if (cfg->be_verbose && do_print){
		w = esl_stopwatch_Create();
		esl_stopwatch_Start(w);
		fprintf(cfg->ofp, "%-40s ... ", "Converting counts to probabilities"); 
		fflush(cfg->ofp);
	}
	PriorifyCM(cm, prior); 

	if( esl_opt_IsOn(go, "--rsearch")) {
		rsearch_CMProbifyEmissions(cm, cfg->fullmat); /* use those probs to set CM probs from cts */
		/*debug_print_cm_params(cm);*/
	}

	if(! esl_opt_GetBoolean(go, "--nodetach")) /* Detach dual inserts where appropriate, if
											   * we get here we've already checked these states */
	{
		cm_find_and_detach_dual_inserts(cm, 
			FALSE, /* Don't check states have 0 counts (they won't due to priors) */
			TRUE); /* Detach the states by setting trans probs into them as 0.0   */
	}

	if(! esl_opt_GetBoolean(go, "--iins")) { 
		/* set all insert emission probabilities equal to the cm->null probabilities */ 
		if((status = flatten_insert_emissions(cm)) != eslOK) ESL_FAIL(status, errbuf, "flatten_insert_emissions() failed");
		/* Note: flatten_insert_emissions() is purposefully a static
		* function local to cmbuild.c b/c once CM files are calibrated no
		* other executable (i.e. cmsearch) should be able to modify the
		* scores of the CM, as that would invalidate the E-value stats */
	}

	CMRenormalize(cm);
	/* don't CMLogoddsify() here, that will come when we configure with cm_Configure() */

	if (cfg->be_verbose && do_print) { 
		fprintf(cfg->ofp, "done.  ");
		esl_stopwatch_Stop(w);
		esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
	}
	if(w != NULL) esl_stopwatch_Destroy(w);

	return eslOK;
}

/* configure_model()
* Configure the model. This determines QDBs and W.
* If niter is 1, we possibly output in verbose mode,
* else we don't.
*/
int
configure_model(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm, int iter)
{
	int status; 
	ESL_STOPWATCH *w = NULL;
	int nstarts, nexits, nd;

	if (iter == 1 && cfg->be_verbose){
		w = esl_stopwatch_Create();
		esl_stopwatch_Start(w);
		fprintf(cfg->ofp, "%-40s ... ", "Configuring model"); 
		fflush(cfg->ofp);
	}

	/* potentially redefine pbegin, pend based on command line options */
	/* (defaults are DEFAULT_PBEGIN, DEFAULT_PEND (set in CreateCMShell()) */
	if(esl_opt_IsUsed(go, "--pbegin")) cm->pbegin = esl_opt_GetReal(go, "--pbegin"); 
	if(esl_opt_IsUsed(go, "--pend"))   cm->pend   = esl_opt_GetReal(go, "--pend"); 

	/* possibly overwrite local begin probs such that all begin points are equiprobable (--pebegin) */
	if(esl_opt_GetBoolean(go, "--pebegin")) {
		nstarts = 0;
		for (nd = 2; nd < cm->nodes; nd++) 
			if (cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd || cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BIF_nd) 
				nstarts++;
		cm->pbegin = 1.- (1./(1+nstarts));
	}
	/* possibly overwrite cm->pend so that local end prob from all legal states is fixed,
	* this is strange in that cm->pend may be placed as a number greater than 1., this number
	* is then divided by nexits in ConfigLocalEnds() to get the prob for each v --> EL transition,
	* this is guaranteed by the way we calculate it to be < 1.,  it's the argument from --pfend */
	if(esl_opt_IsOn(go, "--pfend")) {
		nexits = 0;
		for (nd = 1; nd < cm->nodes; nd++) {
			if ((cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd ||
				cm->ndtype[nd] == MATR_nd || cm->ndtype[nd] == BEGL_nd ||
				cm->ndtype[nd] == BEGR_nd) && 
				cm->ndtype[nd+1] != END_nd)
				nexits++;
		}
		cm->pend = nexits * esl_opt_GetReal(go, "--pfend");
	}

	/* we must calculate QDBs so we can write them to the CM file */
	cm->config_opts |= CM_CONFIG_QDB;   

	/* if --refine, we have to set additional flags and configuration options
	* before calling cm_Configure().
	*/
	if(esl_opt_IsUsed(go, "--refine")) { 
		cm->tau = esl_opt_GetReal(go, "--tau");  /* this will be DEFAULT_TAU unless changed at command line */

		/* update cm->align->opts */
		if     (esl_opt_GetBoolean(go, "--gibbs"))     { cm->align_opts |= CM_ALIGN_SAMPLE; }
		else if(esl_opt_GetBoolean(go, "--cyk"))       { cm->align_opts |= CM_ALIGN_CYK;    }
		else                                           { cm->align_opts |= CM_ALIGN_OPTACC; }
		if     (  esl_opt_GetBoolean(go, "--sub"))     { cm->align_opts |= CM_ALIGN_SUB;    }
		else if(! esl_opt_GetBoolean(go, "--notrunc")) { cm->align_opts |= CM_ALIGN_TRUNC;  }

		if(esl_opt_GetBoolean(go, "--nonbanded"))   { 
			cm->align_opts |=  CM_ALIGN_SMALL; 
			cm->align_opts |=  CM_ALIGN_NONBANDED; 
			cm->align_opts |=  CM_ALIGN_CYK;
			cm->align_opts &= ~CM_ALIGN_OPTACC; /* turn optimal accuracy OFF */
		}
		else cm->align_opts  |= CM_ALIGN_HBANDED;

		if(esl_opt_GetBoolean(go, "--fins")) cm->align_opts  |= CM_ALIGN_FLUSHINSERTS;

		/* update cm->config_opts */
		if     (  esl_opt_GetBoolean(go, "--sub"))     { cm->config_opts |= CM_CONFIG_SUB; }
		else if(! esl_opt_GetBoolean(go, "--notrunc")) { cm->config_opts |= CM_CONFIG_TRUNC; }
		if(esl_opt_GetBoolean(go, "--nonbanded")) cm->config_opts |= CM_CONFIG_NONBANDEDMX;

		if(esl_opt_GetBoolean(go, "-l")) { 
			cm->config_opts |= CM_CONFIG_LOCAL;
			cm->config_opts |= CM_CONFIG_HMMLOCAL;
			cm->config_opts |= CM_CONFIG_HMMEL;
		}
	}

	/* finally, configure the model */
	if((status = cm_Configure(cm, errbuf, -1)) != eslOK) return status;

	if (iter == 1 && cfg->be_verbose) { 
		fprintf(cfg->ofp, "done.  ");
		esl_stopwatch_Stop(w);
		esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
	}

	if(w != NULL) esl_stopwatch_Destroy(w);

	return eslOK;
}


/* set_consensus()
* Set CM consensus using cm->cmcons. We have to do this
* after configuring the model because bit scores must
* be valid.
*/
int
set_consensus(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
	int status = eslOK;
	ESL_STOPWATCH *w = NULL;

	if (cfg->be_verbose) {
		w = esl_stopwatch_Create();
		esl_stopwatch_Start(w);
		fprintf(cfg->ofp, "%-40s ... ", "Setting CM consensus");
		fflush(cfg->ofp);
	}

	if(! (cm->flags & CMH_BITS)) ESL_XFAIL(eslEINVAL, errbuf, "Trying to set cm->consensus before bit scores are valid");

	if ((status = cm_SetConsensus  (cm, cm->cmcons, NULL)) != eslOK) ESL_XFAIL(status, errbuf, "Failed to calculate consensus sequence");

	if (cfg->be_verbose) { 
		fprintf(cfg->ofp, "done.  ");
		esl_stopwatch_Stop(w);
		esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
	}

	if(w != NULL) esl_stopwatch_Destroy(w);
	return eslOK;

ERROR:
	if (cfg->be_verbose) { 
		fprintf(cfg->ofp, "FAILED.  ");
		esl_stopwatch_Stop(w);
		esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
	}
	if(w != NULL) esl_stopwatch_Destroy(w);
	return status;
}

/* build_and_calibrate_p7_filter()
* Build and calibrate the additional p7 HMM filter (differs from the
* ML p7 HMM filter which was built in configure_model()).  
* If <use_mlp7_as_filter>, do just that. This will be true
* if --p7ml was used OR if model has zero basepairs and 
* --noh3pri was NOT used.
*/
int
build_and_calibrate_p7_filter(const ESL_GETOPTS *go, const struct cfg_s *cfg, char *errbuf, ESL_MSA *msa, CM_t *cm, int use_mlp7_as_filter,char *hmmStatsToCalibrate)
{
	int status; 
	CM_t    *acm = NULL; 
	P7_HMM *fhmm = NULL;
	int lmsvL, lvitL, lfwdL, gfwdL;
	int lmsvN, lvitN, lfwdN, gfwdN;
	double agfmu, agflambda;
	float lftailp, gftailp;
	int k, apos, cpos;
	ESL_MSA *amsa = NULL;
	double mlp7_re, fhmm_re;
	double neff;
	ESL_STOPWATCH *w = NULL;

	if (cfg->be_verbose){
		w = esl_stopwatch_Create();
		esl_stopwatch_Start(w);
		fprintf(cfg->ofp, "%-40s ... ", "Calibrating p7 HMM filter"); 
		fflush(cfg->ofp);
	}

	/* calibrate p7 HMM */

	/* Define relevant parameters: */
	/* first, set all to default */
	lmsvL = lvitL = 200;
	lfwdL = 100;
	gfwdL = ESL_MAX(100, 2.*cm->clen);
	lmsvN = esl_opt_GetInteger(go, "--EmN");
	lvitN = esl_opt_GetInteger(go, "--EvN");
	lfwdN = esl_opt_GetInteger(go, "--ElfN");
	gfwdN = esl_opt_GetInteger(go, "--EgfN");
	lftailp = esl_opt_GetReal(go, "--Elftp");
	gftailp = esl_opt_GetReal(go, "--Egftp");

	/* now, modify if nec based on command-line options */
	if(esl_opt_IsUsed(go, "--ElL")) { 
		lmsvL = lvitL = lfwdL = esl_opt_GetInteger(go, "--ElL");
	}
	else if(esl_opt_IsUsed(go, "--Elcmult")) { 
		lmsvL = lvitL = lfwdL = esl_opt_GetReal(go, "--Elcmult") * cm->clen;
	}

	if(esl_opt_IsUsed(go, "--EgL")) { 
		gfwdL = esl_opt_GetInteger(go, "--EgL");
	}
	else if(esl_opt_IsUsed(go, "--Egcmult")) { 
		gfwdL = esl_opt_GetReal(go, "--Egcmult") * cm->clen;
	}

	if(esl_opt_GetBoolean(go, "--Efitlam")) { 
		lmsvN = esl_opt_IsUsed(go, "--EmN") ? esl_opt_GetInteger(go, "--EmN") : 10000;
		lvitN = esl_opt_IsUsed(go, "--EvN") ? esl_opt_GetInteger(go, "--EvN") : 10000;
		lfwdN = esl_opt_IsUsed(go, "--ElfN") ? esl_opt_GetInteger(go, "--ElfN") : 10000;
		gfwdN = esl_opt_IsUsed(go, "--EgfN") ? esl_opt_GetInteger(go, "--EgfN") : 10000;
		lftailp = esl_opt_IsUsed(go, "--Elftp") ? esl_opt_GetReal(go, "--Elftp") : 0.01;
		gftailp = esl_opt_IsUsed(go, "--Egftp") ? esl_opt_GetReal(go, "--Egftp") : 0.01;
	}

	/* Build the HMM filter (cm->fp7) */
	if(use_mlp7_as_filter) { 
		/* use the ML p7 HMM as the HMM filter */
		fhmm = cm->mlp7;
	}
	else { 
		/* Default strategy: 
		* 
		* Build a p7 HMM <fhmm> with p7_builder with entropy weighting
		* and rel ent target of <x> from (--p7ere <x>). Then _overwrite_
		* it's emission probabilities with marginalized emission
		* probabilities from a temporary CM built such that it's ML p7
		* HMM has mean match state entropy the same as <fhmm>.
		*
		* The motivation for this is that rmark3 benchmarking (xref:
		* ~nawrockie/notebook/11_0226_inf_p7s_with_cp9_transitions/00LOG)
		* revealed that HMMs built with H3 transitions and emissions from
		* a marginalized CM are the best for filtering.
		* 
		* The --p7hemit option makes it so HMMER emissions are used instead
		* of the marginalized CM emissions.
		*
		* NOTE: We could avoid this if we had a way or using Infernal's
		* emission priors (there are different prior for base pairs and
		* singlets) in p7_Builder(), but that would require parsing the
		* secondary structure and getting an Infernal function into
		* hmmer). For now, we build a temporary CM and copy it's
		* emissions.
		*
		* First, annotate the msa with RF annotation corresponding to the
		* definition of match columns used by the CM using the cm->map.
		* This will guarantree that the HMM has the same number of match
		* columns as the CM, which is important because we copy
		* marginalized ml emissions from a CM onto the HMM. Note that we
		* also set cfg->fp7_bld->arch_strategy as p7_ARCH_HAND.
		*/
		amsa = esl_msa_Clone(msa);
		if(amsa->rf != NULL) free(amsa->rf);
		ESL_ALLOC(amsa->rf, sizeof(char) * (amsa->alen+1));
		if(! (cm->flags & CMH_MAP)) { cm_Fail("Unable to create additional p7 HMM, CM has no map, this shouldn't happen"); }
		/* init to all inserts, then set match states based on cm->map */
		for (apos = 0; apos <  amsa->alen; apos++) amsa->rf[apos] = '.';
		for (cpos = 1; cpos <= cm->clen;   cpos++) amsa->rf[cm->map[cpos]-1] = 'x'; /* note off by one */
		cfg->fp7_bld->arch_strategy = p7_ARCH_HAND;

		if ((status = p7_Builder(cfg->fp7_bld, amsa, cfg->fp7_bg, &fhmm, NULL, NULL, NULL, NULL)) != eslOK) { strcpy(errbuf, cfg->fp7_bld->errbuf); return status; }
		/* remove the RF annotation, it only exists because we created amsa->rf above */
		if(fhmm->rf != NULL) { 
			free(fhmm->rf);
			fhmm->rf = NULL;
			fhmm->flags &= ~p7H_RF;
		}      
		if(cm->flags & CMH_RF && cm->rf != NULL) { /* copy CM's rf annotation to fhmm, remember they have same # consensus columns */
			ESL_ALLOC(fhmm->rf, sizeof(char) * (cm->clen+2));
			strcpy(fhmm->rf, cm->rf);
			fhmm->flags |= p7H_RF;
		}
		/* overwrite the HMM consensus structure annotation with the CM's it'll be in full WUSS format */
		if(! (fhmm->flags & p7H_CS)) { cm_Fail("additional p7 HMM unexpectedly does not have consensus structure annotation"); }
		fhmm->cs[0] = ' ';
		strcpy(fhmm->cs+1, cm->cmcons->cstr); /* careful: off-by-one */
		fhmm->cs[cm->clen+1] = '\0';

		esl_msa_Destroy(amsa); 

		if(! esl_opt_GetBoolean(go, "--p7hemit")) { 
			/* Overwrite the emission probabilities of the HMM with
			* emissions from a ML HMM built from a CM.  First, build the
			* CM, it will have the same HMM mean match state entropy as the
			* fhmm we just built.
			*/
			if ((status =  build_model(go, cfg, errbuf, FALSE, msa, &acm, NULL, NULL)) != eslOK) return status;
			fhmm_re = p7_MeanMatchRelativeEntropy(fhmm, cfg->fp7_bg);
			status = cm_EntropyWeight(acm, cfg->pri, fhmm_re, esl_opt_GetReal(go, "--eminseq"), TRUE, &mlp7_re, &neff); /* TRUE says: pretend model is an HMM for entropy weighting */
			if      (status == eslEMEM) ESL_FAIL(status, errbuf, "memory allocation failed");
			else if (status != eslOK)   ESL_FAIL(status, errbuf, "internal failure in entropy weighting algorithm");
			acm->eff_nseq = neff;
			cm_Rescale(acm, acm->eff_nseq / (float) msa->nseq);
			if((status = parameterize   (go, cfg, errbuf, FALSE, acm, cfg->pri, msa->nseq)) != eslOK) return status;
			/* We have to configure the model to get cm->W, which gets 
			* copied to cm->mlp7->max_length. Alternatively we could 
			* use p7_Builder_MaxLength() but anecdotally that gives 
			* lengths >> W (more than 2*W commonly).
			* configure_model() will build the mlp7 HMM.
			*/
			if((status = configure_model(go, cfg, errbuf, acm, 2)) != eslOK) return status;

			/* copy the ML p7 emission probs from the CM we just built */
			/* match emissions: copy, then normalize (norm should be unnec actually) */
			for (k = 1; k <= fhmm->M; k++) esl_vec_FCopy(acm->mlp7->mat[k], fhmm->abc->K, fhmm->mat[k]);
			for (k = 1; k <= fhmm->M; k++) esl_vec_FNorm(fhmm->mat[k], fhmm->abc->K);
			/* special case */
			esl_vec_FSet(fhmm->mat[0], fhmm->abc->K, 0.);
			fhmm->mat[0][0] = 1.0;

			/* insert emissions: copy, then normalize (norm should be unnec actually) */
			for (k = 0; k <= fhmm->M; k++) esl_vec_FCopy(acm->mlp7->ins[k], fhmm->abc->K, fhmm->ins[k]);
			for (k = 0; k <= fhmm->M; k++) esl_vec_FNorm(fhmm->ins[k], fhmm->abc->K);
			/* reset HMM composition */
			if ((status = p7_hmm_SetComposition(fhmm)) != eslOK) goto ERROR;
			fhmm->eff_nseq = acm->eff_nseq;

			FreeCM(acm);
		}
	}

	/* calibrate the HMM filter */
	if((status = cm_p7_Calibrate(fhmm, errbuf, 
		lmsvL, lvitL, lfwdL, gfwdL,                 /* length of sequences to search for local (lL) and glocal (gL) modes */    
		lmsvN, lvitN, lfwdN, gfwdN,                 /* number of seqs to search for each alg */
		lftailp,                                    /* fraction of tail mass to fit for local Fwd */
		gftailp,                                    /* fraction of tail mass to fit for glocal Fwd */
		&agfmu, &agflambda,
		hmmStatsToCalibrate))
		!= eslOK) ESL_FAIL(status, errbuf, "Error calibrating additional p7 HMM");

	if((status = cm_p7_hmm_SetConsensus(fhmm)) != eslOK) ESL_FAIL(status, errbuf, "Unable to set the HMM filter consensus annotation");
	if((status = cm_SetFilterHMM(cm, fhmm, agfmu, agflambda))       != eslOK) ESL_FAIL(status, errbuf, "Unable to set the HMM filter for the CM");
	if((status = p7_hmm_AppendComlog (cm->fp7, go->argc, go->argv)) != eslOK) ESL_FAIL(status, errbuf, "Failed to record command log for filter HMM");

	if (cfg->be_verbose) { 
		fprintf(cfg->ofp, "done.  ");
		esl_stopwatch_Stop(w);
		esl_stopwatch_Display(cfg->ofp, w, "CPU time: ");
	}
	if(w != NULL) esl_stopwatch_Destroy(w);

	return eslOK;

ERROR: 
	ESL_FAIL(status, errbuf, "out of memory");
	return status; /* never reached */
}

int
print_countvectors(const struct cfg_s *cfg, char *errbuf, CM_t *cm)
{
	int   v,x;

	if(cfg->cfp == NULL) ESL_FAIL(eslEINCOMPAT, errbuf, "save_countvectors(), but cfg->cfp is NULL, shouldn't happen.");

	/* Print emission counts */
	for (v = 0; v < cm->M; v++) {
		if (cm->sttype[v] == MP_st || cm->sttype[v] == ML_st || cm->sttype[v] == MR_st) { 
			fprintf(cfg->cfp, "E\t%-7s ", UniqueStatetype(cm->stid[v]));
			if (cm->sttype[v] == MP_st) {
				for (x = 0; x < cm->abc->K*cm->abc->K; x++)
					fprintf(cfg->cfp, "%8.3f ", cm->e[v][x]);
			} else {
				for (x = 0; x < cm->abc->K; x++)
					fprintf(cfg->cfp, "%8.3f ", cm->e[v][x]);
			}
			fprintf(cfg->cfp, "\n");
		}
	}

	/* Print transition counts */
	for (v = 0; v < cm->M; v++) {
		if(cm->sttype[v] != B_st && cm->sttype[v] != E_st) {
			fprintf(cfg->cfp, "T\t%-7s : %-2d", UniqueStatetype(cm->stid[v]), cm->ndtype[(cm->ndidx[v] + 1)]);
			for (x = 0; x < cm->cnum[v]; x++) {
				fprintf(cfg->cfp, "%8.3f ", cm->t[v][x]);
			}
			fprintf(cfg->cfp, "\n");
		}
	}
	fprintf(cfg->cfp, "//\n");
	return eslOK;
}

double
set_target_relent(const ESL_GETOPTS *go, const ESL_ALPHABET *abc, int clen, int nbps)
{
	double etarget;
	double re_target;
	double esigma = esl_opt_GetReal(go, "--esigma"); /* default Infernal/HMMER3 sigma is 45.0 */

	if(esl_opt_IsOn(go, "--ere")) { 
		re_target = esl_opt_GetReal(go, "--ere");
	}
	else {
		if(abc->type != eslRNA) cm_Fail("ERROR, alphabet not RNA, user needs to specify target entropy with --ere");
		/* set target differently if we have 0 basepairs or not */
		re_target = (nbps > 0) ? DEFAULT_ETARGET : DEFAULT_ETARGET_HMMFILTER;
	}
	/* the defn of etarget below is identical to how hmmer3 does it in hmmer/src/p7_builder.c as of svn rev 3986 (04.16.12) */
	etarget = (esigma - eslCONST_LOG2R * log( 2.0 / ((double) clen * (double) (clen+1)))) / (double) clen; /* HMMER3.0 default, xref J5/36. */
	etarget = ESL_MAX(etarget, re_target);

	return etarget;
}

double
version_1p0_default_target_relent(const ESL_ALPHABET *abc, int clen, double eX)
{
	double etarget;

	/* HMMER3 default eX = 6.0 as of Tue Jul 10 2007
	*/
	etarget = 6.* (eX + log((double) ((clen * (clen+1)) / 2)) / log(2.))    / (double)(2*clen + 4);

	switch (abc->type) {
case eslRNA:    if (etarget < DEFAULT_ETARGET)   etarget = DEFAULT_ETARGET;   break;
default:        cm_Fail("ERROR in default_target_relent(), alphabet not RNA!\n");
	}
	return etarget;
}

int
flatten_insert_emissions(CM_t *cm)
{
	int v;

	/* Contract check */
	if(cm->abc  == NULL) cm_Fail("flatten_insert_emissions(), cm->abc is NULL.\n");
	if(cm->null == NULL) cm_Fail("flatten_insert_emissions(), cm->null is NULL.\n");

	esl_vec_FNorm(cm->null, cm->abc->K);
	for (v = 0; v < cm->M; v++) {
		if(cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
			esl_vec_FSet(cm->e[v], (cm->abc->K * cm->abc->K), 0.); /* zero them out */
			esl_vec_FCopy(cm->null, cm->abc->K, cm->e[v]); /* overwrite first cm->abc->K values (rest are irrelevant for non-MP states) with cm->null */
		}
	}
	return eslOK;
}

P7_PRIOR *
p7_prior_Read(FILE *fp) 
{
	P7_PRIOR *pri = NULL;
	int        status;
	ESL_FILEPARSER *efp = NULL;

	ESL_ALLOC(pri, sizeof(P7_PRIOR));
	pri->tm = pri->ti = pri->td = pri->em = pri->ei = NULL;

	if ((efp = esl_fileparser_Create(fp)) == NULL) goto ERROR;
	esl_fileparser_SetCommentChar(efp, '#');

	/* Transition section: 3 sets of transitions out of match, out of insert, and out of delete */
	if (esl_mixdchlet_Read(efp, &(pri->tm)) != eslOK) goto ERROR;
	if (esl_mixdchlet_Read(efp, &(pri->ti)) != eslOK) goto ERROR;
	if (esl_mixdchlet_Read(efp, &(pri->td)) != eslOK) goto ERROR;

	/* Emission section: match emissions, then insert emissions */
	if (esl_mixdchlet_Read(efp, &(pri->em)) != eslOK) goto ERROR;
	if (esl_mixdchlet_Read(efp, &(pri->ei)) != eslOK) goto ERROR;

	esl_fileparser_Destroy(efp);

	return pri;

ERROR: 
	if(efp != NULL) esl_fileparser_Destroy(efp);
	if(pri != NULL) p7_prior_Destroy(pri);
	return NULL;
}



struct cfg_s InitCmbuildCfg (ESL_GETOPTS *go)
{
	struct cfg_s cfg;
	cfg.ofp        = NULL;	           
	cfg.fmt        = eslMSAFILE_UNKNOWN;     /* possibly reset below */
	cfg.afp        = NULL;	           /* created in init_cfg() */
	cfg.abc        = NULL;	           /* created in init_cfg() */
	cfg.cmoutfp    = NULL;	           /* opened in init_cfg() */
	cfg.postmsafile= esl_opt_GetString(go, "-O"); /* NULL by default */
	cfg.postmsafp  = NULL;                  
	cfg.null       = NULL;	           /* created in init_cfg() */
	cfg.pri        = NULL;                   /* created in init_cfg() */
	cfg.fullmat    = NULL;                   /* read (possibly) in init_cfg() */
	cfg.r          = NULL;	           /* created (possibly) in init_cfg() */
	cfg.fp7_bg     = NULL;                   /* created (possibly) in init_cfg() */
	cfg.fp7_bld    = NULL;                   /* created (possibly) in init_cfg() */
	/* optional output files, opened in init_cfg(), if at all */
	cfg.cfp        = NULL;
	cfg.escfp      = NULL;
	cfg.tblfp      = NULL;
	cfg.efp        = NULL;
	cfg.gfp        = NULL;
	cfg.gtblfp     = NULL;
	cfg.tfp        = NULL;
	cfg.cdfp       = NULL;
	cfg.refinefp   = NULL;
	cfg.rdfp       = NULL;
	cfg.alifile    = NULL;
	cfg.cmfile     = NULL;
	cfg.pri_zerobp = NULL;
	return cfg;
}

int init_cfg(const ESL_GETOPTS *go, struct cfg_s *cfg, char *errbuf)
{
	int status;

	/* Set the msafile alphabet as RNA, if it's DNA we're fine. 
	* If it's not RNA nor DNA, we can't deal with it anyway,
	* so we're hardcoded to RNA.
	*/
	if((cfg->abc = esl_alphabet_Create(eslRNA)) == NULL) ESL_FAIL(eslEMEM, errbuf, "Failed to create alphabet for sequence file");

	if (esl_opt_IsUsed(go, "-o")) { 
		cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w");
		if (cfg->ofp == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
	} 
	else cfg->ofp = stdout;

	if (cfg->postmsafile) { 
		cfg->postmsafp = fopen(cfg->postmsafile, "w");
		if (cfg->postmsafp == NULL) ESL_FAIL(eslFAIL, errbuf, "Failed to MSA resave file %s for writing", cfg->postmsafile);
	} 
	else cfg->postmsafp = NULL;

	/* Set up the priors */
	if (esl_opt_GetString(go, "--prior") != NULL) { 
		FILE *pfp;
		if ((pfp = fopen(esl_opt_GetString(go, "--prior"), "r")) == NULL) cm_Fail("Failed to open prior file %s\n", esl_opt_GetString(go, "--prior"));
		if ((cfg->pri = Prior_Read(pfp)) == NULL)       	                 cm_Fail("Failed to parse prior file %s\n", esl_opt_GetString(go, "--prior"));
		fclose(pfp);
		cfg->pri_zerobp = NULL;
	}
	else if(esl_opt_GetBoolean(go, "--p56") || esl_opt_GetBoolean(go, "--v1p0") || esl_opt_GetBoolean(go, "--noh3pri")) { 
		cfg->pri        = Prior_Default_v0p56_through_v1p02();
		cfg->pri_zerobp = NULL;
	}
	else { 
		cfg->pri        = Prior_Default(FALSE);
		cfg->pri_zerobp = Prior_Default(TRUE);
	}

	/* Set up the null/random seq model */
	if(esl_opt_GetString(go, "--null") != NULL) /* read freqs from a file and overwrite bg->f */
	{
		if((status = CMReadNullModel(cfg->abc, esl_opt_GetString(go, "--null"), &(cfg->null))) != eslOK)
			cm_Fail("Failure reading the null model, code: %d", status);
	}       
	else /* set up the default null model */
	{
		status = DefaultNullModel(cfg->abc, &(cfg->null)); /* default values, A,C,G,U = 0.25  */
		if(status != eslOK) cm_Fail("Failure creating the null model, code: %d", status);
	}

	/* if --rsearch was enabled, set up RIBOSUM matrix */
	if(esl_opt_GetString(go, "--rsearch") != NULL)
	{
		FILE *matfp;
		if ((matfp = MatFileOpen (esl_opt_GetString(go, "--rsearch"))) == NULL)
			cm_Fail("Failed to open matrix file %s\n", esl_opt_GetString(go, "--rsearch"));
		if (! (cfg->fullmat = ReadMatrix(cfg->abc, matfp)))
			cm_Fail("Failed to read matrix file %s\n", esl_opt_GetString(go, "--rsearch"));
		ribosum_calc_targets(cfg->fullmat); /* overwrite score matrix scores w/target probs */
		fclose(matfp);
	}

	/* if --corig enabled, make sure either --cmaxid, --ctarget, or --call also enabled */
	if (esl_opt_GetBoolean(go, "--corig"))
		if((! esl_opt_IsOn(go, "--ctarget")) && (! esl_opt_IsOn(go, "--cmaxid")) && (! esl_opt_IsOn(go, "--call")))
			cm_Fail("--corig only makes sense in combination with --ctarget, --cmaxid, OR --call");

	/* if --gibbs enabled, open output file for refined MSAs, and seed RNG */
	if(esl_opt_GetBoolean(go, "--gibbs"))
	{
		/* create RNG */
		cfg->r = esl_randomness_Create(esl_opt_GetInteger(go, "--seed"));
		if (cfg->r == NULL) ESL_FAIL(eslEINVAL, errbuf, "Failed to create random number generator: probably out of memory");
	}

	/* set up objects for building additional p7 models to filter with, if nec */
	cfg->fp7_bg = p7_bg_Create(cfg->abc);
	/* create the P7_BUILDER, pass NULL as the <go> argument, this sets all parameters to default */
	cfg->fp7_bld = p7_builder_Create(NULL, cfg->abc);
	cfg->fp7_bld->w_len = -1;
	cfg->fp7_bld->w_beta = p7_DEFAULT_WINDOW_BETA;
	if(esl_opt_IsUsed(go, "--p7prior")) {
		FILE *pfp;
		if (cfg->fp7_bld->prior != NULL) p7_prior_Destroy(cfg->fp7_bld->prior);
		if ((pfp = fopen(esl_opt_GetString(go, "--p7prior"), "r")) == NULL) cm_Fail("Failed to open p7 prior file %s\n", esl_opt_GetString(go, "--p7prior"));
		if((cfg->fp7_bld->prior = p7_prior_Read(pfp)) == NULL) {
			cm_Fail("Failed to parse p7 prior file %s\n", esl_opt_GetString(go, "--p7prior"));
		}
		fclose(pfp);
	}
	else if(! esl_opt_GetBoolean(go, "--p7hprior")) { 
		/* create the default Infernal p7 prior */
		if (cfg->fp7_bld->prior != NULL) p7_prior_Destroy(cfg->fp7_bld->prior);
		cfg->fp7_bld->prior = cm_p7_prior_CreateNucleic();
	}
	cfg->fp7_bld->re_target = esl_opt_IsOn(go, "--p7ere") ?  esl_opt_GetReal(go, "--p7ere") : DEFAULT_ETARGET_HMMFILTER;

	/* open output files */
	/* optionally, open count vector file */
	if (esl_opt_GetString(go, "--cfile") != NULL) {
		if ((cfg->cfp = fopen(esl_opt_GetString(go, "--cfile"), "w")) == NULL) 
			ESL_FAIL(eslFAIL, errbuf, "Failed to open --cfile output file %s\n", esl_opt_GetString(go, "--cfile"));
	}
	/* optionally, open base pair info file */
	if (esl_opt_GetString(go, "--efile") != NULL) {
		if ((cfg->escfp = fopen(esl_opt_GetString(go, "--efile"), "w")) == NULL) 
			ESL_FAIL(eslFAIL, errbuf, "Failed to open --efile output file %s\n", esl_opt_GetString(go, "--efile"));
	}
	/* optionally, open CM tabular file */
	if (esl_opt_GetString(go, "--cmtbl") != NULL) {
		if ((cfg->tblfp = fopen(esl_opt_GetString(go, "--cmtbl"), "w")) == NULL) 
			ESL_FAIL(eslFAIL, errbuf, "Failed to open --cmtbl output file %s\n", esl_opt_GetString(go, "--cmtbl"));
	}
	/* optionally, open emit map file */
	if (esl_opt_GetString(go, "--emap") != NULL) {
		if ((cfg->efp = fopen(esl_opt_GetString(go, "--emap"), "w")) == NULL) 
			ESL_FAIL(eslFAIL, errbuf, "Failed to open --emap output file %s\n", esl_opt_GetString(go, "--emap"));
	}
	/* optionally, open guide tree file */
	if (esl_opt_GetString(go, "--gtree") != NULL) {
		if ((cfg->gfp = fopen(esl_opt_GetString(go, "--gtree"), "w")) == NULL) 
			ESL_FAIL(eslFAIL, errbuf, "Failed to open --gtree output file %s\n", esl_opt_GetString(go, "--gtree"));
	}
	/* optionally, open master tree file */
	if (esl_opt_GetString(go, "--gtbl") != NULL) {
		if ((cfg->gtblfp = fopen(esl_opt_GetString(go, "--gtbl"), "w")) == NULL) 
			ESL_FAIL(eslFAIL, errbuf, "Failed to open --gtbl output file %s\n", esl_opt_GetString(go, "--gtbl"));
	}
	/* optionally, open trace file */
	if (esl_opt_GetString(go, "--tfile") != NULL) {
		if ((cfg->tfp = fopen(esl_opt_GetString(go, "--tfile"), "w")) == NULL) 
			ESL_FAIL(eslFAIL, errbuf, "Failed to open --tfile output file %s\n", esl_opt_GetString(go, "--tfile"));
	}
	/* if --refine enabled, open output file for refined MSAs */
	if (esl_opt_GetString(go, "--refine") != NULL)
	{
		if ((cfg->refinefp = fopen(esl_opt_GetString(go, "--refine"), "w")) == NULL)
			ESL_FAIL(eslFAIL, errbuf, "Failed to open output file %s for writing MSAs from --refine to", esl_opt_GetString(go, "--refine"));
	}
	/* optionally, open --rdump output alignment file */
	if (esl_opt_GetString(go, "--rdump") != NULL) {
		if ((cfg->rdfp = fopen(esl_opt_GetString(go, "--rdump"), "w")) == NULL) 
			ESL_FAIL(eslFAIL, errbuf, "Failed to open --rdump output file %s\n", esl_opt_GetString(go, "--rdump"));
	}
	/* if --cdump enabled, open output file for cluster MSAs */
	if (esl_opt_GetString(go, "--cdump") != NULL)
	{
		/* check to make sure there's a reason for this option, --cmaxid, --ctarget or --call MUST also be enabled */
		if((! esl_opt_IsOn(go, "--ctarget")) && (!esl_opt_IsOn(go, "--cmaxid")) && (!esl_opt_IsOn(go, "--call")))
			cm_Fail("--cdump only makes sense in combination with --ctarget, --cmaxid, OR --call");
		if ((cfg->cdfp = fopen(esl_opt_GetString(go, "--cdump"), "w")) == NULL)
			cm_Fail("Failed to open output file %s for writing MSAs to", esl_opt_GetString(go, "--cdump"));
	}

	if (cfg->pri     == NULL) ESL_FAIL(eslEINVAL, errbuf, "alphabet initialization failed");
	if (cfg->null    == NULL) ESL_FAIL(eslEINVAL, errbuf, "null model initialization failed");

	cfg->nali = 0;
	cfg->ncm_total = 0;
	return eslOK;
}

void destroy_cfg (struct cfg_s *cfg)
{
   if (cfg->postmsafp != NULL) {
     fclose(cfg->postmsafp); 
   }
   if (cfg->cfp != NULL) {
     fclose(cfg->cfp); 
   }
   if (cfg->escfp != NULL) {
     fclose(cfg->escfp); 
   }
   if (cfg->tblfp != NULL) {
     fclose(cfg->tblfp); 
   }
   if (cfg->efp != NULL) {
     fclose(cfg->efp); 
   }
   if (cfg->gfp != NULL) {
     fclose(cfg->gfp); 
   }
   if (cfg->gtblfp != NULL) {
     fclose(cfg->gtblfp); 
   }
   if (cfg->tfp != NULL) {
     fclose(cfg->tfp); 
   }
   if (cfg->cdfp != NULL) {
     fclose(cfg->cdfp); 
   }
   if (cfg->refinefp != NULL) {
     fclose(cfg->refinefp); 
   }
   if (cfg->rdfp != NULL) {
     fclose(cfg->rdfp); 
   }
   if (cfg->afp        != NULL) eslx_msafile_Close(cfg->afp);
   if (cfg->abc        != NULL) esl_alphabet_Destroy(cfg->abc);
   if (cfg->cmoutfp    != NULL) fclose(cfg->cmoutfp);
   if (cfg->pri        != NULL) Prior_Destroy(cfg->pri);
   if (cfg->pri_zerobp != NULL) Prior_Destroy(cfg->pri_zerobp);
   if (cfg->null       != NULL) {
	   free(cfg->null);
   }
   if (cfg->r          != NULL) esl_randomness_Destroy(cfg->r);
   if (cfg->fp7_bg     != NULL) p7_bg_Destroy(cfg->fp7_bg);
   if (cfg->fp7_bld    != NULL) p7_builder_Destroy(cfg->fp7_bld);
}



void cmbuild (const char *cmbuild_commandline_flags,char *alifile,ESL_MSA *msa,CM_t **ret_cm,char *hmmStatsToCalibrate)
{
	char     errbuf[eslERRBUFSIZE];
	int status;

	ESL_GETOPTS *go     = NULL;
	if ((go = esl_getopts_Create(cmbuild_options))   == NULL) { cm_Fail("Internal failure creating options object");}
	/*if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { esl_fatal("Failed to parse command line: %s\n", go->errbuf); }*/
	if (esl_opt_ProcessSpoof(go,cmbuild_commandline_flags) != eslOK) { esl_fatal("Failed to parse command line: %s\n", go->errbuf); }
	if (esl_opt_VerifyConfig(go)               != eslOK)  { esl_fatal("Failed to parse command line: %s\n", go->errbuf); }

	struct cfg_s cfg_var=InitCmbuildCfg(go);
	struct cfg_s *cfg=&cfg_var;
	cfg->be_verbose = esl_opt_GetBoolean(go, "--verbose");

	cfg->alifile=alifile;
	cfg->abc=msa->abc;

	if ((status = init_cfg(go, cfg, errbuf)) != eslOK) cm_Fail(errbuf);

	cfg->nali=1; /* the caller is responsible for ensuring that the input 'msa' is the only one in the file (BTW, if we fail to set this, then set_msa_name complains */

#if 0
	cfg->fp7_bg = p7_bg_Create(cfg->abc);
	/* create the P7_BUILDER, pass NULL as the <go> argument, this sets all parameters to default */
	cfg->fp7_bld = p7_builder_Create(NULL, cfg->abc);
	cfg->fp7_bld->w_len = -1;
	cfg->fp7_bld->w_beta = p7_DEFAULT_WINDOW_BETA;
	if (cfg->fp7_bld->prior != NULL) p7_prior_Destroy(cfg->fp7_bld->prior);
	cfg->fp7_bld->prior = cm_p7_prior_CreateNucleic();
	cfg->fp7_bld->re_target = DEFAULT_ETARGET_HMMFILTER;
#endif

	/* the following code is cut&pasted from the 'master' function in cmbuild.c
	* removed: code having to do with the clustering functions, and the --refine flag, and saving the CM to a file
	* added: returning the CM in ret_cm and only allowing 1 alignment in the file
	*/
	if(set_msa_name(go, cfg, errbuf, msa) != eslOK) cm_Fail(errbuf);
	if(msa->name == NULL)                           cm_Fail("Error naming MSA");
	int ncm = 1;     /* default: only build 1 CM for each MSA in alignment file */

	int c;
	for(c = 0; c < ncm; c++)
	{
          int failed=0;
		cfg->ncm_total++;

		/* if being verbose, print some stuff about what we're about to do.
		*/
		if (cfg->be_verbose) {
			fprintf(cfg->ofp, "Alignment:           %s\n",           msa->name);
			fprintf(cfg->ofp, "Number of sequences: %d\n",           msa->nseq);
			fprintf(cfg->ofp, "Number of columns:   %" PRId64 "\n",  msa->alen);
			if(esl_opt_GetString(go, "--rsearch") != NULL)
				printf ("RIBOSUM Matrix:      %s\n",  cfg->fullmat->name);
			fputs("", cfg->ofp);
			fflush(cfg->ofp);
		}

		/* msa -> cm */
		CM_t    *cm = NULL;
		Parsetree_t  *mtr=NULL;
		Parsetree_t **tr=NULL;
		int i = 0;
		if ((status = process_build_workunit(go, cfg, errbuf, msa, &cm, &mtr, &tr, hmmStatsToCalibrate)) != eslOK) {
			if (status==eslEINCOMPAT) {
				fprintf(stderr,"CM building failed, but I'm assuming it's benign because we just gave a silly alignment (and not because of a bug, or an error that shouldn't happen like out of memory).  This follows from the code in 'cm_from_guide' in src/cm_modelmaker.c in Infernal.  ERROR: %s\n",errbuf);
				*ret_cm=NULL;
				failed=1;
			}
			else {
				cm_Fail(errbuf);
			}
		}

		if (!failed) {
			if ((status = cm_Validate(cm, 0.0001, errbuf)) != eslOK) { esl_fatal("cm_Validate failed"); }

			if(cfg->be_verbose) { 
				fprintf(cfg->ofp, "\n");
				SummarizeCM(cfg->ofp, cm);  
				fprintf(cfg->ofp, "//\n");
			}

			*ret_cm  = cm;
		}

		if (tr != NULL) {
			for (i = 0; i < msa->nseq; i++) FreeParsetree(tr[i]);
			free(tr);
		}
		if (mtr != NULL) FreeParsetree(mtr);

		if (failed) {
		  break;
		}
	}
	
	destroy_cfg(cfg);
	esl_getopts_Destroy(go);
}
