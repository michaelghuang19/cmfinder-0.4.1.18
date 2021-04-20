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

// for debugging which ESL_XFAIL triggers
//#undef ESL_XFAIL
//#define ESL_XFAIL(A,B,C,...) printf("line %d.",__LINE__);

void *MallocOrDie (size_t n)
{
	void *p=malloc(n);
	if (p==NULL) {
		esl_fatal("malloc failed");
	}
	return p;
}

int TriIndex(int i,int j) {
	int coor = ((i < j)? ( j * (j-1)/2 + i) : ( i * (i-1)/2 + j));
	return coor;
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

int**
IntAlloc2DArray(int acol)
{
	int **mat;
	int j;
	mat = (int **) malloc (acol * sizeof (int *));
	for (j = 0; j < acol; j++){
		mat[j] = (int *) malloc ((j+1) * sizeof(int));
		memset(mat[j], 0, sizeof(int) * (j+1));
	}
	return mat;
}

double**
DoubleAlloc2DArray(int acol)
{
	double **mat;
	int j;
	mat = (double **) malloc (acol * sizeof (double *));
	for (j = 0; j < acol; j++){
		mat[j] = (double *) malloc ((j+1) * sizeof(double));
		memset(mat[j], 0, sizeof(double) * (j+1));
	}
	return mat;
}

void
Free2DArray(void **p, int dim1)
{
	int i;

	if (p != NULL) {
		for (i = 0; i < dim1; i++) {
			if (p[i] != NULL) {
				free(p[i]);
			}
		}
		free(p);
	}
}

int IsCanonicalDigitalBpAndNotDegen (int nuc1,int nuc2)
{
	assert(esl_abc_XIsValid(abc,nuc1) && esl_abc_XIsValid(abc,nuc2)); /* should be valid digital nucs */
	assert(esl_abc_DigitizeSymbol(abc,'A')==0 && esl_abc_DigitizeSymbol(abc,'C')==1 && esl_abc_DigitizeSymbol(abc,'G')==2 && esl_abc_DigitizeSymbol(abc,'U')==3);
	
	switch (nuc1) {
	case 0:
		return nuc2==3;
	case 1:
		return nuc2==2;
	case 2:
		return nuc2==1 || nuc2==3;
	case 3:
		return nuc2==0 || nuc2==2;
	default:
		/* we don't consider degenerate nucs to be canonical */
		/* and gaps are obviously not */
		return 0;
	}
}

int* GetPairtable(char* ss)
{
	int  i,j;
	int  sp=0;  
	int  len = (int)(strlen(ss));
	int* stack = malloc(sizeof(int) * len);
	int* pt = malloc(sizeof(int) * len);
	for(i=0; i < len; i++) {
		pt[i] = -1;    
	}

	for(i=0; i < len; i++) {
		if (pair_left(ss[i])) {
			stack[sp++] = i;      
		}
		if (pair_right(ss[i])){
			--sp;
			if (sp < 0) {
				esl_fatal("Structure %s unbalanced base pair at pos %d", ss, i);
			}
			j = stack[sp];
			pt[j] = i;
			pt[i] = j;      
		}	    
	}
	return pt;  
}

/* Remove gaps in a sequence */
char*                         /* The sequence without a gap */
remove_gap(char* seq,                  /* The original sequence */
	   int** ret_idx_map)        /* The mapping of the old sequence index to new sequence index */
{
	int   len = (int)(strlen(seq));
	char* s = malloc(sizeof(char) * (len + 1));
	int   *idx_map = malloc(sizeof(int) * (len + 1));
	char  *s1, *s2;

	s1 = seq;
	s2 = s;
	while (*s1 != '\0') {
		if (esl_abc_CIsGap(abc,*s1)) {
			idx_map[s1 - seq ] = -1;
			s1++;
		}
		else {
			idx_map[s1 - seq ] = s2 - s ;
			*s2 = toupper((int)(*s1));
			s2++;
			s1++;
		}
	}
	*s2 = '\0';
	*ret_idx_map = idx_map;
	return s;
}

ESL_MSA *MakeTextMsaClone(ESL_MSA *digital_msa)
{
	ESL_MSA *text_msa=esl_msa_Clone(digital_msa);
	if (esl_msa_Textize(text_msa)!=eslOK) {
		esl_fatal("problem textizing MSA");
	}
	return text_msa;
}

void SaveMsa (const char *filename,ESL_MSA *msa)
{
	FILE *fp;
	int status;
	fp=fopen(filename, "w");
	if (fp==NULL) {
		esl_fatal("failed to open %s for writing",filename);
	}
	status = eslx_msafile_Write(fp,msa,eslMSAFILE_PFAM);
	if      (status == eslEMEM) esl_fatal("Memory error when outputting alignment to file %s\n",filename);
	else if (status != eslOK)   esl_fatal("Writing alignment to file %s failed with error %d\n",filename,status);
	fclose(fp);
}

void SaveCM(const char *fileName,CM_t *cm)
{
	int status;
	FILE *f=fopen(fileName,"wt");
	if (f==NULL) {
		esl_fatal("could not open file %s",fileName);
	}
	if ((status = cm_file_WriteASCII(f, -1, cm)) != eslOK) esl_fatal("CM save failed");
	fclose(f);
	printf("saved CM to %s\n",fileName);
}

void SaveSeqCacheAsFasta (const char *fileName,ESL_SQCACHE *seqs)
{
	int s,p;
	FILE *f=fopen(fileName,"wt");
	if (f==NULL) {
		esl_fatal("could not open file %s",fileName);
	}
	for (s=0; s<seqs->seq_count; s++) {
		ESL_SQ *sq=&(seqs->sq_list[s]);
		fprintf(f,">%d\n",s);
		for (p=0; p<sq->n; p++) {
			fprintf(f,"%c",abc->sym[sq->dsq[p+1]]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
	printf("saved ESL_SQCACHE to file %s\n",fileName);
}

void InitCmdlineBuilder (CmdlineBuilder *cmdline,char *buf)
{
	cmdline->buf=buf;
	cmdline->buf[0]='\0';
	cmdline->pos=0;
	cmdline->maxpos=COMMANDLINE_FLAGS_MAX;
}
void PrintfCmdlineBuilder (CmdlineBuilder *cmdline,const char *format, ...)
{
	va_list argp;
	int printed;
	int size_remaining=cmdline->maxpos - cmdline->pos-1;
	va_start(argp, format);
	printed=vsnprintf(cmdline->buf+cmdline->pos,size_remaining,format, argp);
	va_end(argp);
	if (printed>=size_remaining || printed<0) {
		esl_fatal("overran buffer or other error in PrintfCmdlineBuilder (%s:%d)",__FILE__,__LINE__);
	}
	cmdline->pos += printed;
}

void cmcalibrate_viterbi (CM_t *cm,CmfinderVars *cmfinderVars,ExpInfo_t ***ret_expA)
{
	int expModesToCalibrate[EXP_NMODES];
	int i;
	char cmcalibrate_commandline_flags[COMMANDLINE_FLAGS_MAX];
	CmdlineBuilder cmdline;
	ESL_STOPWATCH *timer_cmcalibrate=esl_stopwatch_Create();
	
	for (i=0; i<EXP_NMODES; i++) {
		expModesToCalibrate[i]=0;
	}
	if (cmfinderVars->do_local) {
		expModesToCalibrate[EXP_CM_LC]=1;
	}
	else {
		expModesToCalibrate[EXP_CM_GC]=1;
	}
	
	InitCmdlineBuilder(&cmdline,cmcalibrate_commandline_flags);
	PrintfCmdlineBuilder(&cmdline,"executable-dummy-token");
#ifdef HMMER_THREADS
	PrintfCmdlineBuilder(&cmdline," --cpu %d",cmfinderVars->numCpu);
#endif

	if (cmfinderVars->cmcalibrate_tailn>0) {
		PrintfCmdlineBuilder(&cmdline," %s %d",cmfinderVars->do_local?"--ltailn":"--gtailn",cmfinderVars->cmcalibrate_tailn);
	}
	
	printf("running internal cmcalibrate with flags: %s\n",cmcalibrate_commandline_flags);
	esl_stopwatch_Start(timer_cmcalibrate);
	cmcalibrate(cmcalibrate_commandline_flags,cm,expModesToCalibrate,ret_expA);
	esl_stopwatch_Stop(timer_cmcalibrate);

	esl_stopwatch_Include(cmfinderVars->timer_cmcalibrate,timer_cmcalibrate);
	esl_stopwatch_Destroy(timer_cmcalibrate);
}

void cmsearch_wrapper (CM_t *cm,ESL_SQCACHE *dbSeqs,CM_TOPHITS **ret_th,CmfinderVars *cmfinderVars,double cm_scoreThreshold,ScanMode scanMode)
{
	double Zmb_for_cmsearch=(double)(cmfinderVars->inputSeqSize)/1e6;
	ESL_STOPWATCH *timer_cmsearch=esl_stopwatch_Create();
	int do_inside=scanMode.force_inside_alg;
	char cmsearch_commandline_flags[COMMANDLINE_FLAGS_MAX];
	const char *devNullRedir="-o /dev/null";
	/* run cmsearch on all seqs in inputSeqCache, so we can enumerate candidates
	 * -o /dev/null : leave the output code in cmsearch for now, just in case it's useful for debugging/diagnostics later, but don't do anything with it
	 * -g : do glocal (later I might try local)
	 * --notrunc : don't use the truncation-aware algs, for now
	 * --toponly : we never search revcomp in CMfinder
	 * --cyk : use only the CYK algs for now, to be similar to original CMfinder  NOTE: if we use inside, then we have to set the --F6 CYK filter threshold P-value, which'll be hard without being able to run cmcallibrate
	 * --acyk : also only use CYK for aligning, for now
	 * --max : don't use filters, for now
	 * -T cm_scoreThreshold : might help to reduce the amount of E-value code I have to remove
	 */
	CmdlineBuilder cmdline;
	
	if (dbSeqs->seq_count==0) {
		*ret_th=cm_tophits_Create();
		return;
	}
	InitCmdlineBuilder(&cmdline,cmsearch_commandline_flags);
	PrintfCmdlineBuilder(&cmdline,"executable-dummy-token");
#ifdef HMMER_THREADS
	PrintfCmdlineBuilder(&cmdline," --cpu %d",cmfinderVars->numCpu);
#endif
	if (do_inside) {
		/* oops, conflicts with --max "--noF6"; */ /* we can't afford to calc EVDs, and I'm not ready to try any hacks (plus for bacteria most of the sequence will be a hit, so there isn't that much benefit.  Not sure if CYK is used to make bands for Inside, but I think not.) */
	}
	else {
		PrintfCmdlineBuilder(&cmdline," --cyk");
	}
	if (scanMode.need_alignments) {
		if (cmfinderVars->use_maxAlignAccuracy) {
			/* nothing to do, this is the default */
		}
		else {
			PrintfCmdlineBuilder(&cmdline," --acyk");
		}
	}
	else {
		PrintfCmdlineBuilder(&cmdline," --noali");
	}
	if (cmfinderVars->use_fragmentary && cmfinderVars->filter_max) {
		esl_fatal("--fragmentary and --max are both set (at least internally).  These are incompatible because cmsearch in Infernal-1.1rc2 does not implement them.");
	}
	if (scanMode.use_fragmentary) {
	}
	else {
		PrintfCmdlineBuilder(&cmdline," --notrunc");
	}
	if (cmfinderVars->filter_max) {
		PrintfCmdlineBuilder(&cmdline," --max");
	}
	else {
		PrintfCmdlineBuilder(&cmdline," --noF6");
		if (scanMode.disable_filters) {
			PrintfCmdlineBuilder(&cmdline," --noF1 --noF2 --noF2b --noF3 --noF3b --noF4 --noF4b --noF5");
		}
		else {
			if (cmfinderVars->filter_noF4F5) {
				PrintfCmdlineBuilder(&cmdline," --noF4 --noF4b --noF5");
			}
		}
	}
	if (cmfinderVars->do_local) {
		/* this is default for cmsearch */
	}
	else {
		PrintfCmdlineBuilder(&cmdline," -g");
	}
	if (cmfinderVars->use_evalue) {
		PrintfCmdlineBuilder(&cmdline," -E %g",cmfinderVars->evalue);
	}
	else {
		PrintfCmdlineBuilder(&cmdline," -T %g",cm_scoreThreshold);
	}
#ifndef _POSIX_VERSION
	devNullRedir="";
#endif
	PrintfCmdlineBuilder(&cmdline," --toponly %s --mxsize %g --smxsize %g -Z %g",
		devNullRedir,cmfinderVars->cmsearch_mxsize,cmfinderVars->cmsearch_smxsize,Zmb_for_cmsearch);
	if (1) {
		printf("running internal-cmsearch %s\n",cmsearch_commandline_flags);
	}
	
	esl_stopwatch_Start(timer_cmsearch);
	cmsearch(cmsearch_commandline_flags,cm,dbSeqs,abc,ret_th,cmfinderVars->use_evalue);
	esl_stopwatch_Stop(timer_cmsearch);
	
	esl_stopwatch_Include(scanMode.force_inside_alg?cmfinderVars->timer_cmsearch_inside:cmfinderVars->timer_cmsearch_viterbi,timer_cmsearch);
	esl_stopwatch_Destroy(timer_cmsearch);
}


typedef struct {
	char *p;
	size_t n,alloc_size;
} file_as_string;

void file_as_string_open (file_as_string *f)
{
	f->n=0;
	f->alloc_size=10; /* make sure it's tested */
	f->p=(char *)MallocOrDie(f->alloc_size+1);  /* alloc +1 for the extra zero-termination, because I feel it makes the code easier */
}
void file_as_string_close (file_as_string *f)
{
	free(f->p);
}
int file_as_string_printf (file_as_string *f,const char *format, ...)
{
	char *new_p;
	int printed;
	va_list argp;

	while (1) {
		int size_remaining=f->alloc_size - f->n;
		va_start(argp, format);
		printed=vsnprintf(f->p+f->n,size_remaining,format, argp);
		va_end(argp);
		if (printed>=size_remaining || printed<0) {  /* printed<0 is apparently on some older systems that this code probably won't work on anyway */
			f->alloc_size *= 2;
			new_p=(char *)(realloc(f->p,f->alloc_size+1));
			if (new_p==NULL) {
				esl_fatal("realloc %d failed",f->alloc_size);
			}
			f->p=new_p;
		}
		else {
			f->n += printed;
			f->p[f->n]=0;
			break;
		}
	}
	return 1;
}
void file_as_string_puts (const char *s,file_as_string *f)
{
	file_as_string_printf(f,"%s",s);
}
void file_as_string_putc (int ch,file_as_string *f)
{
	file_as_string_printf(f,"%c",ch);
}
size_t file_as_string_fwrite (char *buf,size_t size,size_t n,file_as_string *f)
{
	char *s=(char *)MallocOrDie(size*n+1);
	strncpy(s,buf,size*n);
	s[size*n]=0;
	file_as_string_printf(f,"%s",s);
	return n;
}

int file_as_string_p7_hmmfile_WriteASCII (file_as_string *f,int format,P7_HMM *hmm)
{
	int status;
	char *ascii_hmm;
	status=p7_hmmfile_WriteToString(&ascii_hmm,-1,hmm);
	file_as_string_printf(f,"%s",ascii_hmm);
	free(ascii_hmm);
	return status;
}

/* the following code is copied from cm_file.c, to allow writing to a buffer
 * 
 * WARNING: code to save the HMM was modified (see comment at end of CM write function)
 * 
 * Otherwise, changes were done using #defines to use the file_as_string functions instead of the FILE * functions.
 * */
#define fprintf file_as_string_printf
#define fputs file_as_string_puts
#define fputc file_as_string_putc
#define fwrite file_as_string_fwrite
#define p7_hmmfile_WriteASCII file_as_string_p7_hmmfile_WriteASCII
static char *
prob2ascii(float p, float null)
{
  static char buffer[32];

  if (p == 0.0) return "*";
  sprintf(buffer, "%.3f", sreLOG2(p/null));
  return buffer;
}
static float
ascii2prob(char *s, float null)
{
  return (*s == '*') ? 0. : exp(atof(s)/1.44269504)*null;
}

/* EPN, Tue Aug  7 15:54:15 2007
 * is_integer() and is_real(), savagely ripped verbatim out
 * of Easel's esl_getopts.c, where they were private.
 */
/* Function: is_integer()
 * 
 * Returns TRUE if <s> points to something that atoi() will parse
 * completely and convert to an integer.
 */
static int
is_integer(char *s)
{
  int hex = 0;

  if (s == NULL) return 0;
  while (isspace((int) (*s))) s++;      /* skip whitespace */
  if (*s == '-' || *s == '+') s++;      /* skip leading sign */
				        /* skip leading conversion signals */
  if ((strncmp(s, "0x", 2) == 0 && (int) strlen(s) > 2) ||
      (strncmp(s, "0X", 2) == 0 && (int) strlen(s) > 2))
    {
      s += 2;
      hex = 1;
    }
  else if (*s == '0' && (int) strlen(s) > 1)
    s++;
				/* examine remainder for garbage chars */
  if (!hex)
    while (*s != '\0')
      {
	if (!isdigit((int) (*s))) return 0;
	s++;
      }
  else
    while (*s != '\0')
      {
	if (!isxdigit((int) (*s))) return 0;
	s++;
      }
  return 1;
}


/* is_real()
 * 
 * Returns TRUE if <s> is a string representation
 * of a valid floating point number, convertable
 * by atof().
 */
static int
is_real(char *s)
{
  int gotdecimal = 0;
  int gotexp     = 0;
  int gotreal    = 0;

  if (s == NULL) return 0;

  while (isspace((int) (*s))) s++; /* skip leading whitespace */
  if (*s == '-' || *s == '+') s++; /* skip leading sign */

  /* Examine remainder for garbage. Allowed one '.' and
   * one 'e' or 'E'; if both '.' and e/E occur, '.'
   * must be first.
   */
  while (*s != '\0')
    {
      if (isdigit((int) (*s))) 	gotreal++;
      else if (*s == '.')
	{
	  if (gotdecimal) return 0; /* can't have two */
	  if (gotexp) return 0;     /* e/E preceded . */
	  else gotdecimal++;
	}
      else if (*s == 'e' || *s == 'E')
	{
	  if (gotexp) return 0;	/* can't have two */
	  else gotexp++;
	}
      else if (isspace((int) (*s)))
	break;
      s++;
    }

  while (isspace((int) (*s))) s++;         /* skip trailing whitespace */
  if (*s == '\0' && gotreal) return 1;
  else return 0;
}
static int
multiline(file_as_string *fp, const char *pfx, char *s)
{
  char *sptr  = s;
  char *end   = NULL;
  int   n     = 0;
  int   nline = 1;

  do {
    end = strchr(sptr, '\n');

    if (end != NULL)                  /* if there's no \n left, end == NULL */
      {
  n = end - sptr;                       /* n chars exclusive of \n */
  if (fprintf(fp, "%s [%d] ", pfx, nline++) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (fwrite(sptr, sizeof(char), n, fp)    != n) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");  /* using fwrite lets us write fixed # of chars   */
  if (fprintf(fp, "\n")                     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");  /* while writing \n w/ printf allows newline conversion */
  sptr += n + 1;                       /* +1 to get past \n */
      } 
    else 
      {
  if (fprintf(fp, "%s [%d] %s\n", pfx, nline++, sptr) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); /* last line */
      }
  } while (end != NULL  && *sptr != '\0');   /* *sptr == 0 if <s> terminates with a \n */
  return eslOK;
}
int
cm_file_WriteASCII_file_as_string(file_as_string *fp, int format, CM_t *cm,CmfinderVars *cmfinderVars)
{
  int x, v, nd, y;

  if((cm->flags & CMH_LOCAL_BEGIN) || (cm->flags & CMH_LOCAL_END)) cm_Fail("cm_file_WriteASCII(): CM is in local mode");

  if (format == -1) format = CM_FILE_1a;

  if   (format == CM_FILE_1a) fprintf(fp, "INFERNAL1/a [%s | %s]\n", INFERNAL_VERSION, INFERNAL_DATE);
  else ESL_EXCEPTION(eslEINVAL, "invalid CM file format code");
  
  fprintf(fp, "NAME     %s\n", cm->name);
  if (cm->acc)  fprintf(fp, "ACC      %s\n", cm->acc);
  if (cm->desc) fprintf(fp, "DESC     %s\n", cm->desc);
  fprintf(fp, "STATES   %d\n", cm->M);
  fprintf(fp, "NODES    %d\n", cm->nodes);
  fprintf(fp, "CLEN     %d\n", cm->clen);
  fprintf(fp, "W        %d\n", cm->W);
  fprintf(fp, "ALPH     %s\n", esl_abc_DecodeType(cm->abc->type));
  fprintf(fp, "RF       %s\n", (cm->flags & CMH_RF)   ? "yes" : "no");
  fprintf(fp, "CONS     %s\n", (cm->flags & CMH_CONS) ? "yes" : "no");
  fprintf(fp, "MAP      %s\n", (cm->flags & CMH_MAP)  ? "yes" : "no");
  if (cm->ctime   != NULL) fprintf  (fp, "DATE     %s\n", cm->ctime);
  if (cm->comlog  != NULL) multiline(fp, "COM     ",     cm->comlog);
  fprintf(fp, "PBEGIN   %g\n", cm->pbegin);
  fprintf(fp, "PEND     %g\n", cm->pend);
  fprintf(fp, "WBETA    %g\n", cm->beta_W);
  fprintf(fp, "QDBBETA1 %g\n", cm->qdbinfo->beta1);
  fprintf(fp, "QDBBETA2 %g\n", cm->qdbinfo->beta2);
  fprintf(fp, "N2OMEGA  %6g\n",cm->null2_omega);
  fprintf(fp, "N3OMEGA  %6g\n",cm->null3_omega);
  fprintf(fp, "ELSELF   %.8f\n",cm->el_selfsc);
  fprintf(fp, "NSEQ     %d\n", cm->nseq);
  fprintf(fp, "EFFN     %f\n",cm->eff_nseq);
  if (cm->flags & CMH_CHKSUM)  fprintf(fp, "CKSUM    %u\n", cm->checksum); /* unsigned 32-bit */
  fputs("NULL    ", fp);
  for (x = 0; x < cm->abc->K; x++) { fprintf(fp, "%6s ", prob2ascii(cm->null[x], 1/(float)(cm->abc->K))); }
  fputc('\n', fp);
  if (cm->flags & CMH_GA)  fprintf(fp, "GA       %.2f\n", cm->ga);
  if (cm->flags & CMH_TC)  fprintf(fp, "TC       %.2f\n", cm->tc);
  if (cm->flags & CMH_NC)  fprintf(fp, "NC       %.2f\n", cm->nc);

  if (cm->flags & CMH_FP7) {
	  double mu=cmfinderVars->hmmStatsToCalibrate[CM_p7_GFMU]?cm->fp7_evparam[CM_p7_GFMU]:0;
	  double lambda=cmfinderVars->hmmStatsToCalibrate[CM_p7_GFLAMBDA]?cm->fp7_evparam[CM_p7_GFLAMBDA]:0;
    fprintf(fp, "EFP7GF   %.4f %.5f\n", mu,  lambda);
  }
  if (cm->flags & CMH_EXPTAIL_STATS)
    {
      /* make sure our dbsize values can be cast to a long reliably,
       * even on 32 bit systems (<= 2 Gb), they certainly should be,
       * max value in cmcalibrate is 160Mb. (This is related to bug
       * i31.)
       */
      if(cm->expA[EXP_CM_LC]->dbsize > (2000. * 1000000.)) ESL_EXCEPTION(eslEINVAL, "invalid dbsize (too big) EXP_CM_LC");
      if(cm->expA[EXP_CM_GC]->dbsize > (2000. * 1000000.)) ESL_EXCEPTION(eslEINVAL, "invalid dbsize (too big) EXP_CM_GC");
      if(cm->expA[EXP_CM_LI]->dbsize > (2000. * 1000000.)) ESL_EXCEPTION(eslEINVAL, "invalid dbsize (too big) EXP_CM_LI");
      if(cm->expA[EXP_CM_GI]->dbsize > (2000. * 1000000.)) ESL_EXCEPTION(eslEINVAL, "invalid dbsize (too big) EXP_CM_GI");

      fprintf(fp, "ECMLC    %.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
	      cm->expA[EXP_CM_LC]->lambda, cm->expA[EXP_CM_LC]->mu_extrap, cm->expA[EXP_CM_LC]->mu_orig, 
	      (long) (cm->expA[EXP_CM_LC]->dbsize + 0.5), cm->expA[EXP_CM_LC]->nrandhits, cm->expA[EXP_CM_LC]->tailp);
      fprintf(fp, "ECMGC    %.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
	      cm->expA[EXP_CM_GC]->lambda, cm->expA[EXP_CM_GC]->mu_extrap, cm->expA[EXP_CM_GC]->mu_orig, 
	      (long) (cm->expA[EXP_CM_GC]->dbsize + 0.5), cm->expA[EXP_CM_GC]->nrandhits, cm->expA[EXP_CM_GC]->tailp);
      fprintf(fp, "ECMLI    %.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
	      cm->expA[EXP_CM_LI]->lambda, cm->expA[EXP_CM_LI]->mu_extrap, cm->expA[EXP_CM_LI]->mu_orig, 
	      (long) (cm->expA[EXP_CM_LI]->dbsize + 0.5), cm->expA[EXP_CM_LI]->nrandhits, cm->expA[EXP_CM_LI]->tailp);
      fprintf(fp, "ECMGI    %.5f  %10.5f  %10.5f  %10ld  %10d  %.6f\n", 
	      cm->expA[EXP_CM_GI]->lambda, cm->expA[EXP_CM_GI]->mu_extrap, cm->expA[EXP_CM_GI]->mu_orig, 
	      (long) (cm->expA[EXP_CM_GI]->dbsize + 0.5), cm->expA[EXP_CM_GI]->nrandhits, cm->expA[EXP_CM_GI]->tailp);
    }

  /* main model section */
  fputs("CM\n", fp);

  /* Create emit map if nec, so we can output map, consensus and rf info appropriately */
  if(cm->emap == NULL) { 
    cm->emap = CreateEmitMap(cm);
    if(cm->emap == NULL) ESL_EXCEPTION(eslEINVAL, "unable to create an emit map");
  }  

  for (v = 0; v < cm->M; v++) { 
    nd = cm->ndidx[v];

    /* Node line. node type and additional per-consensus position annotation */
    if (cm->nodemap[nd] == v) { 
      fprintf(fp, "%45s[ %-4s %4d ]", "", Nodetype(cm->ndtype[nd]), nd);

      /* additional annotation */
      /* map (optional) */
      if(cm->flags & CMH_MAP) { 
	if     (cm->ndtype[nd] == MATP_nd) fprintf(fp, " %6d %6d", cm->map[cm->emap->lpos[nd]], cm->map[cm->emap->rpos[nd]]); 
	else if(cm->ndtype[nd] == MATL_nd) fprintf(fp, " %6d %6s", cm->map[cm->emap->lpos[nd]], "-");
	else if(cm->ndtype[nd] == MATR_nd) fprintf(fp, " %6s %6d", "-", cm->map[cm->emap->rpos[nd]]);
	else 	                           fprintf(fp, " %6s %6s", "-", "-");
      }
      else { /* no map annotation */
	fprintf(fp, " %6s %6s", "-", "-");
      }
      /* consensus sequence (mandatory) */
      if     (cm->ndtype[nd] == MATP_nd) fprintf(fp, " %c %c", cm->consensus[cm->emap->lpos[nd]], cm->consensus[cm->emap->rpos[nd]]); 
      else if(cm->ndtype[nd] == MATL_nd) fprintf(fp, " %c %c", cm->consensus[cm->emap->lpos[nd]], '-');
      else if(cm->ndtype[nd] == MATR_nd) fprintf(fp, " %c %c", '-', cm->consensus[cm->emap->rpos[nd]]);
      else 	                         fprintf(fp, " %c %c", '-', '-');
      /* RF (optional) */
      if(cm->flags & CMH_RF) { 
	if     (cm->ndtype[nd] == MATP_nd) fprintf(fp, " %c %c", cm->rf[cm->emap->lpos[nd]], cm->rf[cm->emap->rpos[nd]]); 
	else if(cm->ndtype[nd] == MATL_nd) fprintf(fp, " %c %c", cm->rf[cm->emap->lpos[nd]], '-');
	else if(cm->ndtype[nd] == MATR_nd) fprintf(fp, " %c %c", '-', cm->rf[cm->emap->rpos[nd]]);
	else 	                           fprintf(fp, " %c %c", '-', '-');
      }
      else { /* no RF annotation */
	fprintf(fp, " %c %c", '-', '-');
      }
      fputs("\n", fp);
    }    

    /* State line, w/ parents, children, dmin2, dmin1, dmax1, dmax2, transitions and emissions */
    fprintf(fp, "    %2s %5d %5d %1d %5d %5d %5d %5d %5d %5d ", 
	    Statetype(cm->sttype[v]), v, 
	    cm->plast[v], cm->pnum[v],
	    cm->cfirst[v], cm->cnum[v], 
	    cm->qdbinfo->dmin2[v], cm->qdbinfo->dmin1[v], 
	    cm->qdbinfo->dmax1[v], cm->qdbinfo->dmax2[v]);

    /* Transitions */
    if (cm->sttype[v] != B_st) { 
      for (x = 0; x < cm->cnum[v]; x++) { 
	fprintf(fp, "%7s ", prob2ascii(cm->t[v][x], 1.));
      }
    }
    else {
      x = 0;
    }
    for (; x < 6; x++) {
      fprintf(fp, "%7s ", "");
    }
      
    /* Emissions */ 
    if (cm->sttype[v] == MP_st) {
      for (x = 0; x < cm->abc->K; x++) { 
	for (y = 0; y < cm->abc->K; y++) {
	  fprintf(fp, "%6s ", prob2ascii(cm->e[v][x*cm->abc->K+y], cm->null[x]*cm->null[y]));
	}
      }
    }
    else if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
      for (x = 0; x < cm->abc->K; x++) {
	fprintf(fp, "%6s ", prob2ascii(cm->e[v][x], cm->null[x]));
      }
    }
    fputs("\n", fp);
  }
  fputs("//\n", fp);

	/* WARNING: here's the only change from the infernal code
	 * the CM reading code assumes that the HMM will be either separate or in the same file, but I
	 * don't think I can set cmfp->hfp to point to the same thing.  So I'll save the CM and the HMM into separate files.
	 */
#if 0
  /* print additional p7 hmms if any */
  if(cm->flags & CMH_FP7 && cm->fp7 != NULL) {
    p7_hmmfile_WriteASCII(fp, -1, cm->fp7);
  }
#endif
  return eslOK;
}
#undef fprintf
#undef p7_hmmfile_WriteASCII

/* copied from cm_file.c
 * I can't figure out how to allow it to work without modifying the code
 * the code within 'if (read_fp7)' is changed from the original
 */
static int
read_asc_1p1_cm__with_RAM(CM_FILE *cmfp, int read_fp7, ESL_ALPHABET **ret_abc, CM_t **opt_cm)
{
  int           status;
  ESL_ALPHABET *abc  = NULL;
  CM_t         *cm   = NULL;
  P7_HMM       *hmm  = NULL;
  char         *tag  = NULL;
  char         *tok1 = NULL;
  char         *tok2 = NULL;
  char         *tok3 = NULL;
  char         *tok4 = NULL;
  char         *tok5 = NULL;
  char         *tok6 = NULL;
  int           alphatype;
  off_t         offset = 0;
  int           v, x, y, nd;            /* counters */
  int           read_fp7_stats = FALSE;
  uint32_t      cm_statstracker = 0; /* for making sure we have all CM E-value stats, if we have any */
  int           exp_mode;   
  int           read_el_selfsc = FALSE; /* set to true when we read ELSELF line */

  /* temporary parameters, for storing values prior to their allocation in the CM */
  float *tmp_null          = NULL;
  float  tmp_fp7_gfmu;
  float  tmp_fp7_gflambda;
  double tmp_qdbbeta1;
  double tmp_qdbbeta2;

  /* temporary per-node annotation, will be converted to per-consensus position once the model
   * architecture is known, after the full model is read */
  char *tmp_rf_left    = NULL;
  char *tmp_rf_right   = NULL;
  char *tmp_cons_left  = NULL;
  char *tmp_cons_right = NULL;
  int  *tmp_map_left   = NULL;
  int  *tmp_map_right  = NULL;

  cmfp->errbuf[0] = '\0';

  if (cmfp->newly_opened)
    {
      offset            = 0;
      cmfp->newly_opened = FALSE;
    }
  else
    {
      /* Record where this CM starts on disk */
      if ((! cmfp->do_stdin) && (! cmfp->do_gzip) && (offset = ftello(cmfp->f)) < 0)   ESL_XEXCEPTION(eslESYS, "ftello() failed");

      /* First line of file: "INFERNAL1/a". Allocate shell for CM annotation information (we don't know M,nodes yet) */
      if ((status = esl_fileparser_NextLine(cmfp->efp))                   != eslOK)  goto ERROR;  /* EOF here is normal; could also be a thrown EMEM */
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tag, NULL)) != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "unexpected absence of tokens on data line");

      if      (cmfp->format == CM_FILE_1a) { if (strcmp(tag, "INFERNAL1/a") != 0)    ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Didn't find INFERNAL1/a tag: bad format or not an INFERNAL save file?"); }
      else                                                                           ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "No such CM file format code: this shouldn't happen");
    }

  if ((cm = CreateCMShell()) == NULL)   ESL_XFAIL(eslEMEM,    cmfp->errbuf, "allocation failure, CM shell");
  cm->offset = offset;

  /* Header section */
  while ((status = esl_fileparser_NextLine(cmfp->efp)) == eslOK)
    {
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tag, NULL))     != eslOK)   ESL_XFAIL(status,    cmfp->errbuf, "Premature end of line");

      if (strcmp(tag, "NAME") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))  != eslOK)   ESL_XFAIL(status,    cmfp->errbuf, "No name found on NAME line");
	cm_SetName(cm, tok1);
      } 

      else if (strcmp(tag, "ACC") == 0)  {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))  != eslOK)   ESL_XFAIL(status,    cmfp->errbuf, "No accession found on ACC line");
	cm_SetAccession(cm, tok1); 
      }  

      else if (strcmp(tag, "DESC") == 0) {
	if ((status = esl_fileparser_GetRemainingLine(cmfp->efp, &tok1))      != eslOK)   ESL_XFAIL(status,    cmfp->errbuf, "No description found on DESC line");
	cm_SetDescription(cm, tok1);
      } 

      else if (strcmp(tag, "STATES") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))  != eslOK)   ESL_XFAIL(status,    cmfp->errbuf, "Number of model states not found on STATES line");
	if ((cm->M = atoi(tok1))                                              == 0)  	  ESL_XFAIL(status,    cmfp->errbuf, "Invalid number of states %s on STATES line", tok1);
      }  

      else if (strcmp(tag, "NODES") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))  != eslOK)   ESL_XFAIL(status,    cmfp->errbuf, "Number of model nodes not found on NODES line");
	if ((cm->nodes = atoi(tok1))                                          == 0)  	  ESL_XFAIL(status,    cmfp->errbuf, "Invalid number of nodes %s on NODES line", tok1);
      }  

      else if (strcmp(tag, "CLEN") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))  != eslOK)   ESL_XFAIL(status,    cmfp->errbuf, "No consensus length found on CLEN line");
	if ((cm->clen = atoi(tok1))                                           == 0)   	  ESL_XFAIL(status,    cmfp->errbuf, "Invalid consensus length %s on CLEN line", tok1);
      }

      else if (strcmp(tag, "W") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))  != eslOK)   ESL_XFAIL(status,    cmfp->errbuf, "No consensus length found on W line");
	if ((cm->W = atoi(tok1))                                              == 0)   	  ESL_XFAIL(status,    cmfp->errbuf, "Invalid consensus length %s on W line", tok1);
      }

      else if (strcmp(tag, "ALPH") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))  != eslOK)   ESL_XFAIL(status,    cmfp->errbuf, "No alphabet type found on ALPH");
	if ((alphatype = esl_abc_EncodeType(tok1))                        == eslUNKNOWN)  ESL_XFAIL(status,    cmfp->errbuf, "Unrecognized alphabet type %s", tok1);
	if (*ret_abc == NULL) {
	  if ((abc = esl_alphabet_Create(alphatype))                        == NULL) 	  ESL_XFAIL(eslEMEM,   cmfp->errbuf, "Failed to create alphabet");        
	} else {
	  if ((*ret_abc)->type != alphatype)	                                          ESL_XFAIL(eslEINCOMPAT,cmfp->errbuf,"Alphabet type mismatch: was %s, but current CM says %s", esl_abc_DecodeType( (*ret_abc)->type), tok1);
	  abc = *ret_abc;
	}	  
      } 

      else if (strcmp(tag, "RF") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,    cmfp->errbuf, "No yes/no found for RF line");
	if      (strcasecmp(tok1, "yes") == 0) cm->flags |= CMH_RF;
	else if (strcasecmp(tok1, "no")  != 0)                                            ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "RF header line must say yes/no, not %s", tok1);
      } 

      else if (strcmp(tag, "CONS") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "No yes/no found for CONS line");
	if      (strcasecmp(tok1, "yes") == 0) cm->flags |= CMH_CONS;
	else if (strcasecmp(tok1, "no")  != 0)                                           ESL_XFAIL(eslEFORMAT,  cmfp->errbuf, "CONS header line must say yes/no, not %s", tok1);
      } 

      else if (strcmp(tag, "MAP") == 0) {	
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "No yes/no found for MAP line");
	if      (strcasecmp(tok1, "yes") == 0) cm->flags |= CMH_MAP;
	else if (strcasecmp(tok1, "no")  != 0)                                            ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "MAP header line must say yes/no, not %s", tok1);
      } 

      else if (strcmp(tag, "DATE") == 0) {
	if ((status = esl_fileparser_GetRemainingLine(cmfp->efp, &tok1))       != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "No date found on DATE line");
	if (esl_strdup(tok1, -1, &(cm->ctime))                                 != eslOK)  ESL_XFAIL(eslEMEM,    cmfp->errbuf, "strdup() failed to set date");
      }

      else if (strcmp(tag, "COM") == 0) {
	/* just skip the first token; it's something like [1], numbering the command lines */
	if ((status = esl_fileparser_GetTokenOnLine  (cmfp->efp, &tok1, NULL)) != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "No command number on COM line");
	if ((status = esl_fileparser_GetRemainingLine(cmfp->efp, &tok1))       != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "No command on COM line");
	if (cm->comlog == NULL) {
	  if (esl_strdup(tok1, -1, &(cm->comlog))                              != eslOK)  ESL_XFAIL(eslEMEM,    cmfp->errbuf, "esl_strdup() failed");
	} else {
	  if (esl_strcat(&(cm->comlog), -1, "\n", -1)                          != eslOK)  ESL_XFAIL(eslEMEM,    cmfp->errbuf, "esl_strcat() failed");
	  if (esl_strcat(&(cm->comlog), -1, tok1,  -1)                         != eslOK)  ESL_XFAIL(eslEMEM,    cmfp->errbuf, "esl_strcat() failed");
	}
      }

      else if (strcmp(tag, "PBEGIN") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows PBEGIN tag");
	if ((cm->pbegin = atof(tok1)) <= 0.0f)                                            ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Invalid beta on PBEGIN line: should be a positive real number, not %s", tok1);
      }

      else if (strcmp(tag, "PEND") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows PEND tag");
	if ((cm->pend = atof(tok1)) <= 0.0f)                                              ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Invalid beta on PEND line: should be a positive real number, not %s", tok1);
      }

      else if (strcmp(tag, "WBETA") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows WBETA tag");
	if ((cm->beta_W = atof(tok1)) <= 0.0f)                                            ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Invalid beta on WBETA line: should be a positive real number, not %s", tok1);
      }

      else if (strcmp(tag, "QDBBETA1") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows QDBBETA tag");
	if ((tmp_qdbbeta1 = atof(tok1)) <= 0.0f)                                          ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Invalid beta on QDBBETA line: should be a positive real number, not %s", tok1);
      }

      else if (strcmp(tag, "QDBBETA2") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows QDBBETA tag");
	if ((tmp_qdbbeta2 = atof(tok1)) <= 0.0f)                                          ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Invalid beta on QDBBETA line: should be a positive real number, not %s", tok1);
      }

      else if (strcmp(tag, "N2OMEGA") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows N2OMEGA tag");
	if ((cm->null2_omega = atof(tok1)) <= 0.0f)                                       ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Invalid omega on N2OMEGA line: should be a positive real number, not %s", tok1);
      }

      else if (strcmp(tag, "N3OMEGA") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows N3OMEGA tag");
	if ((cm->null3_omega = atof(tok1)) <= 0.0f)                                       ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Invalid omega on N3OMEGA line: should be a positive real number, not %s", tok1);
      }

      else if (strcmp(tag, "ELSELF") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows ELSELF tag");
	cm->el_selfsc = atof(tok1);
	read_el_selfsc = TRUE;
      }

      else if (strcmp(tag, "NSEQ") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows NSEQ tag");
	if ((cm->nseq = atoi(tok1)) == 0)                                                 ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Invalid nseq on NSEQ line: should be integer, not %s", tok1);
      }

      else if (strcmp(tag, "EFFN") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows EFFN tag");
	if ((cm->eff_nseq = atof(tok1)) < (-1.*eslSMALLX1))                               ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Invalid eff_nseq on EFFN line: should be a positive real number, not %s", tok1);
      }

      else if (strcmp(tag, "CKSUM") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Nothing follows CKSUM tag");
	cm->checksum = atoll(tok1); /* if atoi(), then you may truncate uint32_t checksums > 2^31-1 */
	cm->flags |= CMH_CHKSUM;
      }

      else if (strcmp(tag, "NULL") == 0) { 
	if(abc == NULL) ESL_XFAIL(status,     cmfp->errbuf, "Read NULL line before ALPH line, ALPH line must come first");
	ESL_ALLOC(tmp_null, sizeof(float) * abc->K);
	for (x = 0; x < abc->K; x++) { 
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on NULL line");
	  tmp_null[x] = ascii2prob(tok1, (1./(float) abc->K));
	}
      }

      else if (strcmp(tag, "GA") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on GA line");
	cm->ga     = atof(tok1);
	cm->flags |= CMH_GA;
      }

      else if (strcmp(tag, "TC") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on TC line");
	cm->tc     = atof(tok1);
	cm->flags |= CMH_TC;
      }

      else if (strcmp(tag, "NC") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on NC line");
	cm->nc     = atof(tok1);
	cm->flags |= CMH_NC;
      }

      else if (strcmp(tag, "EFP7GF") == 0) {
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on EFP7GF line"); /* tau */
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok2, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on EFP7GF line"); /* lambda   */
	tmp_fp7_gfmu     = atof(tok1);
	tmp_fp7_gflambda = atof(tok2);
	read_fp7_stats = TRUE;
      }

      else if (strncmp(tag, "ECM", 3) == 0) { /* one of 4 possible CM E-value lines */
	/* determine which one */
	if      (strncmp(tag+3, "LC", 2) == 0) { exp_mode = EXP_CM_LC; cm_statstracker += 1; }
	else if (strncmp(tag+3, "GC", 2) == 0) { exp_mode = EXP_CM_GC; cm_statstracker += 2; }
	else if (strncmp(tag+3, "LI", 2) == 0) { exp_mode = EXP_CM_LI; cm_statstracker += 4; }
	else if (strncmp(tag+3, "GI", 2) == 0) { exp_mode = EXP_CM_GI; cm_statstracker += 8; }
	else                                   { ESL_XFAIL(status, cmfp->errbuf, "Invalid tag beginning with ECM"); }
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on ECM.. line"); /* lambda    */
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok2, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on ECM.. line"); /* mu_extrap */
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok3, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on ECM.. line"); /* mu_orig   */
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok4, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on ECM.. line"); /* dbsize    */
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok5, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on ECM.. line"); /* nrandhits */
	if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok6, NULL))   != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on ECM.. line"); /* tailp     */
	if (cm->expA == NULL) { 
	  ESL_ALLOC(cm->expA, sizeof(ExpInfo_t *) * EXP_NMODES);
	  for(x = 0; x < EXP_NMODES; x++) { cm->expA[x] = CreateExpInfo(); }
	}
	cm->expA[exp_mode]->lambda    = atof(tok1);
	cm->expA[exp_mode]->mu_extrap = atof(tok2);
	cm->expA[exp_mode]->mu_orig   = atof(tok3);
	cm->expA[exp_mode]->dbsize    = atof(tok4); /* store as double, even though it was written as a long */
	cm->expA[exp_mode]->nrandhits = atoi(tok5);
	cm->expA[exp_mode]->tailp     = atof(tok6);
	cm->expA[exp_mode]->is_valid  = TRUE;
      }
      else if (strcmp(tag, "CM") == 0) {  
	/* skip the remainder of this line */
	if ((status = esl_fileparser_NextLine(cmfp->efp)) != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Premature end of data before main model section");
	break;
      }
    } /* end, loop over possible header tags */
  if (status != eslOK) goto ERROR;

  /* Done reading the header information.
   * Check that everything is ok and mandatory info is present before moving on.
   */
  if (cm->M     < 1)      ESL_XFAIL(status, cmfp->errbuf, "Failed to read STATES line in header section");
  if (cm->nodes < 1)      ESL_XFAIL(status, cmfp->errbuf, "Failed to read NODES line in header section");
  if (cm->clen  < 1)      ESL_XFAIL(status, cmfp->errbuf, "Failed to read CLEN line in header section");
  if (cm->W     < 1)      ESL_XFAIL(status, cmfp->errbuf, "Failed to read W line in header section");
  if (! read_el_selfsc)   ESL_XFAIL(status, cmfp->errbuf, "Failed to read ELSELF line in header section");
  if (! read_fp7_stats)   ESL_XFAIL(status, cmfp->errbuf, "Failed to read EFP7GF line in header section");
  if (cm->name == NULL)   ESL_XFAIL(status, cmfp->errbuf, "Failed to read NAME line in header section");
  if (abc      == NULL)   ESL_XFAIL(status, cmfp->errbuf, "Failed to read ALPH line in header section");
  if (tmp_null == NULL)   ESL_XFAIL(status, cmfp->errbuf, "Failed to read NULL line in header section");

  /* Check to make sure we parsed CM E-value stats correctly. 
   */
  if (cm->expA != NULL) { 
    if      (cm_statstracker == 15) cm->flags |= CMH_EXPTAIL_STATS;
    else if (cm_statstracker != 0)  ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Missing one or more ECM.. parameter lines");
  }

  /* Allocate body of CM now that # states (M) and # nodes (nnodes) are known */
  CreateCMBody(cm, cm->nodes, cm->M, cm->clen, abc);

  /* Copy values we stored in temp parameters awaiting CM allocation in CreateCMBody() */
  esl_vec_FCopy(tmp_null, abc->K, cm->null); /* cm->null allocated in CreateCMBody */
  cm->qdbinfo->beta1 = tmp_qdbbeta1;
  cm->qdbinfo->beta2 = tmp_qdbbeta2;
    
  /* Allocate and initialize temporary parameters for values we can't store 
   * until after we've  read the full model and are able to construct an emit map */
  ESL_ALLOC(tmp_map_left,  sizeof(int) * cm->nodes);
  ESL_ALLOC(tmp_map_right, sizeof(int) * cm->nodes);
  esl_vec_ISet(tmp_map_left,  cm->nodes, -1);
  esl_vec_ISet(tmp_map_right, cm->nodes, -1);

  ESL_ALLOC(tmp_rf_left,  sizeof(char) * (cm->nodes+1));
  ESL_ALLOC(tmp_rf_right, sizeof(char) * (cm->nodes+1));
  tmp_rf_left[cm->nodes]  = '\0';
  tmp_rf_right[cm->nodes] = '\0';

  ESL_ALLOC(tmp_cons_left,  sizeof(char) * (cm->nodes+1));
  ESL_ALLOC(tmp_cons_right, sizeof(char) * (cm->nodes+1));
  tmp_cons_left[cm->nodes]  = '\0';
  tmp_cons_right[cm->nodes] = '\0';

  nd = -1;
  cm->clen = 0;

  for (v = 0; v < cm->M; v++)
    {
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))     != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Premature end of data before main model section");
      
      /* Ah, a node line. Process it and get the following line.
       */
      if (*tok1 == '[') 
	{
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on node line: expected %d, got %d", 10, 1);
	  if ((x = NodeCode(tok1)) == -1)                                                ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Invalid node type %s", tok1);             
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on node line: expected %d, got %d", 10, 2);
	  if (!is_integer(tok1))                                                         ESL_XFAIL(status,     cmfp->errbuf, "Invalid node index on node line: should be integer >= 0, not %s", tok1);
	  nd = atoi(tok1);     
	  if (nd <  0)                                                                   ESL_XFAIL(status,     cmfp->errbuf, "Invalid node index on node line: should be integer >= 0, not %s", tok1);
	  if (nd >= cm->nodes)                                                           ESL_XFAIL(status,     cmfp->errbuf, "Invalid node index on node line: should not exceed %d, read %s", cm->nodes, tok1);
	  cm->ndtype[nd]  = x;
	  if     (cm->ndtype[nd] == MATP_nd) cm->clen+=2;
	  else if(cm->ndtype[nd] == MATL_nd) cm->clen++;
	  else if(cm->ndtype[nd] == MATR_nd) cm->clen++;
	  cm->nodemap[nd] = v;

	  /* chew up ']' */
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on node line: expected %d, got %d", 10, 3);

	  /* read annotation: MAP, consensus sequence and RF. Proper format depends on node type. 
	   */
	  /* MAP (optional: CMH_MAP? yes, else no */
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on node line: expected %d, got %d", 10, 4);
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok2, NULL)) != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on node line: expected %d, got %d", 10, 5);
	  if      ((cm->flags & CMH_MAP) && cm->ndtype[nd] == MATP_nd) { 
	    if (!is_integer(tok1))                                                       ESL_XFAIL(status,     cmfp->errbuf, "Invalid 1st MAP value on MATP node line: should be positive integer, not %s", tok1);
	    if (!is_integer(tok2))                                                       ESL_XFAIL(status,     cmfp->errbuf, "Invalid 2nd MAP value on MATP node line: should be positive integer, not %s", tok2);
	    tmp_map_left[nd]  = atoi(tok1);
	    tmp_map_right[nd] = atoi(tok2);
	  }
	  else if((cm->flags & CMH_MAP) && cm->ndtype[nd] == MATL_nd) { 
	    if (!is_integer(tok1))                                                       ESL_XFAIL(status,     cmfp->errbuf, "Invalid 1st MAP value on MATL node line: should be positive integer, not %s", tok1);
	    if (*tok2 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 2nd MAP value on MATL node line: should be '-', not %s", tok2);
	    tmp_map_left[nd]  = atoi(tok1);
	    tmp_map_right[nd] = -1;
	  }
	  else if((cm->flags & CMH_MAP) && cm->ndtype[nd] == MATR_nd) { 
	    if (*tok1 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 1st MAP value on MATR node line: should be '-', not %s", tok1);
	    if (!is_integer(tok2))                                                       ESL_XFAIL(status,     cmfp->errbuf, "Invalid 2nd MAP value on MATR node line: should be positive integer, not %s", tok2);
	    tmp_map_left[nd]  = -1;
	    tmp_map_right[nd] = atoi(tok2);
	  }
	  else { /* either (! (cm->flags & CMH_MAP)) or ndtype is not MATP, MATL nor MATR */
	    if (*tok1 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 1st MAP value on node line: should be '-', not %s", tok1);
	    if (*tok2 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 2nd MAP value on node line: should be '-', not %s", tok2);
	    tmp_map_left[nd]  = -1;
	    tmp_map_right[nd] = -1;
	  }

	  /* consensus sequence (optional: CMH_CONS? yes, else no */
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on node line: expected %d, got %d", 10, 6);
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok2, NULL)) != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on node line: expected %d, got %d", 10, 7);
	  if     ((cm->flags & CMH_CONS) && (cm->ndtype[nd] == MATP_nd)) { 
	    tmp_cons_left[nd]  = *tok1; 
	    tmp_cons_right[nd] = *tok2; 
	  }
	  else if((cm->flags & CMH_CONS) && cm->ndtype[nd] == MATL_nd) { 
	    if (*tok2 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 2nd consensus character on MATL node line: should be '-', not %s", tok2);
	    tmp_cons_left[nd]  = *tok1; 
	    tmp_cons_right[nd] = *tok2; 
	  }
	  else if((cm->flags & CMH_CONS) && cm->ndtype[nd] == MATR_nd) { 
	    if (*tok1 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 1st consensus character on MATR node line: should be '-', not %s", tok1);
	    tmp_cons_left[nd]  = *tok1; 
	    tmp_cons_right[nd] = *tok2; 
	  }
	  else { /* either (! (cm->flags & CMH_CONS) or ndtype is not MATP, MATL nor MATR */
	    if (*tok1 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 1st consensus character on node line: should be '-', not %s", tok1);
	    if (*tok2 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 2nd consensus character on node line: should be '-', not %s", tok2);
	    tmp_cons_left[nd]  = *tok1; 
	    tmp_cons_right[nd] = *tok2; 
	  }
	  
	  /* RF (optional: CMH_RF? yes, else no */
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on node line: expected %d, got %d", 10, 8);
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok2, NULL)) != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on node line: expected %d, got %d", 10, 9);
	  if      ((cm->flags & CMH_RF) && cm->ndtype[nd] == MATP_nd) { 
	    tmp_rf_left[nd]  = *tok1; 
	    tmp_rf_right[nd] = *tok2; 
	  }
	  else if((cm->flags & CMH_RF) && cm->ndtype[nd] == MATL_nd) { 
	    if (*tok2 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 2nd RF character on MATL node line: should be '-', not %s", tok2);
	    tmp_rf_left[nd]  = *tok1; 
	    tmp_rf_right[nd] = *tok2; 
	  }
	  else if((cm->flags & CMH_RF) && cm->ndtype[nd] == MATR_nd) { 
	    if (*tok1 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 1st RF character on MATR node line: should be '-', not %s", tok1);
	    tmp_rf_left[nd]  = *tok1; 
	    tmp_rf_right[nd] = *tok2; 
	  }
	  else { /* either (! (cm->flags & CMH_RF)) or ndtype is not MATP, MATL nor MATR */
	    if (*tok1 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 1st RF character on node line: should be '-', not %s", tok1);
	    if (*tok2 != '-')                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid 2nd RF character on node line: should be '-', not %s", tok2);
	    tmp_rf_left[nd]  = *tok1; 
	    tmp_rf_right[nd] = *tok2; 
	  }
	  if ((status = esl_fileparser_NextLine(cmfp->efp))                    != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Premature end of data in main model: no state %d line", v);
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state line: expected at least %d, got %d", 8, 0);
	} /* done with node line */
      
      /* Process state line.
       * <statecode> <v> <plast> <pnum> <cfirst> <cnum> <dmin2> <dmin1> <dmax1> <dmax2> <transition probs (variable number)> <emission probs (variable number)>
       */

      /* <statecode> */
      if((cm->sttype[v] = StateCode(tok1)) == -1)                                        ESL_XFAIL(status,     cmfp->errbuf, "Invalid state type %s\n", tok1);
      
      /* <v> */
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected at least %d, got %d", v, 8, 1);
      if (! is_integer(tok1))                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid state index on state line: should be integer >= 0, not %s", tok1);
      if (atoi(tok1) != v)                                                               ESL_XFAIL(status,     cmfp->errbuf, "Invalid state index on state line: should be %d, not %s", v, tok1);
      
      /* <plast> */
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected at least %d, got %d", v, 8, 2);
      if (! is_integer(tok1))                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid plast value on state line: should be integer, not %s", tok1);
      cm->plast[v] = atoi(tok1);
      
      /* <pnum> */
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected at least %d, got %d", v, 8, 3);
      if (! is_integer(tok1))                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid pnum value on state line: should be integer, not %s", tok1);
      cm->pnum[v] = atoi(tok1);
      
      /* <cfirst> */
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected at least %d, got %d", v, 8, 4);
      if (! is_integer(tok1))                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid cfirst value on state line: should be integer, not %s", tok1);
      cm->cfirst[v] = atoi(tok1);
      
      /* <cnum> */
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected at least %d, got %d", v, 8, 5);
      if (! is_integer(tok1))                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid cnum value on state line: should be integer, not %s", tok1);
      cm->cnum[v] = atoi(tok1);
      
      /* <dmin2> */
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected at least %d, got %d", v, 8, 6);
      if (! is_integer(tok1))                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid dmin2 value on state line: should be integer, not %s", tok1);
      cm->qdbinfo->dmin2[v] = atoi(tok1);

      /* <dmin1> */
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected at least %d, got %d", v, 8, 6);
      if (! is_integer(tok1))                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid dmin1 value on state line: should be integer, not %s", tok1);
      cm->qdbinfo->dmin1[v] = atoi(tok1);
      
      /* <dmax1> */
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected at least %d, got %d", v, 8, 7);
      if (! is_integer(tok1))                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid dmax1 value on state line: should be integer, not %s", tok1);
      cm->qdbinfo->dmax1[v] = atoi(tok1);

      /* <dmax2> */
      if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected at least %d, got %d", v, 8, 7);
      if (! is_integer(tok1))                                                            ESL_XFAIL(status,     cmfp->errbuf, "Invalid dmax2 value on state line: should be integer, not %s", tok1);
      cm->qdbinfo->dmax2[v] = atoi(tok1);

      /* Transition probabilities. */
      if (cm->sttype[v] != B_st) {
	for (x = 0; x < cm->cnum[v]; x++) {
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected at least %d, got %d", v, 8+cm->cnum[v], 8+x);
	  if ((! is_real(tok1) && (*tok1 != '*')))                                           ESL_XFAIL(status,     cmfp->errbuf, "Invalid transition score %d on state line: should be real number or '*', not %s", x+1, tok1);
	  cm->t[v][x] = ascii2prob(tok1, 1.);
	}
      }
      /* Emission probabilities. */
      if (cm->sttype[v] == ML_st || cm->sttype[v] == MR_st || cm->sttype[v] == IL_st || cm->sttype[v] == IR_st) {
	for (x = 0; x < cm->abc->K; x++) { 
	  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)     ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected %d, got %d", v, 8+cm->cnum[v]+cm->abc->K, 8+cm->cnum[v]+x);
	  if ((! is_real(tok1) && (*tok1 != '*')))                                           ESL_XFAIL(status,     cmfp->errbuf, "Invalid emission score %d on state line: should be real number or '*', not %s", x+1, tok1);
	  cm->e[v][x] = ascii2prob(tok1, cm->null[x]);
	}
      }
      else if (cm->sttype[v] == MP_st) {
	for (x = 0; x < cm->abc->K; x++) {
	  for (y = 0; y < cm->abc->K; y++) { 
	    if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL)) != eslOK)   ESL_XFAIL(status,     cmfp->errbuf, "Too few fields on state %d line: expected %d, got %d", v, 8+cm->cnum[v]+(cm->abc->K*cm->abc->K), 8+cm->cnum[v]+(x*cm->abc->K)+y);
	    if ((! is_real(tok1) && (*tok1 != '*')))                                         ESL_XFAIL(status,     cmfp->errbuf, "Invalid emission score %d on state line: should be real number or '*', not %s", (x*cm->abc->K)+y+1, tok1);
	    cm->e[v][x*cm->abc->K+y] = ascii2prob(tok1, cm->null[x]*cm->null[y]);
	  } 
	}
      }
      cm->ndidx[v] = nd;
      cm->stid[v]  = DeriveUniqueStateCode(cm->ndtype[nd], cm->sttype[v]);
      if ((status = esl_fileparser_NextLine(cmfp->efp))                            != eslOK) ESL_XFAIL(status,     cmfp->errbuf, "Premature end of data in main model: no state %d line", v+1);
    } /* end of loop over states */
  /* The closing // */
  if ((status = esl_fileparser_GetTokenOnLine(cmfp->efp, &tok1, NULL))        != eslOK)  ESL_XFAIL(status,     cmfp->errbuf, "Premature end of data: missing //?");
  if (strcmp(tok1, "//")                                                      != 0)      ESL_XFAIL(eslEFORMAT, cmfp->errbuf, "Expected closing //; found %s instead", tok1);

  /* Finally, read the filter HMM for this CM, unless we're explicitly told not to */
  if (read_fp7) {
	  
	/* deleted some code that's committed to files */

	/* calling p7_hmmfile_Read is really all that 'cm_p7_hmmfile_Read' (in cm_file.c) does anyway */
	if ((status = p7_hmmfile_Read(cmfp->hfp, &abc, &hmm))!=eslOK) goto ERROR;

	/* this line is the same as original code */
	if((status = cm_SetFilterHMM(cm, hmm, tmp_fp7_gfmu, tmp_fp7_gflambda)) != eslOK) ESL_XFAIL(status, cmfp->errbuf, "Unable to set filter HMM for CM");

	/* and deleted more code that's commmited to files */
  } 

  CMRenormalize(cm);
  cm->qdbinfo->setby = CM_QDBINFO_SETBY_CMFILE;
  cm->W_setby        = CM_W_SETBY_CMFILE;

  /* Create emit map now that we know the model architecture */
  cm->emap = CreateEmitMap(cm);
  if(cm->emap == NULL) ESL_XFAIL(eslEINVAL, cmfp->errbuf, "After reading complete model, failed to create an emit map");

  /* Use emit map to map the per-node annotation to per-consensus position */
  if (cm->flags & CMH_RF) { 
    cm->rf[0] = ' ';
    for(nd = 0; nd < cm->nodes; nd++) { 
      if(cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd) cm->rf[cm->emap->lpos[nd]] = tmp_rf_left[nd];
      if(cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATR_nd) cm->rf[cm->emap->rpos[nd]] = tmp_rf_right[nd];
    }
    cm->rf[cm->clen+1] = '\0';
  }
  if (cm->flags & CMH_CONS) { 
    cm->consensus[0] = ' ';
    for(nd = 0; nd < cm->nodes; nd++) { 
      if(cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd) cm->consensus[cm->emap->lpos[nd]] = tmp_cons_left[nd];
      if(cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATR_nd) cm->consensus[cm->emap->rpos[nd]] = tmp_cons_right[nd];
    }
    cm->consensus[cm->clen+1] = '\0';
  }
  if (cm->flags & CMH_MAP) { 
    cm->map[0] = 0;
    for(nd = 0; nd < cm->nodes; nd++) { 
      if(cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATL_nd) cm->map[cm->emap->lpos[nd]] = tmp_map_left[nd];
      if(cm->ndtype[nd] == MATP_nd || cm->ndtype[nd] == MATR_nd) cm->map[cm->emap->rpos[nd]] = tmp_map_right[nd];
    }
  }

  if(tmp_null  != NULL) free(tmp_null);

  /* these get allocated regardless of flag status, free them */
  free(tmp_rf_left);
  free(tmp_rf_right);
  free(tmp_cons_left);
  free(tmp_cons_right);
  free(tmp_map_left);
  free(tmp_map_right);

  if (*ret_abc == NULL) *ret_abc = abc;
  if ( opt_cm != NULL)  *opt_cm = cm; else FreeCM(cm);
  return eslOK;

 ERROR:
  if (*ret_abc == NULL && abc != NULL) esl_alphabet_Destroy(abc);
  if (cm       != NULL) FreeCM(cm);
  if (opt_cm   != NULL) *opt_cm = NULL;
  if      (status == eslEMEM || status == eslESYS) return status; 
  else if (status == eslEOF)                       return status;
  else if (status == eslEINCOMPAT)                 return status;
  else                                             return eslEFORMAT;	/* anything else is a format error: includes premature EOF, EOL, EOD  */
}

static int cm_param_to_hmm_param[CM_p7_NEVPARAM]={p7_MMU,p7_MLAMBDA,p7_VMU,p7_VLAMBDA,p7_FTAU,p7_FLAMBDA,-1,-1};

CM_t * /*new CM copy*/ use_cm_Clone (CM_t *in_cm)
{
	int status;
	char errbuf[eslERRBUFSIZE];
	CM_t *out_cm;
    if((status = cm_Clone(in_cm, errbuf, &out_cm)) != eslOK) cm_Fail("unable to clone CM");
	return out_cm;
}

CM_t * /*new CM copy*/ save_and_load_CM_via_memory_copying (CM_t *in_cm,CmfinderVars *cmfinderVars)
{
	double uninitialized;
	int i;
	int format=-1;
	CM_FILE *cmfp;
	file_as_string f,hmm_f;
	CM_t *out_cm;

#if 0
	int status;
	FILE *cmoutfp=NULL;
	if ((cmoutfp = fopen("normal-save.cm", "w")) == NULL) esl_fatal("Failed to open CM file for writing");
	if ((status = cm_file_WriteASCII(cmoutfp, -1, in_cm)) != eslOK) esl_fatal("CM save failed");
	fclose(cmoutfp);
#endif

	/* for valgrind's benefit, give values to un-calibrated hmm params, and then afterwards set them to undefined values */
	for (i=0; i<CM_p7_NEVPARAM; i++) {
		if (!cmfinderVars->hmmStatsToCalibrate[i]) {
			in_cm->fp7_evparam[i]=0; /* 0 is a totally arbitrary value */
			if (cm_param_to_hmm_param[i]!=-1) {
				in_cm->fp7->evparam[cm_param_to_hmm_param[i]]=0;
			}
		}
	}
	
	file_as_string_open(&f);
	file_as_string_open(&hmm_f);
	if (cm_file_WriteASCII_file_as_string(&f,format,in_cm,cmfinderVars)!=eslOK) {
		esl_fatal("problem in save_and_load_CM_via_memory: cm_file_WriteASCII_file_as_string");
	}
	if(in_cm->flags & CMH_FP7 && in_cm->fp7 != NULL) {
		file_as_string_p7_hmmfile_WriteASCII (&hmm_f,format,in_cm->fp7);
	}

#if 0
	FILE *memfp=NULL;
	memfp=fopen("mem-save.cm","w");
	if (memfp==NULL) {esl_fatal("fopen failed");}
	fprintf(memfp,"%s",f.p);
	fprintf(memfp,"%s",hmm_f.p);
	fclose(memfp);
#endif

	if (cm_file_OpenBuffer(f.p,f.n,0,&cmfp)!=eslOK) {
		esl_fatal("problem in save_and_load_CM_via_memory: cm_file_OpenBuffer");
	}
	if (p7_hmmfile_OpenBuffer(hmm_f.p,hmm_f.n,&(cmfp->hfp))!=eslOK) {
		esl_fatal("problem in save_and_load_CM_via_memory: p7_hmmfile_OpenBuffer");
	}
	if (read_asc_1p1_cm__with_RAM(cmfp,TRUE,&abc,&out_cm)!=eslOK) {  /* use the internal function that we modified */
		esl_fatal("problem in save_and_load_CM_via_memory: cm_file_Read");
	}
	cm_file_Close(cmfp);

	file_as_string_close(&f);
	file_as_string_close(&hmm_f);

	if (out_cm==(CM_t *)1) { /* won't happen */
		uninitialized=20; /* this code stops the compiler from complaining that 'uninitialized' is uninitialized, but it will always be uninitialized */
	}
#if 0 /* positive control that valgrind flags this */
	if (uninitialized>50) {
		printf("something just to make valgrind upset\n");
	}
#endif
	for (i=0; i<CM_p7_NEVPARAM; i++) {
		if (!cmfinderVars->hmmStatsToCalibrate[i]) {
			out_cm->fp7_evparam[i]=uninitialized;
			if (cm_param_to_hmm_param[i]!=-1) {
				in_cm->fp7->evparam[cm_param_to_hmm_param[i]]=uninitialized;
			}
		}
	}
	
	return out_cm;
}

void save_and_load_CM_via_memory_freeing_old (CM_t **ref_cm,CmfinderVars *cmfinderVars)
{
	CM_t *in_cm,*out_cm;
	in_cm=*ref_cm;
	out_cm=save_and_load_CM_via_memory_copying (in_cm,cmfinderVars);
	FreeCM(in_cm);
	*ref_cm=out_cm;
}

void AllocSSconsOrRf(char **p,ESL_MSA *msa)
{
	/* don't free(*p), because at least when called from AddColumnsToMsa for msa->pp, it was an invalid allocation */
	*p=(char *)MallocOrDie(sizeof(char) * (msa->alen+1));
	(*p)[msa->alen] = '\0';
}

void AllocSSconsOrRfIfNotThere(char **p,ESL_MSA *msa)
{
	if (*p==NULL) {
		AllocSSconsOrRf(p,msa);
	}
}

void SetTextColumnarDataToGaps(ESL_MSA *msa,char *s)
{
	int c;
	if (s!=NULL) {
		for (c=0; c<msa->alen; c++) {
			s[c]='.';
		}
	}
}
void CopyTextColumnarDataWithColumnMap(ESL_MSA *msa,ESL_MSA *new_msa,const char *src,char *dest,int *oldToNewColumnMap)
{
	if (src!=NULL) {
		int oldcol;
		for (oldcol=0; oldcol<msa->alen; oldcol++) {
			int newcol=oldToNewColumnMap[oldcol];
			dest[newcol]=src[oldcol];
		}
	}
}
void AddColumnsToMsa (ESL_MSA **get_new_msa,ESL_MSA *msa,ColumnsToAdd *columnsToAdd,int numColumnsToAddEntries,int **get_oldToNewColumnMap,int *get_columnsAdded)
{
	int status;
	int i,j;
	int numColumnsToAdd;
	int e,col,seqnum,newcol;
	int *oldToNewColumnMap;
	
	if ((msa->flags & eslMSA_DIGITAL)==0) {
		esl_fatal("msa should be digital");
	}
	
	numColumnsToAdd=0;
	for (e=0; e<numColumnsToAddEntries; e++) {
		numColumnsToAdd += columnsToAdd[e].numCols;
	}
	*get_columnsAdded=numColumnsToAdd;

	if (numColumnsToAdd==0) {
		*get_new_msa=NULL;
		return;
	}
	
	/* create map */
	oldToNewColumnMap=(int *)MallocOrDie(sizeof(int)*(msa->alen+1));
	e=0;
	newcol=0;
	for (col=0; col<msa->alen+1; col++) {
		if (e<numColumnsToAddEntries) {
			if (e>0) { assert(columnsToAdd[e].colAfterInserts > columnsToAdd[e-1].colAfterInserts); } /* we assume they're sorted by colAfterInserts to make this work.  If it's a hassle for the caller to provide this, we should qsort it ourselves. */
			if (col==columnsToAdd[e].colAfterInserts) {
				newcol += columnsToAdd[e].numCols;
				e++;
			}
		}
		oldToNewColumnMap[col]=newcol;
		newcol++;
	}

	/* create new alignment, then for each field, initialize it and copy over old data */
	int new_alen=msa->alen+numColumnsToAdd;
	ESL_MSA *new_msa=esl_msa_CreateDigital(abc,msa->nseq,new_alen);
	if (new_msa==NULL) {
		esl_fatal("esl_msa_CreateDigital failed");
	}
	/* stuff copied/adapted from esl_msa_Copy follows.  We can't use esl_msa_Copy exactly, because it needs to extend the size of column-specific fields */
	/* also, while we're at it, initialize everything to gaps */
	for (i = 0; i < msa->nseq; i++) {
		esl_strdup(msa->sqname[i], -1, &(new_msa->sqname[i]));
		new_msa->wgt[i] = msa->wgt[i];
	}
	new_msa->flags=msa->flags;
	esl_strdup(msa->name,    -1, &(new_msa->name));
	esl_strdup(msa->desc,    -1, &(new_msa->desc));
	esl_strdup(msa->acc,     -1, &(new_msa->acc));
	esl_strdup(msa->au,      -1, &(new_msa->au));
	if (msa->ss_cons!=NULL) {
		AllocSSconsOrRf(&(new_msa->ss_cons),new_msa);
		SetTextColumnarDataToGaps(new_msa,new_msa->ss_cons);
		CopyTextColumnarDataWithColumnMap(msa,new_msa,msa->ss_cons,new_msa->ss_cons,oldToNewColumnMap);
	}
	if (msa->pp_cons!=NULL) {
		AllocSSconsOrRf(&(new_msa->pp_cons),new_msa);
		SetTextColumnarDataToGaps(new_msa,new_msa->pp_cons);
		CopyTextColumnarDataWithColumnMap(msa,new_msa,msa->pp_cons,new_msa->pp_cons,oldToNewColumnMap);
	}
	if (msa->rf!=NULL) {
		AllocSSconsOrRf(&(new_msa->rf),new_msa);
		SetTextColumnarDataToGaps(new_msa,new_msa->rf);
		CopyTextColumnarDataWithColumnMap(msa,new_msa,msa->rf,new_msa->rf,oldToNewColumnMap);
	}
	if (msa->sqacc != NULL) {
		ESL_ALLOC(new_msa->sqacc, sizeof(char **) * msa->nseq);
		for (i = 0; i < msa->nseq; i++)
			esl_strdup(msa->sqacc[i], -1, &(new_msa->sqacc[i]));
	}
	if (msa->sqdesc != NULL) {
		ESL_ALLOC(new_msa->sqdesc, sizeof(char **) * msa->nseq);
		for (i = 0; i < msa->nseq; i++)
			esl_strdup(msa->sqdesc[i], -1, &(new_msa->sqdesc[i]));
	}
	if (msa->ss != NULL) {
		ESL_ALLOC(new_msa->ss, sizeof(char **) * msa->nseq);
		for (i = 0; i < msa->nseq; i++) {
			AllocSSconsOrRf(&(new_msa->ss[i]),new_msa);
			SetTextColumnarDataToGaps(new_msa,new_msa->ss[i]);
			CopyTextColumnarDataWithColumnMap(msa,new_msa,msa->ss[i],new_msa->ss[i],oldToNewColumnMap);
		}
	}
	if (msa->sa != NULL) {
		esl_fatal("msa->sa copying not implemented");
	}
	if (msa->pp != NULL) {
		ESL_ALLOC(new_msa->pp, sizeof(char **) * msa->nseq);
		for (i = 0; i < msa->nseq; i++) {
			AllocSSconsOrRf(&(new_msa->pp[i]),new_msa);
			SetTextColumnarDataToGaps(new_msa,new_msa->pp[i]);
			CopyTextColumnarDataWithColumnMap(msa,new_msa,msa->pp[i],new_msa->pp[i],oldToNewColumnMap);
		}
	}
	if (msa->ngf > 0) {
		ESL_ALLOC(new_msa->gf_tag, sizeof(char **) * msa->ngf);
		ESL_ALLOC(new_msa->gf,     sizeof(char **) * msa->ngf);
		new_msa->ngf       = msa->ngf;
		new_msa->alloc_ngf = msa->ngf;
		for (i = 0; i < msa->ngf; i++) {
			esl_strdup(msa->gf_tag[i], -1, &(new_msa->gf_tag[i]));
			esl_strdup(msa->gf[i],     -1, &(new_msa->gf[i]));
		}
	}
	if (msa->ngs > 0) {
		ESL_ALLOC(new_msa->gs_tag, sizeof(char **)  * msa->ngs);
		ESL_ALLOC(new_msa->gs,     sizeof(char ***) * msa->ngs);
		new_msa->ngs       = msa->ngs;
		for (i = 0; i < msa->ngs; i++) {
			ESL_ALLOC(new_msa->gs[i], sizeof(char **) * msa->nseq);
			esl_strdup(msa->gs_tag[i], -1, &(new_msa->gs_tag[i]));
			for (j = 0; j < msa->nseq; j++)
				esl_strdup(msa->gs[i][j],  -1, &(new_msa->gs[i][j]));
		}
	}
	if (msa->ngc > 0) {
		esl_fatal("extra #=GC not implemented");
	}
	if (msa->ngr > 0) {
		esl_fatal("extra #=GR not implemented");
	}

	/* now set up the actual MSA (in digital mode) */
	/* set everything in MSA to gaps */
	for (col=0; col<new_msa->alen; col++) {
		for (seqnum=0; seqnum<new_msa->nseq; seqnum++) {
			new_msa->ax[seqnum][col+1]=esl_abc_XGetGap(abc);
		}
	}
	/* copy over MSA */
	for (col=0; col<msa->alen; col++) {
		newcol=oldToNewColumnMap[col];
		for (seqnum=0; seqnum<msa->nseq; seqnum++) {
			new_msa->ax[seqnum][newcol+1]=msa->ax[seqnum][col+1];
		}
	}

	*get_new_msa=new_msa;
	*get_oldToNewColumnMap=oldToNewColumnMap;
	
	return;
	
ERROR:
	esl_fatal("error");
}
