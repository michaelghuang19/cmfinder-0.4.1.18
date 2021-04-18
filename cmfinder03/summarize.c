#include "funcs.h"
#include "squid.h"
#include "structs.h"
#include "msa.h"
#include "fold.h"
#include <string.h>
#include "global.h"


static int min_block_size = 4;
static double gap_threshold = 0.6;
static double cons_threshold = 0.7;
static int min_seq = 4;
static double **bp_pr=NULL;


double entropy(char **dsq, char* r, int num, int alen, double* prior, float* weight, float tot_weight)
{
  int    i, idx;  
  double total_e, e;  
  double f[Alphabet_size];        /* singlet frequency vector            */
  int    sym, c;

  total_e = 0;  
  for (i = 0; i < alen; i++) {    
    for (sym = 0; sym <Alphabet_size; sym++)
      f[sym] = 0;    
    
    for (idx = 0; idx < num; idx++){
      if (dsq[idx][i+1] == DIGITAL_GAP){
	for (sym = 0; sym < Alphabet_size; sym++)
	  f[sym] += prior[sym] * weight[idx];	
      }
      else{	
	sym = dsq[idx][i+1];
	f[sym] += weight[idx];
      }      
    }    
    e = 0;    
    for (sym = 0; sym < Alphabet_size; sym++){
      f[sym] /= tot_weight;            
      e -=  prior[sym]*log(prior[sym]);
      if (f[sym] > 0) 
	e +=  f[sym] * log(f[sym]);            
    }    
    total_e += e *  1.44269504;        
  }
  return total_e;  
}



double mxy(char* cs, char ** dsq, int num, int alen, double* prior, float* weight, float tot_weight)
{
  int     i, j, idx; 
  double  fx[Alphabet_size];        /* singlet frequency vector            */
  double  fy[Alphabet_size];	/* another singlet frequency vector    */
  double  fxy[Alphabet_size][Alphabet_size]; /* pairwise frequency 2D array    */
  int     symi, symj;		/* counters for symbols                */
  float   pairs;		/* counter for pairs in which there are no gaps */  
  double  total_mxy=0;
  double  mxy;  
  int* pt = GetPairtable(cs);  
  
  
  for (i = 0; i < alen; i++) {
    if (pt[i] <= i) continue;
    j = pt[i];
    /* zero counter array */
    for (symj = 0; symj < Alphabet_size; symj++)
      {
	fx[symj] = fy[symj] = 0.0;
	for (symi = 0; symi < Alphabet_size; symi++)
	  fxy[symj][symi] = 0.0;
      }
    /* count symbols in a column */
    for (idx = 0; idx < num; idx++)
      {
	if (dsq[idx][i+1]== DIGITAL_GAP){
	  if (dsq[idx][j+1] == 	 DIGITAL_GAP ){
	      for (symi = 0; symi < Alphabet_size; symi++){		
		for (symj = 0; symj < Alphabet_size; symj++)
		  fxy[symj][symi] += prior[symj] * prior[symi] * weight[idx];			      				
		fx[symi] += prior[symi] * weight[idx];
		fy[symi] += prior[symi] * weight[idx];		
	      }	      
	    }
	    else{     
	      symj = dsq[idx][j+1];	
	      fy[symj] += weight[idx];			      
	      for (symi = 0; symi < Alphabet_size; symi++){		
		fxy[symj][symi] += prior[symi] * weight[idx];		
		fx[symi] += prior[symi] *weight[idx];		
	      }	      
	    }
	}
	else{	  	    	
	  symi = dsq[idx][i+1];
	  fx[symi] += weight[idx];			      
	  if (dsq[idx][j+1] ==  DIGITAL_GAP ){	      
	    for (symj = 0; symj < Alphabet_size; symj++){
	      fxy[symj][symi] += prior[symj] * weight[idx];			      
	      fy[symj] += prior[symj] * weight[idx];		
	    }
	  }
	  else{     
	    symj = dsq[idx][j+1];	
	    fxy[symj][symi] += weight[idx];		
	    fy[symj] += weight[idx];		
	  }	  	  	  
	}	
      }    
    /* convert to frequencies */
    for (symi = 0; symi < Alphabet_size; symi++)
      {
	fx[symi] /=  tot_weight;
	fy[symi] /=  tot_weight;	
	for (symj = 0; symj < Alphabet_size; symj++)
	  fxy[symj][symi] /=  tot_weight;
      } 	       	  
    /* calculate mxy. 144.269504 is a conversion of ln's into
     * bits * 100: i.e. 100 * (1/log(2)) 
     */
    mxy = 0;
    for (symi = 0; symi < Alphabet_size; symi++)
      for (symj = 0; symj < Alphabet_size; symj++)
	{
	  if (fxy[symj][symi] > 0.0)
	    mxy +=  1.44269504 * fxy[symj][symi] *
	      log((fxy[symj][symi] / (fx[symi] * fy[symj])));	    	      
	}	
    if (mxy < -0.00001){
      Die("Error ! Column %d  %d mxy = %f", i, j, mxy);
    }    
    total_mxy += mxy;    
  }  
  return total_mxy;  
}

  

void  bad_base_pair(MSA* msa, double* ret_conflict_bp, double* ret_del_bp)
{
  int i,j;
  int* pt = GetPairtable(msa->ss_cons);  
  double conflict_bp = 0;  
  double del_bp = 0;  
  
  if (bp_pr== NULL){
    Avg_bppr(msa->aseq,msa->nseq, msa->alen, NULL, msa->wgt, &bp_pr);
  }
  for(i=0; i < msa->alen; i++) {
    if (pt[i] > i){
      for(j=0; j < msa->nseq; j++){
	char lc = msa->aseq[j][i];
	char rc = msa->aseq[j][pt[i]];	
	if(!IsBasePair(lc,rc)){
	  if (isgap(lc) && isgap(rc)) {
	    del_bp += msa->wgt[j] * bp_pr[pt[i]+1][i+1]; 
	  }
	  else{
	    conflict_bp += msa->wgt[j] * bp_pr[pt[i]+1][i+1]; 
	  }
	}
      }
    }
  }
  *ret_conflict_bp=conflict_bp;
  *ret_del_bp=del_bp; 
}

double  weighted_base_pair(MSA* msa)
{
  int    i;
  if (bp_pr== NULL){
    Avg_bppr(msa->aseq,msa->nseq, msa->alen, NULL, msa->wgt, &bp_pr);
  }
  int* pt = GetPairtable(msa->ss_cons);  
  double total_bp = 0;  
  for(i=0; i < msa->alen; i++) {
    if (pt[i] > i){
      total_bp += bp_pr[pt[i]+1][i+1];
    }
  }
  free(pt);
  return total_bp;
}



double  average_base_pair(MSA* msa)
{
  int    i;
  int* pt = GetPairtable(msa->ss_cons);  
  double total_bp = 0;  
  for(i=0; i < msa->alen; i++) {
    if (pt[i] > i){
      total_bp ++;
    }
  }
  free(pt);
  return total_bp;
}


double average_score(MSA *msa, float* weight, float tot_weight)
{
  float total_score = 0;
  float score;
  int    i;  
  int    start, end;  
  for(i=0; i < msa->nseq; i++) {
    if (msa->sqdesc && msa->sqdesc[i]){
      sscanf(msa->sqdesc[i], "%d..%d\t%f", &start, &end, &score);    
      total_score += score * weight[i];		    
    }    
  }
  return total_score / tot_weight ;  
}

double average_GC(char **aseq, int num, int alen, float* weight, float tot_weight)
{
  double total_GC = 0;
  int    i,j;
  int    GC_count;
  int    count;  
  double GC;  
  for(i=0; i < num; i++) {
    GC_count=0;
    count = 0;
    for(j=0; j < alen; j++) {      
      if(isgap(aseq[i][j])) continue;
      if(aseq[i][j] == 'G' || aseq[i][j] == 'C'){
	GC_count++;	
      }
      count++;      	
    }
    GC =  ((double) GC_count)/ count;        
    total_GC += GC * weight[i];       
  }  
  return (total_GC / tot_weight);    
}


double average_seq_id(char **aseq, int num, int alen, float* weight, float tot_weight)
{
  int i, j, k;
  double total_id = 0;  
  double id = 0;  
  int    len;
  double pair = 0;  
  for(i=1; i < num; i++) {
    for(j=0; j < i; j++) {
      id = 0;
      len=0;
      for(k=0; k < alen; k++) {
	if (isgap(aseq[i][k]) && isgap(aseq[j][k])){
	  continue;
	}
	len++;
	if (aseq[i][k] == aseq[j][k]) 
	  id++;	
      } 
      total_id += id / len * weight[i] * weight[j];      
      pair+= weight[i] * weight[j];      
    }    
  }   
  return total_id / pair;  
}

double average_seq_len(char **aseq, int nseq, float* weight, float tot_weight)
{
  double total_len=0;
  int    i;
  for(i=0; i < nseq; i++) {
    int length=0;  
    char *sp;    
    for(sp = aseq[i]; *sp ; sp++)
      if (!isgap(*sp)) length++;          
    total_len += length * weight[i];    
  }
  return total_len / tot_weight;  
}

double conserved_position(char **dsq, int nseq, int alen, float* weight, float tot_weight)
{
  int i,j,k;
  int block_count=0;
  int block_size=0;
  int accumulated_block_size = 0;
  double freq[4]={0,0,0,0};
  for(i=1; i <= alen; i++){
    int conserved = 0;
    int  gap_count = 0;
    for(j=0; j < 4; j++){
      freq[j]=0;
    }
    for(j=0; j < nseq; j++){      
      if (dsq[j][i] == DIGITAL_GAP){
	gap_count+= weight[j];
      }
      else{
	if (dsq[j][i] >= 4) Die("Invalid base at sequence %d pos %d", j, i);
	freq[dsq[j][i]]+= weight[j];
      }
    }
    if (gap_count / tot_weight > gap_threshold){
      continue;
    }
    for(j=0; j < 4; j++){
      if (freq[j] / tot_weight > cons_threshold && freq[j] >= 3){
	block_size ++;
	conserved = 1;
      }
    }
    if (!conserved) {
      if (block_size >= min_block_size){
	accumulated_block_size += block_size;
      }
      block_size = 0;
    }
  }
  if (block_size >= min_block_size){
    accumulated_block_size += block_size;
  }
  return accumulated_block_size;
}

double average_energy(MSA  *msa, float* weight, float tot_weight)
{
  char* ss;
  char* seq;
  char  *ssp, *sp1, *sp2;  
  int i,j;  
  double energy;  
  double tot_energy=0;  
  seq = (char *)malloc(sizeof(char) * (msa->alen + 1));     
  ss = (char *)malloc(sizeof(char) * (msa->alen + 1));     
  for(i=0; i < msa->nseq; i++) {        
    sp1 = seq;    
    sp2 = ss;    
    for(j=0; j < msa->alen; j++) {
      if (!isgap(msa->aseq[i][j])) {
	*(sp1++)= msa->aseq[i][j];	
	*sp2= msa->ss[i][j];	
	if (*sp2 == '<') *sp2= '(';
	if (*sp2 == '>') *sp2= ')';
	sp2++;	
      }      
    }
    *sp1=*sp2='\0';
    initialize_fold(strlen(seq));
    energy = energy_of_struct(seq, ss);
    /*printf("\n%d\t%s\n", i, seq);
      printf("\t%s\t",   ss);    
      printf("%f", energy);    */
    // Prevent the outliers 
    if (energy > 10) energy = 10;
    tot_energy += energy * weight[i];    
  }  
  free(seq);
  free(ss);  
  return (tot_energy / tot_weight);  
}

void ss2cons_ss(MSA* msa){
  int i, k;
  if (msa->ss_cons == NULL){
    msa->ss_cons = MallocOrDie(sizeof(char)* (msa->alen+1));
    for(i=0; i < msa->alen; i++){
      msa->ss_cons[i]= '.';
    }
    msa->ss_cons[i]='\0';
    for(i=0; i < msa->alen; i++){
      for(k=0; k < msa->nseq;k++){
	if (pair_left(msa->ss[k][i])){
	  msa->ss_cons[i]='<';
	  break;
	}
	if (pair_right(msa->ss[k][i])){
	  msa->ss_cons[i]='>';
	  break;
	}
      }
    }
  }  
}

void cons_ss2ss(MSA* msa){
  int i,j;
  if (msa->ss == NULL){
    int* pt = GetPairtable(msa->ss_cons);  
    msa->ss = MallocOrDie(sizeof(char *) * msa->nseq);

    for(i=0; i < msa->nseq;i++){
      msa->ss[i] = MallocOrDie(sizeof(char) * (msa->alen+1));
      strcpy(msa->ss[i], msa->ss_cons);
      for(j=0; j < msa->alen; j++){
	if (pt[j] > j){
	  if(!IsBasePair(msa->aseq[i][j],msa->aseq[i][pt[j]])) {
	    msa->ss[i][j] = msa->ss[i][pt[j]] = '.';
	  }
	}
      }      
    }
  }
}

static struct opt_s OPTIONS[] = {
  { "-g", TRUE, sqdARG_FLOAT},
  { "-w", TRUE, sqdARG_NONE},
}
;
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


static char usage[]  = "\
Usage: summarize [-options] <alignment> \n\
where options are:\n\
     -g <GC>  : Background GC content \n\
     -w       : Do GSC weighting \n";



int
main(int argc, char **argv)
{
  
  MSA       *msa;                /* alignment info                            */
  int       idx;
  
  double  e, m, avg_bp, avg_bppr, avg_score, avg_seq_id, avg_seq_len, avg_energy, avg_GC;  
  double  prior[4] = {0.25, 0.25, 0.25, 0.25};
  double  gc = 0.5;  

  char    *optname;                /* name of option found by Getopt()        */
  char    *optarg;                 /* argument found by Getopt()              */
  int     optind;                 /* index in argv[]                         */	

  float*      weight;
  float       tot_weight;  

  char        *aseqfile;          /* file contain the initial alignment of selected cand */
  MSAFILE     *afp = NULL;
  int          format= MSAFILE_STOCKHOLM;
  char        **msa_dsq;
  int         i,k;
  int         cons_pos;
  int         do_gsc_weight=0;
  double      conf_bp=0;
  double      del_bp=0;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if      (strcmp(optname, "-g") == 0) {
      gc  = atof(optarg);
      prior[2]=prior[3] = gc/ 2;
      prior[0]=prior[4] = (1-gc)/2;            
    }  
    else if (strcmp(optname, "-w")== 0){
      /*  GSC sequence weighting */
      do_gsc_weight = 1;      
    }    
  }
  if (argc - optind < 1)
    Die("%s\n", usage);
  aseqfile = argv[argc-1];
  
  
  /* file handle of initial alignment          */
  if ((afp = MSAFileOpen(aseqfile, format, NULL)) == NULL)
    Die("Alignment file %s could not be opened for reading", aseqfile);
  
  if ((msa = MSAFileRead(afp)) == NULL)
    {
      Die("Fail to read alignment from %s", aseqfile);      
    }
  MSAFileClose(afp);        

  if (!msa->ss_cons){
    ss2cons_ss(msa);
  }
  if (!msa->ss){
    cons_ss2ss(msa);
  }

  if (do_gsc_weight){
    float* gsc_weight =MallocOrDie(sizeof(float)* msa->nseq);    
    GSCWeights(msa->aseq, msa->nseq, msa->alen, gsc_weight);
    for(i=0; i < msa->nseq; i++){
      msa->wgt[i] *= gsc_weight[i];      
    }
    free(gsc_weight);
  }
  
  tot_weight=0;  
  for (i = 0; i < msa->nseq; i++){    
          PrepareSequence(msa->aseq[i],0);
    tot_weight += msa->wgt[i];    
  }
  msa_dsq = DigitizeAlignment(msa->aseq, msa->nseq, msa->alen);
  
  e= entropy( msa_dsq, msa->rf, msa->nseq, msa->alen, prior, msa->wgt, tot_weight);  
  m= mxy( msa->ss_cons, msa_dsq,  msa->nseq, msa->alen, prior,  msa->wgt, tot_weight);  
  avg_bppr= weighted_base_pair(msa);
  avg_bp= average_base_pair(msa); 
  avg_score= average_score(msa,  msa->wgt, tot_weight);  
  avg_seq_id = average_seq_id(msa->aseq, msa->nseq, msa->alen,  msa->wgt, tot_weight);  
  avg_seq_len = average_seq_len(msa->aseq, msa->nseq,  msa->wgt, tot_weight);  
  avg_energy = average_energy(msa, msa->wgt, tot_weight);  
  avg_GC = average_GC(msa->aseq, msa->nseq, msa->alen,  msa->wgt, tot_weight);  
  cons_pos = conserved_position(msa_dsq, msa-> nseq, msa-> alen, msa->wgt, tot_weight);
  bad_base_pair(msa, &conf_bp, &del_bp);
  printf("Num=%d\t Weight=%.2f\t Len=%.1f\t Score=%.2f\t Entropy=%.2f\t MI=%.2f\t BP=%.2f\t BP.org=%.2f\t Seq_id=%.2f\t Energy=%.2f\t GC=%.2f\t Conserved_pos=%d\t Conf_bp=%.2f\t Del_bp=%.2f\n", msa->nseq, tot_weight, avg_seq_len, avg_score, e, m, avg_bppr, avg_bp, avg_seq_id, avg_energy, avg_GC, cons_pos, conf_bp, del_bp); 
  MSAFree(msa);
  if (bp_pr){
    Free2DArray((void **)bp_pr, msa->alen+1);
  }
}

  
    
