#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <string.h>
#include "cand.h"
#include "energy_const.h"
#include "fold.h"
#include "part_func.h"
#include "fold_vars.h"
#include "utils.h"
#include "funcs.h"
#include "global.h"

#define MIN_ENERGY   0
#define LOCAL        3


static struct opt_s OPTIONS[] = {
  { "-c", TRUE, sqdARG_INT}, 
  { "-m", TRUE, sqdARG_INT}, 
  { "-M", TRUE, sqdARG_INT}, 
  { "-s", TRUE, sqdARG_INT}, 
  { "-S", TRUE, sqdARG_INT}, 
  { "-o", TRUE, sqdARG_STRING}, 
  { "-r", TRUE, sqdARG_STRING}, 
  { "-t", TRUE, sqdARG_STRING},
  { "-v", TRUE, sqdARG_NONE},
  { "--span-region-file", FALSE, sqdARG_STRING}  /* file has two integers separated by a comma, the first is the maximum value of the start of the candidate, the second number is the minimum value of the end.  the number pairs must be in the same order as the input sequences */
};

typedef struct {
  int maxStart,minEnd;
} SpanRegion;


#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

char usage[] = "\
usage:\n\
candf [-c max_cand] [-o output_file] [-m min_span] [-M max_span] [-s min_hairpin] [-S max_hairpin] [-r range_file] [-t timer-file] <seqfile>\n";

int verbose=0;

char* substr(char* sub, char* seq, int start, int end)
{
  int len = end - start + 1;  
  if (sub == NULL) {
    sub = (char*) space( len + 1);    
  }
  strncpy(sub, seq + start, len);
  sub[len] = '\0';  
  return sub;  
}

int isMaximal(int i, int j, int* tri_indx, int* energy, int range, int total_len)
{
  int k, l;
  int idx = tri_indx[j]  + i;
  int idx1;  
  for(k=i; k >= i-range && k >=0; k--) 
    for(l=j; l <= j+range && l <total_len; l++) {
      if (k==i && l == j) continue;
      idx1 = tri_indx[l] + k;
      if (energy[idx1] < energy[idx]){
	//printf("%d %d energy %d not maximal: %d %d %d\n", i, j, energy[idx], k, l, energy[idx]);
	return 0;
      }      
    }
  for(k=i; k < i+range &&  k < j; k++)
    for(l=j; l > j-range && l > k; l--) {
      if (k==i && l == j) continue;
      idx1 = tri_indx[l] + k;
      if (energy[idx1] <= energy[idx]){
	//printf("%d %d energy %d not maximal: %d %d %d\n", i, j, energy[idx], k, l, energy[idx]);
	return 0;
      }      
    }
  //printf("%d %d energy %d Maximal\n", i, j, energy[idx]);  
  return 1;  
}




void CompCand(char* whole_seq, int start, int stop, int min_hairpin, int max_hairpin, Cand** ret_cand, int*  ret_ncand, SpanRegion *spanRegion)
{
  int    i, j, length, idx, diag;
  int*   tri_indx;  
  Cand*  cand;
  int    ncand;
  Cand** sort_cand;
  int    cand_idx;
  char*  structure=NULL;  
  char   backtrack_type;  
  int    hairpin_num = 0;  
  char*  sp;  
  char*  seq;
  __int64_t i64_length,i64_cand_alloc_size,i64_max_cand_alloc_size;
  int cand_alloc_size;

  if (verbose) {
    printf("CompCand(seq=%s, start=%d, stop=%d, min_hairpin=%d, max_hairpin=%d, [ret_cand], [ret_ncand]) MIN_ENERGY=%d\n",whole_seq,start,stop,min_hairpin,max_hairpin,MIN_ENERGY);
  }

  length  = stop - start + 1;
  seq=substr(NULL, whole_seq, start, stop);
  PrepareSequence(seq,0);
  structure = (char *) MallocOrDie((unsigned) length+1);
  
  tri_indx = (int*) MallocOrDie( sizeof(int) * (length+1));   
  for(j=1; j <= length; j++) 
    tri_indx[j] = j * (j-1)>>1;
  
  initialize_fold(length);
  comp_energy(seq);  

  i64_length=length;
  i64_cand_alloc_size=(__int64_t)(sizeof(Cand)) * i64_length * (i64_length - 1) / 2;
  i64_max_cand_alloc_size=1<<30;
  if (i64_cand_alloc_size>i64_max_cand_alloc_size) {
	  i64_cand_alloc_size=i64_max_cand_alloc_size;
  }
  cand_alloc_size=(int)i64_cand_alloc_size;
  cand = (Cand*) space(cand_alloc_size);
  memset(cand, 0, cand_alloc_size);
  
  ncand = 0;   
  for( diag = 2; diag <= 2 * length - 1; diag++) {
	  if (ncand>=cand_alloc_size) {
		  fprintf(stderr,"didn't allocate enough space for candidates\n");
		  fflush(stdout);
		  fflush(stderr);
		  exit(1);
	  }
    for(j = diag / 2 + 1; j < MIN(length, diag); j++) {
      i = diag - j;
      if ( j - i + 1 <  MINSPAN)
	continue;
      if ( j - i +1 >= MAXSPAN)
	continue;
      if (spanRegion!=NULL) {
        int doPrint=0;
        if (i>spanRegion->maxStart || j<spanRegion->minEnd) {
          if (doPrint) printf("reject since spanRegion: %d,%d vs span %d,%d\n",i,j,spanRegion->maxStart,spanRegion->minEnd);
          continue;
        }
        else {
          if (doPrint) printf("accept since spanRegion: %d,%d vs span %d,%d\n",i,j,spanRegion->maxStart,spanRegion->minEnd);
        }
      }
      idx = tri_indx[j]  + i;
#if 0
      printf("\t\t%d,%d : %d ? %d\n",i,j,bp_energy[idx],MIN_ENERGY);
#endif
      /* Not stable */
      if ( (bp_energy[idx] > MIN_ENERGY && min_hairpin == 1) || 
	   (max_hairpin > 1 && ml_energy[idx] > MIN_ENERGY) )
	continue;

      if( max_hairpin ==1  || bp_energy[idx] <= ml_energy[idx]) {	/* i, j are basepaired */	
	/* not maximal stack */
	if (!isMaximal(i, j, tri_indx, bp_energy, LOCAL, length))
	    continue;

	/*stack of only 1 bp */
	if ( j > i + 2  && bp_energy[tri_indx[j-1] + i+1] >= INF)     
	    continue;

	cand[ncand].start = i;
	cand[ncand].stop = j;
	
	cand[ncand++].energy =  bp_energy[idx];
	      }            
      else if (max_hairpin > 1 ){	
	/* not maximal stack */
	if (!isMaximal(i, j, tri_indx, ml_energy, LOCAL, length)){	  
	  continue;
	}
	cand[ncand].start = i;
	cand[ncand].stop = j;
	cand[ncand++].energy =  ml_energy[idx];	
      }            
    }    
  }

  if (verbose) {
    printf("\tncand=%d\n",ncand);
  }

  cand_idx = 0;  
  *ret_cand = MallocOrDie(sizeof(Cand) * MAXCAND);  
  sort_cand = SortCand(cand, ncand, CompCandByEnergy);      
  for(i=0; i<ncand && cand_idx < MAXCAND; i++) {         
    int contained = 0;     
    int l1 = sort_cand[i]->stop - sort_cand[i]->start + 1;
        
    for(j=0; j < i; j++) {
      int l2 = sort_cand[j]->stop - sort_cand[j]->start + 1;
      int diff = MAX(abs(sort_cand[i]->stop - sort_cand[j]->stop), abs(sort_cand[i]->start -  sort_cand[j]->start));      
      if ( diff <= 0.15 * MIN(l1, l2)  && diff <= 10){	
	contained = 1;	      
	break;
      }      
    }        
    
    if (contained) continue;  /* Too much overlap with in existing candidate */         
    sort_cand[i]->cand_id = cand_idx;    
    if (sort_cand[i]->energy == bp_energy[ tri_indx[  sort_cand[i]->stop ] + sort_cand[i]->start ])
      backtrack_type = 'C';
    else
      backtrack_type = 'M';  
  
    get_structure( sort_cand[i]->start, sort_cand[i]->stop, seq, structure, backtrack_type);    
    hairpin_num = countHairpin( structure);    
    if (hairpin_num < min_hairpin || hairpin_num > max_hairpin)
      continue;    

    //printf("Cand %d_%d energy %.2f \n", sort_cand[i]->start, sort_cand[i]->stop, sort_cand[i]->energy);
    /* remove the dangling bases for multiloop */
    //printf("Str\t%s\n", structure);    
    sp = structure + sort_cand[i]->stop - sort_cand[i]->start;
    while(*sp == '.') {
      sp--;
      sort_cand[i]->stop--;
    }    
    *(sp+1) = '\0';
    sp = structure;
    while(*sp =='.') {
      sp++;
      sort_cand[i]->start++;
    }
    //printf("\t%s\n", sp);
    
    strcpy( sort_cand[i]->ss, sp);                             
    substr(sort_cand[i]->seq, seq, sort_cand[i]->start - 1, sort_cand[i]->stop - 1);
    memcpy( &(*ret_cand)[cand_idx],sort_cand[i], sizeof(Cand));    
    cand_idx++;    
  }
  for(i=0; i < cand_idx; i++){
    (*ret_cand)[i].start += start;
    (*ret_cand)[i].stop += start;
  }

  if (verbose) {
    printf("\tret_ncand=%d\n",cand_idx);
  }
  
  *ret_ncand = cand_idx;  
  free(structure);   
  free(tri_indx);   
  free_arrays();   
  free(cand);
  free(seq);
  
}

void LoadSpanRegionFile(const char *spanRegionFile,SpanRegion *spanRegionPerSeq,SQINFO *sqinfo,int nseq)
{
  int i;
  FILE *f;
  f=fopen(spanRegionFile,"rt");
  if (f==NULL) {
    Die("Failed to open --span-region-file %s",spanRegionFile);
  }
  for (i=0; i<nseq; i++) {
    int numFields=fscanf(f,"%d,%d\n",&(spanRegionPerSeq[i].maxStart),&(spanRegionPerSeq[i].minEnd));
    if (numFields!=2) {
      Die("format problem in zero-based line %d of --span-region-file %s",i,spanRegionFile);
    }
    if (spanRegionPerSeq[i].maxStart>spanRegionPerSeq[i].minEnd) {
      Die("maxStart>minEnd is impossible in zero-based line %d of --span-region-file %s",i,spanRegionFile);
    }
    if (spanRegionPerSeq[i].maxStart<0) {
      Die("maxStart<0 is impossible in zero-based line %d of --span-region-file %s",i,spanRegionFile);
    }
    if (spanRegionPerSeq[i].minEnd>=sqinfo[i].len) {
      Die("minEnd>=%d (sequence len) is impossible in zero-based line %d of --span-region-file %s",sqinfo[i].len,i,spanRegionFile);
    }
  }
  fclose(f);
}

int main(int argc, char* argv[])
{
  char*   seqfile=NULL;
  int     nseq;
  char**  rseqs;
  SQINFO* sqinfo;
  int     format;

  Cand**  cand;
  int*    ncand;
  char*   cand_file=NULL;
  char*   cluster_file=NULL;
  Range* range=NULL;
  int     i,j,k;

  int     min_hairpin = 1;
  int     max_hairpin = 1;
  char  *optname;                /* name of option found by Getopt()        */
  char  *optarg;                 /* argument found by Getopt()              */
  int    optind;                 /* index in argv[]                         */

  char  *spanRegionFile=NULL;
  SpanRegion *spanRegionPerSeq=NULL;

  FILE *timerFile=NULL;
  Stopwatch_t *timer=StopwatchCreate();
  StopwatchStart(timer);

  /*Parse command line */	

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
#endif

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {

    if      (strcmp(optname, "-c") == 0)    MAXCAND       = atoi(optarg);
    else if (strcmp(optname, "-o") == 0)    cand_file     = optarg;
    else if (strcmp(optname, "-r") == 0)    cluster_file  = optarg;
    else if (strcmp(optname, "-m") == 0)    MINSPAN       = atoi(optarg);
    else if (strcmp(optname, "-M") == 0)    MAXSPAN       = atoi(optarg);
    else if (strcmp(optname, "-s") == 0)    min_hairpin   = atoi(optarg);
    else if (strcmp(optname, "-S") == 0)    max_hairpin   = atoi(optarg);
    else if (strcmp(optname, "-v") == 0)    verbose       = 1;
    else if (strcmp(optname, "--span-region-file") == 0) spanRegionFile=optarg;
    else if (strcmp(optname, "-t")==0) {
      timerFile=fopen(optarg,"at");
      if (timerFile==NULL) {
        Die("cannot open file %s\n", optarg);
      }
    }
    else{
      Die("Invalid Option!\n %s", usage);      
    }    
  }
   
  if (max_hairpin < min_hairpin) {
    Die("Max hairpin should be bigger than Min hairpin \n %s\n", usage);    
  }
  
  if (argc - optind < 1)
    Die("%s\n", usage);  
   
   seqfile = argv[argc - 1];      

   /* RNAfold parameter */
   do_backtrack = 1; 
   noLonelyPairs=1;  
  
   /* Read Sequence file */
   format = SQFILE_FASTA;   
   
   if(! ReadMultipleRseqs(seqfile, format, &rseqs, &sqinfo, &nseq))
     Die("Failed to read squences from file %s", seqfile);

   cand = MallocOrDie(sizeof(Cand*) * nseq);  
   ncand = MallocOrDie(sizeof(int) * nseq);

   if (spanRegionFile!=NULL) {
     spanRegionPerSeq=(SpanRegion *)MallocOrDie(sizeof(SpanRegion)*nseq);
     LoadSpanRegionFile(spanRegionFile,spanRegionPerSeq,sqinfo,nseq);
   }
   
   if (cluster_file){
     if (spanRegionPerSeq!=NULL) {
       Die("not implemented");
     }
     range = ReadCluster(cluster_file, sqinfo, nseq);
     for(i=0; i < nseq; i++){
       if (range[i].start == -1){
	 ncand[i]=0;
	 cand[i] = NULL;
       }
       else{
	 CompCand(rseqs[i], range[i].start-1, range[i].stop-1, min_hairpin, max_hairpin, &cand[i], &ncand[i],NULL);
	 for(k=0; k < ncand[i]; k++){
	   cand[i][k].seq_id = i;
	 }
       }
     }
     Write2DCand(cand_file, nseq, cand, ncand);
   }
   else{   
     for(i=0; i < nseq;  i++) {
       SpanRegion *spanRegion=NULL;
       if (spanRegionPerSeq!=NULL) {
         spanRegion=&(spanRegionPerSeq[i]);
       }
       CompCand(rseqs[i], 0, sqinfo[i].len-1, min_hairpin, max_hairpin, &cand[i], &ncand[i] , spanRegion);
       for(j =0; j < ncand[i]; j++) 
	 cand[i][j].seq_id = i;          
     }     
     Write2DCand(cand_file, nseq, cand, ncand);          
   }


   StopwatchStop(timer);
   if (timerFile!=NULL) {
     fprintf(timerFile,"candf\t%lg\t%lg\n",timer->user+timer->sys,timer->elapsed);
    fclose(timerFile);
   }
   StopwatchFree(timer);

   for(i = 0; i < nseq; i++){
     if (cand[i])
       free(cand[i]);
   }
   free(cand);
   free(ncand);  
   for (i = 0; i < nseq; i++)
     FreeSequence(rseqs[i], &(sqinfo[i]));
   if (spanRegionPerSeq!=NULL) {
     free(spanRegionPerSeq);
   }
   free(sqinfo);
   return 0;   
}
