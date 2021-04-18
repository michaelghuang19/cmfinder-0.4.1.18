#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "cand.h"
#include "squid.h"
#include "treedist.h"
#include "edit_cost.h"
#include "dist_vars.h"

static struct opt_s OPTIONS[] = {
  { "-t", TRUE, sqdARG_STRING}
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))
char usage[] = "\
Usage:\n\
canda [-t <timer-append-file>]  <cand_file> <seq_file> <out_file>\n";

int BaseCode(char c)
{
  switch(c) {
  case 'A':
  case 'a':
    return 1;
  case 'C':
  case 'c':
    return 2;
  case 'G':
  case 'g':
    return 3;
  case 'T':
  case 't':
  case 'U':
  case 'u':
    return 4;        
  default:
    return 0;    
  }  
}

MSA* MultipleAlignment(SQINFO* sqinfo, Cand* cand, int ncand, int** map)
{
  int    i, j, k, len, align_len, pos;
  int    *align_pos;    
  int    gap;
  int    temp;  
  MSA    *msa;
  
  //Count the length of the alignment.
  len =  cand[0].len; 
  align_pos=(int *)MallocOrDie(sizeof(int) * (len + 2));
  align_pos[0]=0;

  for(i = 1; i <= len; i++) {
    gap = 1;    
    for(j=1; j < ncand; j++) {
      if (map[j][i] >=0) {	
	temp =0;
	if (i==1 && map[j][i] >= 0) temp = map[j][i];
	else{
	  k= i-1;
	  while(map[j][k] < 0 && k >=0) k--;
	  if (k >= 0)  temp = map[j][i] - map[j][k];
	}		
	if (temp > gap) 
	  gap = temp;	      
      }    
    }
    align_pos[i] = gap + align_pos[i-1];   
    //printf("%d: pos %d\n", i, align_pos[i]);    
  }
  gap = 0;  
  for(j=1; j < ncand; j++) {
    k= len;
    while(map[j][k] < 0 && k >=0) k--;
    temp= cand[j].len - map[j][k];
    if (temp > gap) gap = temp;    
  }
  align_len = align_pos[len+1]=align_pos[len] + gap;    

  //Create the multiple alignment
  msa = MSAAlloc(ncand, align_len);
  msa->nseq = ncand;
  msa->alen = align_len;
  msa->ss   = MallocOrDie(sizeof(char *) * ncand);  

  for (i = 0; i < ncand; i++)
    {
      /* Initialize the aseq with all pads '.' (in insert cols) 
       * and deletes '-' (in match cols).
       */
      msa->ss[i] = MallocOrDie(sizeof(char) * (align_len + 1));      
      for (j = 0; j < align_len; j++){	
	msa->aseq[i][j] = '.';
	msa->ss[i][j] = '.';
      }                 
      msa->aseq[i][align_len] = '\0';
      msa->ss[i][align_len] = '\0';
    }
  for(i=1; i <= len; i++){
    msa->aseq[0][align_pos[i]-1] = cand[0].seq[i-1];        
    msa->ss[0][align_pos[i] -1 ] = cand[0].ss[i-1];                
  }
  for(j=1; j < ncand; j++) {  
    for(i=1; i <= len; i++) {
	  int pos;
	  int orig_pos;
	  int gap_len;
      if (map[j][i] < 0) continue;
      pos = align_pos[i] - 1;
      orig_pos = map[j][i]-1;
      gap_len = 1;
      k= i-1;
      while(map[j][k] < 0 && k >=0) k--;
      if (k >= 0)
	gap_len = map[j][i] - map[j][k];      
      while(gap_len >0){
	msa->aseq[j][pos] = cand[j].seq[orig_pos];	
	msa->ss[j][pos--] = cand[j].ss[orig_pos--];
	gap_len--;
      }      
    }
    k= len;
    while(map[j][k] < 0 && k >=0) k--;
    pos = align_pos[len];
    for(temp= map[j][k]+1; temp <= cand[j].len; temp++){
      msa->aseq[j][pos] = cand[j].seq[temp-1];	
      msa->ss[j][pos++] = cand[j].ss[temp-1];
    }
  }
    
  /* "author" info */
  msa->au   = MallocOrDie(sizeof(char) * 10);
  sprintf(msa->au, "canda ");
  
  for (i = 0; i < ncand; i++)
    {
      int seq_id = cand[i].seq_id;      
      msa->sqname[i] = sre_strdup(sqinfo[seq_id].name, -1);
      msa->sqlen[i]  = sqinfo[seq_id].len;
      if (sqinfo[seq_id].flags & SQINFO_ACC)
        MSASetSeqAccession(msa, i, sqinfo[seq_id].acc);
      msa->wgt[i] = 1.0;
    }
  msa->ss_cons = (char*)MallocOrDie(align_len + 1);
  strcpy(msa->ss_cons, msa->ss[0]);
  return msa;
}


int* TransformAlignment(char* align[2], int* len)
{
  int max_len= strlen(align[0]);  
  int *map = (int*) MallocOrDie(sizeof(int) * max_len);
  int *stack1 = (int*) MallocOrDie(sizeof(int) * max_len);
  int *stack2 = (int*) MallocOrDie(sizeof(int) * max_len);
  int sp1=0;
  int sp2=0;
  int pos1 =1;
  int pos2 =1;    
  int i;
  for(i = 1; i < max_len; i++) {
    map[i] = -1;
  }
  map[0]=0;
  //printf("%s\n",align[0]);
  //printf("%s\n",align[1]);
  for(i = 1; i < max_len; i++) {  //Ignore the beginning '(' and ending ')'
    //printf("i %d, %c %c (%d %d): ", i, align[0][i], align[1][i],pos1,pos2);
    if (align[0][i] == '('){
      stack1[sp1++] = pos1;
      if (!BaseCode(align[0][i+1])) pos1++;
    }
    if (align[1][i] == '('){
      stack2[sp2++] = pos2;
      if (!BaseCode(align[1][i+1])) pos2++;
    }
    if (align[0][i] == ')') sp1--;      
    if (align[1][i] == ')') sp2--;      

    //pair position
    if (BaseCode(align[0][i]) && BaseCode(align[0][i+1]) || 
	BaseCode(align[1][i]) && BaseCode(align[1][i+1])){
	  char right1;
	  char right2;
      char left1 = align[0][i];
      char left2 = align[1][i];
      if (!isgap(left1)){
	if (!isgap(left2)){
	  map[stack1[sp1-1]] = stack2[sp2-1];      
	  //printf("(pos1 %d pos2 %d)", stack1[sp1-1],stack2[sp2-1]);
	}
      }
      i++;
      right1 = align[0][i];
      right2 = align[1][i];
      if (!isgap(right1) && !isgap(right2)){
	map[pos1]=pos2;
	//printf("(pos1 %d pos2 %d)", pos1, pos2);
      }
      if (!isgap(right1)){
	pos1++;
      }
      if (!isgap(right2)){
	pos2++;
      }
    }    
    else if (BaseCode(align[0][i]) || BaseCode(align[1][i])){
      if (!isgap(align[1][i]) && !isgap(align[0][i])){
	map[pos1] = pos2;
	//printf("(pos1 %d pos2 %d)", pos1, pos2);
      }
      if (!isgap(align[1][i])){
	pos2++;
      }
      if (!isgap(align[0][i]))
	pos1++;
    }      
    //printf("\n");
  }
  *len = pos1 -1;
  free(stack1);
  free(stack2);
  return map;  
}


int main(int argc, char* argv[])
{
  
  char*   seqfile; 
  int     nseq;  
  char**  rseqs;
  SQINFO* sqinfo;
  int     format;  
  char*   align_file;
  FILE*   align_fout;  

  char*  candfile;
  int    ncand;
  Cand*  cand;  
  int    i,j;
  
  Tree**   trees;
  double   dist;
  char*    struc; 
  int**    map;
    
  MSA     *msa;

  char  *optname;                /* name of option found by Getopt()        */
  char  *optarg;                 /* argument found by Getopt()              */
  int    optind;                 /* index in argv[]                         */

  FILE *timerFile=NULL;
  Stopwatch_t *timer=StopwatchCreate();
  StopwatchStart(timer);

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {

    if (strcmp(optname, "-t")==0) {
      timerFile=fopen(optarg,"at");
      if (timerFile==NULL) {
        Die("cannot open file %s\n", optarg);
      }
    }
    else{
       Die("Invalid Option %s\n %s", optname, usage);      
     } 
  }

  if (argc-optind != 3) {
    Die("%s\n",usage);
  }
    
  candfile = argv[argc-3];
  seqfile = argv[argc-2];
  align_file = argv[argc-1];  
 
  /*Search Motifs in each sequence. 
   * Produce allCands */
  nseq = 0;  
  cand = Read1DCand(candfile, &ncand);
    
  trees = (Tree**) MallocOrDie(sizeof(Tree*) * ncand);
  for(i =0 ; i < ncand; i++) {
    struc = ExpandFull(cand[i].ss, cand[i].seq);      
    trees[i] = make_tree(struc);      
    free(struc);
  }
  map = (int**)MallocOrDie(sizeof(int*) * ncand);
  
  i=0;  
  for(j=i+1; j < ncand; j++) {
    int temp;
    dist = tree_edit_distance(trees[i], trees[j]) ;    
    map[j]=TransformAlignment(aligned_line,&temp);
  }  

  format = SQFILE_FASTA;
  
  /* read the training seqs from file */
  if (! ReadMultipleRseqs(seqfile, format, &rseqs, &sqinfo, &nseq))
    Die("Failed to read sequences from file %s", seqfile);

  msa = MultipleAlignment(sqinfo, cand, ncand, map);
  
  if (align_file != NULL && (align_fout = fopen(align_file, "w")) != NULL) 
    {
      WriteStockholm(align_fout, msa);
      printf("Alignment saved in file %s\n", align_file);
      fclose(align_fout);
    }
  else
    WriteStockholm(stdout, msa);


  StopwatchStop(timer);
  if (timerFile!=NULL) {
    fprintf(timerFile,"canda\t%lg\t%lg\n",timer->user+timer->sys,timer->elapsed);
    fclose(timerFile);
  }
  StopwatchFree(timer);
  
  for(j = 0; j < ncand; j++) {
    free_tree(trees[j]);
    if (j > 0) free(map[j]);
  }
  free(trees);
  free(map);
  free(cand);  
  
  for (i = 0; i < nseq; i++) 
    {
      FreeSequence(rseqs[i], &(sqinfo[i]));
    }
  free(sqinfo);  
  return 0;
}
