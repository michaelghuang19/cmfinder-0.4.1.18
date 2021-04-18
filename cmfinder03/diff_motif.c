#include "squid.h"
#include "msa.h"
#include "sqfuncs.h"
#include <string.h>
#include <ctype.h>


#define MIN_ALIGN_LENGTH 20
#define MAX_GAP   60
#define MIN_OVERLAP 2
#define max(x,y) (x<y? y: x)
#define min(x,y) (x<y? x: y)
#define pair_left(c) (c == '(' || c == '<' || c == '{')
#define pair_right(c) (c == ')' || c == '>' || c == '}')
#define IsBasePair(l, r) ((l=='A'||l=='G')&&r=='U'||l=='U'&&(r=='A'||r=='G')||l=='G'&&r=='C'||l=='C'&&r=='G')  
int   do_weight=0;

static char usage[]  = "\
Usage: diff_motif [-options] <motif1> <motif2> \n\
where options are:\n\
     -o      : Show overlapped sequences \n\
     -m     : Show sequences in motif1 that are not in motif2\n\
     -f     : Show sequences in motif2 that are not in motif1\n\
     -h      : print short help and version info\n\
";


static struct opt_s OPTIONS[] = {
  { "-o", TRUE, sqdARG_NONE},
  { "-m", TRUE, sqdARG_NONE},
  { "-f", TRUE, sqdARG_NONE},
  { "-v", TRUE, sqdARG_NONE},
};

struct motif_coor_s {		
  char*          acc;      
  int            start;	   
  int            end;	   
  int            strand;
};
typedef struct motif_coor_s MotifCoor;


struct structure_diff_s{
  int  motif1;
  int  motif2;
  int  overlap_len;
  int  TP;
  int  FP;
  int  FN;
  int  FM;
  int  W_FP;
  int  W_FN;
};
typedef struct structure_diff_s StructureDiff;



int no_gap_length(char* s)
{
  int i=0;
  while(*s){
    if (!isgap(*s)) i++;
    s++;
  }
  return i;
}

void remove_gap(char* seq, char* ss)
{
  int i,j,len;
  i=j=0;
  len = strlen(seq);
  while(j < len){
    if (!isgap(seq[j])){
      seq[i] = seq[j];
      ss[i] = ss[j];
      i++;
    }
    j++;
  }
  seq[i]='\0';
  ss[i]='\0';
}

MotifCoor* parse_motif_coor(MSA* msa)
{
  int i;  
  MotifCoor* motif_coor = (MotifCoor*) MallocOrDie(sizeof(MotifCoor) * msa->nseq);  
  for(i=0; i < msa->nseq; i++){
    char* temp = (char*)MallocOrDie(sizeof(char) * (strlen(msa->sqname[i])+1));    
    char* acc ;
    int   seq_start=0;
    int   seq_end=0;
    int   motif_start=0;
    int   motif_end=0;
    int   strand;
    char* acc_index;
    strcpy(temp, msa->sqname[i]);
    acc_index= strchr(temp, '/');        
    if (acc_index){
      *acc_index='\0';
      acc_index++;
      sscanf(acc_index, "%d-%d", &seq_start, &seq_end);
      strand =  (seq_start <= seq_end) ? 1: -1;
    }
    else{
      seq_start=1;
      strand = 1;
    }
    acc = temp;    
    if (msa->sqdesc && sscanf(msa->sqdesc[i], "%d..%d", &motif_start, &motif_end)>0){      
      if (abs(seq_start - seq_end) - abs(motif_start - motif_end) >= 0){
	if (strand == 1){
	  motif_start += seq_start -1;
	  motif_end += seq_start -1;
	}
	else{
	  motif_start = seq_start - motif_start + 1;
	  motif_end = seq_start - motif_end + 1;      
	}
      }
    }
    else{
      motif_start = seq_start;
      motif_end = seq_end;
    }
    motif_coor[i].acc = acc;
    motif_coor[i].start = motif_start;
    motif_coor[i].end = motif_end;    
    motif_coor[i].strand = strand;
    //printf("%s %d %d\n", acc, motif_start, motif_end);
  }  
  return motif_coor;
}


/* Remove gaps in a sequence */
int**  alignment_map(MSA* msa, MotifCoor* motif_coor)
{
  int   i, j, k;
  int** align_map = (int**)MallocOrDie(sizeof(int*) * msa->nseq);
  for(k=0; k < msa->nseq; k++){        
    char* seq=msa->aseq[k];
    align_map[k] = (int*) MallocOrDie(sizeof(int) * msa->alen);    
    j = motif_coor[k].start;
    for(i=0; i < msa->alen; i++){
      if (!isgap(seq[i])){
	align_map[k][i] = j;
	if (motif_coor[k].strand == 1){
	  j++;
	}
	else{
	  j--;
	}		
      }	    
      else{
	align_map[k][i]=-1;
      }      
    }      
  }
  return align_map;
}


int* pair_table(char* ss, char* seq, int len)
{
  int  i,j;
  int  sp=0;
  int* stack = (int*)MallocOrDie(sizeof(int) * len);
  int* pt = (int*)MallocOrDie(sizeof(int) * len);
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
        Die("unbalanced base pair at pos %d", i);
      }
      j = stack[sp];
      if (seq == NULL || IsBasePair(seq[i], seq[j])){
	pt[j] = i;
	pt[i] = j;
      }
    }
  }
  return pt;
}

StructureDiff * no_overlap(int alen1, int* pt1, int alen2, int* pt2)
{
  int i;
  StructureDiff * result = (StructureDiff*)MallocOrDie(sizeof(StructureDiff));  
  result->overlap_len = 0;
  result->TP = 0;
  result->FP = 0;
  result->FN = 0;
  result->FM = 0;
  result->W_FN = result->W_FP = 0;
  if (pt1){
    for(i=0; i < alen1; i++)
      if (pt1[i] >= 0) result->W_FP++;
  }
  if (pt2){
    for(i=0; i < alen2; i++)
      if (pt2[i] >= 0) result->W_FN++;
  }
  return result;    
}


StructureDiff * overlap(MotifCoor * coor1, MotifCoor * coor2, int alen1, int alen2, int * align_map1, int * align_map2, int * pt1, int * pt2){
  int i,j,i1,j1,k;
  int TP,FP, FN, TN, FM,W_FN, W_FP;
  TP= FP= FN= TN= FM= W_FN = W_FP= 0;  
  int olap_start, olap_end, all_start, all_end, olap;
  int strand;


  
  if (strcmp(coor1->acc, coor2->acc) != 0 && coor1->strand != coor2->strand) 
    return NULL;

  strand = coor1->strand;  
  if (coor1->strand == 1){
    olap_start = max(coor1->start, coor2->start);
    olap_end =  min(coor1->end, coor2->end);    
    all_start =  min(coor1->start, coor2->start);
    all_end =  max(coor1->end, coor2->end);    
  }
  else{
    olap_start = min(coor2->start, coor1->start);
    olap_end =  max(coor2->end, coor1->end);
    all_start =  max(coor1->start, coor2->start);
    all_end =  min(coor1->end, coor2->end);    
  }
  if ( (olap_start - olap_end)*strand >= 0)
    return NULL;
  
  olap = abs(olap_start - olap_end) + 1;
  
  //printf("%d %d\n", olap_start, olap_end);
  i=j=0;

  for(; align_map1[i] != olap_start; i++)
    if(align_map1[i] > 0 && pt1[i] >=0)   W_FP++;
  
  for(; align_map2[j] != olap_start; j++)
    if(align_map2[j] > 0 && pt2[j] >=0)   W_FN++;
  
  for(k=olap_start; k!=olap_end; k += strand){
    for(; align_map1[i] != k; i++);
    for(; align_map2[j] != k; j++);
    if (align_map1[i]!= k || align_map2[j]!=k) break;    
    
    if (pt1[i] >= 0 && pt2[j] >= 0){
      i1 = pt1[i];
      j1 = pt2[j];
      if (align_map1[i1] == align_map2[j1]){
	TP ++;
      }
      else{
	FM ++;
      }
    }  
    else if (pt1[i] >= 0){
      FP ++;
    }
    else if (pt2[j] >= 0) {
      FN++;
    }
  }

  for(; i < alen1; i++)
    if(pt1[i] >=0)     W_FP++;
  
  for(; j < alen2; j++)
    if(pt2[j] >=0)     W_FN++;

  StructureDiff * result = (StructureDiff*)MallocOrDie(sizeof(StructureDiff));    
  result->overlap_len = olap;
  result->TP = TP;
  result->FP = FP;
  result->FN = FN;
  result->FM = FM;
  result->W_FN = W_FN;
  result->W_FP = W_FP;
  return result;
}

void preproc(MSA* msa)
{  
  int i,j;
  int* pt;
  if (msa->ss == NULL && msa->ss_cons){
    msa->ss = (char**)MallocOrDie(sizeof(char*) * msa->nseq);
    for(i=0; i < msa->nseq; i++){ 
      msa->ss[i] = (char*)MallocOrDie(sizeof(char) * (strlen(msa->ss_cons) +1));      
      strcpy(msa->ss[i],msa->ss_cons);      
    }
  }
  for(i=0; i < msa->nseq; i++){
    pt = pair_table(msa->ss[i], msa->aseq[i], msa->alen);
    for(j=0; j < msa->alen; j++){
      if (pt[j] > j){
	char l = toupper(msa->aseq[i][j]);
	char r = toupper(msa->aseq[i][pt[j]]);
	if (!IsBasePair(l,r)){
	  msa->ss[i][j]=msa->ss[i][pt[j]]='.';	  
	}	  
      }
    } 
    free(pt);
  }
  if (!do_weight){
    for(i=0; i < msa->nseq; i++){
      msa->wgt[i] = 1;
    }
  }
}






#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

int
main(int argc, char **argv)
{
  int        format   =MSAFILE_STOCKHOLM;  
  char      *sto_file1=NULL;      /* training sequence file                    */
  char      *sto_file2=NULL;      /* training sequence file                    */
  MSA       *msa1,*msa2;
  MotifCoor *motif_coor1  = NULL;
  MotifCoor *motif_coor2  = NULL;
  int      **align_map1, **align_map2;
  int      **pt1 = NULL;
  int      **pt2 = NULL;
  MSAFILE   *afp      = NULL;     /* file handle of stockholm alignment        */    
  double     tot_weight1,tot_weight2,len1,len2;
  char       *seq_buf,*ss_buf;
  StructureDiff ** best1=NULL;
  StructureDiff ** best2=NULL;
  
  int  i,j;
  double TP1,TP2, FP, FN, FM1,FM2, W_FP, W_FN, sens,spec;
  
  char  *optname;                /* name of option found by Getopt()        */
  char  *optarg;                 /* argument found by Getopt()              */
  int    optind;                 /* index in argv[]                         */

  int   show_overlap=0;
  int   show_nolap1=0;
  int   show_nolap2=0;


  int  *olap1;
  int  *olap2;
  char temp1[100], temp2[100];

  if (argc < 2){
      puts(usage);
      exit(EXIT_SUCCESS);
  }
  
  /*Parse command line */
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg)) {

    if      (strcmp(optname, "-o") == 0)  show_overlap=1;      
    else if (strcmp(optname, "-m") == 0)  show_nolap1=1;      
    else if (strcmp(optname, "-f") == 0)  show_nolap2=1;      
    else if (strcmp(optname, "-w") == 0)  do_weight = 1;
    else if (strcmp(optname, "-h") == 0) {
      puts(usage);
      exit(EXIT_SUCCESS);
    }
  }
  
  sto_file1= argv[optind++];
  sto_file2= argv[optind++];

  if ((afp = MSAFileOpen(sto_file1, format, NULL)) == NULL)
    Die("Alignment file %s could not be opened for reading", sto_file1);
  if ((msa1 = MSAFileRead(afp)) == NULL)
    Die("Invalid stockholm format");
  MSAFileClose(afp);
  preproc(msa1);
  
  
  if ((afp = MSAFileOpen(sto_file2, format, NULL)) == NULL)
    Die("Alignment file %s could not be opened for reading", sto_file2); 
  if ((msa2 = MSAFileRead(afp)) == NULL)
    Die("Invalid stockholm format");
  MSAFileClose(afp); 
  preproc(msa2);

  
  motif_coor1= parse_motif_coor(msa1);
  align_map1=alignment_map(msa1,motif_coor1);
  motif_coor2= parse_motif_coor(msa2);
  align_map2=alignment_map(msa2,motif_coor2);
  
  pt1 = (int**)MallocOrDie(sizeof(int*) * msa1->nseq);
  pt2 = (int**)MallocOrDie(sizeof(int*) * msa2->nseq);
  olap1 = (int*)MallocOrDie(sizeof(int) * msa1->nseq);
  olap2 = (int*)MallocOrDie(sizeof(int) * msa2->nseq);
  best1 = (StructureDiff**)MallocOrDie(sizeof(StructureDiff*) * msa1->nseq);
  best2 = (StructureDiff**)MallocOrDie(sizeof(StructureDiff*) * msa2->nseq);  
  memset(olap1, 0, sizeof(int) * msa1->nseq);
  memset(olap2, 0, sizeof(int) * msa2->nseq);
    
  tot_weight1=tot_weight2=len1=len2=0;
  for(i=0; i < msa1->nseq; i++){
    pt1[i]= pair_table(msa1->ss[i], msa1->aseq[i], msa1->alen);
    tot_weight1 += msa1->wgt[i];    
    len1 += (abs(motif_coor1[i].start - motif_coor1[i].end) + 1) * msa1->wgt[i];
  }
  for(i=0; i < msa2->nseq; i++){
    pt2[i]= pair_table(msa2->ss[i], msa2->aseq[i], msa2->alen);
    tot_weight2 += msa2->wgt[i];
    len2 += (abs(motif_coor2[i].start - motif_coor2[i].end) + 1) * msa2->wgt[i];
  }
  len1 /= tot_weight1;
  len2 /= tot_weight2;

  double olap_seq1, olap_len1, olap_seq2, olap_len2;
  olap_seq1= olap_len1=  olap_seq2= olap_len2= 0;
  TP1=TP2=FP=FN=FM1=FM2=W_FP=W_FN=sens=spec=0;
  
  int l=msa1->alen > msa2->alen ? msa1->alen +1 : msa2->alen+1;
  seq_buf = MallocOrDie(l *sizeof(char));
  ss_buf = MallocOrDie(l *sizeof(char));
  

  /* Here we are dealing with the cases in which the best match relatioinship is not mutual. */
  for(i=0; i < msa1->nseq; i++){    
    best1[i]= no_overlap(msa1->alen, pt1[i], msa2->alen, NULL);    
    best1[i]->motif1 = i;
    best1[i]->motif2 = -1;
    for(j=0; j < msa2->nseq; j++){
      if (strcmp(motif_coor1[i].acc, motif_coor2[j].acc)==0 ){
	StructureDiff * temp = overlap(&motif_coor1[i], &motif_coor2[j], msa1->alen, msa2->alen, align_map1[i], align_map2[j], pt1[i], pt2[j]);
	if (temp && temp->overlap_len > best1[i]->overlap_len){
	  free(best1[i]);	  
	  best1[i] = temp;
	  best1[i]->motif2 = j;
	}	
      }      
    }
  }

  for(j=0; j < msa2->nseq; j++){    
    best2[j]= no_overlap(msa1->alen, NULL, msa2->alen, pt2[j]);    
    best2[j]->motif1 = -1;
    best2[j]->motif1 = j;
    for(i=0; i < msa1->nseq; i++){
      if (strcmp(motif_coor1[i].acc, motif_coor2[j].acc)==0 ){
	StructureDiff * temp = overlap(&motif_coor1[i], &motif_coor2[j], msa1->alen, msa2->alen, align_map1[i], align_map2[j], pt1[i], pt2[j]);
	if (temp && temp->overlap_len > best2[j]->overlap_len){
	  free(best2[j]);
	  best2[j] = temp;
	  best2[j]->motif1 = i;
	}	
      }      
    }
  }
  
  for(i=0; i < msa1->nseq; i++){    
    StructureDiff* diff = best1[i];        
    sprintf(temp1, "%s/%d-%d", motif_coor1[i].acc, motif_coor1[i].start, motif_coor1[i].end);	    
    if (diff->overlap_len > 0){
      j = diff->motif2;
      sprintf(temp2, "%s/%d-%d", motif_coor2[j].acc, motif_coor2[j].start, motif_coor2[j].end);	    
      if (show_overlap){
	strcpy(seq_buf, msa1->aseq[i]);
	strcpy(ss_buf, msa1->ss[i]);
	remove_gap(seq_buf, ss_buf);

	printf("%-30s %s\n", temp1, seq_buf, "\n");
	printf("%-30s %s\n", "", ss_buf);
	
	strcpy(seq_buf, msa2->aseq[j]);
	strcpy(ss_buf, msa2->ss[j]);
	remove_gap(seq_buf, ss_buf);
	printf("%-30s %s\n", temp2, seq_buf, "\n");
	printf("%-30s %s\n", "",ss_buf);
	
	printf("%s len=%d TP=%d FP=%d FN=%d FM=%d W_FP=%d W_FN=%d\n", motif_coor1[i].acc, diff->overlap_len, diff->TP, 
	       diff->FP, diff->FN, diff->FM,diff->W_FP, diff->W_FN);
	printf("\n");
      }
      olap_seq1 ++;
      olap_len1 += diff->overlap_len;
      TP1 +=  diff->TP;
      FP +=  diff->FP;
      FM1 +=  diff->FM;	      
    }
    if (diff->overlap_len == 0 && show_nolap1){      
      strcpy(seq_buf, msa1->aseq[i]);
      strcpy(ss_buf, msa1->ss[i]);
      remove_gap(seq_buf, ss_buf);
      printf("%-30s %s\n", temp1, seq_buf, "\n");
      printf("%-30s %s\n", "", ss_buf);
    }
    W_FP += diff->W_FP;
  } 

  for(j=0; j < msa2->nseq; j++){    
    StructureDiff* diff = best2[j];    
    if (diff->overlap_len > 0){
      olap_seq2 ++;
      olap_len2 += diff->overlap_len;
      TP2 +=  diff->TP;
      FN +=  diff->FN;
      FM2 +=  diff->FM;	      
    }   
    if (diff->overlap_len == 0 && show_nolap2){      
      sprintf(temp2, "%s/%d-%d", motif_coor2[j].acc, motif_coor2[j].start, motif_coor2[j].end);	    
      strcpy(seq_buf, msa2->aseq[j]);
      strcpy(ss_buf, msa2->ss[j]);
      remove_gap(seq_buf, ss_buf);      
      printf("%-30s %s\n", temp2, seq_buf, "\n");
      printf("%-30s %s\n", "", ss_buf);
    }
    W_FN += diff->W_FN;
  } 
  
  free(seq_buf);
  free(ss_buf);
  free(olap1);
  free(olap2);  
  
  if (olap_seq1 > 0){
    olap_len1 /= olap_seq1;
    if (TP1 > 0) spec= ((double)TP1)/(TP1 + FP + FM1 + W_FP);    
    TP1 /= olap_seq1;
    FP /= olap_seq1;
    FM1 /= olap_seq1;
  }  
  W_FP /= tot_weight1;
  if (olap_seq2 > 0){
    olap_len2 /= olap_seq2;
    if (TP2 > 0) sens= ((double)TP2)/(TP2 + FN + FM2 + W_FN);
    TP2 /= olap_seq2;
    FN  /= olap_seq2;
    FM2 /= olap_seq2;
  }
  W_FN /= tot_weight2;
  printf("%s\t %s\t olap_seq=%.1f\t nolap_seq1=%.1f\t nolap_seq2=%.1f\t olap_len=%.1f\t nolap_len1=%.1f\t nolap_len2=%.1f\t TP=%.1f\t FP=%.1f\t FN=%.1f\t FM=%.1f\t W_FP=%.1f\t W_FN=%.1f\t Sens=%.2f\t Spec=%.2f\n",
	 sto_file1,sto_file2,olap_seq1,tot_weight1-olap_seq1, tot_weight2-olap_seq2, 
	 olap_len1, len1 - olap_len1, len2 - olap_len2, 
	 TP1, FP,FN,FM1,W_FP, W_FN, sens,spec);
  
  
  for(i=0; i < msa1->nseq; i++){
    free(pt1[i]);
    free(motif_coor1[i].acc);
  }
  free(motif_coor1);
  free(pt1);
  for(i=0; i < msa2->nseq; i++){
    free(pt2[i]);
    free(motif_coor2[i].acc);
  }
  free(pt2);
  free(motif_coor2);


  MSAFree(msa1);
  MSAFree(msa2);
}
