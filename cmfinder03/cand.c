#include "cand.h"
#include "edit_cost.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int MAXCAND = 100;
int MINSPAN = 20;
int MAXSPAN = MAXLENGTH;


void replace(char* s, char old, char new)
{
  char* p = s;
  while(*p) {
    if (*p == old)
      *p = new;
    p++;  
  }
}


int compDouble(const void* a, const void* b)
{
  double* atemp = (double*)a;
  double* btemp = (double*)b;
  if (*atemp > *btemp)
    return -1;
  if (*atemp < *btemp)
    return 1;
  return 0;
}

int CompCandByScore(const void* a, const void* b)
{
  Cand** atemp = (Cand**) a;
  Cand** btemp = (Cand**) b;
  if ((*btemp)->score < (*atemp)->score )
    return -1;

  if ((*btemp)->score > (*atemp)->score )
    return 1;
  
  return 0;  
}

int CompCandByEnergy(const void* a, const void* b)
{
  Cand** atemp = (Cand**) a;
  Cand** btemp = (Cand**) b;
  int  l1 = (*atemp)->stop - (*atemp)->start + 1;
  int  l2 = (*btemp)->stop - (*btemp)->start + 1;
  double diff = (*atemp)->energy / l1 - (*btemp)->energy / l2 ; 
  return diff < 0 ? -1 : (diff > 0 ? 1 : 0); 
}





Cand** SortCand(Cand* cand, int ncand, int (*criterion)(const void* a, const void* b))
{
  Cand** sort_cand;
  int i;  
  sort_cand = (Cand**) malloc(sizeof(Cand*) * ncand);
  for(i=0; i< ncand; i++)
    sort_cand[i] = &cand[i];
  qsort(sort_cand, ncand, sizeof(Cand*), criterion);
  return sort_cand;  
}


void Reverse(char* seq)
{
  int len;
  char* start;
  char* end;  
  char temp;  
  len = strlen(seq);
  for(start = seq, end = seq + len - 1; start < end; start++, end--) {
    temp = *start;
    *start = *end;
    *end = temp;    
  }  
}


int hairpin(int* pair_table, int i, int j)
{
  int k;
  for(k=i; k<=j; k++) {
    if (pair_table[k] >= 0) 
      return 0;
  }
  return 1;     
}



int multiloop(int* pair_table, int i, int j)
{
  int k;
  int prev_left = i;
  int prev_right = j;
  
  for(k= i+1; k < j; k++) {
    if (pair_table[k] > k) {      
      if (k > prev_right){	
	//rintf("New branch %d %d %d %d\n", k, pair_table[k], pair_table[prev_right], prev_right);	
	return 1;      
      }
      
      prev_left = k;
    }
    if (pair_table[k] >= 0 && pair_table[k] < k) {
      prev_right = k;
    }
  }
  return 0;  
}


int isHairpin(char* ss)
{
  int  len;  
  int* pair_table;
  len = strlen(ss);
  pair_table= GetPairtable(ss);
  return hairpin(pair_table, 0, len - 1); 
}

int isMultiloop(char* ss)
{
  int  len;  
  int* pair_table;
  len = strlen(ss);
  pair_table= GetPairtable(ss);
  return multiloop(pair_table, 0, len - 1); 
}

int countHairpin(char* ss)
{
  int hairpin = 0;
  char* s = ss;
  int begin = 0;
  while(*s) {
    begin = 0;
    while(*s && *s != ')'){      
      s++;      
      if (*s == '(') 
	begin = 1;      
    }
    if (begin) 
      hairpin++;    
    while(*s && *s != '(')
      s++;        
  }
  return hairpin;
}


void Write2DCand(char* cand_file, int nseq, Cand** cand, int* ncand)
{
  int i,j;
  FILE * cand_out;
  
  if (cand_file){    
    if ( (cand_out = fopen(cand_file, "w")) == NULL) {
      Die("Fail to open file %s to write Candidates", cand_file);
    }
  }
  else{
    cand_out = stdout;
  }
    
  for(i=0;  i< nseq; i++) {
    if (cand[i] == NULL) continue;    
    for(j=0; j < ncand[i]; j++) {      
      fprintf(cand_out, ">%d_%d_%d_%d %f\n", cand[i][j].seq_id, cand[i][j].cand_id, cand[i][j].start, cand[i][j].stop, cand[i][j].energy);	
      fprintf(cand_out, "%s\n", cand[i][j].seq);
      fprintf(cand_out, "%s\n", cand[i][j].ss);      
    }
  }    
  if (cand_file) {
    fclose(cand_out);    
  }
  
}

void Write1DCand(char* cand_file, Cand** cand, int ncand)
{
  int i;
  FILE * cand_out;
  
  if ( (cand_out = fopen(cand_file, "w")) == NULL) {
    Die("Fail to open file %s to write Candidates", cand_file);
  }  
  for(i=0;  i< ncand; i++) {
    if (cand[i] == NULL) continue;    
    fprintf(cand_out, ">%d_%d_%d_%d %f\n", cand[i]->seq_id, cand[i]->cand_id, cand[i]->start, cand[i]->stop, cand[i]->energy);	
    fprintf(cand_out, "%s\n", cand[i]->seq);
    fprintf(cand_out, "%s\n", cand[i]->ss);      
  }      
}

/*** Assume that at least one cand for each sequence */
Cand** Read2DCand(char* cand_file,  int nseq, int ** ret_ncand, int * ret_max_cand)
{
  FILE * cand_fin;
  char buffer[MAXLINE]; 
  char seq[MAXLINE];
  char anno[MAXLINE];
  int i,j;  
  int start;
  int stop;  
  int seq_idx, cand_idx;    
  float energy;
  Cand** cand;
  int* ncand;  
  int  num_seq;
  int max_cand = 0;
  
  
  /* Count the number of seqs */
  if ( (cand_fin = fopen(cand_file, "r") ) == NULL) 
    Die("Fail to read file \"%s\"", cand_file);
  num_seq = 0;  
  while(fgets(buffer,MAXLINE, cand_fin) > 0) {    
    if (strstr(buffer, ">")){      
      sscanf(buffer+1, "%d_%d_%d_%d", &seq_idx, &cand_idx, &start, &stop);  
      if (seq_idx + 1 > num_seq) 
	num_seq = seq_idx + 1;      
    }    
  }
  fclose(cand_fin);
  if (num_seq > nseq) 
    Die("Sequence index %d is greater than the number of input sequences %d", num_seq, nseq);  
        
  ncand = (int*) malloc(sizeof(int) * nseq);  
  memset(ncand, 0, sizeof(int) * nseq);
  
  cand = (Cand **) malloc(sizeof(Cand*) * nseq);  
  
  for(i=0; i < nseq; i++) {
    cand[i] =  (Cand *) malloc(sizeof(Cand) * MAXCAND);
  }


  if ( (cand_fin = fopen(cand_file, "r") ) == NULL) 
    Die("Fail to read file \"%s\"", cand_file);

  while(fgets(buffer,MAXLINE, cand_fin) > 0) {    
    if (strstr(buffer, ">")){            
      sscanf(buffer+1, "%d_%d_%d_%d %f", &seq_idx, &cand_idx, &start, &stop, &energy);      
      if (ncand[seq_idx] >= MAXCAND)
	continue;      

      if (fgets(buffer, MAXLINE, cand_fin) <= 0) break;      
      if (strlen(buffer) > MAXSPAN)
	continue;

      sscanf(buffer, "%s", seq);
      cand[seq_idx][ncand[seq_idx]].seq_id  = seq_idx;
      cand[seq_idx][ncand[seq_idx]].cand_id = cand_idx;
      cand[seq_idx][ncand[seq_idx]].start   = start;
      cand[seq_idx][ncand[seq_idx]].stop    = stop;      
      cand[seq_idx][ncand[seq_idx]].len     = abs(stop - start) + 1;      
      cand[seq_idx][ncand[seq_idx]].energy  = energy;
      
      strcpy(cand[seq_idx][ncand[seq_idx]].seq, seq);

      if (fgets(buffer, MAXLINE, cand_fin) <= 0) continue;      
      sscanf(buffer, "%s", anno);
      strcpy(cand[seq_idx][ncand[seq_idx]].ss, anno);

      ncand[seq_idx] ++;
      if (ncand[seq_idx] > max_cand)
	max_cand = ncand[seq_idx];      
    }
    
  }

  /*
  for(i=0; i < nseq; i++){
    for(j=0; j < ncand[i]; j++){
      printf("Check block %d %d (%d %d) %d_%d\n", i,j, cand[i][j].seq_id, cand[i][j].cand_id, cand[i][j].start, cand[i][j].stop);    }
  }
  */

  *ret_ncand = ncand;  
  *ret_max_cand = max_cand;  
  return(cand);
}


/*** Assume that at least one cand for each sequence */
Cand* Read1DCand(char* cand_file, int * ret_ncand)
{
  FILE * cand_fin;
  char buffer[MAXLINE]; 
  char seq[MAXLINE];
  char anno[MAXLINE];

  int i;  
  int start;
  int stop;  
  int seq_idx, cand_idx;    
  float energy;
  Cand* cand;
  int ncand;  
  
  /* Count the number of cand */
  if ( (cand_fin = fopen(cand_file, "r") ) == NULL) 
    Die("Fail to read file \"%s\"", cand_file);
  ncand = 0;
  while(fgets(buffer,MAXLINE, cand_fin) > 0) {    
    if (strstr(buffer, ">")){      
      ncand++;
    }
  }
  fclose(cand_fin);
        
  cand = (Cand *) malloc(sizeof(Cand) * ncand);

  
  if ( (cand_fin = fopen(cand_file, "r") ) == NULL) 
    Die("Fail to read file \"%s\"", cand_file);

  i = 0;
  while(fgets(buffer,MAXLINE, cand_fin) > 0) {    
    if (strstr(buffer, ">")){            
      sscanf(buffer+1, "%d_%d_%d_%d %f", &seq_idx, &cand_idx, &start, &stop, &energy);      
      if (fgets(buffer, MAXLINE, cand_fin) <= 0) break;      
      if (strlen(seq) > MAXSPAN)
	continue;

      sscanf(buffer, "%s", seq);
      cand[i].seq_id  = seq_idx;
      cand[i].cand_id = cand_idx;
      cand[i].start   = start;
      cand[i].stop    = stop;      
      cand[i].energy  = energy;
      cand[i].weight  = 1;
      cand[i].len     = abs(stop - start) + 1;      
      strcpy(cand[i].seq, seq);
      if (fgets(buffer, MAXLINE, cand_fin) <= 0) continue;      
      sscanf(buffer, "%s", anno);
      strcpy(cand[i].ss, anno);

      i++;
    }    
  }

  *ret_ncand = ncand;  
  return(cand);
}

/* transform the format to Vienna package */
char* ExpandFull(char* anno, char* seq)
{   
  int i,j, pp;    
  int pos[1000];
  /* stack, keep the positions of unresolved '(' */  
  int len;      
  char* struc;  
  char c1, c2;
    
  len= strlen(anno);  
  struc = (char *)malloc( len * 5);    
  
  pp = 0;
  struc[0] = '(';
  j = 1;  
  for(i = 0; i < len; i++) {    
    switch( anno[i] ) {      
    case '(':
      pos[pp++] = i;            
      struc[j++] = anno[i];            
      break;      
    case ')':
      c1 = seq[ pos[--pp] ];      
      c2 = seq[i];      
      j += PairCodeVienna(c1, c2, struc+j);      
      struc[j++] = anno[i];
      break;      
    default:
      struc[j++] = '(';      
      struc[j++] = seq[i];
      struc[j++] = ')';      
    }        
  }  
  struc[j] = '\0';  
  /* add a virtual root */
  strcat(struc, "Null)");  
  return(struc);
}


SQINFO* Cand2Sqinfo(Cand** cand, int ncand, SQINFO* seq_sqinfo)		   
{
  int i;  
  SQINFO* sqinfo;
  
  int flags = SQINFO_NAME | SQINFO_START | SQINFO_STOP | SQINFO_SS | SQINFO_TYPE |SQINFO_DESC| SQINFO_LEN;    
  if ((sqinfo =  (SQINFO*)malloc(sizeof(SQINFO) * ncand)) == NULL) 
    Die("Memory failure at Cand2Sqinfo for sqinfo");
  memset(sqinfo, 0, sizeof(SQINFO)* ncand);
  
  for(i=0; i < ncand; i++){    
    if (cand[i] == NULL)
      continue;
    sqinfo[i].len = strlen(cand[i]->seq);    
    strcpy(sqinfo[i].name, seq_sqinfo[cand[i]->seq_id].name);
    sqinfo[i].start = cand[i]->start;
    sqinfo[i].stop = cand[i]->stop;
    sqinfo[i].ss   = cand[i]->ss;    
    replace(sqinfo[i].ss, '(', '<');
    replace(sqinfo[i].ss, ')', '>');    
    sprintf(sqinfo[i].desc, "%3d..%3d\t%f", cand[i]->start, cand[i]->stop, cand[i]->score);      
    sqinfo[i].type =  kRNA;
    sqinfo[i].flags = flags;    
  }
  return sqinfo;
}


int Overlap(Cand* c1, Cand* c2, int* olap_min, int* olap_max, int* min, int* max)
{
  if (c1->start <= c2->start){
    *min = c1->start;
    if (c2->start <= c1 ->stop){
      *olap_min = c2->start;      
      if (c2->stop <= c1->stop){
	*olap_max = c2->stop;
	*max = c1->stop;
      }
      else{
	*olap_max = c1->stop;
	*max = c2->stop;
      }	
    }
    else{
      return 0;
    }
  }
  else{
    *min = c2->start;
    if (c1->start <= c2 ->stop){
      *olap_min = c1->start;
      if (c1->stop <= c2-> stop){
	*olap_max = c1->stop;
	*max = c2->stop;
      }
      else{
	*olap_max = c2->stop;
	*max = c1->stop;
      }
    }
    else{
      return 0;
    }
  }
  return *olap_max - *olap_min + 1;
}
