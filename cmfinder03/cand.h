#ifndef CAND_H
#define CAND_H

#define MAXLENGTH 500

#include "squid.h"
#include "global.h"

extern int MAXCAND;
extern int MAXSPAN;
extern int MINSPAN;


typedef struct
{
  int seq_id;
  int cand_id;  
  int start;  
  int stop;
  int len;  
  double score;  
  double weight;    
  double energy;  
  char seq[MAXLENGTH];  
  char ss[MAXLENGTH];  
}Cand;


/* I/O */
void Write2DCand(char* cand_file, int nseq, Cand** cand, int* ncand);
void Write1DCand(char* cand_file, Cand** cand, int ncand);
Cand** Read2DCand(char* cand_file,  int nseq, int ** ret_ncand, int * ret_max_cand);
Cand* Read1DCand(char* cand_file,  int * ret_ncand);
char* ExpandFull(char* anno, char* seq);
int Overlap(Cand* c1, Cand* c2, int* olap_min, int* olap_max, int* min, int* max);

/* Transform to Squid format */
SQINFO* Cand2Sqinfo(Cand** cand, int ncand, SQINFO* seq_sqinfo);

int compDouble(const void* a, const void* b);
int CompCandByScore(const void* a, const void* b);
int CompCandByEnergy(const void* a, const void* b);
Cand** SortCand(Cand* cand, int ncand, int (*f) (const void* a, const void* b));
int* GetPairtable(char* ss);
int isHairpin(char* ss);
int isMultiloop(char* ss);
int countHairpin(char* ss);
#endif
