#ifndef  MATCH_CONSTR_H
#define MATCH_CONSTR_H
#include "global.h"

typedef struct 
{
  int       seq_id1;
  int       seq_id2;
  int       start1;
  int       stop1;
  int       start2;
  int       stop2;    
  double    e_val;  
  int       valid;
}
MatchRecord;

struct MatchNode;
typedef struct MatchNode* MatchPtr;

struct MatchNode
{
  MatchRecord node;
  MatchPtr     next;  
}
; 


MatchPtr** ReadMatchConstr(char* matchfile, int nseq, SQINFO* sqinfo, double threshold, Range* range);
int CheckMatch(int start1, int stop1, int start2, int stop2, MatchPtr mc, double threshold);
double CheckConserve(int seq_id, int start, int stop, MatchPtr mc);

#endif


   
