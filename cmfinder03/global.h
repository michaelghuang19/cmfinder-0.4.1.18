#ifndef GLOBAL_H
#define GLOBAL_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "squid.h"
#include "structs.h"
#include "hmmband.h"         
#include "cm_postprob.h"

#include "part_func.h"
#include "utils.h"
#include "fold_vars.h"
#include "global.h"

#define MAXLINE  1000
#define NEGINFINITY  -9999999
#define POSINFINITY  9999999
#define FLT_MAX      1000000
#define FLANK        100
#define MIN_WEIGHT   2.5

#define pair_left(c) (c == '(' || c == '<' || c == '{')
#define pair_right(c) (c == ')' || c == '>' || c == '}')
#ifndef max
#define max(x,y) (x < y ? y : x)
#define min(x,y) (x > y ? y : x)
#endif

#define IsBasePairDigital(l, r) (l + r== 3 || l+r == 5)

extern int IsBasePair(char l, char r);
extern int* GetPairtable(char* ss);
extern int  TriIndex(int i, int j);
extern void Avg_bppr(char    **aseq,    
		int       nseq,		
		int       alen,
		double  **bp_pr,
		float*    weights,    
		double ***ret_bppr);
extern char*  remove_gap(char* seq,  int** ret_idx_map);
extern  double* bppr_seq(char* seq);

int** IntAlloc2DArray( int alen);
double** DoubleAlloc2DArray( int alen);

extern int   Alphabet_type;
extern int   Alphabet_size;
extern int   Alphabet_iupac;
extern char *Alphabet;

typedef struct{
  int   start;
  int   stop;
}
Range;

extern int PrepareSequence(char *seq,int use_fragmentary);
Range*   ReadCluster(char* cluster_file, SQINFO* sqinfo, int nseq);
void     WriteCluster(Range* range, char* cluster_file, SQINFO* sqinfo, int nseq);

typedef struct cplan9_s  HMM;
typedef CP9Map_t   HMM_CM_Map;    
typedef struct cp9_dpmatrix_s HMM_Matrix;
typedef CP9Bands_t HMM_Band;



#define NO_BP -1
#define UNDEF_BP -2


#endif  
