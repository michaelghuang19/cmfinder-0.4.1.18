#include "global.h"
#include <ctype.h>

int IsBasePair(char l, char r){
  if (l=='A' && r == 'U' ||
      l=='U' && r == 'A' ||
      l=='G' && r == 'C' ||
      l=='C' && r == 'G' ||
      l=='G' && r == 'U' ||
      l=='U' && r == 'G'){
    return 1;
  }
  return 0;
}

int TriIndex(int i,int j){
  int coor = ((i < j)? ( j * (j-1)/2 + i) : ( i * (i-1)/2 + j));
  return coor;
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
	Die("Structure %s unbalanced base pair at pos %d", ss, i);
      }
      j = stack[sp];
      pt[j] = i;
      pt[i] = j;      
    }	    
  }
  return pt;  
}


Range* ReadCluster(char* cluster_file, SQINFO* sqinfo, int nseq)
{
  FILE *fin;
  int  i;
  char name[MAXLINE];
  char buffer[MAXLINE];
  int  cluster_count=-1;
  int  alloc_cluster=10;
  Range* range = (Range*)MallocOrDie(sizeof(Range) * nseq);  

  if ( (fin = fopen(cluster_file, "r") ) == NULL)
    Die("Fail to read file %s", cluster_file);

  for(i=0; i < nseq; i++){
    range[i].start = -1;
    range[i].stop = -1;
  }
    
  while(fgets(buffer,MAXLINE, fin) > 0) {
    if (isspace(buffer[0])){
      continue;
    }
    else{
      int start;
      int stop;
      memset(name, 0, sizeof(char)*MAXLINE);
      sscanf(buffer, "%s %d %d", name, &start, &stop);    
      for(i=0; i < nseq; i++) {
	if(strncmp(sqinfo[i].name, name, strlen(sqinfo[i].name))==0)
	  break;
      }
      if (i== nseq){
	Warn("Unknow sequence name %s in file %s", name, cluster_file);
	continue;
      }
      range[i].start = start;
      range[i].stop = stop;    
      if (stop > sqinfo[i].len){
	range[i].stop = sqinfo[i].len;
      }
    }
  }
  return range;
}




void WriteCluster(Range* range, char* cluster_file, SQINFO* sqinfo, int nseq)
{
  int i;
  FILE* fout = fopen (cluster_file, "w");
  if (fout == NULL) Die("Can't open file %s for writing ", cluster_file);
  for(i=0; i < nseq; i++){
    if (range[i].start==-1) continue;
    fprintf(fout, "%s %d %d\n", sqinfo[i].name, range[i].start, range[i].stop);    
  } 
  fclose(fout);
}

int
PrepareSequence(char *seq,int use_fragmentary)
{
  char *sym;

  if (use_fragmentary) {
	  Die("implement me");
  }

  for (sym = seq; *sym != '\0'; sym++)
    {
      *sym = toupper((int)*sym);
      if (*sym == 'T') *sym= 'U';
      /* sym is in alphabet, or a gap? ok, go to next one */
      if (isgap(*sym)) continue;

      /* then it's a degenerate symbol.
       * According to alphabet, choose a single symbol to represent it.
       * watch out for too-clever scheme for random choice: "ABC"[random() % 3]
       */
      switch (*sym) {
      case 'B': *sym = "CGU"[CHOOSE(3)]; break;
      case 'D': *sym = "AGU"[CHOOSE(3)]; break;
      case 'H': *sym = "ACU"[CHOOSE(3)]; break;
      case 'K': *sym = "GU" [CHOOSE(2)]; break;
      case 'M': *sym = "AC" [CHOOSE(2)]; break;
      case 'R': *sym = "AG" [CHOOSE(2)]; break;
      case 'S': *sym = "CG" [CHOOSE(2)]; break;
      case 'T': *sym = 'U';                 break;
      case 'V': *sym = "ACG"[CHOOSE(3)]; break;
      case 'W': *sym = "AU" [CHOOSE(2)]; break;
      case 'Y': *sym = "CU" [CHOOSE(2)]; break;
      case 'N': *sym = Alphabet[CHOOSE(4)];
        break;
      }
    }
  return 1;
}


void
Avg_bppr(char    **aseq,            /* array of aligned sequences, flushed right  */
	 int       nseq,		/* number of aligned sequences */
	 int       alen,		/* length of each sequence (all the same)  */
	 double  **bp_pr,
	 float*    weights,    
	 double ***ret_bppr)        /* RETURN: bppr array           */
{
  double **bppr;
  char   *nogap_seq;
  int    *idx_map;
  int    k,i,j, i1, j1;
  double tot_weight=0;
  bppr = DoubleAlloc2DArray(alen+1);

  for(j=1; j < alen+1; j++){
    for(i=0; i <j ; i++)
      bppr[j][i] = 0;
  }
  for(k=0;k< nseq; k++) tot_weight+= weights[k];

  if (bp_pr==NULL){
    bp_pr = (double**)MallocOrDie(sizeof(double*) * nseq);
    memset(bp_pr, 0, sizeof(double**) * nseq); 
  }
  for (k = 0; k < nseq; k++) {
    nogap_seq = remove_gap(aseq[k], &idx_map);
    if (bp_pr[k]==NULL){
      bp_pr[k]= bppr_seq(nogap_seq);
    }
    if (bp_pr[k] == NULL) continue;    
    for(j = 1; j < alen; j++)
      for(i = 0; i < j; i++){
	i1 = idx_map[i];
	j1 = idx_map[j];
	if (i1 >=0 && j1 >= 0){
	  bppr[j+1][i+1] += bp_pr[k][TriIndex(i1,j1)] * weights[k]/tot_weight;
	}
      }
    free(nogap_seq);
    free(idx_map);
  }
  *ret_bppr = bppr;
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
    if (*s1 == '.' || *s1 == '-' || *s1=='_') {
      idx_map[s1 - seq ] = -1;
      s1++;
    }
    else{
      idx_map[s1 - seq ] = s2 - s ;
      *s2 ++ = toupper(*s1++);
    }
  }
  *s2 = '\0';
  *ret_idx_map = idx_map;
  return s;
}


double* bppr_seq(char* seq){
  int i,j;
  int size;
  double* bp_pr;
  char* structure;
  int len = (int)(strlen(seq));
  if (len == 0) return NULL;
  size = TriIndex(len,len-1) - 1;
  bp_pr = (double*) MallocOrDie(sizeof(double)*size);
  memset(bp_pr, 0, sizeof(double) * size);
  structure = (char *) malloc( sizeof(char) * (len + 1));
  init_pf_fold(len);
  pf_fold(seq, structure);
  for(j = 1; j < len; j++)
    for(i = 0; i < j; i++){
      bp_pr[TriIndex(i,j)] += pr[iindx[i + 1] - (j + 1)];
    }
  free(structure);
  free_pf_arrays();
  return (bp_pr);
}
