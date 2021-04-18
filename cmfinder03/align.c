#include <stdio.h>
#include "squid.h"
#include "sqfuncs.h"
#include "global.h"

#define DIM 4
int    DNA_mat[DIM][DIM];
int    gap_open;
int    gap_ext;

typedef struct{
  int  row;
  int  col;
}Coor;



typedef struct Match Match;
struct Match{
  Coor   start;
  Coor   stop;
  int    score;
  double eval;
  Match  *next;
};


static struct opt_s OPTIONS[] = {
  { "-m", TRUE, sqdARG_STRING}, 
  { "-e", TRUE, sqdARG_FLOAT}, 
  { "-s", TRUE, sqdARG_INT}, 
  { "-w", TRUE, sqdARG_INT}, 
  { "-f", TRUE, sqdARG_INT}, 
  { "-o", TRUE, sqdARG_STRING},
  { "-v", TRUE, sqdARG_NONE},
  { "--minh", FALSE, sqdARG_INT},  
  { "--maxh", FALSE, sqdARG_INT},  
  { "-t", TRUE, sqdARG_STRING}
};


#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

char usage[] = "\
usage:\n\
align [-m score_matrix] <seqfile>\n";


int** Alloc2DArray(int len1, int len2 )
{
  int **mat;
  int i,j;
  mat = (int **) MallocOrDie (len1 * sizeof (int *));
  for (j = 0; j < len1; j++){
    mat[j] = (int *) MallocOrDie ( len2 * sizeof(int));
    for(i=0; i < len2; i++) mat[j][i]=0;
  }
  return mat;
}

/*
void Free2DArray(int** mat, int len )
{
  int i;
  for(i=0; i < len; i++){
    if (mat[i]) free(mat[i]);
  }
  free(mat);
}
*/

void read_matrix(char* file, int matrix[DIM][DIM], int* ret_gap_open, int* ret_gap_ext)
{
  char buffer[MAXLINE];
  int i,j;
  char* go_tag="gap_open:";
  char* ge_tag="gap_ext:";
  char* sp=NULL;
  int   gap_open, gap_ext;
  FILE*  fin;
  if ( (fin = fopen(file, "r")) == NULL){
    Die("Can't open file %s for reading", file);
  }
  for(i = 0; i < DIM; i++){
    for(j=0; j < DIM; j++){
      fscanf(fin, "%d", &DNA_mat[i][j]);
    }
  }
  
  while(fgets(buffer, MAXLINE, fin)){
    sp = buffer;
    if (strstr(sp, go_tag)) {
      sp += strlen(go_tag);      
      sscanf(sp, "%d", &gap_open);
    }
    if (strstr(sp, ge_tag)) {
      sp += strlen(ge_tag);
      sscanf(sp, "%d", &gap_ext);
    }
  }
  *ret_gap_open = gap_open;
  *ret_gap_ext = gap_ext;
}



void add_match(Match*** ret_matches, int* allocated_match, int* n_match, Match* m)
{
  Match** matches = *ret_matches;
  static int block_size = 100;
  if(*n_match + 1 > *allocated_match){
    *allocated_match += block_size;
    matches = realloc(matches, sizeof(Match*) * (*allocated_match));	    
    *ret_matches = matches;
  }
  matches[*n_match] =m;
  (*n_match)++;
}

double Evalue(int score, int n0, int  n1)
{
  static double Lambda= 0.192;
  static double K=0.177;
  static double H=0.36;
  static double min= 1e-100;
  double p_val;
  double mp, np, a_n0, a_n0f, a_n1, a_n1f;
  int    z_size = (int) 2 * 1000000 / (n0 + n1);  
  a_n0 = (double)n0;
  a_n1 = (double)n1;
  if (H > 0.0) {
    a_n0f = log(a_n0)/H;
    a_n1f = log(a_n1)/H;
  }
  else a_n0f = a_n1f = 0.0;
  mp = a_n0 - a_n0f - a_n1f;
  np = a_n1 - a_n0f - a_n1f;
  if (np < 1.0) np = 1.0;
  if (mp < 1.0) mp = 1.0;
  p_val = K * np * mp * exp ( - Lambda * score);
  if (p_val > 0.01) p_val = 1.0 - exp(-p_val);
  if (p_val < min) p_val = min;  
  return p_val * (double)z_size;
}


int CompMatchByScore(const void* a, const void* b)
{
  Match** atemp = (Match**) a;
  Match** btemp = (Match**) b;
  if ((*btemp)->score > (*atemp)->score )
    return 1;  
  if ((*btemp)->score < (*atemp)->score )
    return -1;  
  return 0;  
}

void copy_coor(Coor* coor1, Coor* coor2)
{
  coor1->row = coor2->row;
  coor1->col = coor2->col;
}



void smith_waterman(char* dsq1, int len1, Match** row_matches, int row_offset,
		    char* dsq2, int len2, Match** col_matches, int col_offset,
		    int score_threshold)
{    
  int* M_mat[2];
  int* Ix_mat[2];
  int* Iy_mat[2];
  Coor*  M_trace[2]; 
  Coor*  Ix_trace[2];
  Coor*  Iy_trace[2];
  int     i,j,k;
  int     cur_row=0;
  int     prev_row=0;
  int     row, col;
  int     score;
  int     n_match = 0;
  Match   m;

  for(i=0; i< 2; i++){
    M_mat[i] =  (int*) MallocOrDie(sizeof(int) * (len2 + 1));
    memset (M_mat[i], 0, sizeof(int) * (len2 + 1));
    Ix_mat[i] = (int*) MallocOrDie(sizeof(int) * (len2 + 1));
    memset (Ix_mat[i], 0, sizeof(int) * (len2 + 1));
    Iy_mat[i] = (int*) MallocOrDie(sizeof(int) * (len2 + 1));
    memset (Iy_mat[i], 0, sizeof(int) * (len2 + 1));

    M_trace[i] = (Coor*)MallocOrDie(sizeof(Coor) * (len2 + 1));
    memset (M_trace[i], 0, sizeof(Coor) * (len2 + 1));
    Ix_trace[i] = (Coor*)MallocOrDie(sizeof(Coor) * (len2 + 1));
    memset (Ix_trace[i], 0, sizeof(Coor) * (len2 + 1));
    Iy_trace[i] = (Coor*)MallocOrDie(sizeof(Coor) * (len2 + 1));
    memset (Iy_trace[i], 0, sizeof(Coor) * (len2 + 1));
  } 

  for(i=0; i <= len1; i++){  
    prev_row = (cur_row + 1) %2;
    for(j=0; j <= len2; j++){
      /* Initialize */
      if (i==0 || j==0) {
	M_mat[cur_row][j] = 0;
	Ix_mat[cur_row][j] = 0;
	Iy_mat[cur_row][j] = 0;
	
      }
      else{
	/*M_mat */
	M_mat[cur_row][j] =0;	  
	M_trace[cur_row][j].row = i;
	M_trace[cur_row][j].col = j;
	
	if (M_mat[cur_row][j] < M_mat[prev_row][j-1]){
	  M_mat[cur_row][j] =  M_mat[prev_row][j-1];
	  memcpy(&M_trace[cur_row][j], &M_trace[prev_row][j-1], sizeof(Coor));
	}
	if (M_mat[cur_row][j] < Ix_mat[prev_row][j-1]){
	  M_mat[cur_row][j] = Ix_mat[prev_row][j-1];
	  memcpy(&M_trace[cur_row][j], &Ix_trace[prev_row][j-1], sizeof(Coor));
	}
	if (M_mat[cur_row][j] < Iy_mat[prev_row][j-1]){
	  M_mat[cur_row][j] = Iy_mat[prev_row][j-1];
	  memcpy(&M_trace[cur_row][j], &Iy_trace[prev_row][j-1],sizeof(Coor));
	}
	M_mat[cur_row][j] += DNA_mat[dsq1[i]][dsq2[j]];	

	/* If a good match*/
	if (M_mat[cur_row][j] > score_threshold && DNA_mat[dsq1[i]][dsq2[j]] > 0 
	    && (i == len1 || j == len2 || DNA_mat[dsq1[i+1]][dsq2[j+1]] < 0)){
	  row = M_trace[cur_row][j].row;
	  col = M_trace[cur_row][j].col;
	  score = M_mat[cur_row][j];	  
	  m.start.row = row + row_offset;
	  m.start.col = col + col_offset;
	  m.stop.row = i + row_offset;
	  m.stop.col=j + col_offset;
	  m.score= score;	  
	  //printf("%d-%d %d-%d score %d\n", row, col, i,j, score);
	  if (!row_matches[row_offset + i]  || score > row_matches[row_offset + i]->score){
	    if (!row_matches[row_offset + i]) {
	      row_matches[row_offset + i]=(Match*)MallocOrDie(sizeof(Match));
	    }
	    memcpy(row_matches[row_offset + i], &m, sizeof(Match));
	  }
	  if (!col_matches[col_offset + j]  || score > col_matches[col_offset + j]->score){
	    if (!col_matches[col_offset + j]) {
	      col_matches[col_offset + j]=(Match*)MallocOrDie(sizeof(Match));
	    }
	    memcpy(col_matches[col_offset + j], &m, sizeof(Match));	    
	  }
	}      
	
	/* Upate Ix_mat */
	if (M_mat[prev_row][j] + gap_open < Ix_mat[prev_row][j] + gap_ext){
	  Ix_mat[cur_row][j] = Ix_mat[prev_row][j] + gap_ext;
	  memcpy(&Ix_trace[cur_row][j], &Ix_trace[prev_row][j], sizeof(Coor));
	}
	else{
	  Ix_mat[cur_row][j] = M_mat[prev_row][j] + gap_open;
	  memcpy(&Ix_trace[cur_row][j], &M_trace[prev_row][j],sizeof(Coor));
	}

	/* Upate Iy_mat */
	if (M_mat[cur_row][j-1] + gap_open < Iy_mat[cur_row][j-1] + gap_ext){
	  Iy_mat[cur_row][j] = Iy_mat[cur_row][j-1] + gap_ext;
	  memcpy(&Iy_trace[cur_row][j], &Iy_trace[cur_row][j-1],sizeof(Coor));
	}
	else{
	  Iy_mat[cur_row][j] = M_mat[cur_row][j-1] + gap_open;
	  memcpy(&Iy_trace[cur_row][j], &M_trace[cur_row][j-1],sizeof(Coor));
	}
      }
    }
    cur_row = (cur_row + 1) %2;      
  }  
  for(i=0; i < 2; i++){
    free(M_mat[i]);
    free(Ix_mat[i]);
    free(Iy_mat[i]);
    free(M_trace[i]);
    free(Ix_trace[i]);
    free(Iy_trace[i]);
  }
}

int is_contained(Match* m1, Match* m2)
{
  static int olap_frac = 0.8;
  int len1_row = m1->stop.row - m1->start.row;
  int len2_row = m2->stop.row - m2->start.row;
  int len1_col = m1->stop.col - m2->start.col;
  int len2_col = m2->stop.col - m2->start.col;
  int olap_start_row = m1->start.row > m2->start.row ? m1->start.row : m2->start.row;
  int olap_start_col = m1->start.col > m2->start.col ? m1->start.col : m2->start.col;
  int olap_stop_row = m1->stop.row < m2->stop.row ? m1->stop.row : m2->stop.row;
  int olap_stop_col = m1->stop.col < m2->stop.col ? m1->stop.col : m2->stop.col;
  int olap_row = olap_stop_row - olap_start_row;
  int olap_col = olap_stop_col - olap_start_col;
      
  if (olap_row > olap_frac * len1_row && olap_col > olap_frac * len1_col){
    return 1;
  }
  if (olap_row > olap_frac * len2_row && olap_col > olap_frac * len2_col){
    return 2;
  }
  return 0;
}

int main(int argc, char* argv[])
{  
  char*   seqfile=NULL; 
  int     nseq;  
  char**  rseqs;
  char**  dsq;
  SQINFO* sqinfo;
  int     format;     
  char*   matrixfile = NULL;
  
  char  *optname;                /* name of option found by Getopt()        */
  char  *optarg;                 /* argument found by Getopt()              */
  int    optind;                 /* index in argv[]                         */	

  int    i0, i1,i2,j0,j1, j2,i,j,k,l,n;
  int    verbose = 0;
  double evalue_threshold = 0.01 ;
  int    min_hits = 30;
  int    max_hits = 200;
  int    window = 200;
  int    offset = 100;
  int    score_threshold = 70;
  char  *outfile=NULL;
  FILE* fout;

  FILE *timerFile=NULL;
  Stopwatch_t *timer=StopwatchCreate();
  StopwatchStart(timer);

  static int local_range = 20;
  

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {    
    if (strcmp(optname, "-m") == 0)    matrixfile     = optarg;    
    if (strcmp(optname, "-e") == 0)    evalue_threshold = atof(optarg);
    if (strcmp(optname, "-s") == 0)    score_threshold = atoi(optarg);
    if (strcmp(optname, "--minh") == 0)    min_hits = atoi(optarg);
    if (strcmp(optname, "--maxh") == 0)    max_hits = atoi(optarg);
    if (strcmp(optname, "-w") == 0)    window = atoi(optarg);
    if (strcmp(optname, "-f") == 0)    offset = atoi(optarg);
    if (strcmp(optname, "-v") == 0)    verbose = 1;
    if (strcmp(optname, "-o") == 0)    outfile= optarg;
    if (strcmp(optname, "-t")==0) {
      timerFile=fopen(optarg,"at");
      if (timerFile==NULL) {
        Die("cannot open file %s\n", optarg);
      }
    }
  }
   
   
  seqfile = argv[optind++];  
  if(! ReadMultipleRseqs(seqfile, SQFILE_FASTA, &rseqs, &sqinfo, &nseq))
    Die("Failed to read squences from file %s", seqfile);
  dsq = (char **) malloc(sizeof(char*) * nseq);
  for (i = 0; i < nseq; i++){
    PrepareSequence(rseqs[i],0);
    dsq[i] = DigitizeSequence(rseqs[i], sqinfo[i].len);
  }         
  read_matrix(matrixfile, DNA_mat, &gap_open, &gap_ext);
  fout = stdout;
  if (outfile){
    fout = fopen(outfile, "w");
    if (fout == NULL) Die ("Can't open file %s for writing", outfile);	  
  }

  for(i=0; i < nseq; i++)
    for(j=i+1; j < nseq; j++){
	  Match** col_matches;
      int len1 = sqinfo[i].len;
      int len2 = sqinfo[j].len;
	  int all_match_count;
	  int n_match;

      Match** all_matches = (Match**) MallocOrDie(sizeof(Match*) * (len1  + len2));      
      Match** row_matches = (Match**)MallocOrDie(sizeof(Match*) * (len1+1));
      memset(row_matches, 0, sizeof(Match*) * (len1+1));
      col_matches = (Match**)MallocOrDie(sizeof(Match*) * (len2+1));
      memset(col_matches, 0, sizeof(Match*) * (len2+1));   

      all_match_count=0;
      //fprintf(stderr, "Align %s %s \n", sqinfo[i].name, sqinfo[j].name);
      for(k = 0; k < len2; k+= offset){
		char* cur_dsq;
        int cur_len = window;
        if (k+ cur_len > len2){
          cur_len = len2 - k;
        }
        cur_dsq = DigitizeSequence(rseqs[j] + k, cur_len);			
	smith_waterman(dsq[i], sqinfo[i].len, row_matches, 0, cur_dsq,cur_len, col_matches, k, score_threshold);
	free(cur_dsq);
      }
      
      /* Check to local optimality */
      for(i0 =1; i0 <= len1; i0 ++){    
		  int i1,i2;
	if (!row_matches[i0]) continue;
	i1 = i0 - local_range;
	if (i1 <= 0) i1 = 1;
	i2 = i0 + local_range;
	if  (i2 > len1 ) i2 = len1;
	for(k = i1 ; k < i2; k++){
	  if (i0 == k) continue;
	  if (row_matches[k] && abs(row_matches[k]->stop.col - row_matches[i0]->stop.col) < local_range){
	    if (row_matches[k]->score < row_matches[i0]->score){
	      free(row_matches[k]);	  
	      row_matches[k] = NULL;
	    }
	    else{
	      free(row_matches[i0]);
	      row_matches[i0] = NULL;
	      break;
	    }
	  }
	}
	if (row_matches[i0]){
	  for(k = i0 + local_range + 1 ; k < i0 + window && k < len1; k++){
	    if (row_matches[k]){
	      if (is_contained(row_matches[i0], row_matches[k])){
		if (row_matches[k]->score <= row_matches[i0]->score){
		  free(row_matches[k]);
		  row_matches[k] = NULL;
		}
		else{
		  free(row_matches[i0]);
		  row_matches[i0] = NULL;
		  break;
		}
	      }
	    }	    
	  }
	}
      }
      for(j0=1; j0 <= len2; j0++){    
		  int j1,j2;
	if (!col_matches[j0]) continue;
	j1 = j0 - local_range;
	if (j1 <= 0) j1 = 1;
	j2 = j0 + local_range;
	if  (j2 > len2 ) j2 = len2;
	for(k = j1 ; k < j2; k++){
	  if (j0 == k) continue;
	  if (col_matches[k] && abs(col_matches[k]->stop.row - col_matches[j0]->stop.row) < local_range){	  
	    if (col_matches[k]->score < col_matches[j0]->score){	  
	      free(col_matches[k]);
	      col_matches[k] = NULL;
	    }
	    else{
	      free(col_matches[j0]);
	      col_matches[j0] = NULL;
	      break;
	    }
	  }
	}
	if (col_matches[j0]){
	  for(k = j0 + local_range + 1 ; k < j0 + window && k < len2; k++){
	    if (col_matches[k]){
	      if (is_contained(col_matches[j0], col_matches[k])){
		if (col_matches[k]->score <= col_matches[j0]->score){
		  free(col_matches[k]);
		  col_matches[k] = NULL;
		}
		else{
		  free(col_matches[j0]);
		  col_matches[j0] = NULL;
		  break;
		}
	      }
	    }	    
	  }
	}
      }
      n_match = 0;
      for(i0=0; i0 <= len1; i0++){
	if (row_matches[i0]){
	  all_matches[n_match++] = row_matches[i0];
	}
      }
      for(j0=0; j0 <= len2 ; j0++){
	if (col_matches[j0]){
	  all_matches[n_match++] = col_matches[j0];            
	}
      }
      free(row_matches);
      free(col_matches);
      
      qsort(all_matches, n_match, sizeof(Match*), CompMatchByScore);      
      all_match_count=0;
      for(n=0; n < n_match;n++){
	Match* m= all_matches[n];
	for (k=0; k < n; k++){
	  if (is_contained(all_matches[k], all_matches[n])) break;
	}
	if (k < n) continue;
	m->eval=Evalue(m->score, len1, len2);
	if (m->eval < evalue_threshold || n < min_hits){
	  fprintf(fout, "%s : %d - %d\t %g\t \t%d\n", 
		  sqinfo[i].name, m->start.row, m->stop.row, m->eval, m->score);      	     
	  fprintf(fout, "%s : %d - %d\n\n", sqinfo[j].name, m->start.col, m->stop.col);
	  all_match_count++;
	}	
	if (all_match_count >= max_hits) break;
      }
      
      for(n=0; n < all_match_count; n++){
	free(all_matches[n]);	     
      }
      free(all_matches);
    }


  StopwatchStop(timer);
  if (timerFile!=NULL) {
    fprintf(timerFile,"cands\t%lg\t%lg\n",timer->user+timer->sys,timer->elapsed);
    fclose(timerFile);
  }
  StopwatchFree(timer);

  if (outfile)fclose(fout);
  for (i = 0; i < nseq; i++)
     FreeSequence(rseqs[i], &(sqinfo[i]));
   free(sqinfo);
   free(dsq);
   free(rseqs);   

   return 0;
}
