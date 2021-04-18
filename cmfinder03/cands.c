#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cand.h"
#include "match_constr.h"
#include "treedist.h"
#include "edit_cost.h"
#include "squid.h"
#include "utils.h"


#define MAX_LENGTH_DIFF  0.3
#define LENGTH_CONSTANT  5
#define SCORE_THRESHOLD  0.2
#define MIN_SCORE        -0.25

static struct opt_s OPTIONS[] = {
  { "-n", TRUE, sqdARG_INT}, 
  { "-f", TRUE, sqdARG_FLOAT}, 
  { "-m", TRUE, sqdARG_STRING}, 
  { "-r", TRUE, sqdARG_STRING},
  { "-l", TRUE, sqdARG_STRING},
  { "-t", TRUE, sqdARG_STRING}
};

#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))


char usage[] = "\
Usage:\n\
cands [-n Seed] [-f fraction]  [-r range_file] [-m match_constrain_file] [-l file that will contain list of cand files created] [-t timer_append_file] <seqfile> <candfile> \n";

int nseed = 0;
int nseq = 0;
Cand**  seeds = NULL;
Cand**  cand = NULL;
int*    ncand = NULL;
Cand*** chosen = NULL;
int*    chosen_size = NULL;
float*   edit_cost;
float*** best_match_cost;
int***   best_match;  
double constr_threshold = 0.001;  
double match_threshold = 0.1;  
long long edit_cost_size;

int compFloat(const void * a, const void* b)
{
  float a1 = *(float*)a;
  float b1 = *(float*)b;
  if (a1 < b1) return 1;
  if (a1 > b1) return -1;
  return 0;
}


long long my_index(int i, int j, int k, int l)
{
  long long idx, idx1, idx2;
  idx1 = (i < k )? i * MAXCAND + j: k * MAXCAND + l;
  idx2 = (i > k )? i * MAXCAND + j: k * MAXCAND + l;      
  idx = idx2 * (idx2 + 1) / 2 + idx1;
  assert(idx>=0);
  return idx;
}


double CheckCandConserve(int start1, int stop1, int start2, int stop2, MatchPtr mc)
{ 
  MatchPtr curr = mc;
  MatchRecord* r;    
  int     mstart=0, mstop=0;
  double  max = 0;
  
  while(curr){
	int ostart1;
	int ostart2;
	int ostop1;
	int ostop2;
	double olap1;
	double olap2;
	double olap;
    r = &(curr->node);
    ostart1 =  start1 > r->start1 ? start1 : r->start1;
    ostart2 =  start2 > r->start2 ? start2 : r->start2;

    ostop1 =  stop1 > r->stop1 ? stop1 : r->stop1;
    ostop2 =  stop2 < r->stop2 ? stop2 : r->stop2;
    
    olap1= ((double)(ostop1 - ostart1))/ (stop1 - start1);
    olap2= ((double)(ostop2 - ostart2))/ (stop2 - start2);
    olap = olap1 < olap2 ? olap1 : olap2;
    if (olap > max){
      //printf("olap %.2f %d / %d\n", olap, ostop - ostop, stop - stop);
      max = olap;
    }
    curr = curr->next;    
  }  
  return max;
}

int choose(Cand* seed, Cand** select, double fraction)
{
  Cand* curr = 0;  
  Cand* best_match= NULL; 
  int   count;  
  int   i, j;  
  short* seq_flags = space(sizeof(short) * nseq);
  long long num_edit_cost_accesses=0;
  memset(seq_flags, 0, sizeof(short) * nseq);
  for (i = 0; i < nseq; i++) 
    for (j = 0; j < ncand[i]; j++) 
      cand[i][j].weight = 0;  
  
  curr = seed;
  count = 1;
  select[curr->seq_id] = curr;
  while (count < nseq) {
    double max_score = INF;
    best_match=NULL;
    seq_flags[ curr-> seq_id] = 1;
    for (i = 0; i < nseq; i++) {
      if (seq_flags[i] ) continue;
      for (j = 0; j < ncand[i]; j++) {
	cand[i][j].weight += edit_cost[ my_index(i, j, curr->seq_id, curr->cand_id)];
        num_edit_cost_accesses++;
	if (cand[i][j].weight > max_score ) {
	  max_score = cand[i][j].weight;
	  best_match = &cand[i][j];
	}
      }
    }
    if (best_match==NULL) {
      break;
    }

    printf("choose: access edit_cost array %Ld times out of %Ld elements\n",num_edit_cost_accesses,edit_cost_size);
    
    curr = best_match;
    printf("select cand in seq %d weight %.2f, count %d, %.2f\n",
	   curr->seq_id, curr->weight, count, curr->weight/count);
    if (curr->weight / count < SCORE_THRESHOLD && count >= nseq * fraction)
      break;
    select[curr->seq_id] = curr;
    count++;
  }
  return count;
}


void SelectCand(double fraction, int max_seed, MatchPtr** match_constr)
{
  Cand**   all_cand;  
  int      total_cand;
  float    score;
  float    dist;  
  int      sflag;  
  int      i, j, k, l, m, n;
  char     *struc;  
  Tree***   trees;
  MatchPtr  curr, prev;
  long long num_edit_cost_accesses=0;

  if (MAXCAND==0) {
    edit_cost_size=0;
    edit_cost=NULL;
  }
  else {
    edit_cost_size=my_index(nseq-2,  MAXCAND -1, nseq-1, MAXCAND -1 ) + 1;
    //printf("edit_cost_size=%Ld (%d,%d,%d,%d)\n",edit_cost_size,nseq-2,  MAXCAND -1, nseq-1, MAXCAND -1);
    edit_cost = (float*) calloc( edit_cost_size,sizeof(float*));
    if (edit_cost==NULL) {
      fprintf(stderr,"could not alloc memory for edit_cost\n");
      exit(1);
    }
  }
  best_match_cost = (float***) space(sizeof(float**) * nseq);
  best_match = (int***) space(sizeof(int**) * nseq);
  trees = (Tree***) space(sizeof(Tree**) * nseq);
  seeds = (Cand**) space(sizeof(Cand*) * max_seed);  

  for(i=0; i< nseq; i++) {
    best_match_cost[i] = (float**) space(sizeof(float*) * ncand[i]);    
    best_match[i] = (int**) space(sizeof(int*) * ncand[i]);    
    trees[i] = (Tree**) space(sizeof(Tree*) * ncand[i]);    
    for(j =0 ; j < ncand[i]; j++) {
      best_match_cost[i][j] = (float*) space(sizeof(float) * nseq);
      best_match[i][j] = (int*) space(sizeof(int) * nseq);      
      struc = ExpandFull(cand[i][j].ss, cand[i][j].seq);      
      trees[i][j] = make_tree(struc);      
      free(struc);
      
      for(k=0; k < nseq; k++) {	
	best_match_cost[i][j][k] = INF;
	best_match[i][j][k] = -1;
      }
    }    
  }
      
  for( i=0; i < nseq; i++) { 
    for( j=0; j < ncand[i]; j++) {
      for(k = i+1; k < nseq; k++) {
	if (k==i) continue;	
	for( l = 0; l < ncand[k]; l++) {
	  double dist;
	  /* If length of the two candidates differ significantly, no need to compare them */
	  int l1 = abs(cand[i][j].stop - cand[i][j].start)+1;
	  int l2 = abs(cand[k][l].stop - cand[k][l].start)+1;
	  int length_diff = abs(l1 - l2);
	  if ( length_diff > MAX_LENGTH_DIFF * l1 ||length_diff > MAX_LENGTH_DIFF  * l2) {
	    edit_cost[my_index(i,j,k,l)] = INF;
            num_edit_cost_accesses++;
	    continue;
	  }
	  if(match_constr && match_constr[i][k]) {
	    if (!CheckMatch(cand[i][j].start, cand[i][j].stop, cand[k][l].start, cand[k][l].stop,
			    match_constr[i][k], constr_threshold)) {
	      edit_cost[my_index(i,j,k,l)] = INF;
              num_edit_cost_accesses++;
	      continue;
	    }
	  }

	  /* The distance of two candidates is the tree-edit distance normalized by length */
	  dist = tree_edit_distance(trees[i][j], trees[k][l]);      
	  dist /= sqrt(l1 * l2) + LENGTH_CONSTANT;	  	  

	  if (match_constr){
	    double olap=CheckCandConserve(cand[i][j].start, cand[i][j].stop, 
					  cand[k][l].start, cand[k][l].stop, 
					  match_constr[i][k]);
	    /*
	    if (olap > 0)	      
	      printf("Cand %d_%d %d-%d : %d_%d %d-%d olap %.2f dist %.2f\n", 
		     i, j, cand[i][j].start, cand[i][j].stop, 
		     k, l, cand[k][l].start, cand[k][l].stop, 
		     olap, dist);	    
	    */
	    dist *= (1 + olap); 
	  }
					 
          long long index=my_index(i, j, k, l);
	  edit_cost[ index ] = dist;	  
          num_edit_cost_accesses++;
	  if (dist > best_match_cost[i][j][k]){
	    best_match_cost[i][j][k] = dist;
	    best_match[i][j][k] = l;
	  }
	  if (dist > best_match_cost[k][l][i]) {
	    best_match_cost[k][l][i] = dist;
	    best_match[k][l][i] = j;
	  }
	}
      }
      
    }
  }    
 
  for(i=0, total_cand = 0; i < nseq; i++) total_cand += ncand[i];  
  all_cand = (Cand **) malloc( sizeof(Cand *) * total_cand);  
  for(i=0, k=0; i < nseq; i++)
    for(j = 0; j < ncand[i]; j++) 
      all_cand[k++] = &cand[i][j];          
  
  //printf("Choose\n");  

  /* Choose seed candidate */
  for( nseed = 0; nseed < max_seed && total_cand > 0; nseed ++) {    /* NOTE: used to be a single Ampersand, I changed it to double, since it looks like it's doing a logical Boolean op */
    for(m=0; m < total_cand; m++){      
	  double* match_cost;
      i = all_cand[m]->seq_id;
      j = all_cand[m]->cand_id;
      match_cost= (double*)MallocOrDie(sizeof(double) * nseq);      
      for(k = 0; k < nseq ; k++) {
	int max_olap=0;	  
	double frac = 1;
	match_cost[k] = INF;
	if (i == k) continue;	  
	if (best_match[i][j][k] != -1){
	  Cand* match = & cand[k][best_match[i][j][k]];	  
	  //Count the overlap with other seed alignment
	  if (match){
	    for (n=0; n < nseed; n++) {
	      int t1,t2,t3,t4;
	      if (chosen[n][k]){
		int olap = Overlap(match, chosen[n][k], &t1,&t2,&t3,&t4);
		if (olap > max_olap){
		  max_olap=olap;
		}
	      }
	    }
	    frac = ( 1- max_olap / (double)(match->len));			
	  }
	}
	match_cost[k] = best_match_cost[i][j][k] * frac;
      }            
      qsort(match_cost, nseq, sizeof(double), compDouble);
      score = 0;
      for(k=0; k < nseq; k++){	
	if (match_cost[k] < SCORE_THRESHOLD && k >= nseq * fraction) 
	  break;	
	if (match_cost[k] < MIN_SCORE){
	  match_cost[k] = MIN_SCORE;
	}
	score += match_cost[k];	
      }
      free(match_cost);
      all_cand[m]->score = score / (nseq -1);                        
      //printf("%d_%d score %.2f %.2f\n", all_cand[m]->seq_id, all_cand[m]->cand_id, score, all_cand[m]->score);
    }
    
    qsort(all_cand, total_cand, sizeof(Cand*), CompCandByScore);         
    if (all_cand[0]->score < SCORE_THRESHOLD){
      break;
    }
    seeds[nseed] = all_cand[0];
    
    printf("Seq_%d_Cand%d_%d_%d  %f\n", seeds[nseed]->seq_id, seeds[nseed]->cand_id, 
	   seeds[nseed]->start, seeds[nseed]->stop, seeds[nseed]->score);
    printf("%s\n", seeds[nseed]->seq);
    printf("%s\n", seeds[nseed]->ss);    


    chosen[nseed] = (Cand**)space(sizeof(Cand*) * nseq);
    memset(chosen[nseed], 0, sizeof(Cand*) * nseq);          
    chosen_size[nseed] = choose(all_cand[0], chosen[nseed], fraction);
    //Remove candidates in the seed alignment
    for(m=0; m < total_cand; m++) {
      i = all_cand[m]-> seq_id;
      if (all_cand[m] == chosen[nseed][i]) {
	all_cand[m] = all_cand[total_cand-1];
	total_cand -- ;
      }
    }
  }
  //printf("Free\n");
  free(all_cand);
  if (match_constr) {
    for(i=0; i< nseq; i++) {
      for(j=0; j < nseq; j++) {
	curr = match_constr[i][j];
	while(curr) {
	  prev = curr;
	  curr = curr->next;
	  free(prev);
	}
      }
      free(match_constr[i]);
    }
    free(match_constr);
  }
  
  for(i=0; i < nseq; i++) {
    
    for(j=0; j < ncand[i]; j++) {
      free(best_match_cost[i][j]);
      free(best_match[i][j]);
      free_tree(trees[i][j]);
    }
    free(best_match_cost[i]);
    free(best_match[i]);
    free(trees[i]);
  }
  
  free(edit_cost);
  free(best_match_cost);
  free(best_match);
  free(trees);

  printf("SelectCand: access edit_cost array %Ld times out of %Ld elements\n",num_edit_cost_accesses,edit_cost_size);
}


int main(int argc, char* argv[])
{
  int     max_seed = 1;  
  char*   seqfile= NULL;
  char*    candfile; 
  char*   matchfile=NULL;
  char*   clusterfile=NULL;
  char*   fileListFileName=NULL;
  FILE*   fileListFile=NULL;
  char**  rseqs;
  SQINFO* sqinfo;
  int     format; 
  MatchPtr** match_constr = NULL;
  MatchPtr** match_block = NULL;
  int   i,j,k;  
  int   max_cand;  
  Range* range=NULL;

  int          optc;
  float        fraction = 1;  

  char  *optname;                /* name of option found by Getopt()        */
  char  *optarg;                 /* argument found by Getopt()              */
  int    optind;                 /* index in argv[]                         */
  FILE *timerFile=NULL;
  Stopwatch_t *timer=StopwatchCreate();
  StopwatchStart(timer);

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {

    if (strcmp(optname, "-m") == 0)  matchfile     = optarg;   
    else if (strcmp(optname, "-r") == 0)  clusterfile     = optarg;     
    else if (strcmp(optname, "-n") == 0)  max_seed         = atoi(optarg);  
    else if (strcmp(optname, "-f") == 0)  fraction      = atof(optarg);
    else if (strcmp(optname, "-l") == 0)  fileListFileName=optarg;
    else if (strcmp(optname, "-t")==0) {
      timerFile=fopen(optarg,"at");
      if (timerFile==NULL) {
        Die("cannot open file %s\n", optarg);
      }
    }
    else{
       Die("Invalid Option %s\n %s", optname, usage);      
     } 
  }

  if (argc - optind < 2)
    Die("%s\n", usage);
    
  seqfile = argv[argc-2];
  candfile=(char *)malloc(strlen(argv[argc-1])+64); /* +64 since we sometimes append a number to it, and +64 is much more than needed.  BTW, we need to alloc memory because we sprintf into candfile */
  if (candfile==NULL) {
    Die("couldn't alloc memory for candfile string");
  }
  strcpy(candfile, argv[argc-1]);

  if (fileListFileName!=NULL) {
          fileListFile=fopen(fileListFileName,"wt");
          if (fileListFile==NULL) {
            Die("cannot open -l file \"%s\" for writing",fileListFileName);
          }
  }
  
   /* Read Sequence file */
  format = SQFILE_FASTA;   
  
  
  if(! ReadMultipleRseqs(seqfile, format, &rseqs, &sqinfo, &nseq))
    Die("Failed to read squences from file %s", seqfile);

  if (clusterfile){
     range = ReadCluster(clusterfile, sqinfo, nseq);     
  }

  /* Read Candidates */ 			 
  cand = Read2DCand(candfile,  nseq,  &ncand, &max_cand);

  if (max_cand < MAXCAND)
    MAXCAND = max_cand;

  /* Read match constraint */
  if (matchfile) {    
    match_constr = ReadMatchConstr(matchfile, nseq, sqinfo, match_threshold,range);
    /*
    for(i=0; i < nseq; i++){
      for(j=0; j < ncand[i]; j++){
	cand[i][j].conserve_num = 1;
	for(k=0; k < nseq; k++){
	  if (i==k) continue;
	  if(match_constr[i][k]){
	    double conserve = CheckConserve(cand[i][j].seq_id, 
					    cand[i][j].start, cand[i][j].stop, match_constr[i][k]);
	    cand[i][j].conserve_num += conserve;	    
	  }
	}
      }
    }
    */
  }
  
  chosen = (Cand***) malloc( sizeof(Cand**) * max_seed );  
  chosen_size = (int*) malloc( sizeof(int) * max_seed );  
  SelectCand(fraction, max_seed, match_constr);
  for(i = 0; i < nseed; i++) {
    sprintf(candfile, "%s_%d", argv[argc-1], i+1);    
    // Put the seed candidate first
    chosen[i][seeds[i]->seq_id] = chosen[i][0];
    chosen[i][0] = seeds[i];    
    Write1DCand(candfile, chosen[i], nseq);
    if (fileListFileName!=NULL) {
            fprintf(fileListFile,"%s\n",candfile);
    }
    free(chosen[i]);        
  }

  StopwatchStop(timer);
  if (timerFile!=NULL) {
    fprintf(timerFile,"cands\t%lg\t%lg\n",timer->user+timer->sys,timer->elapsed);
    fclose(timerFile);
  }
  StopwatchFree(timer);

  if (fileListFile!=NULL) {
          fclose(fileListFile);
  }
  
  free(chosen);  
  free(seeds);  
  free(chosen_size);
  for(i=0; i <nseq; i++) {    
    free(cand[i]);
  }
  
  free(cand);  
  free(ncand);  
  return 0;
}
