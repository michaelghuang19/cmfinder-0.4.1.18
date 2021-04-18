#include "squid.h"
#include "global.h"
#include <string.h>

#define MAXSEQ   100
#define DIRECT 0
#define INDIRECT 1
#define NEIGHBORHOOD 2

#define SCALE_EVAL 100
#define MAX_SCORE  200
#define MIN_SCORE  0
#define SCORE_THRESHOLD 50
#define RANGE      300
#define MIN_SEQ    3

typedef struct Match Match;
typedef struct MatchNode MatchNode;
typedef struct Segment Segment;
typedef struct SegmentNode SegmentNode;

struct Segment{
  int    start;
  int    stop;  
  int    len;
  int    seq_id;
  int    segment_idx;
  int    value;       //reserved for algorithm. 
  int    selected_len;
  int    accu_value;       //reserved for algorithm. 
  int    selected;
  SegmentNode* neighbor;
  MatchNode* match[MAXSEQ];  //matches to all other sequences
};

int overlap(Segment* s1, Segment* s2, int* ret_range)
{
  int min_start;
  int max_start;
  int max_stop;
  int min_stop;
  int olap;
  if  (s1->seq_id != s2->seq_id){
    return 0;
  }
  min_start = min(s1->start, s2->start);
  max_start = max(s1->start, s2->start);
  max_stop =  max(s1->stop, s2->stop);
  min_stop = min(s1->stop, s2->stop);

  olap = min_stop - max_start + 1;
  if (olap < 0) olap = 0;
  if (ret_range){
    *ret_range = max_stop - min_start + 1;
  }
  return olap;
}

struct SegmentNode
{
  Segment* segment;
  SegmentNode* next;
};

struct Match
{  
  Segment* s1;
  Segment* s2;
  int      weight;
  int      type;
  int      neighborhood_weight;
};

struct  MatchNode
{
  Match* match;
  MatchNode* next;
};

SegmentNode* all_segments[MAXSEQ];
char*        all_seqnames[MAXSEQ];
int          nseq;
int          nsegments[MAXSEQ];	

void insert_segment_simple(Segment* s, SegmentNode** head, SegmentNode** tail)
{
  SegmentNode* n = (SegmentNode*) MallocOrDie(sizeof(SegmentNode));
  n->segment = s;
  n->next = NULL;
  if (*head == NULL){
    *head = n;
    *tail = n;
    return;
  }  
  (*tail)->next = n;
  *tail = n;
  return;
}


void insert_segment(Segment* s, SegmentNode** head)
{
  SegmentNode* cur;
  SegmentNode* prev;
  SegmentNode* n = (SegmentNode*) MallocOrDie(sizeof(SegmentNode));
  n->segment = s;
  n->next = NULL;
  if (*head == NULL){
    *head = n;
    return;
  }
  cur = *head;
  prev = NULL;
  while(cur){
    if (cur->segment->start > s->start) break;
    prev = cur;
    cur= cur->next;
  }
  if (prev){
    prev->next = n;    
  }
  else{
    *head = n;
  }  
  n->next = cur;
}

Segment* create_segment(int seq_id, int start, int stop)
{
  Segment* seg = (Segment*) MallocOrDie(sizeof(Segment));
  seg->seq_id = seq_id;
  seg->start = start;
  seg->len = abs(stop - start) + 1;
  seg->stop = stop;
  seg->neighbor= NULL;
  seg->selected_len = 0;
  seg->selected = 0;
  insert_segment(seg, &(seg->neighbor));
  memset(seg->match, 0, sizeof(MatchNode*)* MAXSEQ);  
  return seg;
}

Match* create_match(Segment* s1, Segment* s2, int weight, int type)
{
  Match* m = (Match*) MallocOrDie(sizeof(Match));
  m->s1 = s1;
  m->s2 = s2;
  m->weight = weight;
  m->type = type;
  m->neighborhood_weight = weight;
  return m;    
}

void insert_match(Match* m, MatchNode** head)
{
  MatchNode* n = (MatchNode*) MallocOrDie(sizeof(MatchNode));
  n->match = m;
  n->next = NULL;
  if (*head ){
    n->next = *head;
  }  
  *head = n;
  while(n){
    n= n->next;
  }
}

int eval_to_weight(double e_val)
{
  int weight = (int) (log(SCALE_EVAL) - log(e_val));  
  if (weight > MAX_SCORE) weight = MAX_SCORE;  
  return weight;
}

Match* get_match(Segment* s1, Segment* s2)
{
  MatchNode* m = s1->match[s2->seq_id];
  Match* found1 = NULL;
  while(m){
    if (m->match->s2 == s2 || m->match->s1 == s2){
      found1 = m->match;
      break;
    }
    m = m->next;
  }
  return found1;
}

void read_match(char* filename)
{
  char buffer[MAXLINE];
  char name1[100];
  char name2[100];
  int start1, start2, stop1, stop2, score;
  double e_val;
  FILE*  fin;
  int    i;
  SegmentNode* n1;
  if ( (fin = fopen(filename, "r") ) == NULL){
    printf("Fail to read file %s", filename);
    exit(-1);
  }  
  while(fgets(buffer,MAXLINE, fin) > 0) {
	int  seq_id1;
	int  seq_id2;
	Segment* s1;
	Segment* s2;
	Match* m;
	SegmentNode* n2;
    //Read the file
    if (strlen(buffer) < 5) continue;
    memset(name1, 0, sizeof(char)*100);
    memset(name2, 0, sizeof(char)*100);    
    sscanf(buffer, "%s :\t %d - %d\t %lf\t %d", name1, &start1, &stop1, &e_val, &score);    
    if (fgets(buffer,MAXLINE, fin) > 0) {
      sscanf(buffer, "%s :\t %d - %d", name2, &start2, &stop2);
    }
    else break;
    if (strcmp(name1, name2)== 0) continue;
    seq_id1=-1;
    seq_id2=-1;
    for(i=0; i < nseq; i++){
      if (strcmp(all_seqnames[i], name1) == 0){
	seq_id1 = i;
	break;
      }
    }
    if (seq_id1 == -1){
      all_seqnames[nseq] = (char*) MallocOrDie(sizeof(char) * (strlen(name1) + 1));      
      strcpy(all_seqnames[nseq],name1);      
      all_segments[nseq] = NULL;
      nsegments[nseq]= 0;
      seq_id1 = nseq;
      nseq++;
    }
    for(i=0; i < nseq; i++){
      if (strcmp(all_seqnames[i], name2) == 0){
	seq_id2 = i;
	break;
      }
    }
    if (seq_id2 == -1){
      all_seqnames[nseq] = (char*) MallocOrDie(sizeof(char) * (strlen(name2) + 1));      
      strcpy(all_seqnames[nseq],name2);      
      all_segments[nseq] = NULL;
      nsegments[nseq]= 0;
      seq_id2 = nseq;
      nseq++;
    }
    s1 = create_segment(seq_id1,start1, stop1);
    s2 = create_segment(seq_id2,start2, stop2);
    if (score > MAX_SCORE) score = MAX_SCORE;
    score += -log10(e_val);
    m = create_match(s1, s2, score, DIRECT);

    n1 = all_segments[seq_id1];
    while(n1){
      int range =0;
      int olap = overlap(n1->segment, s1, &range);
      if (range < RANGE){
	insert_segment(n1->segment, &(s1->neighbor));
	insert_segment(s1, &(n1->segment->neighbor));
      }
      n1 = n1->next;
    }

    n2 = all_segments[seq_id2];
    while(n2){
      int range = 0;
      int olap = overlap(n2->segment, s2, &range);
      if (range < RANGE){
	insert_segment(n2->segment, &(s2->neighbor));
	insert_segment(s2, &(n2->segment->neighbor));
      }
      n2 = n2->next;
    }
    insert_segment(s1, &(all_segments[seq_id1]));
    insert_segment(s2, &(all_segments[seq_id2]));    
    insert_match(m, &(s1->match[seq_id2]));
    insert_match(m, &(s2->match[seq_id1]));
    nsegments[seq_id1]++;
    nsegments[seq_id2]++;
  }
}



int choose(Segment* seed, Segment** select)
{
  Segment* curr = 0;  
  Segment* best_match; 
  SegmentNode* n;
  int   count;  
  int   i, j;  
  short* seq_flags = (short*)MallocOrDie(sizeof(short) * nseq);
  memset(seq_flags, 0, sizeof(short) * nseq);
  for( i = 0; i < nseq; i++) {
    select[i] = NULL;
    n = all_segments[i];
    while(n){
      n->segment->value = 0;
      n = n->next;
    }
  }

  curr = seed;
  count = 1;
  select[curr->seq_id] = curr;    
  while(count < nseq) {
    double max_score = NEGINFINITY;
    best_match=NULL;
    seq_flags[ curr-> seq_id] = 1;
    for( i = 0; i < nseq; i++) {
      if (seq_flags[i] ) continue;      
      for(n = all_segments[i]; n; n = n->next){
	Match* m = get_match(curr, n->segment);
	if (m){
	  n->segment->value += m-> neighborhood_weight;
	  if (m && m->neighborhood_weight > 0){
	    /*
	    fprintf(stderr, "Match %s %d-%d: %s %d-%d value %d accu %d\n", 
		    all_seqnames[m->s1->seq_id], m->s1->start, m->s1->stop, 
		    all_seqnames[m->s2->seq_id], m->s2->start, m->s2->stop, 
		    m->neighborhood_weight, n->segment->value);
	    */
	  }
	}
	else{
	  n->segment->value += MIN_SCORE;
	}
	if (n->segment->value > max_score ){
	  max_score = n->segment->value;
	  best_match = n->segment;
	}
      }      
    }
    if (best_match==NULL || 
	best_match->value/count < SCORE_THRESHOLD && best_match->value < SCORE_THRESHOLD * 3 ) 
      break;
    curr = best_match;
    
    //fprintf(stderr, "Select %s %d-%d value %d\n", all_seqnames[curr->seq_id], curr->start, curr->stop, curr->value/count);
	    

    select[curr->seq_id] = curr;    
    count++;
  } 
  
  return count;  
}


void write_match(FILE* fout, Match* m)
{
  fprintf(fout, "%s:%d-%d\t %s:%d-%d\t score %d (type %d) neighbor_score %d\n", 
	  all_seqnames[m->s1->seq_id], m->s1->start, m->s1->stop,
	  all_seqnames[m->s2->seq_id], m->s2->start, m->s2->stop,
	  m->weight, m->type, m->neighborhood_weight);	  	  
}

void write_all_match(FILE* fout)
{
  int i,j;
  for(i=0; i < nseq; i++){
    SegmentNode* n;
    for(n = all_segments[i];n; n= n->next){
      fprintf(fout, "Segment %s:%d-%d\n", all_seqnames[n->segment->seq_id], n->segment->start, n->segment->stop);
      for(j=0; j < nseq; j++){
	MatchNode* mn;
	for(mn = n->segment->match[j];mn; mn= mn->next){
	  Match* m = mn->match;
	  write_match(fout, m);
	}	
      }
    }
  }
}

void neighbor_score()
{
  int i,j;
  for(i=0; i < nseq; i++){    
    SegmentNode* n;
    for(n = all_segments[i]; n; n = n->next){
      MatchNode* mn;
      for(j=0; j < nseq; j++){
	for(mn = n->segment->match[j]; mn!= NULL; mn = mn->next){
	  Match* m = mn->match;	  
	  if (m->type== DIRECT){
	    SegmentNode *n1;		
	    SegmentNode *n2;
	    for(n1 = m->s1->neighbor; n1 != NULL ; n1 = n1->next)
	      for(n2 = m->s2->neighbor; n2 != NULL; n2= n2->next){
		float  frac1;
		float  frac2;
		float  new_weight;
		Match* old_m = get_match(n1->segment, n2->segment);
		int    range1, range2, olap1, olap2;
		olap1 = overlap(n1->segment, m->s1, &range1);
		olap2 = overlap(n2->segment, m->s2, &range2);
		//Overlap weight
		frac1 = ((double)(olap1 + olap2))/(m->s1->len + m->s2->len);
		//Nonoverlap weight
		frac2 = ((double)(m->s1->len + m->s2->len))/(range1 + range2);
		new_weight = m->weight * (frac1 * 0.8 + frac2 * 0.2);
		if (old_m == NULL){
		  Match* new_m = create_match(n1->segment, n2->segment, new_weight, NEIGHBORHOOD);
		  insert_match(new_m, &(n1->segment->match[n2->segment->seq_id]));
		  insert_match(new_m, &(n2->segment->match[n1->segment->seq_id]));
		}
		else if (old_m->neighborhood_weight < new_weight){
		  old_m->neighborhood_weight = new_weight;
		}
	      }
	  }
	}
      }
    }
  }
}


void find_best_match()
{
  Segment* best_segment = NULL;
  int i,j;
  for(i=0; i < nseq; i++){    
    SegmentNode* n;
    for(n = all_segments[i]; n; n = n->next){
      Segment* s = n->segment;
      int     sum_weight = 0;
      for(j=0; j < nseq; j++){
	MatchNode* best_match=NULL;
	MatchNode* best_match_prev= NULL;
	MatchNode* prev=NULL;
	MatchNode* mn;	
	for(mn = s->match[j]; mn!= NULL; prev= mn, mn = mn->next){
	  if (best_match == NULL || mn->match->neighborhood_weight > best_match->match->neighborhood_weight){
	    best_match = mn;
	    best_match_prev = prev;
	  }
	}
	//Move the best match to the head of the list
	if (best_match && best_match_prev){
	  best_match_prev->next = best_match->next;
	  best_match->next = s->match[j];	  
	  s->match[j] = best_match;
	}
      }
    }
  }
}

Segment* select_seed(Segment*** selected,  int nseed)
{
  Segment* best_segment = NULL;
  int i,j,k;
  for(i=0; i < nseq; i++){    
    SegmentNode* n;    
	int tmp_range;
    int sum_weight;
    for(n = all_segments[i]; n; n = n->next){
      Segment* s = n->segment;
      for (k=0; k < nseed; k++){
	if (!selected[k][i]) continue;
	if (selected[k][i] == s) break;	
	overlap(selected[k][i], s, &tmp_range);
	if (tmp_range < RANGE) break;
      }
      if (k < nseed) continue;
      sum_weight = 0;
      for(j=0; j < nseq; j++){
	      int match_weight;
		  MatchNode* best_match;
		  Segment* s2;
	if (i==j) continue;
	best_match= s->match[j];
	s2= NULL;	
	if (best_match) {
	  Match* match = best_match->match;
	  match_weight = match->neighborhood_weight;
	  if (match->s1 == s) s2 = match->s2;
	  else s2 = match->s1;
	  match_weight *= ((double)( s2->len - s2->selected_len))/ s2->len;
	}
	else{
	  match_weight = MIN_SCORE;
	}
	sum_weight += match_weight;	    
      }
      s->value = sum_weight;
      if (!s->selected && (best_segment == NULL || s->value > best_segment->value)) {	
	best_segment = s;	
      }
    }
  }
  for(j=0; j < nseq; j++){
	MatchNode* best_match;
    if (best_segment->seq_id==j) continue;
    best_match= best_segment->match[j];
    if (best_match) {
      Match* m = best_match->match;
      /*
      fprintf(stderr, "Match %s %d-%d: %s %d-%d value %d\n", 
	      all_seqnames[m->s1->seq_id], m->s1->start, m->s1->stop, 
	      all_seqnames[m->s2->seq_id], m->s2->start, m->s2->stop, 
	      m->neighborhood_weight);
      */
    }
  }
  return best_segment;
}

Segment** Extend_selected(Segment** selected )
{
  int k,j;
  Segment** selected_final = (Segment**)MallocOrDie(sizeof(Segment*) * nseq);  
  for(j=0; j < nseq; j++){
    if (selected[j]){
	  Segment* s;
      int start = selected[j]->start ;
      int stop =  selected[j]->stop ;
      SegmentNode* neighbor = selected[j]->neighbor;
      while(neighbor){
	if (neighbor->segment->start < start || neighbor->segment->stop > stop){
	  Match* found_match = NULL;
	  int skip_match = 0;
	  for(k=0; k < nseq; k++){
		SegmentNode* neighbor1;
	    if (j == k|| !selected[k]) continue;
	      neighbor1 = selected[k]->neighbor;
	      while(neighbor1){
		Match* m =  get_match(neighbor->segment, neighbor1->segment);
		if (m && m->type == DIRECT){
		  found_match = m;
		  break;
		}
		neighbor1 = neighbor1->next;
	      }
	  }
	  //Remove segments that are at the boundary of the selected seed (their matches are 
	  //outside range of the corresponding seed 
	  if (found_match){
		int included;
	    write_match(stdout, found_match);
	    included = 1;
	    for(k=0; k < nseq; k++){
	      int range;
		  MatchNode* mn;
	      if (!selected[k]) continue;
	      mn = neighbor->segment->match[k];		
	      while(mn &&  mn->match->weight > 200){
		Segment* s1 = mn->match->s1;
		if (s1->seq_id != k) s1 = mn->match->s2;
		  range = 0;
		  overlap(selected[k], s1, &range);
		  if (range > RANGE + 50 && range < 2 * RANGE){
		    included = 0;
		    break;
		  }
		  mn = mn->next;
	      }
	    }
	    if (included){
	      if (neighbor->segment->start < start){
		start =neighbor->segment->start;	
	      }
	      if (neighbor->segment->stop > stop){
		stop =neighbor->segment->stop;	
	      }
	    }
	  }
	}
	neighbor= neighbor->next;
      }
      if (stop - start < RANGE){
	int diff = RANGE - (stop - start +1);
	start -= diff/2;
	stop +=  diff/2;
	if (start < 1) start =1;
      }
      s=create_segment(selected[j]->seq_id, start, stop);
      selected_final[j] = s;
    }      
  }  
  return selected_final;
}


Segment** Extend_selected2(Segment** selected )
{
  int i,j,k;
  Segment* best_ext_seg= NULL;
  int best_ext_sum = 0;
  SegmentNode* candidates = NULL;
  SegmentNode* candidates_tail = NULL;
  Segment** selected_final = (Segment**)MallocOrDie(sizeof(Segment*) * nseq); 
  memset(selected_final, 0, sizeof(Segment*) * nseq);   
  for(j=0; j < nseq; j++){
    if (selected[j]){
      selected_final[j] = (Segment*) MallocOrDie(sizeof(Segment));      
      memcpy(selected_final[j], selected[j], sizeof(Segment));
    }    
  }

  for(i=0; i < nseq; i++){    
    SegmentNode* n;    
    for(n = all_segments[i]; n; n = n->next){
      n->segment->value = 0;
    }
  }

  //Put all segments that have match to the selected into the extension candidates list
  for(j=0; j < nseq; j++){
    if (selected_final[j] ){
      for(k=0; k < nseq; k++){
	MatchNode* mn;      
	for(mn= selected_final[j]->match[k]; mn; mn=mn->next){	  
	  double frac;
	  Segment* s = mn->match->s1;
	  int range=0;
	  if (s->seq_id == j) s = mn->match->s2;
	  if (s->value == -1) continue;
	  if (selected_final[s->seq_id]){
	    overlap(s, selected_final[s->seq_id], &range);
	    if (range <= s->stop - s->start + 1 || range > 2 * RANGE){
	      s->value = -1;
	      continue;	
	    }
	  }
	  if (s->value == 0) insert_segment_simple(s,&candidates, &candidates_tail);	
	  s->value += mn->match->weight;
	  frac = 1 - ( (double)s->selected_len) / s->len;	  
	  if (s->value * frac > best_ext_sum){
	    best_ext_sum = s->value * frac;
	    best_ext_seg = s;
	  }
	}
      }
    }
  }
  
  while(best_ext_seg){
    SegmentNode *n;
    int seq_id = best_ext_seg->seq_id;
    int range=0;
    if (selected_final[seq_id]){
      overlap(best_ext_seg, selected_final[seq_id], &range);      
      if (best_ext_seg->start < selected_final[seq_id]->start){
	selected_final[seq_id]->start = best_ext_seg->start;
      }	    
      if (best_ext_seg->stop > selected_final[seq_id]->stop){
	selected_final[seq_id]->stop = best_ext_seg->stop;
      }	
      selected_final[seq_id]->len = selected_final[seq_id]->stop - selected_final[seq_id]->start + 1;
    }
    else{
      selected_final[seq_id] = (Segment*) MallocOrDie(sizeof(Segment));      
      memcpy(selected_final[seq_id], best_ext_seg, sizeof(Segment));
    }
    //printf("seq %d  [%d-%d] range %d\n", seq_id, selected_final[seq_id]->start, selected_final[seq_id]->stop, range);
    //Update the candidate values
    best_ext_seg->value = -1;    
    for(j=0; j < nseq; j++){
      MatchNode* mn;      
      for(mn= best_ext_seg->match[j]; mn; mn=mn->next){
	Segment* s = mn->match->s1;
	int range;
	if (s->seq_id != j) s = mn->match->s2;
	if (s->value == -1)continue;
	if (s->value == 0) {	
	  if (selected_final[s->seq_id]){
	    overlap(s, selected_final[s->seq_id], &range);
	    if (range > 2 * RANGE) {
	      s->value = -1;
	      continue;
	    }	    
	  }
	  insert_segment_simple(s, &candidates, &candidates_tail);	
	}
	s->value += mn->match->weight;	
      }
    }
    //Select the new segment to extend
    best_ext_seg= NULL;
    best_ext_sum = 0;
    for(n=candidates; n; n=n->next){
      double frac;
	  Segment* s = n->segment;
      if (s->value == -1) continue;
      if (selected_final[s->seq_id]){
	int range;
	overlap(s, selected_final[s->seq_id], &range);      
	if ( range > 2 * RANGE || selected_final[s->seq_id]->stop - selected_final[s->seq_id]->start >=  RANGE){	
	  s->value = -1;
	  continue;
	}
      }
      frac = 1 - ( (double)s->selected_len) / s->len;	  
      if (s->value * frac > best_ext_sum){
	best_ext_sum = s->value * frac;
	best_ext_seg = s;
      }
    }
  }
  return selected_final;
}

int main(int argc, char* argv[])
{
  int i,j,k;
  char* match_file = argv[1];
  int   nseed = atof(argv[2]);  
  char* outfile = argv[3];
  FILE* fout;   
  int   trial = 0;
  
  Segment*** selected= (Segment***)MallocOrDie(sizeof(Segment**) * nseed);  
  Segment*** selected_final = (Segment***)MallocOrDie(sizeof(Segment**) * nseed);  
  
  read_match(match_file);
  for(i=0; i < nseed; i++){
    selected[i] = (Segment**)MallocOrDie(sizeof(Segment*) * nseq);
    selected_final[i] = (Segment**)MallocOrDie(sizeof(Segment*) * nseq);
    memset(selected[i], 0, sizeof(Segment*) * nseq);
    memset(selected_final[i], 0, sizeof(Segment*) * nseq);
  }
  fprintf(stderr, "Neigbhor Score\n");
  neighbor_score();  
  //write_match(stdout);  
  fprintf(stderr, "Find best_match\n");
  find_best_match();

  for(i=0; i < nseed && trial < nseed * 10 ; trial++){
	int total;
	double overlap_frac;
    char temp[100];
    Segment* seed = select_seed(selected_final, i);    
    if (seed == NULL) continue;
    seed->selected = 1;
    if (choose(seed, selected[i]) < MIN_SEQ) break;
    fprintf(stderr, "\nSeed %s:%d-%d value %d\n", all_seqnames[seed->seq_id], seed->start, seed->stop, seed->value);      
    selected_final[i] = Extend_selected2(selected[i]);

    //Check overlap
    total = 0;
    overlap_frac=0;
    for(j=0; j < nseq; j++){
      if (selected_final[i][j]){
	int max_olap = 0;
	for(k=0; k < i; k++){
	  if (selected_final[k][j]){
	    int olap= overlap(selected_final[i][j], selected_final[k][j], NULL);
	    if (olap > max_olap){
	      max_olap = olap;
	    }
	  }
	}
	overlap_frac += ((double) max_olap)/ selected_final[i][j]->len;	    
	total++;
      }
    }
    if (overlap_frac > total * 0.7 ) {
      continue;
    }
    sprintf(temp, "%s.%d",outfile, i+1);
    printf("Write %s\n", temp);
    if ( (fout = fopen(temp, "w") ) == NULL){
      printf("Fail to write file %s", temp);
    }
    for(j=0; j < nseq; j++){
      if (selected_final[i][j]){
	SegmentNode* n;    
	int start = selected_final[i][j]->start;
	int stop = selected_final[i][j]->stop;      
	start -= 100;
	stop +=  100;
	if (stop - start < RANGE){
	  int diff = RANGE - (stop - start +1);
	  start -= diff/2;
	  stop +=  diff/2;
	}
	if (start < 1) start =1;
	fprintf(fout, "%s %d %d\n", 
		all_seqnames[selected_final[i][j]->seq_id], start, stop);

	for(n = all_segments[j]; n; n = n->next){
	  Segment* s = n->segment;      
	  int tmp_olap=overlap(s, selected_final[i][j], NULL);
	  if (tmp_olap > s->selected_len){
	    s->selected_len = tmp_olap;
	  }
	}	
      }
    }      
    fclose(fout);
    i++;
  }  
}

