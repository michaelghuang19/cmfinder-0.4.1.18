#define MAX_CAND_LENGTH 500
#define MAX_CAND_FILE_LINE_LENGTH 1000

extern int maxNumCand;
extern int maxCandSpan;
extern int minCandSpan;

typedef struct {
  int seq_id;
  int cand_id;
  int start;
  int stop;
  int len;
  double evalue; /* -1 if not set */
  double score;
  double weight;
  double energy;
  char seq[MAX_CAND_LENGTH];
  char ss[MAX_CAND_LENGTH];
  CM_ALIDISPLAY *ad;
} Cand;

void InitCand (Cand ***ret_cand,int **ret_ncand,int nseq);
void DeleteCand(Cand **cand,int *ncand,int nseq);
void DeleteOneCand (Cand *cand);
int GetTotalNumCand (int nseq,int *ncand,Cand **cand);
int GetLongestCandLen (int nseq,int *ncand,Cand **cand);

/* I/O */
void Write2DCand(char* cand_file, int nseq, Cand** cand, int* ncand);
void Write1DCand(char* cand_file, Cand** cand, int ncand);
Cand** Read2DCand(char* cand_file,  int nseq, int ** ret_ncand, int * ret_max_cand);
Cand* Read1DCand(char* cand_file,  int * ret_ncand);
char* ExpandFull(char* anno, char* seq);
int Overlap(Cand* c1, Cand* c2, int* olap_min, int* olap_max, int* min, int* max);

/* Transform to Squid format */
/* SQINFO* Cand2Sqinfo(Cand** cand, int ncand, SQINFO* seq_sqinfo); */

int compDouble(const void* a, const void* b);
int CompCandByScore(const void* a, const void* b);
int CompCandByEnergy(const void* a, const void* b);
Cand** SortCand(Cand* cand, int ncand, int (*f) (const void* a, const void* b));
int* GetPairtable(char* ss);
int isHairpin(char* ss);
int isMultiloop(char* ss);
int countHairpin(char* ss);
