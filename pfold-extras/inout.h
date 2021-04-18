#include "grammar.h"
#include "newcolprob.h"
#include "llist.h"
#include "edouble.h"

/* Inside variables are defined as a pointer to a list (an entry for
   each nonterminal) of pointers to a list of pointers to
   probabilities:

Inside->prob
      |
      |
      v
  +--------+
  | No 0   |
  +--------+        +--------+
  | No 1   | -----> | No 0   | -----> List of w   Probs (for seqlen 0)
  +--------+        +--------+
  .        .        | No 1   | -----> List of w-1 Probs (for seqlen 1)
  .        .        +--------+
  +--------+        .        .
  | No n-1 |        .        .
  +--------+        +--------+
                    | No w-1 | -----> List of 1   Prob  (for seqlen w-1)
 n = number of      +--------+
   nonterminals
                         w = window length
*/

/* wlen is window length. An Inside structure belongs to a specific
   grammar. */

typedef struct tagInside
{
  int wlen;
  Grammar *grammar;
  Edouble ***prob;
} Inside;

typedef struct tagOutside
{
  int wlen;
  Grammar *grammar;
  Edouble ***prob;
} Outside;

typedef struct tagPairprob
{
  int wlen;
  Grammar *grammar;
  Edouble **dbl;
  Edouble *sgl;
} Pairprob;

typedef struct tagPairaffinity
{
  int wlen;
  Grammar *grammar;
  Edouble **dbl;
  Edouble *sgl;
} Pairaffinity;

typedef struct tagTrace
{
  int type;       /* Rule type: 0: ruls, 1: ruln, 2: rulsn,
		     3: rulns, 4: ruldnd, 5: rulnn 6: ruldd */
  int wlen;        /* Wlen for first split in rulnn */
  void *rule;     /* Pointer to rule */
} Trace;

typedef struct tagCyk
{
  int wlen;
  Grammar *grammar;
  Edouble ***prob;
  Trace ***trace;
  int *structure;
  int numpair;
} Cyk;

typedef struct tagFCyk
{
  int wlen;
  Grammar *grammar;
  Fdouble ***prob;
  Trace ***trace;
  int *structure;
  int numpair;
} FCyk;

typedef struct tagPpCyk
{
  int wlen;
  Grammar *grammar;
  Edouble **prob;
  Trace **trace;
  int *structure;
} PpCyk;

typedef struct tagPaircnt
{
  int wlen;
  Grammar *grammar;
  Edouble ***cnt;
} Paircnt;

Inside *MakeInside(Grammar *grammar, int wlen);
void FreeInside(Grammar *grammar, Inside *inside);
void InitInside(Inside *e, Colprob *colprob);
void MoveInside(Inside *e, Colprob *colprob, int pos);

Outside *MakeOutside(Grammar *grammar, int wlen);
void InitOutside(Inside *e, Outside *f, Colprob *colprob);

Pairprob *MakePairprob(int wlen);
void FreePairprob(Pairprob *pairprob);
void InitPairprob(Pairprob *pp, Inside *e, Outside *f, Colprob *colprob);
void PpmPairprob(FILE *fp, Pairprob *pp);
void PrintPairprob(FILE *fp, Pairprob *pp, Pairprob *pp2, Alignstr *alignstr);
void PsPairprob(FILE *fp, Pairprob *pp, int *map);

void InitPairaff(Pairprob *pp, Inside *e, Outside *f, Colprob *colprob, double prior);

Cyk *MakeCyk(Grammar *grammar, int wlen);
void InitCyk(Cyk *e, Colprob *colprob);
void MoveCyk(Cyk *e, Colprob *colprob, int pos);
void ModifyCyk(Cyk *e, Pairprob *pp);
void PrintCyk(FILE *fp, Cyk *e);
int *CykStr(Cyk *cyk);

void AddStr(Entry *entry, Alignstr *alignstr, Pairprob *pp);

FCyk *MakeFCyk(Grammar *grammar, int wlen);
void InitFCyk(FCyk *e, FColprob *fcolprob);
void MoveFCyk(FCyk *e, FColprob *fcolprob, int pos);
int *FCykStr(FCyk *cyk);

PpCyk *MakePpCyk(int wlen);
void InitPpCyk(PpCyk *e, Pairprob *pairprob);
int *PpCykStr(PpCyk *cyk);
void PrintPpCyk(FILE *fp, PpCyk *e);

Paircnt *MakePaircnt(Grammar *grammar, int wlen);
void FreePaircnt(Grammar *grammar, Paircnt *paircnt);
void InitPaircnt(Inside *e, Paircnt *paircnt, Colprob *colprob);
void MovePaircnt(Inside *e, Paircnt *paircnt, Colprob *colprob, int pos);
