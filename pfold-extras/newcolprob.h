#include "grammar.h"
#include <stdio.h>
#include "align.h"
#include "phyl.h"
#include "matrix.h"
#include "col.h"
#include "llist.h"
#include "edouble.h"
#include "search.h"

#ifndef __colprob_h__   /* Only define the following once */
#define __colprob_h__

#define SINGLE -2
#define DOUBLE -3

typedef struct tagAlignstr {
  Align *align;
  int *str;
  int *map;
} Alignstr;

typedef struct tagColprob
{
  int wlen;
  Grammar *grammar;
  Edouble **sgl;
  Edouble ***dbl;
} Colprob;

typedef struct tagFColprob
{
  int wlen;
  Grammar *grammar;
  Fdouble **sgl;
  Fdouble ***dbl;
} FColprob;

Colprob *MakeColprob(Grammar *grammar, int wlen);
void FreeColprob(Grammar *grammar, Colprob *colprob);
void InitColprob(Colprob *colprob, Alignstr *alignstr);
void MoveColprob(Colprob *colprob, Alignstr *alignstr, int pos);

FColprob *MakeFColprob(Grammar *grammar, int wlen);
void FreeFColprob(Grammar *grammar, FColprob *colprob);
void InitFColprob(FColprob *colprob, Alignstr *alignstr);
void MoveFColprob(FColprob *colprob, Alignstr *alignstr, int pos);

void InitCol(Phyl *phyl, Grammar *grammar, int wlen,
	     Alignstr *alignstr, double logp_limit);
void FinishCol(Phyl *phyl, Grammar *grammar);

Alignstr *Removegap(Grammar *grammar, Alignstr *alignstr, double limit);

Edouble ScolProb(Alignstr *alignstr, int sgrpno, int pos);
Edouble DcolProb(Alignstr *alignstr, int dgrpno, int pos1, int dist);

Fdouble FScolProb(Alignstr *alignstr, int sgrpno, int pos);
Fdouble FDcolProb(Alignstr *alignstr, int dgrpno, int pos1, int dist);

#endif
