#include "inout.h"

typedef struct tagBacktrack
{
  int nont;
  int pos;
  int wlen;
} Backtrack;

/*------------------------------------------------------------------*/

/* Allocates memory for Inside variable for a window of length
   wlen. */

Inside *MakeInside(Grammar *grammar, int wlen)
{
  Inside* inside = (Inside *)malloc(sizeof(Inside));
  int n = grammar->nont->count;
  int i, j;

  inside->wlen = wlen;
  inside->grammar = grammar;
  inside->prob = (Edouble ***)malloc(n * sizeof(Edouble **));

  for (j = 0; j < n; j++)
    {
      inside->prob[j] = (Edouble **)malloc(wlen * sizeof(Edouble *));
      for (i = 0; i < wlen; i++)
	inside->prob[j][i] = (Edouble *)malloc( (wlen - i) * sizeof(Edouble));
    }

  return inside;
}

/*------------------------------------------------------------------*/

void FreeInside(Grammar *grammar, Inside *inside)
{
  int n = grammar->nont->count;
  int i, j;
  int wlen;

  wlen = inside->wlen;

  for (j = 0; j < n; j++) {
    for (i = 0; i < wlen; i++)
      free(inside->prob[j][i]);
    free(inside->prob[j]);
  }

  free(inside->prob);
  free(inside);
}

/*------------------------------------------------------------------*/

/* This initializes the inside variables according to the
 probabilities in colprob. */

void InitInside(Inside *e, Colprob *colprob)
{
  int n;   /* nonterminal iterator */
  int p;   /* position iterator */
  int w;   /* string length iterator */
  int w2;   /* string length iterator for left string in N -> NN */
  int wl = e->wlen;
  LListCounter *lcount = MakeCounter(e->grammar->ruls, FIRST);
  Ruls *ruls;
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Ruldd *ruldd;
  Rulnn *rulnn;
  double a, b, c;

  /* Set Inside variables to zero */

  for (n = 0; n < e->grammar->nont->count; n++)
    for (p = 0; p < wl; p++)
      for (w = 0; w + p < wl; w++)
	e->prob[n][p][w] = Dbl2Edbl(0.);

  /* stringlength 0 initialization */

  InitCounter(lcount, e->grammar->ruls, FIRST);  /* N -> s rules */
  while((ruls = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      AddEdouble(&e->prob[ruls->lnt][p][0],
	      ProdEdouble(ruls->prob,
		      colprob->sgl[ruls->rsg][p]));

  /* ---------- NOTICE ------------- */
  /* ERROR: SORT Ruln, CORRECT LATER */
  /* ---------- NOTICE ------------- */

  InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
  while((ruln = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      AddEdouble(&e->prob[ruln->lnt][p][0],
	      ProdEdouble(ruln->prob,
		      e->prob[ruln->rnt][p][0]));

  /* Use of the rest of the rules */

  for (w = 1; w < wl; w++) {

    InitCounter(lcount, e->grammar->rulsn, FIRST);  /* N -> sN rules */
    while((rulsn = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
      AddEdouble(&e->prob[rulsn->lnt][p][w],
	      ProdEdouble(rulsn->prob,
	      ProdEdouble(colprob->sgl[rulsn->rsg][p],
			 e->prob[rulsn->rnt][p+1][w-1])));
	
    InitCounter(lcount, e->grammar->rulns, FIRST);  /* N -> Ns rules */
    while((rulns = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
      AddEdouble(&e->prob[rulns->lnt][p][w],
	      ProdEdouble(rulns->prob,
              ProdEdouble(colprob->sgl[rulns->rsg][p+w],
			 e->prob[rulns->rnt][p][w-1])));
	
    if (w >= 2) {  /* minimum stringlength of 3 */
      InitCounter(lcount, e->grammar->ruldnd, FIRST);  /* N -> dNd rules */
      while((ruldnd = Next(lcount)) != NULL)
	for (p = 0; w + p < wl; p++)
	  AddEdouble(&e->prob[ruldnd->lnt][p][w],
		  ProdEdouble(ruldnd->prob,
		  ProdEdouble(colprob->dbl[ruldnd->rdb][p][w],
			  e->prob[ruldnd->rnt][p+1][w-2])));
    }

    if (w == 1) {  /* stringlength of 2 */
      InitCounter(lcount, e->grammar->ruldd, FIRST);  /* N -> dd rules */
      while((ruldd = Next(lcount)) != NULL)
	for (p = 0; w + p < wl; p++)
	  AddEdouble(&e->prob[ruldd->lnt][p][w],
		     ProdEdouble(ruldd->prob, colprob->dbl[ruldd->rdb][p][w]));
    }

    InitCounter(lcount, e->grammar->rulnn, FIRST);  /* N -> NN rules */
    while((rulnn = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	for (w2 = 0; w2 < w; w2++)
	  AddEdouble(&e->prob[rulnn->lnt][p][w],
		  ProdEdouble(rulnn->prob,
		  ProdEdouble(e->prob[rulnn->rnt1][p][w2],
			  e->prob[rulnn->rnt2][p+w2+1][w-w2-1])));

    InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
    while((ruln = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
      AddEdouble(&e->prob[ruln->lnt][p][w],
	      ProdEdouble(ruln->prob,
		      e->prob[ruln->rnt][p][w]));

    if (w % 10 == 0) {
      a = b = -Edbl2Dbl(LogEdouble(
		   e->prob[FindNont('S', e->grammar->nont)][0][w]));
      for (p = 1; w + p < wl; p++) {
	c = -Edbl2Dbl(LogEdouble(
		 e->prob[FindNont('S', e->grammar->nont)][p][w]));
	if (a > c)
	  a = c;
	if (b < c)
	  b = c;
      }
      //      printf(" %5d  %8.4f  %8.4f\n", w, a, b);
    }
  }

  free(lcount);
}
      
/*------------------------------------------------------------------*/

/* This moves the inside variables according to the probabilities in
   colprob. */

void MoveInside(Inside *e, Colprob *colprob, int pos)
{
  int n;   /* nonterminal iterator */
  int w;   /* string length iterator */
  int w2;   /* string length iterator for left string in N -> NN */
  LListCounter *lcount = MakeCounter(e->grammar->ruls, FIRST);
  int wl = e->wlen;
  Ruls *ruls;
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Ruldd *ruldd;
  Rulnn *rulnn;

  /* Set Inside variables to zero */

  for (n = 0; n < e->grammar->nont->count; n++)
    for (w = 0; w < wl; w++)
      e->prob[n][(pos-w)%(wl-w)][w] = Dbl2Edbl(0.);

  /* stringlength 0 initialization */

  InitCounter(lcount, e->grammar->ruls, FIRST);  /* N -> s rules */
  while((ruls = Next(lcount)) != NULL)
    AddEdouble(&e->prob[ruls->lnt][pos % wl][0],
	    ProdEdouble(ruls->prob,
		    colprob->sgl[ruls->rsg][pos % wl]));

  /* ---------- NOTICE ------------- */
  /* ERROR: SORT Ruln, CORRECT LATER */
  /* ---------- NOTICE ------------- */

  InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
  while((ruln = Next(lcount)) != NULL)
    AddEdouble(&e->prob[ruln->lnt][pos % wl][0],
	    ProdEdouble(ruln->prob,
		    e->prob[ruln->rnt][pos % wl][0]));

  /* Use of the rest of the rules */

  for (w = 1; w < wl; w++) {

    InitCounter(lcount, e->grammar->rulsn, FIRST);  /* N -> sN rules */
    while((rulsn = Next(lcount)) != NULL)
                                   /* below: get positive remainder */
      AddEdouble(&e->prob[rulsn->lnt][(pos-w)%(wl-w)][w],
	      ProdEdouble(rulsn->prob,
	 ProdEdouble(colprob->sgl[rulsn->rsg][(pos-w) % wl],
		      e->prob[rulsn->rnt][(pos-w+1)%(wl-(w-1))][w-1])));

    InitCounter(lcount, e->grammar->rulns, FIRST);  /* N -> Ns rules */
    while((rulns = Next(lcount)) != NULL)
      AddEdouble(&e->prob[rulns->lnt][(pos-w)%(wl-w)][w],
	      ProdEdouble(rulns->prob,
		 ProdEdouble(colprob->sgl[rulns->rsg][pos % wl],
		      e->prob[rulns->rnt]
		      [(pos-w) % (wl-(w-1))][w-1])));

    if (w >= 2) {  /* minimum stringlength of 3 */
      InitCounter(lcount, e->grammar->ruldnd, FIRST);  /* N -> dNd rules */
      while((ruldnd = Next(lcount)) != NULL)
	  AddEdouble(&e->prob[ruldnd->lnt][(pos-w)%(wl-w)][w],
		  ProdEdouble(ruldnd->prob,
              ProdEdouble(colprob->dbl[ruldnd->rdb][(pos-w)%(wl-w)][w],
			  e->prob[ruldnd->rnt][(pos-w+1)%(wl-(w-2))][w-2])));
    }

    if (w == 2) {  /* stringlength of 2 */
      InitCounter(lcount, e->grammar->ruldd, FIRST);  /* N -> dd rules */
      while((ruldd = Next(lcount)) != NULL)
	  AddEdouble(&e->prob[ruldd->lnt][(pos-w)%(wl-w)][w],
		  ProdEdouble(ruldd->prob,
			      colprob->dbl[ruldd->rdb][(pos-w)%(wl-w)][w]));
    }

    InitCounter(lcount, e->grammar->rulnn, FIRST);  /* N -> NN rules */
    while((rulnn = Next(lcount)) != NULL)
	for (w2 = 0; w2 < w; w2++)
	  AddEdouble(&e->prob[rulnn->lnt][(pos-w)%(wl-w)][w],
		  ProdEdouble(rulnn->prob,
		  ProdEdouble(e->prob[rulnn->rnt1][(pos-w)%(wl-w2)][w2],
			  e->prob[rulnn->rnt2]\
			  [(pos-w+w2+1)%(wl-(w-w2-1))][w-w2-1])));

    InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
    while((ruln = Next(lcount)) != NULL)
      AddEdouble(&e->prob[ruln->lnt][(pos-w)%(wl-w)][w],
	      ProdEdouble(ruln->prob,
		      e->prob[ruln->rnt][(pos-w)%(wl-w)][w]));
  }

  free(lcount);
}
      
/*------------------------------------------------------------------*/

/* Allocates memory for Outside variables for a window of length
   wlen. */

Outside *MakeOutside(Grammar *grammar, int wlen)
{
  Outside* outside = (Outside *)malloc(sizeof(Outside));
  int n = grammar->nont->count;
  int i, j;

  outside->wlen = wlen;
  outside->grammar = grammar;
  outside->prob = (Edouble ***)malloc(n * sizeof(Edouble **));

  for (j = 0; j < n; j++)
    {
      outside->prob[j] = (Edouble **)malloc(wlen * sizeof(Edouble *));
      for (i = 0; i < wlen; i++)
	outside->prob[j][i] = (Edouble *)malloc( (wlen - i) * sizeof(Edouble));
    }

  return outside;
}

/*------------------------------------------------------------------*/

/* This initializes the outside variables according to the
 probabilities in colprob. */

void InitOutside(Inside *e, Outside *f, Colprob *colprob)
{
  int n;   /* nonterminal iterator */
  int p;   /* position iterator */
  int w;   /* string length iterator */
  int w2;   /* string length iterator for left string in N -> NN */
  int wl = f->wlen;
  int s = FindNont('S', f->grammar->nont);  /* the number of 'S' */
  LListCounter *lcount = MakeCounter(f->grammar->ruls, FIRST);
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Rulnn *rulnn;

  /* Set Outside variables to zero */

  for (n = 0; n < f->grammar->nont->count; n++)
    for (p = 0; p < wl; p++)
      for (w = 0; w + p < wl; w++)
	f->prob[n][p][w] = Dbl2Edbl(0.);

  /* stringlength wl initialization */

  f->prob[s][0][wl-1] = Dbl2Edbl(1.);  /* S -> -> S: probability 1 */

  /* ---------- NOTICE ------------- */
  /* ERROR: SORT Ruln, CORRECT LATER */
  /* ---------- NOTICE ------------- */

  InitCounter(lcount, f->grammar->ruln, LAST);  /* N -> N rules */
  while((ruln = Prev(lcount)) != NULL)
    AddEdouble(&f->prob[ruln->rnt][0][wl-1],
	    ProdEdouble(ruln->prob,
		    f->prob[ruln->lnt][0][wl-1]));

  /* Use of the rest of the rules */

  for (w = wl-2; w >= 0; w--) {

    InitCounter(lcount, f->grammar->rulsn, FIRST);  /* N -> sN rules */
    while((rulsn = Next(lcount)) != NULL)
      for (p = 1; w + p < wl; p++)
      AddEdouble(&f->prob[rulsn->rnt][p][w],
	      ProdEdouble(rulsn->prob,
	      ProdEdouble(colprob->sgl[rulsn->rsg][p-1],
		      f->prob[rulsn->lnt][p-1][w+1])));
	
    InitCounter(lcount, f->grammar->rulns, FIRST);  /* N -> Ns rules */
    while((rulns = Next(lcount)) != NULL)
      for (p = 0; w + p+1 < wl; p++)
      AddEdouble(&f->prob[rulns->rnt][p][w],
	      ProdEdouble(rulns->prob,
	      ProdEdouble(colprob->sgl[rulns->rsg][p+w+1],
		      f->prob[rulns->lnt][p][w+1])));
	
    InitCounter(lcount, f->grammar->ruldnd, FIRST);  /* N -> dNd rules */
    while((ruldnd = Next(lcount)) != NULL)
      for (p = 1; w + p+1 < wl; p++) {
	/* printf("%d %d %d\n", ruldnd->lnt, ruldnd->rdb, ruldnd->rnt);*/
	AddEdouble(&f->prob[ruldnd->rnt][p][w],
		ProdEdouble(ruldnd->prob,
			    ProdEdouble(colprob->dbl[ruldnd->rdb][p-1][w+2],
				f->prob[ruldnd->lnt][p-1][w+2])));
      }

    InitCounter(lcount, f->grammar->rulnn, FIRST);  /* N -> NN rules */
    while((rulnn = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++) {
	for (w2 = 0; w2+p+w+1 < wl; w2++)
	  AddEdouble(&f->prob[rulnn->rnt1][p][w],
		  ProdEdouble(rulnn->prob,
		  ProdEdouble(f->prob[rulnn->lnt][p][w+w2+1],
			  e->prob[rulnn->rnt2][p+w+1][w2])));
	for (w2 = 0; w2 < p; w2++)
	  AddEdouble(&f->prob[rulnn->rnt2][p][w],
		  ProdEdouble(rulnn->prob,
		  ProdEdouble(f->prob[rulnn->lnt][p-1-w2][w+w2+1],
			  e->prob[rulnn->rnt1][p-1-w2][w2])));
      }

    InitCounter(lcount, f->grammar->ruln, FIRST);  /* N -> N rules */
    while((ruln = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
      AddEdouble(&f->prob[ruln->rnt][p][w],
	      ProdEdouble(ruln->prob,
		      f->prob[ruln->lnt][p][w]));
  }

  free(lcount);
}
      
/*------------------------------------------------------------------*/

Pairprob *MakePairprob(int wlen)
{
  Pairprob* pairprob = (Pairprob *)malloc(sizeof(Pairprob));
  int i;

  pairprob->wlen = wlen;
  pairprob->dbl = (Edouble **)malloc(wlen * sizeof(Edouble *));

  for (i = 0; i < wlen; i++)
    pairprob->dbl[i] = (Edouble *)malloc( (wlen - i) * sizeof(Edouble));

  pairprob->sgl = (Edouble *)malloc(wlen * sizeof(Edouble));

  return pairprob;
}

/*------------------------------------------------------------------*/

void FreePairprob(Pairprob *pairprob)
{
  int i;
  int wlen;

  wlen = pairprob->wlen;

  free(pairprob->sgl);

  for (i = 0; i < wlen; i++)
    free(pairprob->dbl[i]);

  free(pairprob->dbl);
  free(pairprob);
}

/*------------------------------------------------------------------*/

void InitPairprob(Pairprob *pp, Inside *e, Outside *f, Colprob *colprob)
{
  int p, w;
  int wl = pp->wlen;
  LListCounter *lcount = MakeCounter(e->grammar->ruldnd, FIRST);
  Ruldnd *ruldnd;
  Ruldd *ruldd;
  Ruls *ruls;
  Rulns *rulns;
  Rulsn *rulsn;
  int s = FindNont('S', f->grammar->nont);  /* the number of 'S' */
  Edouble totprob = e->prob[s][0][wl-1];

  for (w = 0; w < wl; w++)
    for (p = 0; p + w < wl; p++)
      pp->dbl[p][w] = Dbl2Edbl(0.);

  while ((ruldnd = Next(lcount)) != NULL)
    for (w = 2; w < wl; w++)
      for (p = 0; p + w < wl; p++)
	AddEdouble(&pp->dbl[p][w], ProdEdouble(ruldnd->prob,
		 ProdEdouble(f->prob[ruldnd->lnt][p][w],
		 ProdEdouble(colprob->dbl[ruldnd->rdb][p][w],
			     e->prob[ruldnd->rnt][p+1][w-2]))));

  lcount = MakeCounter(e->grammar->ruldd, FIRST);
  while ((ruldd = Next(lcount)) != NULL) {
    w = 1;
    for (p = 0; p + w < wl; p++)
      AddEdouble(&pp->dbl[p][w], ProdEdouble(ruldd->prob,
		 ProdEdouble(f->prob[ruldd->lnt][p][w],
			     colprob->dbl[ruldd->rdb][p][w])));
  }

  for (w = 0; w < wl; w++)
    for (p = 0; p + w < wl; p++)
      DivEdouble(&pp->dbl[p][w], totprob);


  for (p = 0; p < wl; p++)
    pp->sgl[p] = Dbl2Edbl(0.);

  lcount = MakeCounter(e->grammar->ruls, FIRST);
  while ((ruls = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      AddEdouble(&pp->sgl[p], ProdEdouble(ruls->prob,
		 ProdEdouble(colprob->sgl[ruls->rsg][p],
			     f->prob[ruls->lnt][p][0])));

  lcount = MakeCounter(e->grammar->rulns, FIRST);
  while ((rulns = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      for (w = 1; w <= p; w++)
	AddEdouble(&pp->sgl[p], ProdEdouble(rulns->prob,
		 ProdEdouble(f->prob[rulns->lnt][p-w][w],
		 ProdEdouble(colprob->sgl[rulns->rsg][p],
			     e->prob[rulns->rnt][p-w][w-1]))));

  lcount = MakeCounter(e->grammar->rulsn, FIRST);
  while ((rulsn = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      for (w = 1; w < wl-p; w++)
	AddEdouble(&pp->sgl[p], ProdEdouble(rulsn->prob,
		 ProdEdouble(f->prob[rulsn->lnt][p][w],
		 ProdEdouble(colprob->sgl[rulsn->rsg][p],
			     e->prob[rulsn->rnt][p+1][w-1]))));

  for (p = 0; p < wl; p++)
    DivEdouble(&pp->sgl[p], totprob);

  free(lcount);
}

/*------------------------------------------------------------------*/

void PrintPairprob(FILE *fp, Pairprob *pp, Pairprob *pp2, Alignstr *alignstr)
{
  int i, j, k, p, w;
  int wl = pp->wlen;
  int *map = alignstr->map;
  int len = map[pp->wlen];

  fprintf(fp, "%d\n", len);

  for (k = 0; k < map[0]-1; k++) {
    for (j = 0; j < len; j++)
      fprintf(fp," 0.000000");
    fprintf(fp, "\n");
  }

  for (i = 0; i < wl; i++) {
    for (k = 0; k < map[0]-1; k++)
      fprintf(fp," 0.000000");

    for (j = 0; j < wl; j++) {
      w = (i > j) ? i - j : j - i;
      p = (i > j) ? j : i;
      if (i > j)
	fprintf(fp," %8.6f", Edbl2Dbl(pp->dbl[p][w]));
      else
	fprintf(fp," %8.6f", Edbl2Dbl(pp2->dbl[p][w]));

      if (j == wl-1)
	for (k = 0; k < map[j+1]-map[j]; k++)
	  fprintf(fp," 0.000000");
      else
	for (k = 0; k < map[j+1]-map[j]-1; k++)
	  fprintf(fp," 0.000000");
    }

    fprintf(fp, "\n");

    if (i == wl-1)
      for (k = 0; k < map[i+1]-map[i]; k++) {
	for (j = 0; j < len; j++)
	  fprintf(fp," 0.000000");
	fprintf(fp, "\n");
      }
    else
      for (k = 0; k < map[i+1]-map[i]-1; k++) {
	for (j = 0; j < len; j++)
	  fprintf(fp," 0.000000");
	fprintf(fp, "\n");
      }
  }
  fprintf(fp, "\n");

  for (k = 0; k < map[0]-1; k++)
    fprintf(fp," 0.000000");

  for (p = 0; p < wl; p++) {
    fprintf(fp," %8.6f", Edbl2Dbl(pp->sgl[p]));
    
    if (p == wl-1)
      for (k = 0; k < map[p+1]-map[p]; k++)
	fprintf(fp," 0.000000");
    else
      for (k = 0; k < map[p+1]-map[p]-1; k++)
	fprintf(fp," 0.000000");
  }
  fprintf(fp, "\n");
}

/*------------------------------------------------------------------*/

void PsPairprob(FILE *fp, Pairprob *pp, int *map)
{
  int i, j, p, w;
  int wl;
  int size;
  double scale;

  wl = pp->wlen;

  size = map[wl];

  scale = 300./size;

  fprintf(fp, "%%!PS-Adobe-2.0\n");
  fprintf(fp, "%%%%BoundingBox: %f %f %f %f\n", -2*scale-5, -2*scale-5, size*scale+5, size*scale+5);
  fprintf(fp, "%%!\n");
  fprintf(fp, "\n");
  fprintf(fp, "%f %f scale\n", scale, scale);
  fprintf(fp, "\n");

  fprintf(fp, "/b {\n");
  fprintf(fp, "  newpath\n");
  fprintf(fp, "  3 1 roll\n");
  fprintf(fp, "  moveto\n");
  fprintf(fp, "  sqrt 0.9 mul dup dup dup\n");
  fprintf(fp, "  2 div 0.5 exch sub dup rmoveto\n");
  fprintf(fp, "  0 exch rlineto\n");
  fprintf(fp, "  0 rlineto\n");
  fprintf(fp, "  0 exch neg rlineto\n");
  fprintf(fp, "  fill\n");
  fprintf(fp, "  stroke\n");
  fprintf(fp, "} def\n");
  /*
  fprintf(fp, "/b {\n");
  fprintf(fp, "  newpath\n");
  fprintf(fp, "  1 exch sub dup dup setrgbcolor  \n");
  fprintf(fp, "  moveto\n");
  fprintf(fp, "  0.05 0.05 rmoveto\n");
  fprintf(fp, "  0 0.9 rlineto\n");
  fprintf(fp, "  0.9 0 rlineto\n");
  fprintf(fp, "  0 -0.9 rlineto\n");
  fprintf(fp, "  fill\n");
  fprintf(fp, "  stroke\n");
  fprintf(fp, "} def\n");
  */
  fprintf(fp, "\n");
  fprintf(fp, "0.1 setlinewidth\n");
  fprintf(fp, "\n");
  fprintf(fp, "newpath\n");
  fprintf(fp, "0 0 moveto\n");
  fprintf(fp, "%d 0 rlineto\n", size);
  fprintf(fp, "0 %d rlineto\n", size);
  fprintf(fp, "-%d 0 rlineto\n", size);
  fprintf(fp, "closepath\n");
  fprintf(fp, "stroke\n");
  fprintf(fp, "\n");
  fprintf(fp, "newpath\n");
  fprintf(fp, "0 0 moveto\n");
  fprintf(fp, "%d %d rlineto\n", size, size);

  for (i = 10; i < size; i+= 10) {
    fprintf(fp, "0 %d moveto\n", i);
    fprintf(fp, "%d 0 rlineto\n", size);
    fprintf(fp, "%d 0 moveto\n", i);
    fprintf(fp, "0 %d rlineto\n", size);
  }

  fprintf(fp, "stroke\n");
  fprintf(fp, "\n");

  for (i = 0; i < wl; i++) {
    for (j = 0; j < wl; j++) {
      w = (i > j) ? i - j : j - i;
      p = (i > j) ? j : i;
      if (Edbl2Dbl(pp->dbl[p][w]) > 0.02)
	fprintf(fp, "%d %d %8.6f b\n", map[i], map[j], Edbl2Dbl(pp->dbl[p][w]));
    }
  }

  fprintf(fp, "\n");

  for (p = 0; p < wl; p++) {
    if (Edbl2Dbl(pp->sgl[p]) > 0.02) {
      fprintf(fp, "-2 %d %8.6f b\n", map[p], Edbl2Dbl(pp->sgl[p]));
      fprintf(fp, "%d -2 %8.6f b\n", map[p], Edbl2Dbl(pp->sgl[p]));
    }
  }

  fprintf(fp, "\nshowpage\n");

}

/*------------------------------------------------------------------*/

void PpmPairprob(FILE *fp, Pairprob *pp)
{
  int i, j, p, w, k, l;
  int wl = pp->wlen;

  fprintf(fp, "P3\n");
  fprintf(fp, "%d %d\n", wl*8, wl*8+8);
  fprintf(fp, "255\n");

  for (i = 0; i < wl; i++) {
    for (j = 0; j < wl; j++) {
      if (j % 10 == 0 || i % 10 == 0)
	fprintf(fp,"    0   0   0");
      else
	fprintf(fp,"  255 255 255");
      for (k = 0; k < 7; k++) {
	if (i % 10 == 0)
	  fprintf(fp,"    0   0   0");
	else
	  fprintf(fp,"  255 255 255");
      }
    }
    fprintf(fp, "\n");
    for (l = 0; l < 7; l++) {
      for (j = 0; j < wl; j++) {
	w = (i > j) ? i - j : j - i;
	p = (i > j) ? j : i;
	if (j % 10 == 0)
	  fprintf(fp,"    0   0   0");
	else
	  fprintf(fp,"  255 255 255");
	for (k = 0; k < 7; k++) {
	  fprintf(fp,"  %3d", (int)(255*(1-Edbl2Dbl(pp->dbl[p][w]))));
	  fprintf(fp," %3d", (int)(255*(1-Edbl2Dbl(pp->dbl[p][w]))));
	  fprintf(fp," %3d", (int)(255*(1-Edbl2Dbl(pp->dbl[p][w]))));
	}
      }
      fprintf(fp, "\n");
    }
  }

  for (l = 0; l < 8; l++) {
    for (p = 0; p < wl; p++) {
      for (k = 0; k < 8; k++) {
	fprintf(fp,"  %3d", (int)(255*(1-Edbl2Dbl(pp->sgl[p]))));
	fprintf(fp," %3d", (int)(255*(1-Edbl2Dbl(pp->sgl[p]))));
	fprintf(fp," %3d", (int)(255*(1-Edbl2Dbl(pp->sgl[p]))));
      }
    }
    fprintf(fp, "\n");
  }
}

/*------------------------------------------------------------------*/

void InitPairaff(Pairprob *pp, Inside *e, Outside *f, Colprob *colprob, double prior)
{
  int p, w, g;
  int wl = pp->wlen;
  Edouble prob;

  for (w = 0; w < wl; w++)
    for (p = 0; p + w < wl; p++)
      pp->dbl[p][w] = Dbl2Edbl(0.);

  for (p = 0; p < wl; p++)
    pp->sgl[p] = Dbl2Edbl(0.);

  for (w = 0; w < wl; w++)
    for (p = 0; p + w < wl; p++) {
      prob = Dbl2Edbl(0.);
      for (g = 0; g < e->grammar->dgrp->count; g++) {
	AddEdouble(&prob, ProdEdouble(Dbl2Edbl(prior), colprob->dbl[g][p][w]));
	AddEdouble(&pp->dbl[p][w], ProdEdouble(Dbl2Edbl(prior), colprob->dbl[g][p][w]));
      }
      for (g = 0; g < e->grammar->sgrp->count; g++)
	AddEdouble(&prob, ProdEdouble(Dbl2Edbl(1-prior), ProdEdouble(colprob->sgl[g][p],
				      colprob->sgl[g][p+w])));
      DivEdouble(&pp->dbl[p][w], prob);
    }
}

/*------------------------------------------------------------------*/

/* Allocates memory for Cyk variable for a window of length
   wlen. */

Cyk *MakeCyk(Grammar *grammar, int wlen)
{
  Cyk *cyk = (Cyk *)malloc(sizeof(Cyk));
  int n = grammar->nont->count;
  int i, j;

  cyk->wlen = wlen;
  cyk->grammar = grammar;
  cyk->prob = (Edouble ***)malloc(n * sizeof(Edouble **));

  for (j = 0; j < n; j++)
    {
      cyk->prob[j] = (Edouble **)malloc(wlen * sizeof(Edouble *));
      for (i = 0; i < wlen; i++)
	cyk->prob[j][i] = (Edouble *)malloc( (wlen - i) * sizeof(Edouble));
    }

  cyk->trace = (Trace ***)malloc(n * sizeof(Trace **));

  for (j = 0; j < n; j++)
    {
      cyk->trace[j] = (Trace **)malloc(wlen * sizeof(Trace *));
      for (i = 0; i < wlen; i++)
	cyk->trace[j][i] = (Trace *)malloc( (wlen - i) * sizeof(Trace));
    }
  
  cyk->structure = (int *)malloc(wlen * sizeof(int));

  return cyk;
}

/*------------------------------------------------------------------*/

/* This initializes the cyk variables according to the
 probabilities in colprob. */
void InitCyk(Cyk *e, Colprob *colprob)
{
  int n;   /* nonterminal iterator */
  int p;   /* position iterator */
  int w;   /* string length iterator */
  int w2;   /* string length iterator for left string in N -> NN */
  int wl = e->wlen;
  LListCounter *lcount = MakeCounter(e->grammar->ruls, FIRST);
  Ruls *ruls;
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Ruldd *ruldd;
  Rulnn *rulnn;

  Backtrack *bt;
  void *rule;
  LList *stack;
  int nont, pos, wlen;
  int stp;

  /* Set Inside variables to zero */

  for (n = 0; n < e->grammar->nont->count; n++)
    for (p = 0; p < wl; p++)
      for (w = 0; w + p < wl; w++)
	e->prob[n][p][w] = Dbl2Edbl(0.);

  /* stringlength 0 initialization */

  InitCounter(lcount, e->grammar->ruls, FIRST);  /* N -> s rules */
  while((ruls = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      if (MaxEdouble(&e->prob[ruls->lnt][p][0],
	      ProdEdouble(ruls->prob,
			  colprob->sgl[ruls->rsg][p]))) {
	e->trace[ruls->lnt][p][0].type = 0;
	e->trace[ruls->lnt][p][0].rule = (Ruls *) ruls;
      }

  /* ---------- NOTICE ------------- */
  /* ERROR: SORT Ruln, CORRECT LATER */
  /* ---------- NOTICE ------------- */

  InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
  while((ruln = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      if (MaxEdouble(&e->prob[ruln->lnt][p][0],
	      ProdEdouble(ruln->prob,
		      e->prob[ruln->rnt][p][0]))) {
	e->trace[ruln->lnt][p][0].type = 1;
	e->trace[ruln->lnt][p][0].rule = (Ruln *) ruln;
      }

  /* Use of the rest of the rules */

  for (w = 1; w < wl; w++) {

    InitCounter(lcount, e->grammar->rulsn, FIRST);  /* N -> sN rules */
    while((rulsn = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	if (MaxEdouble(&e->prob[rulsn->lnt][p][w],
	      ProdEdouble(rulsn->prob,
	      ProdEdouble(colprob->sgl[rulsn->rsg][p],
			 e->prob[rulsn->rnt][p+1][w-1])))) {
	  e->trace[rulsn->lnt][p][w].type = 2;
	  e->trace[rulsn->lnt][p][w].rule = (Rulsn *) rulsn;
	}
	
    InitCounter(lcount, e->grammar->rulns, FIRST);  /* N -> Ns rules */
    while((rulns = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	if (MaxEdouble(&e->prob[rulns->lnt][p][w],
	      ProdEdouble(rulns->prob,
              ProdEdouble(colprob->sgl[rulns->rsg][p+w],
			 e->prob[rulns->rnt][p][w-1])))) {
	  e->trace[rulns->lnt][p][w].type = 3;
	  e->trace[rulns->lnt][p][w].rule = (Rulns *) rulns;
	}
	
    if (w >= 2) {  /* minimum stringlength of 3 */
      InitCounter(lcount, e->grammar->ruldnd, FIRST);  /* N -> dNd rules */
      while((ruldnd = Next(lcount)) != NULL)
	for (p = 0; w + p < wl; p++)
	  if (MaxEdouble(&e->prob[ruldnd->lnt][p][w],
		  ProdEdouble(ruldnd->prob,
		  ProdEdouble(colprob->dbl[ruldnd->rdb][p][w],
			  e->prob[ruldnd->rnt][p+1][w-2])))) {
	    e->trace[ruldnd->lnt][p][w].type = 4;
	    e->trace[ruldnd->lnt][p][w].rule = (Ruldnd *) ruldnd;
	  }

    }

    if (w == 1) {  /* stringlength of 2 */
      InitCounter(lcount, e->grammar->ruldd, FIRST);  /* N -> dd rules */
      while((ruldd = Next(lcount)) != NULL)
	for (p = 0; w + p < wl; p++)
	  if (MaxEdouble(&e->prob[ruldd->lnt][p][w],
		  ProdEdouble(ruldd->prob,
			      colprob->dbl[ruldd->rdb][p][w]))) {
	    e->trace[ruldd->lnt][p][w].type = 6;
	    e->trace[ruldd->lnt][p][w].rule = (Ruldd *) ruldd;
	}
    }

    InitCounter(lcount, e->grammar->rulnn, FIRST);  /* N -> NN rules */
    while((rulnn = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	for (w2 = 0; w2 < w; w2++)
	  if (MaxEdouble(&e->prob[rulnn->lnt][p][w],
		  ProdEdouble(rulnn->prob,
		  ProdEdouble(e->prob[rulnn->rnt1][p][w2],
			  e->prob[rulnn->rnt2][p+w2+1][w-w2-1])))) {
	  e->trace[rulnn->lnt][p][w].type = 5;
	  e->trace[rulnn->lnt][p][w].wlen = w2;
	  e->trace[rulnn->lnt][p][w].rule = (Rulnn *) rulnn;
	}

    InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
    while((ruln = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	if (MaxEdouble(&e->prob[ruln->lnt][p][w],
	      ProdEdouble(ruln->prob,
		      e->prob[ruln->rnt][p][w]))) {
	  e->trace[ruln->lnt][p][w].type = 1;
	  e->trace[ruln->lnt][p][w].rule = (Ruln *) ruln;
	}
  }

  free(lcount);

  e->numpair = 0;

  stack = MakeLList();

  nont = FindNont('S', e->grammar->nont);
  pos = 0;
  wlen = wl-1;
  stp = 0;

  for (;;) {
    rule = e->trace[nont][pos][wlen].rule;
    switch (e->trace[nont][pos][wlen].type) {
    case 0:  /* Ruls */
      e->structure[pos] = -1;
      if (stack->count == 0) {
	stp = 1;
	break; }
      bt = (Backtrack *)Pop(stack);
      nont = bt->nont;
      pos = bt->pos;
      wlen = bt->wlen;
      free(bt);
      break;
    case 1: /* Ruln */
      nont = ((Ruln *)rule)->rnt;
      break;
    case 2: /* Rulsn */
      e->structure[pos] = -1;
      pos++;
      wlen--;
      nont = ((Rulsn *)rule)->rnt;
      break;
    case 3: /* Rulns */
      e->structure[pos+wlen] = -1;
      wlen--;
      nont = ((Rulns *)rule)->rnt;
      break;
    case 4: /* Ruldnd */
      e->structure[pos] = pos+wlen;
      e->structure[pos+wlen] = pos;
      pos++;
      wlen -= 2;
      nont = ((Ruldnd *)rule)->rnt;
      e->numpair++;
      break;
    case 5: /* Rulnn */
      bt = (Backtrack *)malloc(sizeof(Backtrack));
      bt->nont = ((Rulnn *)rule)->rnt2;
      bt->pos = pos+e->trace[nont][pos][wlen].wlen+1;
      bt->wlen = wlen-e->trace[nont][pos][wlen].wlen-1;
      Push(stack, (void *)bt);
      wlen = e->trace[nont][pos][wlen].wlen;
      nont = ((Rulnn *)rule)->rnt1;
      break;
    case 6: /* Ruldd */
      e->structure[pos] = pos+wlen;
      e->structure[pos+wlen] = pos;
      e->numpair++;
      if (stack->count == 0) {
	stp = 1;
	break; }
      bt = (Backtrack *)Pop(stack);
      nont = bt->nont;
      pos = bt->pos;
      wlen = bt->wlen;
      free(bt);
      break;
    }
    if (stp == 1)
      break;
  }
}
      
/*------------------------------------------------------------------*/

/* This initializes the cyk variables according to the
 probabilities in colprob. */
void MoveCyk(Cyk *e, Colprob *colprob, int pos)
{
  int n;   /* nonterminal iterator */
  int w;   /* string length iterator */
  int w2;   /* string length iterator for left string in N -> NN */
  int wl = e->wlen;
  LListCounter *lcount = MakeCounter(e->grammar->ruls, FIRST);
  Ruls *ruls;
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Ruldd *ruldd;
  Rulnn *rulnn;

  Backtrack *bt;
  void *rule;
  LList *stack;
  int nont, p, wlen;
  int stp;

  /* Set Inside variables to zero */

  for (n = 0; n < e->grammar->nont->count; n++)
      for (w = 0; w < wl; w++)
	e->prob[n][(pos-w)%(wl-w)][w] = Dbl2Edbl(0.);

  /* stringlength 0 initialization */

  InitCounter(lcount, e->grammar->ruls, FIRST);  /* N -> s rules */
  while((ruls = Next(lcount)) != NULL)
    if (MaxEdouble(&e->prob[ruls->lnt][pos % wl][0],
		   ProdEdouble(ruls->prob,
			       colprob->sgl[ruls->rsg][pos % wl]))) {
      e->trace[ruls->lnt][pos % wl][0].type = 0;
      e->trace[ruls->lnt][pos % wl][0].rule = (Ruls *) ruls;
    }

  /* ---------- NOTICE ------------- */
  /* ERROR: SORT Ruln, CORRECT LATER */
  /* ---------- NOTICE ------------- */

  InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
  while((ruln = Next(lcount)) != NULL)
    if (MaxEdouble(&e->prob[ruln->lnt][pos % wl][0],
		   ProdEdouble(ruln->prob,
			       e->prob[ruln->rnt][pos % wl][0]))) {
      e->trace[ruln->lnt][pos % wl][0].type = 1;
      e->trace[ruln->lnt][pos % wl][0].rule = (Ruln *) ruln;
      }

  /* Use of the rest of the rules */

  for (w = 1; w < wl; w++) {

    InitCounter(lcount, e->grammar->rulsn, FIRST);  /* N -> sN rules */
    while((rulsn = Next(lcount)) != NULL)
      if (MaxEdouble(&e->prob[rulsn->lnt][(pos-w)%(wl-w)][w],
	      ProdEdouble(rulsn->prob,
	      ProdEdouble(colprob->sgl[rulsn->rsg][(pos-w) % wl],
			 e->prob[rulsn->rnt][(pos-w+1)%(wl-(w-1))][w-1])))) {
	  e->trace[rulsn->lnt][(pos-w)%(wl-w)][w].type = 2;
	  e->trace[rulsn->lnt][(pos-w)%(wl-w)][w].rule = (Rulsn *) rulsn;
	}
	
    InitCounter(lcount, e->grammar->rulns, FIRST);  /* N -> Ns rules */
    while((rulns = Next(lcount)) != NULL)
      if (MaxEdouble(&e->prob[rulns->lnt][(pos-w)%(wl-w)][w],
	      ProdEdouble(rulns->prob,
              ProdEdouble(colprob->sgl[rulns->rsg][pos % wl],
			 e->prob[rulns->rnt]
			  [(pos-w) % (wl-(w-1))][w-1])))) {
	  e->trace[rulns->lnt][(pos-w)%(wl-w)][w].type = 3;
	  e->trace[rulns->lnt][(pos-w)%(wl-w)][w].rule = (Rulns *) rulns;
      }
	
    if (w >= 2) {  /* minimum stringlength of 3 */
      InitCounter(lcount, e->grammar->ruldnd, FIRST);  /* N -> dNd rules */
      while((ruldnd = Next(lcount)) != NULL)
	if (MaxEdouble(&e->prob[ruldnd->lnt][(pos-w)%(wl-w)][w],
		  ProdEdouble(ruldnd->prob,
		  ProdEdouble(colprob->dbl[ruldnd->rdb][(pos-w)%(wl-w)][w],
			  e->prob[ruldnd->rnt][(pos-w+1)%(wl-(w-2))][w-2])))) {
	  e->trace[ruldnd->lnt][(pos-w)%(wl-w)][w].type = 4;
	  e->trace[ruldnd->lnt][(pos-w)%(wl-w)][w].rule = (Ruldnd *) ruldnd;
	}

    }

    if (w == 1) {  /* stringlength of 2 */
      InitCounter(lcount, e->grammar->ruldd, FIRST);  /* N -> dd rules */
      while((ruldd = Next(lcount)) != NULL)
	if (MaxEdouble(&e->prob[ruldd->lnt][(pos-w)%(wl-w)][w],
		  ProdEdouble(ruldd->prob,
			      colprob->dbl[ruldd->rdb][(pos-w)%(wl-w)][w]))) {
	  e->trace[ruldd->lnt][(pos-w)%(wl-w)][w].type = 6;
	  e->trace[ruldd->lnt][(pos-w)%(wl-w)][w].rule = (Ruldd *) ruldd;
	}

    }

    InitCounter(lcount, e->grammar->rulnn, FIRST);  /* N -> NN rules */
    while((rulnn = Next(lcount)) != NULL)
	for (w2 = 0; w2 < w; w2++)
	  if (MaxEdouble(&e->prob[rulnn->lnt][(pos-w)%(wl-w)][w],
		  ProdEdouble(rulnn->prob,
		  ProdEdouble(e->prob[rulnn->rnt1][(pos-w)%(wl-w2)][w2],
			  e->prob[rulnn->rnt2]
			      [(pos-w+w2+1)%(wl-(w-w2-1))][w-w2-1])))) {
	  e->trace[rulnn->lnt][(pos-w)%(wl-w)][w].type = 5;
	  e->trace[rulnn->lnt][(pos-w)%(wl-w)][w].wlen = w2;
	  e->trace[rulnn->lnt][(pos-w)%(wl-w)][w].rule = (Rulnn *) rulnn;
	}

    InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
    while((ruln = Next(lcount)) != NULL)
	if (MaxEdouble(&e->prob[ruln->lnt][(pos-w)%(wl-w)][w],
	      ProdEdouble(ruln->prob,
		      e->prob[ruln->rnt][(pos-w)%(wl-w)][w]))) {
	  e->trace[ruln->lnt][(pos-w)%(wl-w)][w].type = 1;
	  e->trace[ruln->lnt][(pos-w)%(wl-w)][w].rule = (Ruln *) ruln;
	}
  }

  free(lcount);

  e->numpair = 0;

  stack = MakeLList();

  nont = 0;
  p = 0;
  wlen = wl-1;
  stp = 0;

  for (;;) {
    rule = e->trace[nont][(pos-wl+p+1)%(wl-wlen)][wlen].rule;
    switch (e->trace[nont][(pos-wl+p+1)%(wl-wlen)][wlen].type) {
    case 0:  /* Ruls */
      e->structure[p] = -1;
      if (stack->count == 0) {
	stp = 1;
	break; }
      bt = (Backtrack *)Pop(stack);
      nont = bt->nont;
      p = bt->pos;
      wlen = bt->wlen;
      free(bt);
      break;
    case 1: /* Ruln */
      nont = ((Ruln *)rule)->rnt;
      break;
    case 2: /* Rulsn */
      e->structure[p] = -1;
      p++;
      wlen--;
      nont = ((Rulsn *)rule)->rnt;
      break;
    case 3: /* Rulns */
      e->structure[p+wlen] = -1;
      wlen--;
      nont = ((Rulsn *)rule)->rnt;
      break;
    case 4: /* Ruldnd */
      e->structure[p] = p+wlen;
      e->structure[p+wlen] = p;
      p++;
      wlen -= 2;
      nont = ((Rulsn *)rule)->rnt;
      e->numpair++;
      break;
    case 5: /* Rulnn */
      bt = (Backtrack *)malloc(sizeof(Backtrack));
      bt->nont = ((Rulnn *)rule)->rnt2;
      bt->pos = p+e->trace[nont][(pos-wl+p+1)%(wl-wlen)][wlen].wlen+1;
      bt->wlen = wlen-e->trace[nont][(pos-wl+p+1)%(wl-wlen)][wlen].wlen-1;
      Push(stack, (void *)bt);
      wlen = e->trace[nont][(pos-wl+p+1)%(wl-wlen)][wlen].wlen;
      nont = ((Rulnn *)rule)->rnt1;
      break;
    case 6: /* Ruldd */
      e->structure[p] = p+wlen;
      e->structure[p+wlen] = p;
      e->numpair++;
      if (stack->count == 0) {
	stp = 1;
	break; }
      bt = (Backtrack *)Pop(stack);
      nont = bt->nont;
      pos = bt->pos;
      wlen = bt->wlen;
      free(bt);
      break;
    }
    if (stp == 1)
      break;
  }
}
      
/*------------------------------------------------------------------*/

void ModifyCyk(Cyk *e, Pairprob *pp)
{
  int i, j;
  int maxj;
  double max;
  
  for (i = 0; i < e->wlen; i++) {
    max = 0;
    maxj = -1;
    for (j = 0; j < e->wlen; j++) {
      if (j > i && Edbl2Dbl(pp->dbl[i][j-i]) > max) {
	max = Edbl2Dbl(pp->dbl[i][j-i]);
	maxj = j;
      }
      else if (i > j && Edbl2Dbl(pp->dbl[j][i-j]) > max) {
	max = Edbl2Dbl(pp->dbl[j][i-j]);
	maxj = j;
      }
    }
    if (Edbl2Dbl(pp->sgl[i]) > max)
      e->structure[i] = -1;
    else
      e->structure[i] = maxj;
  }
  
  for (i = 0; i < e->wlen; i++)
    if (e->structure[i] != -1 &&
	e->structure[e->structure[i]] != i)
      if (e->structure[e->structure[i]] == -1) {
	if (2*Edbl2Dbl(pp->dbl[i][
	      e->structure[i]>i?e->structure[i]-i:i-e->structure[i]
                ]) <
	    Edbl2Dbl(pp->sgl[i])+Edbl2Dbl(pp->sgl[e->structure[i]])) {
	  e->structure[e->structure[i]] = -1;
	  e->structure[i] = -1;
	}
	else {
	  e->structure[e->structure[i]] = i;
	}
      }
      else
	e->structure[i] = -1;
}

/*------------------------------------------------------------------*/

void PrintCyk(FILE *fp, Cyk *e)
{
  int wl = e->wlen;
  int i;

  for (i = 0; i < wl; i++) {
    /*
    for (w = 0; w < wl; w++)
      if (w+i < wl)
	fprintf(fp, " %5.2f", -Edbl2Dbl(LogEdouble(e->prob[0][i][w])));
      else
	fprintf(fp, "      ");
    */
    if (e->structure[i] == -1)
      fprintf(fp, "     .\n");
    else
      fprintf(fp, " %5d\n", e->structure[i]);
  }
}


int *CykStr(Cyk *cyk)
{
  return cyk->structure;
}

/*------------------------------------------------------------------*/

void AddStr(Entry *entry, Alignstr *alignstr, Pairprob *pp)
{
  int wl = alignstr->align->len;
  int i;
  int pos;
  int align_bp_col, alignpos_col, seq_bp_col, seqpos_col;
  int cert_col;
  int *map;

  map = alignstr->map;

  align_bp_col = ReadColno(entry, "align_bp");
  alignpos_col = ReadColno(entry, "alignpos");
  seq_bp_col = ReadColno(entry, "seq_bp");
  seqpos_col = ReadColno(entry, "seqpos");

  if (align_bp_col == 0 && seq_bp_col == 0) {
    if (alignpos_col != 0)
      align_bp_col = EnsureCol(entry, "align_bp", "    .");
    else if (seqpos_col != 0)
      seq_bp_col = EnsureCol(entry, "seq_bp", "    .");
    else
      return;
  }

  if (pp != NULL)
    cert_col = EnsureCol(entry, "certainty", "0.0000");

  pos = 1;
  for (i = 0; i < wl; i++) {
    for (;pos < map[i]; pos++) {
      SetPair(entry, pos, 0,
	      align_bp_col, alignpos_col, seq_bp_col, seqpos_col);
      if (pp != NULL)
	ChgField(entry, pos, cert_col, "0.0000");
    }
    if (alignstr->str[i] >= 0) {
      SetPair(entry, pos, map[alignstr->str[i]],
	      align_bp_col, alignpos_col, seq_bp_col, seqpos_col);
      if (pp != NULL) {
	if (i > alignstr->str[i])
	  ChgField(entry, pos, cert_col, "%6.4f",
		   Edbl2Dbl(pp->dbl[alignstr->str[i]]
			    [i - alignstr->str[i]]));
	else
	  ChgField(entry, pos, cert_col, "%6.4f",
		   Edbl2Dbl(pp->dbl[i][alignstr->str[i] - i]));
      }
    }
    else {
      SetPair(entry, pos, 0,
	      align_bp_col, alignpos_col, seq_bp_col, seqpos_col);
      if (pp != NULL)
	ChgField(entry, pos, cert_col, "%6.4f", Edbl2Dbl(pp->sgl[i]));
    }
    pos++;
  }
}

/*------------------------------------------------------------------*/

/* Allocates memory for FCyk variable for a window of length
   wlen. */

FCyk *MakeFCyk(Grammar *grammar, int wlen)
{
  FCyk *cyk = (FCyk *)malloc(sizeof(FCyk));
  int n = grammar->nont->count;
  int i, j;

  cyk->wlen = wlen;
  cyk->grammar = grammar;
  cyk->prob = (Fdouble ***)malloc(n * sizeof(Fdouble **));

  for (j = 0; j < n; j++)
    {
      cyk->prob[j] = (Fdouble **)malloc(wlen * sizeof(Fdouble *));
      for (i = 0; i < wlen; i++)
	cyk->prob[j][i] = (Fdouble *)malloc( (wlen - i) * sizeof(Fdouble));
    }

  cyk->trace = (Trace ***)malloc(n * sizeof(Trace **));

  for (j = 0; j < n; j++)
    {
      cyk->trace[j] = (Trace **)malloc(wlen * sizeof(Trace *));
      for (i = 0; i < wlen; i++)
	cyk->trace[j][i] = (Trace *)malloc( (wlen - i) * sizeof(Trace));
    }
  
  cyk->structure = (int *)malloc(wlen * sizeof(int));

  return cyk;
}

/*------------------------------------------------------------------*/

/* This initializes the cyk variables according to the
 probabilities in colprob. */
void InitFCyk(FCyk *e, FColprob *fcolprob)
{
  int n;   /* nonterminal iterator */
  int p;   /* position iterator */
  int w;   /* string length iterator */
  int w2;   /* string length iterator for left string in N -> NN */
  int wl = e->wlen;
  LListCounter *lcount = MakeCounter(e->grammar->ruls, FIRST);
  Ruls *ruls;
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Rulnn *rulnn;

  Backtrack *bt;
  void *rule;
  LList *stack;
  int nont, pos, wlen;
  int stp;

  /* Set Inside variables to zero */

  for (n = 0; n < e->grammar->nont->count; n++)
    for (p = 0; p < wl; p++)
      for (w = 0; w + p < wl; w++)
	e->prob[n][p][w] = Dbl2Fdbl(0.);

  /* stringlength 0 initialization */

  InitCounter(lcount, e->grammar->ruls, FIRST);  /* N -> s rules */
  while((ruls = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      if (MaxFdouble(&e->prob[ruls->lnt][p][0],
	      ProdFdouble(ruls->fprob,
			  fcolprob->sgl[ruls->rsg][p]))) {
	e->trace[ruls->lnt][p][0].type = 0;
	e->trace[ruls->lnt][p][0].rule = (Ruls *) ruls;
      }

  /* ---------- NOTICE ------------- */
  /* ERROR: SORT Ruln, CORRECT LATER */
  /* ---------- NOTICE ------------- */

  InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
  while((ruln = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      if (MaxFdouble(&e->prob[ruln->lnt][p][0],
	      ProdFdouble(ruln->fprob,
		      e->prob[ruln->rnt][p][0]))) {
	e->trace[ruln->lnt][p][0].type = 1;
	e->trace[ruln->lnt][p][0].rule = (Ruln *) ruln;
      }

  /* Use of the rest of the rules */

  for (w = 1; w < wl; w++) {

    InitCounter(lcount, e->grammar->rulsn, FIRST);  /* N -> sN rules */
    while((rulsn = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	if (MaxFdouble(&e->prob[rulsn->lnt][p][w],
	      ProdFdouble(rulsn->fprob,
	      ProdFdouble(fcolprob->sgl[rulsn->rsg][p],
			 e->prob[rulsn->rnt][p+1][w-1])))) {
	  e->trace[rulsn->lnt][p][w].type = 2;
	  e->trace[rulsn->lnt][p][w].rule = (Rulsn *) rulsn;
	}
	
    InitCounter(lcount, e->grammar->rulns, FIRST);  /* N -> Ns rules */
    while((rulns = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	if (MaxFdouble(&e->prob[rulns->lnt][p][w],
	      ProdFdouble(rulns->fprob,
              ProdFdouble(fcolprob->sgl[rulns->rsg][p+w],
			 e->prob[rulns->rnt][p][w-1])))) {
	  e->trace[rulns->lnt][p][w].type = 3;
	  e->trace[rulns->lnt][p][w].rule = (Rulns *) rulns;
	}
	
    if (w >= 2) {  /* minimum stringlength of 3 */
      InitCounter(lcount, e->grammar->ruldnd, FIRST);  /* N -> dNd rules */
      while((ruldnd = Next(lcount)) != NULL)
	for (p = 0; w + p < wl; p++)
	  if (MaxFdouble(&e->prob[ruldnd->lnt][p][w],
		  ProdFdouble(ruldnd->fprob,
		  ProdFdouble(fcolprob->dbl[ruldnd->rdb][p][w],
			  e->prob[ruldnd->rnt][p+1][w-2])))) {
	  e->trace[ruldnd->lnt][p][w].type = 4;
	  e->trace[ruldnd->lnt][p][w].rule = (Ruldnd *) ruldnd;
	}

    }

    InitCounter(lcount, e->grammar->rulnn, FIRST);  /* N -> NN rules */
    while((rulnn = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	for (w2 = 0; w2 < w; w2++)
	  if (MaxFdouble(&e->prob[rulnn->lnt][p][w],
		  ProdFdouble(rulnn->fprob,
		  ProdFdouble(e->prob[rulnn->rnt1][p][w2],
			  e->prob[rulnn->rnt2][p+w2+1][w-w2-1])))) {
	  e->trace[rulnn->lnt][p][w].type = 5;
	  e->trace[rulnn->lnt][p][w].wlen = w2;
	  e->trace[rulnn->lnt][p][w].rule = (Rulnn *) rulnn;
	}

    InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
    while((ruln = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	if (MaxFdouble(&e->prob[ruln->lnt][p][w],
	      ProdFdouble(ruln->fprob,
		      e->prob[ruln->rnt][p][w]))) {
	  e->trace[ruln->lnt][p][w].type = 1;
	  e->trace[ruln->lnt][p][w].rule = (Ruln *) ruln;
	}
  }

  free(lcount);

  e->numpair = 0;

  stack = MakeLList();

  nont = 0;
  pos = 0;
  wlen = wl-1;
  stp = 0;

  for (;;) {
    rule = e->trace[nont][pos][wlen].rule;
    switch (e->trace[nont][pos][wlen].type) {
    case 0:  /* Ruls */
      e->structure[pos] = -1;
      if (stack->count == 0) {
	stp = 1;
	break; }
      bt = (Backtrack *)Pop(stack);
      nont = bt->nont;
      pos = bt->pos;
      wlen = bt->wlen;
      free(bt);
      break;
    case 1: /* Ruln */
      nont = ((Ruln *)rule)->rnt;
      break;
    case 2: /* Rulsn */
      e->structure[pos] = -1;
      pos++;
      wlen--;
      nont = ((Rulsn *)rule)->rnt;
      break;
    case 3: /* Rulns */
      e->structure[pos+wlen] = -1;
      wlen--;
      nont = ((Rulsn *)rule)->rnt;
      break;
    case 4: /* Ruldnd */
      e->structure[pos] = pos+wlen;
      e->structure[pos+wlen] = pos;
      pos++;
      wlen -= 2;
      nont = ((Rulsn *)rule)->rnt;
      e->numpair++;
      break;
    case 5: /* Rulnn */
      bt = (Backtrack *)malloc(sizeof(Backtrack));
      bt->nont = ((Rulnn *)rule)->rnt2;
      bt->pos = pos+e->trace[nont][pos][wlen].wlen+1;
      bt->wlen = wlen-e->trace[nont][pos][wlen].wlen-1;
      Push(stack, (void *)bt);
      wlen = e->trace[nont][pos][wlen].wlen;
      nont = ((Rulnn *)rule)->rnt1;
      break;
    }
    if (stp == 1)
      break;
  }
}
      
/*------------------------------------------------------------------*/

/* This initializes the cyk variables according to the
 probabilities in colprob. */
void MoveFCyk(FCyk *e, FColprob *fcolprob, int pos)
{
  int n;   /* nonterminal iterator */
  int w;   /* string length iterator */
  int w2;   /* string length iterator for left string in N -> NN */
  int wl = e->wlen;
  LListCounter *lcount = MakeCounter(e->grammar->ruls, FIRST);
  Ruls *ruls;
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Rulnn *rulnn;

  Backtrack *bt;
  void *rule;
  LList *stack;
  int nont, p, wlen;
  int stp;

  /* Set Inside variables to zero */

  for (n = 0; n < e->grammar->nont->count; n++)
      for (w = 0; w < wl; w++)
	e->prob[n][(pos-w)%(wl-w)][w] = Dbl2Fdbl(0.);

  /* stringlength 0 initialization */

  InitCounter(lcount, e->grammar->ruls, FIRST);  /* N -> s rules */
  while((ruls = Next(lcount)) != NULL)
    if (MaxFdouble(&e->prob[ruls->lnt][pos % wl][0],
		   ProdFdouble(ruls->fprob,
			       fcolprob->sgl[ruls->rsg][pos % wl]))) {
      e->trace[ruls->lnt][pos % wl][0].type = 0;
      e->trace[ruls->lnt][pos % wl][0].rule = (Ruls *) ruls;
    }

  /* ---------- NOTICE ------------- */
  /* ERROR: SORT Ruln, CORRECT LATER */
  /* ---------- NOTICE ------------- */

  InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
  while((ruln = Next(lcount)) != NULL)
    if (MaxFdouble(&e->prob[ruln->lnt][pos % wl][0],
		   ProdFdouble(ruln->fprob,
			       e->prob[ruln->rnt][pos % wl][0]))) {
      e->trace[ruln->lnt][pos % wl][0].type = 1;
      e->trace[ruln->lnt][pos % wl][0].rule = (Ruln *) ruln;
      }

  /* Use of the rest of the rules */

  for (w = 1; w < wl; w++) {

    InitCounter(lcount, e->grammar->rulsn, FIRST);  /* N -> sN rules */
    while((rulsn = Next(lcount)) != NULL)
      if (MaxFdouble(&e->prob[rulsn->lnt][(pos-w)%(wl-w)][w],
	      ProdFdouble(rulsn->fprob,
	      ProdFdouble(fcolprob->sgl[rulsn->rsg][(pos-w) % wl],
			 e->prob[rulsn->rnt][(pos-w+1)%(wl-(w-1))][w-1])))) {
	  e->trace[rulsn->lnt][(pos-w)%(wl-w)][w].type = 2;
	  e->trace[rulsn->lnt][(pos-w)%(wl-w)][w].rule = (Rulsn *) rulsn;
	}
	
    InitCounter(lcount, e->grammar->rulns, FIRST);  /* N -> Ns rules */
    while((rulns = Next(lcount)) != NULL)
      if (MaxFdouble(&e->prob[rulns->lnt][(pos-w)%(wl-w)][w],
	      ProdFdouble(rulns->fprob,
              ProdFdouble(fcolprob->sgl[rulns->rsg][pos % wl],
			 e->prob[rulns->rnt]
			  [(pos-w) % (wl-(w-1))][w-1])))) {
	  e->trace[rulns->lnt][(pos-w)%(wl-w)][w].type = 3;
	  e->trace[rulns->lnt][(pos-w)%(wl-w)][w].rule = (Rulns *) rulns;
      }
	
    if (w >= 2) {  /* minimum stringlength of 3 */
      InitCounter(lcount, e->grammar->ruldnd, FIRST);  /* N -> dNd rules */
      while((ruldnd = Next(lcount)) != NULL)
	if (MaxFdouble(&e->prob[ruldnd->lnt][(pos-w)%(wl-w)][w],
		  ProdFdouble(ruldnd->fprob,
		  ProdFdouble(fcolprob->dbl[ruldnd->rdb][(pos-w)%(wl-w)][w],
			  e->prob[ruldnd->rnt][(pos-w+1)%(wl-(w-2))][w-2])))) {
	  e->trace[ruldnd->lnt][(pos-w)%(wl-w)][w].type = 4;
	  e->trace[ruldnd->lnt][(pos-w)%(wl-w)][w].rule = (Ruldnd *) ruldnd;
	}

    }

    InitCounter(lcount, e->grammar->rulnn, FIRST);  /* N -> NN rules */
    while((rulnn = Next(lcount)) != NULL)
	for (w2 = 0; w2 < w; w2++)
	  if (MaxFdouble(&e->prob[rulnn->lnt][(pos-w)%(wl-w)][w],
		  ProdFdouble(rulnn->fprob,
		  ProdFdouble(e->prob[rulnn->rnt1][(pos-w)%(wl-w2)][w2],
			  e->prob[rulnn->rnt2]
			      [(pos-w+w2+1)%(wl-(w-w2-1))][w-w2-1])))) {
	  e->trace[rulnn->lnt][(pos-w)%(wl-w)][w].type = 5;
	  e->trace[rulnn->lnt][(pos-w)%(wl-w)][w].wlen = w2;
	  e->trace[rulnn->lnt][(pos-w)%(wl-w)][w].rule = (Rulnn *) rulnn;
	}

    InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
    while((ruln = Next(lcount)) != NULL)
	if (MaxFdouble(&e->prob[ruln->lnt][(pos-w)%(wl-w)][w],
	      ProdFdouble(ruln->fprob,
		      e->prob[ruln->rnt][(pos-w)%(wl-w)][w]))) {
	  e->trace[ruln->lnt][(pos-w)%(wl-w)][w].type = 1;
	  e->trace[ruln->lnt][(pos-w)%(wl-w)][w].rule = (Ruln *) ruln;
	}
  }

  free(lcount);

  e->numpair = 0;

  stack = MakeLList();

  nont = 0;
  p = 0;
  wlen = wl-1;
  stp = 0;

  for (;;) {
    rule = e->trace[nont][(pos-wl+p+1)%(wl-wlen)][wlen].rule;
    switch (e->trace[nont][(pos-wl+p+1)%(wl-wlen)][wlen].type) {
    case 0:  /* Ruls */
      e->structure[p] = -1;
      if (stack->count == 0) {
	stp = 1;
	break; }
      bt = (Backtrack *)Pop(stack);
      nont = bt->nont;
      p = bt->pos;
      wlen = bt->wlen;
      free(bt);
      break;
    case 1: /* Ruln */
      nont = ((Ruln *)rule)->rnt;
      break;
    case 2: /* Rulsn */
      e->structure[p] = -1;
      p++;
      wlen--;
      nont = ((Rulsn *)rule)->rnt;
      break;
    case 3: /* Rulns */
      e->structure[p+wlen] = -1;
      wlen--;
      nont = ((Rulsn *)rule)->rnt;
      break;
    case 4: /* Ruldnd */
      e->structure[p] = p+wlen;
      e->structure[p+wlen] = p;
      p++;
      wlen -= 2;
      nont = ((Rulsn *)rule)->rnt;
      e->numpair++;
      break;
    case 5: /* Rulnn */
      bt = (Backtrack *)malloc(sizeof(Backtrack));
      bt->nont = ((Rulnn *)rule)->rnt2;
      bt->pos = p+e->trace[nont][(pos-wl+p+1)%(wl-wlen)][wlen].wlen+1;
      bt->wlen = wlen-e->trace[nont][(pos-wl+p+1)%(wl-wlen)][wlen].wlen-1;
      Push(stack, (void *)bt);
      wlen = e->trace[nont][(pos-wl+p+1)%(wl-wlen)][wlen].wlen;
      nont = ((Rulnn *)rule)->rnt1;
      break;
    }
    if (stp == 1)
      break;
  }
}
      
/*------------------------------------------------------------------*/

void PrintFCyk(FILE *fp, FCyk *e)
{
  int wl = e->wlen;
  int i;

  for (i = 0; i < wl; i++) {
    /*
    for (w = 0; w < wl; w++)
      if (w+i < wl)
	fprintf(fp, " %5d", (int) (e->prob[0][i][w]));
      else
	fprintf(fp, "      ");
    */
    if (e->structure[i] == -1)
      fprintf(fp, "     .\n");
    else
      fprintf(fp, " %5d\n", e->structure[i]);
  }
}

int *FCykStr(FCyk *cyk)
{
  return cyk->structure;
}
/*------------------------------------------------------------------*/

/* Allocates memory for Cyk variable for a window of length
   wlen. */

PpCyk *MakePpCyk(int wlen)
{
  PpCyk *cyk = (PpCyk *)malloc(sizeof(PpCyk));
  int i;

  cyk->wlen = wlen;

  cyk->prob = (Edouble **)malloc(wlen * sizeof(Edouble *));
  for (i = 0; i < wlen; i++)
    cyk->prob[i] = (Edouble *)malloc( (wlen - i) * sizeof(Edouble));

  cyk->trace = (Trace **)malloc(wlen * sizeof(Trace *));
  for (i = 0; i < wlen; i++)
    cyk->trace[i] = (Trace *)malloc( (wlen - i) * sizeof(Trace));
  
  cyk->structure = (int *)malloc(wlen * sizeof(int));

  return cyk;
}

/*------------------------------------------------------------------*/

/* This initializes the cyk variables according to the
 probabilities in colprob. */
void InitPpCyk(PpCyk *e, Pairprob *pp)
{
  int p;   /* position iterator */
  int w;   /* string length iterator */
  int w2;   /* string length iterator for left string in N -> NN */
  int wl = e->wlen;

  Backtrack *bt;
  LList *stack;
  int pos, wlen;
  int stp;

  if (wl == 0)
    return;

  /* Set Inside variables to zero */

  for (p = 0; p < wl; p++)
     for (w = 0; w + p < wl; w++)
      e->prob[p][w] = Dbl2Edbl(0.);
  
  /* stringlength 0 initialization */

  for (p = 0; p < wl; p++) {
    e->prob[p][0] = pp->sgl[p];
    e->trace[p][0].type = 0;
  }

  /* Use of the rest of the rules */

  for (w = 1; w < wl; w++) {
    for (p = 0; w + p < wl; p++)
      if (MaxEdouble(&e->prob[p][w],
		     SumEdouble(pp->sgl[p],
				e->prob[p+1][w-1])))
	e->trace[p][w].type = 2;

    for (p = 0; w + p < wl; p++)
      if (MaxEdouble(&e->prob[p][w],
		     SumEdouble(pp->sgl[p+w],
				e->prob[p][w-1])))
	e->trace[p][w].type = 3;
	
    if (w >= 2) {  /* minimum stringlength of 3 */
      for (p = 0; w + p < wl; p++) {
	if (MaxEdouble(&e->prob[p][w],
		       SumEdouble(pp->dbl[p][w],
		       SumEdouble(pp->dbl[p][w],
				  e->prob[p+1][w-2])))) {
	  e->trace[p][w].type = 4;
	}
      }
    }

    if (w == 1) {  /* stringlength of 2 */
      for (p = 0; w + p < wl; p++) {
	if (MaxEdouble(&e->prob[p][w],
		       SumEdouble(pp->dbl[p][w],
				  pp->dbl[p][w]))) {
	  e->trace[p][w].type = 6;
	}
      }
    }

    for (p = 0; w + p < wl; p++)
      for (w2 = 0; w2 < w; w2++)
	if (MaxEdouble(&e->prob[p][w],
		       SumEdouble(e->prob[p][w2],
				  e->prob[p+w2+1][w-w2-1]))) {
	  e->trace[p][w].type = 5;
	  e->trace[p][w].wlen = w2;
	}
  }

  stack = MakeLList();

  pos = 0;
  wlen = wl-1;
  stp = 0;

  for (;;) {
    switch (e->trace[pos][wlen].type) {
    case 0:  /* Ruls */
      e->structure[pos] = -1;
      if (stack->count == 0) {
	stp = 1;
	break; }
      bt = (Backtrack *)Pop(stack);
      pos = bt->pos;
      wlen = bt->wlen;
      free(bt);
      break;
    case 2: /* Rulsn */
      e->structure[pos] = -1;
      pos++;
      wlen--;
      break;
    case 3: /* Rulns */
      e->structure[pos+wlen] = -1;
      wlen--;
      break;
    case 4: /* Ruldnd */
      e->structure[pos] = pos+wlen;
      e->structure[pos+wlen] = pos;
      pos++;
      wlen -= 2;
      break;
    case 5: /* Rulnn */
      bt = (Backtrack *)malloc(sizeof(Backtrack));
      bt->pos = pos+e->trace[pos][wlen].wlen+1;
      bt->wlen = wlen-e->trace[pos][wlen].wlen-1;
      Push(stack, (void *)bt);
      wlen = e->trace[pos][wlen].wlen;
      break;
    case 6: /* Ruldd */
      e->structure[pos] = pos+wlen;
      e->structure[pos+wlen] = pos;
      if (stack->count == 0) {
	stp = 1;
	break; }
      bt = (Backtrack *)Pop(stack);
      pos = bt->pos;
      wlen = bt->wlen;
      free(bt);
      break;
    }
    if (stp == 1)
      break;
  }
}

int *PpCykStr(PpCyk *cyk)
{
  return cyk->structure;
}

void PrintPpCyk(FILE *fp, PpCyk *e)
{
  int wl = e->wlen;
  int i, w;

  for (i = 0; i < wl; i++) {
    for (w = 0; w < wl; w++)
      if (w+i < wl)
	fprintf(fp, " %5.2f", Edbl2Dbl(e->prob[i][w]));
      else
	fprintf(fp, "      ");
    if (e->structure[i] == -1)
      fprintf(fp, "     .\n");
    else
      fprintf(fp, " %5d\n", e->structure[i]);
  }
}

/*------------------------------------------------------------------*/

/* Allocates memory for Pair counter for a window of length wlen. */

Paircnt *MakePaircnt(Grammar *grammar, int wlen)
{
  Paircnt* paircnt = (Paircnt *)malloc(sizeof(Paircnt));
  int n = grammar->nont->count;
  int i, j;

  paircnt->wlen = wlen;
  paircnt->grammar = grammar;
  paircnt->cnt = (Edouble ***)malloc(n * sizeof(Edouble **));

  for (j = 0; j < n; j++)
    {
      paircnt->cnt[j] = (Edouble **)malloc(wlen * sizeof(Edouble *));
      for (i = 0; i < wlen; i++)
	paircnt->cnt[j][i] = (Edouble *)malloc( (wlen - i) * sizeof(Edouble));
    }

  return paircnt;
}

/*------------------------------------------------------------------*/

void FreePaircnt(Grammar *grammar, Paircnt *paircnt)
{
  int n = grammar->nont->count;
  int i, j;
  int wlen;

  wlen = paircnt->wlen;

  for (j = 0; j < n; j++) {
    for (i = 0; i < wlen; i++)
      free(paircnt->cnt[j][i]);
    free(paircnt->cnt[j]);
  }

  free(paircnt->cnt);
  free(paircnt);
}

/*------------------------------------------------------------------*/

/* This initializes the inside variables according to the
 probabilities in colprob. */

void InitPaircnt(Inside *e, Paircnt *paircnt, Colprob *colprob)
{
  int n;   /* nonterminal iterator */
  int p;   /* position iterator */
  int w;   /* string length iterator */
  int w2;   /* string length iterator for left string in N -> NN */
  int wl = e->wlen;
  LListCounter *lcount = MakeCounter(e->grammar->ruls, FIRST);
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Ruldd *ruldd;
  Rulnn *rulnn;

  /* Set counters to zero */

  for (n = 0; n < e->grammar->nont->count; n++)
    for (p = 0; p < wl; p++)
      for (w = 0; w + p < wl; w++)
	paircnt->cnt[n][p][w] = Dbl2Edbl(0.);

  /* Use of the rules */

  for (w = 1; w < wl; w++) {

    InitCounter(lcount, e->grammar->rulsn, FIRST);  /* N -> sN rules */
    while((rulsn = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	AddEdouble(&paircnt->cnt[rulsn->lnt][p][w],
		ProdEdouble(paircnt->cnt[rulsn->rnt][p+1][w-1],
		ProdEdouble(rulsn->prob,
		ProdEdouble(colprob->sgl[rulsn->rsg][p],
			e->prob[rulsn->rnt][p+1][w-1]))));
	
    InitCounter(lcount, e->grammar->rulns, FIRST);  /* N -> Ns rules */
    while((rulns = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	AddEdouble(&paircnt->cnt[rulns->lnt][p][w],
		ProdEdouble(paircnt->cnt[rulns->rnt][p][w-1],
		ProdEdouble(rulns->prob,
		ProdEdouble(colprob->sgl[rulns->rsg][p+w],
			e->prob[rulns->rnt][p][w-1]))));
	
    if (w >= 2) {  /* minimum stringlength of 3 */
    InitCounter(lcount, e->grammar->ruldnd, FIRST);  /* N -> dNd rules */
    while((ruldnd = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	AddEdouble(&paircnt->cnt[ruldnd->lnt][p][w],
		ProdEdouble(SumEdouble(paircnt->cnt[ruldnd->rnt][p+1][w-2],
				      Dbl2Edbl(1.)),
                ProdEdouble(ruldnd->prob,
       		ProdEdouble(colprob->dbl[ruldnd->rdb][p][w],
			e->prob[ruldnd->rnt][p+1][w-2]))));
    }

    if (w == 1) {  /* minimum stringlength of 2 */
    InitCounter(lcount, e->grammar->ruldd, FIRST);  /* N -> dd rules */
    while((ruldd = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	AddEdouble(&paircnt->cnt[ruldd->lnt][p][w],
		   ProdEdouble(ruldd->prob,
			       colprob->dbl[ruldd->rdb][p][w]));
    }

    InitCounter(lcount, e->grammar->rulnn, FIRST);  /* N -> NN rules */
    while((rulnn = Next(lcount)) != NULL)
      for (p = 0; w + p < wl; p++)
	for (w2 = 0; w2 < w; w2++)
 	  AddEdouble(&paircnt->cnt[rulnn->lnt][p][w], ProdEdouble(
		SumEdouble(paircnt->cnt[rulnn->rnt1][p][w2],
			paircnt->cnt[rulnn->rnt2][p+w2+1][w-w2-1]),
	        ProdEdouble(rulnn->prob,
		ProdEdouble(e->prob[rulnn->rnt1][p][w2],
			e->prob[rulnn->rnt2][p+w2+1][w-w2-1]))));

    InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
    while((ruln = Next(lcount)) != NULL) {
      for (p = 0; w + p < wl; p++)
	AddEdouble(&paircnt->cnt[ruln->lnt][p][w],
		ProdEdouble(paircnt->cnt[ruln->rnt][p][w],
			ruln->prob));
    }
    for (n = 0; n < e->grammar->nont->count; n++)
      for (p = 0; w + p < wl; p++) {
	if (!IsZeroEdouble(paircnt->cnt[n][p][w]))
	  DivEdouble(&paircnt->cnt[n][p][w], e->prob[n][p][w]);
      }
  }

  free(lcount);
}
      
/*------------------------------------------------------------------*/

/* This moves the inside variables according to the probabilities in
   colprob. */

void MovePaircnt(Inside *e, Paircnt *paircnt, Colprob *colprob, int pos)
{
  int n;   /* nonterminal iterator */
  int w;   /* string length iterator */
  int w2;   /* string length iterator for left string in N -> NN */
  LListCounter *lcount = MakeCounter(e->grammar->ruls, FIRST);
  int wl = e->wlen;
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Ruldd *ruldd;
  Rulnn *rulnn;

  /* Set Inside variables and counters to zero */

  for (n = 0; n < e->grammar->nont->count; n++)
    for (w = 0; w < wl; w++)
      paircnt->cnt[n][(pos-w)%(wl-w)][w] = Dbl2Edbl(0.);

  /* Use of the rules */

  for (w = 1; w < wl; w++) {

    InitCounter(lcount, e->grammar->rulsn, FIRST);  /* N -> sN rules */
    while((rulsn = Next(lcount)) != NULL)
      AddEdouble(&paircnt->cnt[rulsn->lnt][(pos-w)%(wl-w)][w],
	      ProdEdouble(paircnt->cnt[rulsn->rnt][(pos-w+1)%(wl-(w-1))][w-1],
	      ProdEdouble(rulsn->prob,
              ProdEdouble(colprob->sgl[rulsn->rsg][(pos-w)],
		      e->prob[rulsn->rnt][(pos-w+1)%(wl-(w-1))][w-1]))));

    InitCounter(lcount, e->grammar->rulns, FIRST);  /* N -> Ns rules */
    while((rulns = Next(lcount)) != NULL)
      AddEdouble(&paircnt->cnt[rulns->lnt][(pos-w)%(wl-w)][w],
	      ProdEdouble(paircnt->cnt[rulns->rnt][(pos-w)%(wl-(w-1))][w-1],
	      ProdEdouble(rulns->prob,
	      ProdEdouble(colprob->sgl[rulns->rsg][pos],
		      e->prob[rulns->rnt][(pos-w)%(wl-(w-1))][w-1]))));

    if (w >= 2) {  /* minimum stringlength of 3 */
    InitCounter(lcount, e->grammar->ruldnd, FIRST);  /* N -> dNd rules */
    while((ruldnd = Next(lcount)) != NULL)
      AddEdouble(&paircnt->cnt[ruldnd->lnt][(pos-w)%(wl-w)][w],
	      ProdEdouble(SumEdouble(
		     paircnt->cnt[ruldnd->rnt][(pos-w+1)%(wl-(w-2))][w-2],
				      Dbl2Edbl(1.)),
              ProdEdouble(ruldnd->prob,
       	      ProdEdouble(colprob->dbl[ruldnd->rdb][(pos-w)%(wl-w)][w],
		      e->prob[ruldnd->rnt][(pos-w+1)%(wl-(w-2))][w-2]))));
    }

    if (w == 1) {  /* stringlength of 2 */
    InitCounter(lcount, e->grammar->ruldd, FIRST);  /* N -> dd rules */
    while((ruldd = Next(lcount)) != NULL)
      AddEdouble(&paircnt->cnt[ruldd->lnt][(pos-w)%(wl-w)][w],
              ProdEdouble(ruldd->prob,
			  colprob->dbl[ruldd->rdb][(pos-w)%(wl-w)][w]));
    }

    InitCounter(lcount, e->grammar->rulnn, FIRST);  /* N -> NN rules */
    while((rulnn = Next(lcount)) != NULL)
	for (w2 = 0; w2 < w; w2++)
 	  AddEdouble(&paircnt->cnt[rulnn->lnt][(pos-w)%(wl-w)][w],
		  ProdEdouble(SumEdouble(
		        paircnt->cnt[rulnn->rnt1][(pos-w)%(wl-w2)][w2],
			paircnt->cnt[rulnn->rnt2]
			       [(pos-w+w2+1)%(wl-(w-w2-1))][w-w2-1]),
		  ProdEdouble(rulnn->prob,
		  ProdEdouble(e->prob[rulnn->rnt1][(pos-w)%(wl-w2)][w2],
			  e->prob[rulnn->rnt2]\
			  [(pos-w+w2+1)%(wl-(w-w2-1))][w-w2-1]))));

    InitCounter(lcount, e->grammar->ruln, FIRST);  /* N -> N rules */
    while((ruln = Next(lcount)) != NULL)
      AddEdouble(&paircnt->cnt[ruln->lnt][(pos-w)%(wl-w)][w],
	      ProdEdouble(paircnt->cnt[ruln->rnt][(pos-w)%(wl-w)][w],
		      ruln->prob));

    for (n = 0; n < e->grammar->nont->count; n++)
      if (!IsZeroEdouble(paircnt->cnt[n][(pos-w)%(wl-w)][w]))
	DivEdouble(&paircnt->cnt[n][(pos-w)%(wl-w)][w],
		   e->prob[n][(pos-w)%(wl-w)][w]);
    
  }

  free(lcount);
}
      
/*------------------------------------------------------------------*/

