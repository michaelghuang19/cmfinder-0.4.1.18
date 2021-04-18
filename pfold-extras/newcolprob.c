#include "newcolprob.h"

double fmul = 0.96;//0.894603;//0.868534;

typedef struct tagSNodeInfo {
  Grammar *grammar;    /* The grammar */
  SuffixTree *stree;   /* Search tree */
  int *col;            /* Place for search string */
  int noseqs;          /* Number of leaf sequences below this node */
  int *seqs;           /* List of leaf sequences below this node */
  Matrix *qmatrix;     /* Matrix for up branch evolution */
  Matrix *pvector;     /* Probability distribution at top of up branch */
  Matrix *nucvector;   /* Nucleotide distribution at node */
} SNodeInfo;

typedef struct tagDNodeInfo {
  Grammar *grammar;    /* The grammar */
  SuffixTree *stree;   /* Search tree */
  int *col;            /* Place for search string */
  int noseqs;          /* Number of leaf sequences below this node */
  int *seqs;           /* List of leaf sequences below this node */
  Matrix *qmatrix;     /* Matrix for up branch evolution */
  Matrix *pvector;     /* Probability distribution at top of up branch */
  Matrix *nucvector;   /* Nucleotide distribution at node */
  SuffixTree *singletree; /* Single column search tree */
} DNodeInfo;

typedef struct Init_sarg {
  Grammar *grammar;
  Sgrp *sgrp;
  Alignstr *alignstr;
} Init_sarg;

typedef struct Init_darg {
  Grammar *grammar;
  Dgrp *dgrp;
  Alignstr *alignstr;
} Init_darg;

static Phyl **sphyl = NULL;     /* Phylogenies for single groups */
static Phyl **dphyl = NULL;     /* Phylogenies for double groups */
static Edouble p_limit;

static Grammar *curr_grammar;
static SuffixTree *sfstree = NULL;
static SuffixTree *dfstree = NULL;

static int *sfcol;
static int *dfcol;

/*------------------------------------------------------------------*/

/* Allocates memory for Colprob for a window of length wlen. */

Colprob *MakeColprob(Grammar *grammar, int wlen)
{
  Colprob* colprob = (Colprob *)malloc(sizeof(Colprob));
  int s = grammar->sgrp->count;
  int d = grammar->dgrp->count;
  int i, j;

  colprob->wlen = wlen;
  colprob->grammar = grammar;
  colprob->sgl = (Edouble **)malloc(s * sizeof(Edouble *));
  colprob->dbl = (Edouble ***)malloc(d * sizeof(Edouble **));

  for (j = 0; j < s; j++)
    colprob->sgl[j] = (Edouble *)malloc(wlen * sizeof(Edouble));

  for (j = 0; j < d; j++)
    {
      colprob->dbl[j] = (Edouble **)malloc(wlen * sizeof(Edouble *));
      for (i = 0; i < wlen; i++)
	colprob->dbl[j][i] = (Edouble *)malloc( (wlen - i) * sizeof(Edouble));
    }

  return colprob;
}

/*------------------------------------------------------------------*/

/* Allocates memory for Colprob for a window of length wlen. */

void FreeColprob(Grammar *grammar, Colprob *colprob)
{
  int s = grammar->sgrp->count;
  int d = grammar->dgrp->count;
  int i, j;
  int wlen;

  wlen = colprob->wlen;

  for (j = 0; j < d; j++) {
    for (i = 0; i < wlen; i++)
      free(colprob->dbl[j][i]);
    free(colprob->dbl[j]);
  }

  for (j = 0; j < s; j++)
    free(colprob->sgl[j]);

  free(colprob->sgl);
  free(colprob->dbl);
  free(colprob);
}

/*------------------------------------------------------------------*/

/* This initializes colprob according to the sequence stored in
   seq. Pos 0 is the start of seq. */

void InitColprob(Colprob *colprob, Alignstr *alignstr)
{
  LListCounter *lcount = MakeCounter(colprob->grammar->sgrp, FIRST);
  Sgrp *sgrp;
  Dgrp *dgrp;
  int wl = colprob->wlen;
  int p;   /* position iterator */
  int w;   /* string length itererator */

  while ((sgrp = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      colprob->sgl[lcount->pos][p] =
	ScolProb(alignstr, lcount->pos, p);

  InitCounter(lcount, colprob->grammar->dgrp, FIRST);
  while ((dgrp = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      for (w = 0; w + p < wl; w++)
	colprob->dbl[lcount->pos][p][w] =
	  DcolProb(alignstr, lcount->pos, p, w);

  /* New stuff */

#ifndef ETYPE
  InitCounter(lcount, colprob->grammar->dgrp, FIRST);
  while ((dgrp = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      for (w = 0; w + p < wl; w++) {
	colprob->dbl[lcount->pos][p][w] /= colprob->sgl[0][p]/fmul;
	colprob->dbl[lcount->pos][p][w] /= colprob->sgl[0][p+w]/fmul;
      }

  InitCounter(lcount, colprob->grammar->dgrp, LAST);
  while ((sgrp = Prev(lcount)) != NULL)
    for (p = 0; p < wl; p++) {
      colprob->sgl[lcount->pos][p] /= colprob->sgl[0][p]/fmul;
    }
#endif

  /* End new stuff */

  free(lcount);
}

/*------------------------------------------------------------------*/

/* This moves colprob one position with window rightmost (new)
   terminal term has postion pos. Pos 0 is the start of seq. */

void MoveColprob(Colprob *colprob, Alignstr *alignstr, int pos)
{
  LListCounter *lcount = MakeCounter(colprob->grammar->sgrp, FIRST);
  Sgrp *sgrp;
  Dgrp *dgrp;
  int wl = colprob->wlen;
  int w;   /* string length itererator */
  
  while ((sgrp = Next(lcount)) != NULL)
    colprob->sgl[lcount->pos][pos % wl]
      = ScolProb(alignstr, lcount->pos, pos);

  InitCounter(lcount, colprob->grammar->dgrp, FIRST);
  while ((dgrp = Next(lcount)) != NULL)
    for (w = 0; w < wl; w++)
      colprob->dbl[lcount->pos][(pos-w)%(wl-w)][w]
	= DcolProb(alignstr, lcount->pos, pos-w, w);

  free(lcount);
}

/*------------------------------------------------------------------*/

/* Allocates memory for Colprob for a window of length wlen. */

FColprob *MakeFColprob(Grammar *grammar, int wlen)
{
  FColprob* colprob = (FColprob *)malloc(sizeof(FColprob));
  int s = grammar->sgrp->count;
  int d = grammar->dgrp->count;
  int i, j;

  colprob->wlen = wlen;
  colprob->grammar = grammar;
  colprob->sgl = (Fdouble **)malloc(s * sizeof(Fdouble *));
  colprob->dbl = (Fdouble ***)malloc(d * sizeof(Fdouble **));

  for (j = 0; j < s; j++)
    colprob->sgl[j] = (Fdouble *)malloc(wlen * sizeof(Fdouble));

  for (j = 0; j < d; j++)
    {
      colprob->dbl[j] = (Fdouble **)malloc(wlen * sizeof(Fdouble *));
      for (i = 0; i < wlen; i++)
	colprob->dbl[j][i] = (Fdouble *)malloc( (wlen - i) * sizeof(Fdouble));
    }

  return colprob;
}

/*------------------------------------------------------------------*/

/* Allocates memory for Colprob for a window of length wlen. */

void FreeFColprob(Grammar *grammar, FColprob *colprob)
{
  int s = grammar->sgrp->count;
  int d = grammar->dgrp->count;
  int i, j;
  int wlen;

  wlen = colprob->wlen;

  for (j = 0; j < d; j++) {
    for (i = 0; i < wlen; i++)
      free(colprob->dbl[j][i]);
    free(colprob->dbl[j]);
  }

  for (j = 0; j < s; j++)
    free(colprob->sgl[j]);

  free(colprob->sgl);
  free(colprob->dbl);
  free(colprob);
}

/*------------------------------------------------------------------*/

/* This initializes colprob according to the sequence stored in
   seq. Pos 0 is the start of seq. */

void InitFColprob(FColprob *colprob, Alignstr *alignstr)
{
  LListCounter *lcount = MakeCounter(colprob->grammar->sgrp, FIRST);
  Sgrp *sgrp;
  Dgrp *dgrp;
  int wl = colprob->wlen;
  int p;   /* position iterator */
  int w;   /* string length itererator */

  while ((sgrp = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      colprob->sgl[lcount->pos][p] =
	FScolProb(alignstr, lcount->pos, p);

  InitCounter(lcount, colprob->grammar->dgrp, FIRST);
  while ((dgrp = Next(lcount)) != NULL)
    for (p = 0; p < wl; p++)
      for (w = 0; w + p < wl; w++)
	colprob->dbl[lcount->pos][p][w] =
	  FDcolProb(alignstr, lcount->pos, p, w);
  free(lcount);
}

/*------------------------------------------------------------------*/

/* This moves colprob one position with window rightmost (new)
   terminal term has postion pos. Pos 0 is the start of seq. */

void MoveFColprob(FColprob *colprob, Alignstr *alignstr, int pos)
{
  LListCounter *lcount = MakeCounter(colprob->grammar->sgrp, FIRST);
  Sgrp *sgrp;
  Dgrp *dgrp;
  int wl = colprob->wlen;
  int w;   /* string length itererator */
  
  while ((sgrp = Next(lcount)) != NULL)
    colprob->sgl[lcount->pos][pos % wl]
      = FScolProb(alignstr, lcount->pos, pos);

  InitCounter(lcount, colprob->grammar->dgrp, FIRST);
  while ((dgrp = Next(lcount)) != NULL)
    for (w = 0; w < wl; w++)
      colprob->dbl[lcount->pos][(pos-w)%(wl-w)][w]
	= FDcolProb(alignstr, lcount->pos, pos-w, w);

  free(lcount);
}

/* ------------------------------------------------------------ */

int SeqNumber(Align *align, char *name);
static void singlepointer(PhylNode *pnode, PhylNode *pnode_single);
static void initsphyl(PhylNode *pnode, void *isa);
static void initdphyl(PhylNode *pnode, void *ida);
static void freesphyl(PhylNode *pnode, void *arg);
static void freedphyl(PhylNode *pnode, void *arg);

/* Initiates phylogenies for the groups */
void InitCol(Phyl *phyl, Grammar *grammar, int wlen,
	     Alignstr *alignstr, double logp_limit)
{
  LListCounter *sgrpcnt;
  LListCounter *dgrpcnt;
  int i, j, k;
  Sgrp *sgrp;
  Dgrp *dgrp;
  Matrix *equilibrium;
  int numterm, numterm2;
  Init_sarg *isa;
  Init_darg *ida;

  p_limit = ExpEdouble(Dbl2Edbl(logp_limit*log(10.)));

  /* Single groups */

  isa = (Init_sarg *)malloc(sizeof(Init_sarg));

  isa->grammar = grammar;
  isa->alignstr = alignstr;

  numterm = grammar->term->count;

  sphyl = (Phyl **)malloc(grammar->sgrp->count * sizeof(Phyl *));

  sgrpcnt = MakeCounter(grammar->sgrp, FIRST);
  for (i = 0; (sgrp = Next(sgrpcnt)) != NULL; i++) {
    /* Initiate phylogeny */
    sphyl[i] = CopyPhyl(phyl);
    isa->sgrp = sgrp;
    PostOrderTraversePhyl(sphyl[i], initsphyl, (void *)isa);

    /* Replace root matrix */
    equilibrium = ZeroMatrix(numterm, numterm);
    for (j = 0; j < numterm; j++)
      equilibrium->entry[j][0] = sgrp->freq->entry[0][j];

    FreeMatrix(((SNodeInfo *)(sphyl[i]->root->elm))->qmatrix);
    ((SNodeInfo *)sphyl[i]->root->elm)->qmatrix = equilibrium;
  }

  free(isa);

  /* Double groups */

  ida = (Init_darg *)malloc(sizeof(Init_darg));

  ida->grammar = grammar;
  ida->alignstr = alignstr;

  numterm2 = numterm * numterm;

  dphyl = (Phyl **)malloc(grammar->dgrp->count * sizeof(Phyl *));

  dgrpcnt = MakeCounter(grammar->dgrp, FIRST);
  for (i = 0; (dgrp = Next(dgrpcnt)) != NULL; i++) {
    /* Initiate phylogeny */
    dphyl[i] = CopyPhyl(phyl);
    ida->dgrp = dgrp;
    PostOrderTraversePhyl(dphyl[i], initdphyl, (void *)ida);

    singlepointer(dphyl[i]->root, sphyl[0]->root);

    /* Replace root matrix */
    equilibrium = ZeroMatrix(numterm2, numterm2);
    for (j = 0; j < numterm2; j++)
      equilibrium->entry[j][0] = dgrp->freq->entry[0][j];

    FreeMatrix(((DNodeInfo *)(dphyl[i]->root->elm))->qmatrix);
    ((DNodeInfo *)dphyl[i]->root->elm)->qmatrix = equilibrium;
  }

  free(ida);

  curr_grammar = grammar;

  sfstree = InitSuffixSearch(grammar->term->count+1);
  dfstree = InitSuffixSearch(grammar->term->count+1);

  sfcol = (int *)malloc((alignstr->align->numseq+1) * sizeof(int));
  dfcol = (int *)malloc((2*alignstr->align->numseq+1) * sizeof(int));
}

void singlepointer(PhylNode *pnode, PhylNode *pnode_single)
{
  PhylNode *child, *child_single;

  ((DNodeInfo *)(pnode->elm))->singletree =
    ((SNodeInfo *)(pnode_single->elm))->stree;

  for (child = Child(pnode), child_single = Child(pnode_single);
       child != NULL;
       child = Brother(child), child_single = Brother(child_single))
    singlepointer(child, child_single);
}

/* Frees phylogenies for the groups etc.*/
void FinishCol(Phyl *phyl, Grammar *grammar)
{
  LListCounter *sgrpcnt;
  LListCounter *dgrpcnt;
  int i, j;

  /* Single groups */

  sgrpcnt = MakeCounter(grammar->sgrp, FIRST);
  for (i = 0; Next(sgrpcnt) != NULL; i++) {
    PostOrderTraversePhyl(sphyl[i], freesphyl, NULL);
  }

  /* Double groups */

  dgrpcnt = MakeCounter(grammar->dgrp, FIRST);
  for (i = 0; Next(dgrpcnt) != NULL; i++) {
    PostOrderTraversePhyl(dphyl[i], freedphyl, NULL);
  }
}

/* Sets names and makes room for SNodeInfo */
static void initsphyl(PhylNode *pnode, void *isa)
{
  Grammar *grammar;
  Alignstr *alignstr;
  Sgrp *sgrp;
  char *name;
  SNodeInfo *pn_info;
  PhylNode *child;
  int i, j;

  /* The sgrp of this phylogeny */
  sgrp = ((Init_sarg *)isa)->sgrp;

  /* The grammar of this phylogeny */
  grammar = ((Init_sarg *)isa)->grammar;

  /* The alignment */
  alignstr = ((Init_sarg *)isa)->alignstr;

  /* Read name */
  name = (char *)pnode->elm;

  /* Make room for node info */
  pnode->elm = malloc(sizeof(SNodeInfo));
  pn_info = (SNodeInfo *)(pnode->elm);
  pn_info->grammar = grammar;

  /* set sequence numbers */
  if (Child(pnode) == NULL) { /* Leaf */
    pn_info->noseqs = 1;
    pn_info->seqs = (int *)malloc(1 * sizeof(int));      
    pn_info->seqs[0] = SeqNumber(alignstr->align, name);
    if (pn_info->seqs[0] == -1) {
      fprintf(stderr, "Sequence not found: '%s'\n", name);
      exit(1);
    }
  }
  else { /* Internal node */
    /* Count the number of sequences below this one */
    pn_info->noseqs = 0;
    for (child = Child(pnode); child != NULL; child = Brother(child))
      pn_info->noseqs += ((SNodeInfo *)(child->elm))->noseqs;

    /* Make room */
    pn_info->seqs = (int *)malloc(pn_info->noseqs * sizeof(int));

    /* Insert them in list */
    for (child = Child(pnode), i = 0; child != NULL; child = Brother(child))
      for (j = 0; j < ((SNodeInfo *)(child->elm))->noseqs; j++, i++)
	pn_info->seqs[i] = ((SNodeInfo *)(child->elm))->seqs[j];
  }

  /* Make room for search tree */
  pn_info->stree = InitSuffixSearch(grammar->term->count+1);

  /* Make room only once for search strings */
  pn_info->col = (int *)malloc((pn_info->noseqs+1) * sizeof(int));

  /* Set matrix for up branch */
  pn_info->qmatrix = 
    TransposeMatrix(ExpMatrix(pnode->uplen, sgrp->eigen,
	      sgrp->diag, sgrp->inveigen));

  /* Make nucleotide and probability matricces */
  pn_info->nucvector = MakeMatrix(1, sgrp->freq->cols);
  pn_info->pvector = NULL;
}

/* Sets names and makes room for DNodeInfo */
static void initdphyl(PhylNode *pnode, void *ida)
{
  Grammar *grammar;
  Alignstr *alignstr;
  Dgrp *dgrp;
  char *name;
  DNodeInfo *pn_info;
  PhylNode *child;
  int i, j;

  /* The dgrp of this phylogeny */
  dgrp = ((Init_darg *)ida)->dgrp;

  /* The grammar of this phylogeny */
  grammar = ((Init_darg *)ida)->grammar;

  /* The alignment */
  alignstr = ((Init_darg *)ida)->alignstr;

  /* Read name */
  name = (char *)pnode->elm;

  /* Make room for node info */
  pnode->elm = malloc(sizeof(DNodeInfo));
  pn_info = (DNodeInfo *)(pnode->elm);
  pn_info->grammar = grammar;

  /* set sequence numbers */
  if (Child(pnode) == NULL) { /* Leaf */
    pn_info->noseqs = 1;
    pn_info->seqs = (int *)malloc(1 * sizeof(int));      
    pn_info->seqs[0] = SeqNumber(alignstr->align, name);
    if (pn_info->seqs[0] == -1) {
      fprintf(stderr, "Sequence not found: '%s'\n", name);
      exit(1);
    }
  }
  else { /* Internal node */
    /* Count the number of sequences below this one */
    pn_info->noseqs = 0;
    for (child = Child(pnode); child != NULL; child = Brother(child))
      pn_info->noseqs += ((DNodeInfo *)(child->elm))->noseqs;

    /* Make room */
    pn_info->seqs = (int *)malloc(pn_info->noseqs * sizeof(int));

    /* Insert them in list */
    for (child = Child(pnode), i = 0; child != NULL; child = Brother(child))
      for (j = 0; j < ((DNodeInfo *)(child->elm))->noseqs; j++, i++)
	pn_info->seqs[i] = ((DNodeInfo *)(child->elm))->seqs[j];
  }

  /* Make room for search tree */
  pn_info->stree = InitSuffixSearch(grammar->term->count+1);

  /* Make room only once for search strings */
  pn_info->col = (int *)malloc((2*pn_info->noseqs+1) * sizeof(int));

  /* Set matrix for up branch */
  pn_info->qmatrix = 
    TransposeMatrix(ExpMatrix(pnode->uplen, dgrp->eigen,
	      dgrp->diag, dgrp->inveigen));

  /* Make nucleotide and probability matricces */
  pn_info->nucvector = MakeMatrix(1, dgrp->freq->cols);
  pn_info->pvector = NULL;
}

static void freesphyl(PhylNode *pnode, void *arg)
{
  SNodeInfo *pn_info;

  pn_info = (SNodeInfo *)(pnode->elm);

  FinishSuffixSearch(pn_info->stree);
}

static void freedphyl(PhylNode *pnode, void *arg)
{
  DNodeInfo *pn_info;

  pn_info = (DNodeInfo *)(pnode->elm);

  FinishSuffixSearch(pn_info->stree);
}

/* ------------------------------------------------------------ */

Alignstr *Removegap(Grammar *grammar, Alignstr *alignstr, double limit)
{
  int i, j;
  Alignstr *newalignstr;
  int newpos;
  int numgap;
  int gap;
  int *map;
  int *localmap;
  int *invlocalmap;

  gap = FindSym('-', grammar);

  newalignstr = (Alignstr *)malloc(sizeof(Alignstr));

  newalignstr->align = MakeAlign(alignstr->align->numseq);

  for (j = 0; j < alignstr->align->numseq; j++) {
    newalignstr->align->name[j] = (char *)malloc(
	   (strlen(alignstr->align->name[j])+1)
	   * sizeof(char));
    strcpy(newalignstr->align->name[j], alignstr->align->name[j]);
  }

  newpos = 0;
  for (i = 0; i < alignstr->align->len; i++) {
    numgap = 0;
    for (j = 0; j < alignstr->align->numseq; j++) {
      if (alignstr->align->seq[j][i] == gap)
	numgap++;
    }
    if (1.-((double) numgap)/alignstr->align->numseq > limit)
      newpos++;
  }

  newalignstr->align->len = newpos;
  map = (int *)malloc((newpos+1) * sizeof(int));
  localmap = (int *)malloc(newpos * sizeof(int));
  invlocalmap = (int *)malloc(alignstr->align->len * sizeof(int));

  map[newpos] = alignstr->map[alignstr->align->len];

  for (j = 0; j < alignstr->align->numseq; j++)
    newalignstr->align->seq[j] = (int *)malloc(newpos * sizeof(int));
  newalignstr->str = (int *)malloc(newpos * sizeof(int));

  newpos = 0;
  for (i = 0; i < alignstr->align->len; i++) {
    invlocalmap[i] = -1;
    numgap = 0;
    for (j = 0; j < alignstr->align->numseq; j++) {
      if (alignstr->align->seq[j][i] == gap)
	numgap++;
    }
    if (1.-((double) numgap)/alignstr->align->numseq > limit) {
      map[newpos] = alignstr->map[i];
      localmap[newpos] = i;
      invlocalmap[i] = newpos;
      for (j = 0; j < alignstr->align->numseq; j++)
	newalignstr->align->seq[j][newpos] = alignstr->align->seq[j][i];
      newpos++;
    }
  }

  for (i = 0; i < newpos; i++) {
    if (alignstr->str[localmap[i]] < 0)
      newalignstr->str[i] = alignstr->str[localmap[i]];
    else
      newalignstr->str[i] = invlocalmap[alignstr->str[localmap[i]]];
  }

  newalignstr->map = map;

  return newalignstr;
}

/* ------------------------------------------------------------ */

void scolrec(PhylNode *pnode, Alignstr *alignstr, int pos);

/* Return column probability for single group */
Edouble ScolProb(Alignstr *alignstr, int sgrpno, int pos)
{
  if (alignstr->str[pos] == DOUBLE)
    return Dbl2Edbl(0.);

  if (alignstr->str[pos] >= 0)
    return Dbl2Edbl(0.);

  scolrec(sphyl[sgrpno]->root, alignstr, pos);

  return ((SNodeInfo *)(sphyl[sgrpno]->root->elm))->pvector->entry[0][0];
}

void scolrec(PhylNode *pnode, Alignstr *alignstr, int pos)
{
  SuffixNode *SNode;
  int i;
  int numterm;
  int do_search;
  PhylNode *child;
  int symno;
  Matrix *nucvector, *mat;
  int size;
  SNodeInfo *pn_info;
  Edouble sum;
  int noseqs;

  pn_info = (SNodeInfo *)(pnode->elm);
  noseqs = pn_info->noseqs;

  numterm = pn_info->grammar->term->count;

  do_search = 1;
  for (i = 0; i < noseqs; i++) {
    pn_info->col[i] = alignstr->align->seq[pn_info->seqs[i]][pos];
    if (pn_info->col[i] >= numterm+1)
      do_search = 0;
  }
  pn_info->col[i] = -1;

  if (do_search == 1) {
    SNode = SuffixSearch(pn_info->stree, pn_info->col);

    if (SNode->elm != NULL) {
      SetMatrix(pn_info->pvector, (Matrix *)SNode->elm);
      return;
    }
  }

  /* The calculation has not been done before */

  /* Do calculations for children */
  for (child = Child(pnode); child != NULL; child = Brother(child))  
    scolrec(child, alignstr, pos);

  nucvector = pn_info->nucvector;

  /* Set nucleotide probability vector */

  if (Child(pnode) != NULL) { /* Internal node */
    for (i = 0; i < numterm; i++)
      nucvector->entry[0][i] = Dbl2Edbl(1.);
  }
  else { /* Leaf */
    symno = alignstr->align->seq[pn_info->seqs[0]][pos];
    for (i = 0; i < numterm; i++)
      nucvector->entry[0][i] = pn_info->grammar->quickdist[symno][i];
  }

  /* Calculate evolution */
  size = nucvector->cols;

  for (child = Child(pnode); child != NULL; child = Brother(child))
    for (i = 0; i < size; i++)
      MulEdouble(&nucvector->entry[0][i],
	((SNodeInfo *)(child->elm))->pvector->entry[0][i]);

  FreeMatrix(pn_info->pvector);

  pn_info->pvector = MulMatrix(nucvector, pn_info->qmatrix);

  /* Insert result in search tree */
  if (do_search == 1) {
    sum = Dbl2Edbl(0.);
    AddEdouble(&sum, p_limit);
    for (i = 0; i < pn_info->pvector->cols; i++) {
      SubEdouble(&sum, pn_info->pvector->entry[0][i]);
    }
    if (!PositiveEdouble(sum)) {
      /* only insert if sufficiently likely */
      mat = CopyMatrix(pn_info->pvector);
      SuffixInsert(pn_info->stree, pn_info->col, (void *)mat);
    }
  }
}

void dcolrec(PhylNode *pnode, Alignstr *alignstr, int pos1, int dist);

/* Return column probability for double group */
Edouble DcolProb(Alignstr *alignstr, int dgrpno, int pos1, int dist)
{
  if (alignstr->str[pos1] == SINGLE || alignstr->str[pos1+dist] == SINGLE) {
    return Dbl2Edbl(0.);
  }

  if ((alignstr->str[pos1] >= 0 || alignstr->str[pos1+dist] >= 0) &&
      alignstr->str[pos1] != pos1+dist)
    return Dbl2Edbl(0.);

  if (dist < curr_grammar->mindist)
    return Dbl2Edbl(0.);

  dcolrec(dphyl[dgrpno]->root, alignstr, pos1, dist);

  return ((DNodeInfo *)(dphyl[dgrpno]->root->elm))->pvector->entry[0][0];
}

void dcolrec(PhylNode *pnode, Alignstr *alignstr, int pos1, int dist)
{
  SuffixNode *SNode;
  int i, j;
  int numterm, numterm2;
  int do_search;
  PhylNode *child;
  int symno1, symno2;
  Matrix *nucvector, *mat;
  int size;
  DNodeInfo *pn_info;
  int temp;
  int noseqs;

  pn_info = (DNodeInfo *)(pnode->elm);
  noseqs = pn_info->noseqs;

  numterm = pn_info->grammar->term->count;
  numterm2 = numterm*numterm;

  do_search = 1;
  for (i = 0; i < noseqs; i++) {
    pn_info->col[i] = alignstr->align->seq[pn_info->seqs[i]][pos1];
    pn_info->col[i+noseqs] = alignstr->align->seq[pn_info->seqs[i]][pos1+dist];
    if (pn_info->col[i] >= numterm+1 || pn_info->col[i+noseqs] >= numterm+1)
      do_search = 0;
  }
  pn_info->col[i*2] = -1;

  if (do_search == 1) {
    SNode = SuffixSearch(pn_info->stree, pn_info->col);

    if (SNode->elm != NULL) {
      SetMatrix(pn_info->pvector, (Matrix *)SNode->elm);
      return;
    }
  }

  /* The calculation has not been done before */

  /* Do calculations for children */
  for (child = Child(pnode); child != NULL; child = Brother(child))  
    dcolrec(child, alignstr, pos1, dist);

  nucvector = pn_info->nucvector;

  /* Set nucleotide probability vector */

  if (Child(pnode) != NULL) { /* Internal node */
    for (i = 0; i < numterm2; i++)
      nucvector->entry[0][i] = Dbl2Edbl(1.);
  }
  else { /* Leaf */
    symno1 = alignstr->align->seq[pn_info->seqs[0]][pos1];
    symno2 = alignstr->align->seq[pn_info->seqs[0]][pos1+dist];
    for (i = 0; i < numterm; i++)
      for (j = 0; j < numterm; j++)
	nucvector->entry[0][i+numterm*j]
	  = ProdEdouble(pn_info->grammar->quickdist[symno1][i],
			pn_info->grammar->quickdist[symno2][j]);
  }

  /* Calculate evolution */
  size = nucvector->cols;

  for (child = Child(pnode); child != NULL; child = Brother(child))
    for (i = 0; i < size; i++)
      MulEdouble(&nucvector->entry[0][i],
	((DNodeInfo *)(child->elm))->pvector->entry[0][i]);

  FreeMatrix(pn_info->pvector);

  pn_info->pvector = MulMatrix(nucvector, pn_info->qmatrix);

  /* Insert result in search tree */
  if (do_search == 1) {
    SNode = SuffixSearch(pn_info->singletree, pn_info->col+noseqs);

    if (SNode->elm != NULL) {
      temp = pn_info->col[noseqs];
      pn_info->col[noseqs] = -1;
      
      SNode = SuffixSearch(pn_info->singletree, pn_info->col);

      pn_info->col[noseqs] = temp;

      if (SNode->elm != NULL) {
	/* Both columns found in single tree, insert in this tree */
	mat = CopyMatrix(pn_info->pvector);
	SuffixInsert(pn_info->stree, pn_info->col, (void *)mat);
      }
    }
  }
}


/* Return column probability for single group */
Fdouble FScolProb(Alignstr *alignstr, int sgrpno, int pos)
{
  SuffixNode *SNode;
  int i;
  int do_search;
  Fdouble fd;
  Fdouble *fp;
  int numterm;

  if (alignstr->str[pos] == DOUBLE)
    return Dbl2Fdbl(0.);

  if (alignstr->str[pos] >= 0)
    return Dbl2Fdbl(0.);

  numterm = curr_grammar->term->count;

  do_search = 1;
  for (i = 0; i < alignstr->align->numseq; i++) {
    sfcol[i] = alignstr->align->seq[i][pos];
    if (sfcol[i] >= numterm+1)
      do_search = 0;
  }
  sfcol[i] = -1;

  if (do_search == 1) {
    SNode = SuffixSearch(sfstree, sfcol);
    if (SNode->elm != NULL)
      return *((Fdouble *)(SNode->elm));
  }

  fd = Edbl2Fdbl(ScolProb(alignstr, sgrpno, pos));

  if (do_search == 1) {
    fp = (Fdouble *)malloc(sizeof(Fdouble));
    *fp = fd;
    SuffixInsert(sfstree, sfcol, (void *)fp);
  }

  return fd;
}

/* Return column probability for double group */
Fdouble FDcolProb(Alignstr *alignstr, int dgrpno, int pos1, int dist)
{
  SuffixNode *SNode;
  int i;
  int do_search;
  Fdouble fd;
  Fdouble *fp;
  int numterm;

  if (alignstr->str[pos1] == SINGLE || alignstr->str[pos1+dist] == SINGLE)
    return Dbl2Fdbl(0.);

  if ((alignstr->str[pos1] >= 0 || alignstr->str[pos1+dist] >= 0) &&
      alignstr->str[pos1] != pos1+dist)
    return Dbl2Fdbl(0.);

  if (dist < curr_grammar->mindist)
    return Dbl2Fdbl(0.);

  numterm = curr_grammar->term->count;

  do_search = 1;
  for (i = 0; i < alignstr->align->numseq; i++) {
    dfcol[i*2] = alignstr->align->seq[i][pos1];
    dfcol[i*2+1] = alignstr->align->seq[i][pos1+dist];
    if (dfcol[i*2] >= numterm+1 || dfcol[i*2+1] >= numterm+1)
      do_search = 0;
  }
  dfcol[i*2] = -1;

  if (do_search == 1) {
    SNode = SuffixSearch(dfstree, dfcol);
    if (SNode->elm != NULL)
      return *((Fdouble *)(SNode->elm));
  }

  fd = Edbl2Fdbl(DcolProb(alignstr, dgrpno, pos1, dist));

  if (do_search == 1) {
    fp = (Fdouble *)malloc(sizeof(Fdouble));
    *fp = fd;
    SuffixInsert(dfstree, dfcol, (void *)fp);
  }

  return fd;
}
