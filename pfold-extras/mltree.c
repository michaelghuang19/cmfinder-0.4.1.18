#include "file.h"
#include "phyl.h"

typedef struct NodeInfo {
  int seqnum;               /* Sequence at this node */
  Matrix *pmatrix;          /* Matrix for up branch evolution */
  Matrix **p_top_up;         /* Probability distribution of tree above
			       top of up branch */
  Matrix **p_top_down;       /* Probability distribution of tree below
			       top of up branch */
  Matrix **p_bot_up;         /* Probability distribution of tree above
			       bottom of up branch */
  Matrix **p_bot_down;       /* Probability distribution of tree below
			       bottom of up branch */
  Grammar *grammar;
} NodeInfo;

void usage(void);
void initphyl(PhylNode *pnode, Grammar *grammar, Align *align);
void setmat(PhylNode *pnode, Grammar *grammar);
void downcalc(PhylNode *pnode, Align *align);
void upcalc(PhylNode *pnode, Align *align);
double optimize(PhylNode *pnode, Align *align, Grammar *grammar, Phyl *phyl);
PhylNode *findnode(PhylNode *pnode, int number);

int main(int argc, char **argv)
{
  FILE *fp, *grammarfp;
  Grammar *grammar;
  Header *header;
  Entry *entry, *entry_tree;
  Entry **entry_list;
  Align *align;
  int read_error;   /* For keeping track of errors in reading entries */
  Phyl *phyl;
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */
  int size;
  int i;
  double nll, oldnll;
  int option_quiet;

  cmdarg = InitArgument(argc, argv);

  option_quiet = 0;

  while ((s = GetArgument(cmdarg)) != NULL) {
    if (strcmp(s, "-quiet") == 0) {
      option_quiet = 1;
    }
    else {
      usage();
      return 1;
    }
  }

  /* Open ratefile */

  if ((s = GetFilename(cmdarg)) == NULL) {
    usage();
    return 1;
  }
  else if ((grammarfp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "mltree: Error in opening file '%s'\n", s);
    return 1; }

  /* Open col file */

  if ((s = GetFilename(cmdarg)) == NULL)
    fp = stdin;
  else if ((fp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "mltree: Error in opening file '%s'\n", s);
    return 1; }

  header = MakeHeader();

  if (ReadHeader(fp, header) != 0)
    return 1;

  PrintHeader(stdout, header);

  entry = MakeEntry();

  phyl = NULL;

  /* Read tree first */
  while ((read_error = ReadEntry(fp, entry)) == 0) {
    if (!ReadType(entry, "TREE")) 
      continue;

    phyl = ReadColPhyl(entry);
    break;

    PrintEntry(stdout, entry);
  }

  entry_tree = entry;

  /* Read the sequences */
  size = 1;

  entry_list = (Entry **)malloc(size*sizeof(Entry *));

  for (i = 0; (entry = MakeEntry()) != NULL &&
	      (read_error = ReadEntry(fp, entry)) == 0; i++) {
    if (i >= size) {
      size *= 2;
      entry_list = (Entry **)realloc(entry_list, size*sizeof(Entry *));
    }
    entry_list[i] = entry;
  }

  if (fp != stdin && fclose(fp) != 0) {
    fprintf(stderr, "mltree: Error in closing col file\n");
    return 1; }

  if (i+1 != size)
    entry_list = (Entry **)realloc(entry_list, (i+1)*sizeof(Entry *));

  entry_list[i] = NULL;

  if (read_error == 1)
    exit(1);

  grammar = ReadGrammar(grammarfp);

  align = Col2Align(entry_list, (int (*)(char, void *))FindSym, (void *)grammar);

  if (option_quiet == 0) fprintf(stdout, "Init\n");

  initphyl(phyl->root, grammar, align);
  
  oldnll = -1;
  nll = -1;
  while (oldnll == -1 || oldnll-nll > 0.1) {
    oldnll = nll;
    if (option_quiet == 0) fprintf(stdout, "Setting matrices\n");
    setmat(phyl->root, grammar);
    
    if (option_quiet == 0) fprintf(stdout, "Down calculations\n");
    downcalc(phyl->root, align);
    
    if (option_quiet == 0) fprintf(stdout, "Up calculations\n");
    upcalc(phyl->root, align);
    
    if (option_quiet == 0) fprintf(stdout, "Optimize\n");
    nll = optimize(phyl->root, align, grammar, phyl);
    if (option_quiet == 0) fprintf(stdout, "  NLL: %f\n", nll);
  }

  UpdatePhylEntry(phyl, entry_tree);
  PrintEntry(stdout, entry_tree);

  for (i = 0; entry_list[i] != NULL; i++)
    PrintEntry(stdout, entry_list[i]);

  return 0;
}

void usage(void)
{
  fprintf(stderr,
	  "usage: mltree <ratefile> [<file>]\n");
}

/*
   Sets alignment names, and numbers.
*/
void initphyl(PhylNode *pnode, Grammar *grammar, Align *align)
{
  PhylNode *child;
  NodeInfo *pn_info;
  char *name;
  int numterm;
  int len;
  int i;

  numterm = grammar->term->count;

  name = (char *)pnode->elm;

  pnode->elm = (NodeInfo *)malloc(sizeof(NodeInfo));

  if (pnode->uplen < 0.001)
    pnode->uplen = 0.001;

  pn_info = (NodeInfo *)pnode->elm;
  len = align->len;

  pn_info->seqnum = SeqNumber(align, name);
  pn_info->grammar = grammar;
  pn_info->pmatrix = NULL;

  pn_info->p_bot_down = (Matrix **)malloc(len * sizeof(Matrix *));
  for (i = 0; i < len; i ++)
    pn_info->p_bot_down[i] = MakeMatrix(1, numterm);

  pn_info->p_top_down = (Matrix **)malloc(len * sizeof(Matrix *));
  for (i = 0; i < len; i ++)
    pn_info->p_top_down[i] = NULL;

  pn_info->p_bot_up = (Matrix **)malloc(len * sizeof(Matrix *));
  for (i = 0; i < len; i ++)
    pn_info->p_bot_up[i] = NULL;

  pn_info->p_top_up = (Matrix **)malloc(len * sizeof(Matrix *));
  for (i = 0; i < len; i ++)
    pn_info->p_top_up[i] = MakeMatrix(1, numterm);

  for (child = pnode->child; child != NULL; child = child->brother)
    initphyl(child, grammar, align);
}

/*
   Calculates evolutionary matrices.
*/
void setmat(PhylNode *pnode, Grammar *grammar)
{
  Sgrp *sgrp;
  LListCounter *lcount;
  PhylNode *child;
  NodeInfo *pn_info;
  Matrix *temp;
  int seqnum;
  int numterm;
  int number;

  number = pnode->number;
  pn_info = ((NodeInfo *)pnode->elm);
  seqnum = pn_info->seqnum;
  numterm = pn_info->grammar->term->count;

  lcount = MakeCounter(grammar->sgrp, FIRST);
  sgrp = Next(lcount);
  free(lcount);

  FreeMatrix(pn_info->pmatrix);
  pn_info->pmatrix = 
    TransposeMatrix(temp = ExpMatrix(pnode->uplen, sgrp->eigen,
				     sgrp->diag, sgrp->inveigen));
  FreeMatrix(temp);
  
  for (child = pnode->child; child != NULL; child = child->brother)
    setmat(child, grammar);
}

/* Calculate probability distributions downwards in the tree */
void downcalc(PhylNode *pnode, Align *align)
{
  PhylNode *child;
  NodeInfo *pn_info;
  int pos;
  int i;
  int numterm;
  int symno;

  pn_info = (NodeInfo *)(pnode->elm);

  numterm = pn_info->grammar->term->count;

  /* Do calculations for children */
  for (child = Child(pnode); child != NULL; child = Brother(child))  
    downcalc(child, align);

  for (pos = 0; pos < align->len; pos++) {
    /* Set nucleotide probability vector */
    if (pn_info->seqnum == -1) { /* Internal node */
      for (i = 0; i < numterm; i++)
	pn_info->p_bot_down[pos]->entry[0][i] = Dbl2Edbl(1.);
    }
    else { /* Sequence here */
      symno = align->seq[pn_info->seqnum][pos];
      for (i = 0; i < numterm; i++)
	pn_info->p_bot_down[pos]->entry[0][i]
	  = pn_info->grammar->quickdist[symno][i];
    }

    /* Calculate evolution */
    for (child = Child(pnode); child != NULL; child = Brother(child))
      for (i = 0; i < numterm; i++)
	MulEdouble(&pn_info->p_bot_down[pos]->entry[0][i],
		   ((NodeInfo *)(child->elm))->p_top_down[pos]->entry[0][i]);
    
    FreeMatrix(pn_info->p_top_down[pos]);

    pn_info->p_top_down[pos]
      = MulMatrix(pn_info->p_bot_down[pos], pn_info->pmatrix);
  }
}

/* Calculate probability distributions upwards in the tree */
void upcalc(PhylNode *pnode, Align *align)
{
  PhylNode *child;
  NodeInfo *pn_info;
  int pos;
  int i;
  int numterm;

  pn_info = (NodeInfo *)(pnode->elm);

  numterm = pn_info->grammar->term->count;

  for (pos = 0; pos < align->len; pos++) {

    /* Set nucleotide probability vector */
    for (i = 0; i < numterm; i++)
      pn_info->p_top_up[pos]->entry[0][i] = Dbl2Edbl(1.);

    /* Calculate evolution */
    if (pnode->parent != NULL) {
      for (i = 0; i < numterm; i++)
	MulEdouble(&pn_info->p_top_up[pos]->entry[0][i],
	     ((NodeInfo *)pnode->parent->elm)->p_bot_up[pos]->entry[0][i]);

      for (child = Child(pnode->parent); child != NULL; child = Brother(child))
	if (child != pnode) {
	  for (i = 0; i < numterm; i++)
	    MulEdouble(&pn_info->p_top_up[pos]->entry[0][i],
		    ((NodeInfo *)(child->elm))->p_top_down[pos]->entry[0][i]);
	}
    }

    FreeMatrix(pn_info->p_bot_up[pos]);

    pn_info->p_bot_up[pos]
      = MulMatrix(pn_info->p_top_up[pos], pn_info->pmatrix);
  }

  /* Do calculations for children */
  for (child = Child(pnode); child != NULL; child = Brother(child))  
    upcalc(child, align);
}

/* Optimize branch lengths */
double optimize(PhylNode *pnode, Align *align, Grammar *grammar, Phyl *phyl)
{
  PhylNode *child;
  NodeInfo *pn_info;
  LListCounter *lcount;
  Sgrp *sgrp;
  Edouble e, sum;
  Matrix *temp, *p_temp, *old, *mat, *new;
  int pos;
  int i;
  int numterm;
  double nll;
  double l, x;
  int error;
  int finish;

  /* Do calculations for children */
  for (child = Child(pnode); child != NULL; child = Brother(child))  
    nll = optimize(child, align, grammar, phyl);

  if (pnode->parent == NULL)
    return nll;

  lcount = MakeCounter(grammar->sgrp, FIRST);
  sgrp = Next(lcount);
  free(lcount);

  pn_info = (NodeInfo *)(pnode->elm);

  numterm = pn_info->grammar->term->count;

  p_temp = NULL;

  x = log(pnode->uplen+0.001);
  l = exp(x);
  finish = 0;
  InitMinimize(x-0.2, x+0.2, 0.001);
  while ((error = Minimize(&x, &nll)) == 0) {
    if (x < -10) {
      x = -10;
      finish = 10;
    }
    else if (x > 6) {
      x = 6;
      finish++;
    }
    if (finish > 4)
      break;
    l = exp(x);
    nll = 0;

    p_temp = 
      TransposeMatrix(temp = ExpMatrix(l, sgrp->eigen,
				       sgrp->diag, sgrp->inveigen));
    FreeMatrix(temp);

    for (pos = 0; pos < align->len; pos++) {
      temp = MulMatrix(pn_info->p_bot_down[pos], p_temp);
      sum = Dbl2Edbl(0.);
      for (i = 0; i < numterm; i++) {
	e = ProdEdouble(sgrp->freq->entry[0][i], 
			ProdEdouble(temp->entry[0][i],
				    pn_info->p_top_up[pos]->entry[0][i]));
	AddEdouble(&sum, e);
      }
      FreeMatrix(temp);
      nll -= Edbl2Dbl(LogEdouble(sum));
    }
    FreeMatrix(p_temp);
  }

  if (error == 1) {
    fprintf(stderr, "problems\n");
    exit(1);
  }

  p_temp = TransposeMatrix(temp = ExpMatrix(l, sgrp->eigen,
					    sgrp->diag, sgrp->inveigen));
  FreeMatrix(temp);

  pnode->uplen = exp(x);
  /* Adjust probabilities for relevant neighbours */

  if (pnode->parent != NULL) {
    for (pos = 0; pos < align->len; pos++) {
      old = pn_info->p_top_down[pos];
      new = MulMatrix(pn_info->p_bot_down[pos], p_temp);
      mat = ((NodeInfo *)pnode->parent->elm)->p_bot_down[pos];
      for (i = 0; i < numterm; i++) {
	if (Edbl2Dbl(old->entry[0][i]) != 0)
	  DivEdouble(&mat->entry[0][i], old->entry[0][i]);
	MulEdouble(&mat->entry[0][i], new->entry[0][i]);
      }

      for (child = Child(pnode->parent); child != NULL; child = Brother(child))
	if (child != pnode) {
	  mat = ((NodeInfo *)child->elm)->p_top_up[pos];
	  for (i = 0; i < numterm; i++) {
	    if (Edbl2Dbl(old->entry[0][i]) != 0)
	      DivEdouble(&mat->entry[0][i], old->entry[0][i]);
	    MulEdouble(&mat->entry[0][i], new->entry[0][i]);
	  }
	}
      FreeMatrix(new);
    }
  }

  FreeMatrix(p_temp);

  return nll;
}

PhylNode *findnode(PhylNode *pnode, int number)
{
  PhylNode *p, *child;

  if (pnode->number == number)
    return pnode;

  for (child = pnode->child; child != NULL; child = child->brother)
    if ((p = findnode(child, number)) != NULL)
      return p;

  return NULL;
}
