#include "grammar.h"

/* This allocates and initializes a Grammar */
Grammar *MakeGrammar()
{
  Grammar *grammar = (Grammar *)malloc(sizeof(Grammar));

  /* Make LLists */
  grammar->nont = MakeLList();
  grammar->term = MakeLList();
  grammar->iprt = MakeLList();
  grammar->trns = MakeLList();
  grammar->sgrp = MakeLList();
  grammar->dgrp = MakeLList();
  grammar->ruls = MakeLList();
  grammar->ruln = MakeLList();
  grammar->rulsn = MakeLList();
  grammar->rulns = MakeLList();
  grammar->ruldnd = MakeLList();
  grammar->ruldd = MakeLList();
  grammar->rulnn = MakeLList();
  grammar->quickdist = NULL;

  return grammar;
}

int readnont(FILE *fp, Grammar *grammar);
int readterm(FILE *fp, Grammar *grammar);
int readiprt(FILE *fp, Grammar *grammar);
int readtrns(FILE *fp, Grammar *grammar);
int readsgrp(FILE *fp, Grammar *grammar);
int readdgrp(FILE *fp, Grammar *grammar);
int readprod(FILE *fp, Grammar *grammar);

Grammar *ReadGrammar(FILE *fp)
{
  Grammar *grammar;
  char *s;      /* For headers */
  LListCounter *lcount;
  Iprt *elm;
  int numsym, numterm;
  int i, j;

  grammar = MakeGrammar();

  /* Get ready for reading connected lines */
  InitConnect();

  while ((s = GetConnect(fp)) != NULL) {
    if (StrCmp(s, "> Nonterminal:\n") == 0)
      if (readnont(fp, grammar) != 0)
	exit(1);
    if (StrCmp(s, "> Terminal:\n") == 0)
      if (readterm(fp, grammar) != 0)
	exit(1);
    if (StrCmp(s, "> Interpret:\n") == 0)
      if (readiprt(fp, grammar) != 0)
	exit(1);
    if (StrCmp(s, "> Translate:\n") == 0)
      if (readtrns(fp, grammar) != 0)
	exit(1);
    if (StrCmp(s, "> Single group:\n") == 0)
      if (readsgrp(fp, grammar) != 0)
	exit(1);
    if (StrCmp(s, "> Double group:\n") == 0)
      if (readdgrp(fp, grammar) != 0)
	exit(1);
    if (StrCmp(s, "> Production:\n") == 0)
      if (readprod(fp, grammar) != 0)
	exit(1);
  }

  numterm = grammar->term->count;
  numsym = numterm + grammar->iprt->count;

  grammar->quickdist = (Edouble **)malloc(numsym * sizeof(Edouble *));

  for (i = 0; i < numsym; i++)
    grammar->quickdist[i] = (Edouble *)malloc(numterm * sizeof(Edouble));

  /* Set terminal distributions */
  for (i = 0; i < numterm; i++) {
    for (j = 0; j < numterm; j++)
      grammar->quickdist[i][j] = Dbl2Edbl(0.);
    grammar->quickdist[i][i] = Dbl2Edbl(1.);
  }

  lcount = MakeCounter(grammar->iprt, FIRST);

  /* Set interpret distributions */
  for (i = numterm; (elm = Next(lcount)) != NULL; i++) {
    for (j = 0; j < numterm; j++)
      grammar->quickdist[i][j] = elm->freq->entry[0][j];
  }
  free(lcount);

  grammar->mindist = 1;

  return grammar;
}

int readnont(FILE *fp, Grammar *grammar)
{
  char *s;    /* Line with nonterminals */
  int ptr;    /* Pointer in line */
  Nont *nont; /* For last read nonterminal */

  if ((s = GetConnect(fp)) == NULL) {
    fprintf(stderr, "Error in reading nonterminals\n");
    return 1; }

  ptr = 0;

  for (;;) {
    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] == '\0')
      break;

    nont = (Nont *)malloc(sizeof(Nont));
    nont->sym = s[ptr++];
    Enqueue(grammar->nont, nont);
  }

  return 0;
}

int readterm(FILE *fp, Grammar *grammar)
{
  char *s;    /* Line with terminals */
  int ptr;    /* Pointer in line */
  Term *term; /* For last read terminal */

  if ((s = GetConnect(fp)) == NULL) {
    fprintf(stderr, "Error in reading terminals\n");
    return 1; }

  ptr = 0;

  for (;;) {
    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] == '\0')
      break;

    term = (Term *)malloc(sizeof(Term));
    term->sym = s[ptr++];
    Enqueue(grammar->term, term);
  }

  return 0;
}

Iprt *MakeIprt(Grammar *grammar);
Sgrp *MakeSgrp(Grammar *grammar);
Dgrp *MakeDgrp(Grammar *grammar);

Iprt *MakeIprt(Grammar *grammar)
{
  Iprt *iprt;
  int numterm;  /* Number of terminals */

  numterm = grammar->term->count;

  iprt = (Iprt *)malloc(sizeof(Iprt));
  iprt->freq = MakeMatrix(1, numterm);

  return iprt;
}
  
Sgrp *MakeSgrp(Grammar *grammar)
{
  Sgrp *sgrp;
  int numterm;  /* Number of terminals */

  numterm = grammar->term->count;

  sgrp = (Sgrp *)malloc(sizeof(Sgrp));
  sgrp->freq = MakeMatrix(1, numterm);
  sgrp->eigen = MakeMatrix(numterm, numterm);
  sgrp->diag = MakeMatrix(numterm, numterm);
  sgrp->inveigen = MakeMatrix(numterm, numterm);

  return sgrp;
}
  
Dgrp *MakeDgrp(Grammar *grammar)
{
  Dgrp *dgrp;
  int numterm2;  /* Number of terminals squared */

  numterm2 = grammar->term->count * grammar->term->count;

  dgrp = (Dgrp *)malloc(sizeof(Dgrp));
  dgrp->freq = MakeMatrix(1, numterm2);
  dgrp->eigen = MakeMatrix(numterm2, numterm2);
  dgrp->diag = MakeMatrix(numterm2, numterm2);
  dgrp->inveigen = MakeMatrix(numterm2, numterm2);

  return dgrp;
}
  

int readiprt(FILE *fp, Grammar *grammar)
{
  char *s;      /* Line with terminals */
  int ptr;      /* Pointer in line */
  Iprt *iprt;   /* For last read terminal */
  int numterm;  /* Number of terminals */
  int len;      /* Length of read numbers */
  int i;        /* For matrix entries */
  double d;     /* For reading doubles */

  numterm = grammar->term->count;

  while ((s = GetConnect(fp)) != NULL && s[0] != '>') {
    ptr = 0;

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] == '\0') {
      fprintf(stderr, "No symbol in reading interpret\n");
      return 1; }

    iprt = MakeIprt(grammar);
    iprt->sym = s[ptr++];
    Enqueue(grammar->iprt, iprt);

    for (i = 0; i < numterm; i++) {
      if (sscanf(s+ptr, "%lf%n", &d, &len) != 1) {
	fprintf(stderr, "Not enough interpret frequencies\n");
	return 1; }
      iprt->freq->entry[0][i] = Dbl2Edbl(d);

      ptr += len;
    }

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "Too many frequencies in reading interpret\n");
      return 1; }
  }

  PutConnect(s);

  return 0;
}

int readtrns(FILE *fp, Grammar *grammar)
{
  char *s;      /* Line with terminals */
  int ptr;      /* Pointer in line */
  Trns *trns;   /* For last read translation */
  char newsym;

  while ((s = GetConnect(fp)) != NULL && s[0] != '>') {
    ptr = 0;

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] == '\0') {
      fprintf(stderr, "No symbol in reading translate\n");
      return 1; }

    newsym = s[ptr++];

    for (;;) {
      while (isspace(s[ptr]))
	ptr++;

      if (s[ptr] == '\0')
	break;

      trns = (Trns *)malloc(sizeof(Trns));
      trns->newsym = newsym;
      trns->sym = s[ptr++];
      Enqueue(grammar->trns, trns);
    }
  }

  PutConnect(s);

  return 0;
}

int readsgrp(FILE *fp, Grammar *grammar)
{
  char *s;      /* Line with terminals */
  int ptr;      /* Pointer in line */
  Sgrp *sgrp;   /* For last read terminal */
  int numterm;  /* Number of terminals */
  int len;      /* Length of read numbers */
  int i, j;     /* For matrix entries */
  double d;     /* For reading doubles */

  numterm = grammar->term->count;

  while ((s = GetConnect(fp)) != NULL && s[0] != '>') {
    ptr = 0;

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] == '\0') {
      fprintf(stderr, "No symbol in reading single group\n");
      return 1; }

    sgrp = MakeSgrp(grammar);
    sgrp->sym = s[ptr++];
    Enqueue(grammar->sgrp, sgrp);

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "Too many symbols in reading single group\n");
      return 1; }

    /* Frequencies */

    if ((s = GetConnect(fp)) == NULL) {
      fprintf(stderr, "Error in reading single group frequencies\n");
      return 1; }

    ptr = 0;

    for (i = 0; i < numterm; i++) {
      if (sscanf(s+ptr, "%lf%n", &d, &len) != 1) {
	fprintf(stderr, "Not enough single group frequencies\n");
	return 1; }
      sgrp->freq->entry[0][i] = Dbl2Edbl(d);

      ptr += len;
    }

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "Too many frequencies in reading single group\n");
      return 1; }
    
    /* Eigenvector matrix */
    
    if ((s = GetConnect(fp)) == NULL) {
      fprintf(stderr, "Error in reading single group eigen matrix\n");
      return 1; }
  
    ptr = 0;

    for (j = 0; j < numterm; j++)
      for (i = 0; i < numterm; i++) {
	if (sscanf(s+ptr, "%lf%n", &d, &len) != 1) {
	  fprintf(stderr, "Not enough single group eigen matrix entries\n");
	  return 1; }
	sgrp->eigen->entry[j][i] = Dbl2Edbl(d);

	ptr += len;
      }
    
    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "Too many eigen matrix entries in reading single group\n");
      return 1; }

    /* Diagonal matrix */

    if ((s = GetConnect(fp)) == NULL) {
      fprintf(stderr, "Error in reading single group diag matrix\n");
      return 1; }

    ptr = 0;

    for (j = 0; j < numterm; j++)
      for (i = 0; i < numterm; i++) {
	if (sscanf(s+ptr, "%lf%n", &d, &len) != 1) {
	  fprintf(stderr, "Not enough single group diag matrix entries\n");
	  return 1; }
	sgrp->diag->entry[j][i] = Dbl2Edbl(d);

	ptr += len;
      }
    
    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "Too many diag matrix entries in reading single group\n");
      return 1; }

    /* Inverse eigenvector matrix */

    if ((s = GetConnect(fp)) == NULL) {
      fprintf(stderr, "Error in reading single group inveigen matrix\n");
      return 1; }

    ptr = 0;

    for (j = 0; j < numterm; j++)
      for (i = 0; i < numterm; i++) {
	if (sscanf(s+ptr, "%lf%n", &d, &len) != 1) {
	  fprintf(stderr, "Not enough single group inveigen matrix entries\n");
	  return 1; }
	sgrp->inveigen->entry[j][i] = Dbl2Edbl(d);

	ptr += len;
      }

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "Too many inveigen matrix entries in reading single group\n");
      return 1; }

  }

  PutConnect(s);
  /*
  PrintMatrix(stdout, sgrp->eigen);
  PrintMatrix(stdout, sgrp->diag);
  PrintMatrix(stdout, sgrp->inveigen);
  PrintMatrix(stdout, MulMatrix(sgrp->eigen, sgrp->inveigen));
  PrintMatrix(stdout, ExpMatrix(0.01, sgrp->eigen, sgrp->diag, sgrp->inveigen));
  a = MakeMatrix(4, 1);
  a->entry[0][0] = Dbl2Edbl(1);
  a->entry[1][0] = Dbl2Edbl(1);
  a->entry[2][0] = Dbl2Edbl(1);
  a->entry[3][0] = Dbl2Edbl(1);
  PrintMatrix(stdout, a);
  PrintMatrix(stdout, MulMatrix(ExpMatrix(0.01, sgrp->eigen, sgrp->diag, sgrp->inveigen), a));
  PrintMatrix(stdout, MulMatrix(ExpMatrix(0.1, sgrp->eigen, sgrp->diag, sgrp->inveigen), a));
  PrintMatrix(stdout, MulMatrix(ExpMatrix(1, sgrp->eigen, sgrp->diag, sgrp->inveigen), a));
  PrintMatrix(stdout, MulMatrix(ExpMatrix(10, sgrp->eigen, sgrp->diag, sgrp->inveigen), a));
  PrintMatrix(stdout, MulMatrix(ExpMatrix(100, sgrp->eigen, sgrp->diag, sgrp->inveigen), a));
  PrintMatrix(stdout, MulMatrix(ExpMatrix(1000, sgrp->eigen, sgrp->diag, sgrp->inveigen), a));
  */
  return 0;
}

int readdgrp(FILE *fp, Grammar *grammar)
{
  char *s;      /* Line with terminals */
  int ptr;      /* Pointer in line */
  Dgrp *dgrp;   /* For last read terminal */
  int numterm2;  /* Number of terminals squared */
  int len;      /* Length of read numbers */
  int i, j;     /* For matrix entries */
  double d;     /* For reading doubles */

  numterm2 = grammar->term->count*grammar->term->count;

  /* First line with symbol */

  while ((s = GetConnect(fp)) != NULL && s[0] != '>') {
    ptr = 0;

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] == '\0') {
      fprintf(stderr, "No symbol in reading double group\n");
      return 1; }

    dgrp = MakeDgrp(grammar);
    dgrp->sym = s[ptr++];
    Enqueue(grammar->dgrp, dgrp);

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "%s", s);
      fprintf(stderr, "Too many symbols in reading double group\n");
      return 1; }

    /* Frequencies */

    if ((s = GetConnect(fp)) == NULL) {
      fprintf(stderr, "Error in reading double group frequencies\n");
      return 1; }

    ptr = 0;

    for (i = 0; i < numterm2; i++) {
      if (sscanf(s+ptr, "%lf%n", &d, &len) != 1) {
	fprintf(stderr, "Not enough double group frequencies\n");
	return 1; }
      dgrp->freq->entry[0][i] = Dbl2Edbl(d);

      ptr += len;
    }

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "Too many frequencies in reading double group\n");
      return 1; }
    
    /* Eigenvector matrix */
    
    if ((s = GetConnect(fp)) == NULL) {
      fprintf(stderr, "Error in reading double group eigen matrix\n");
      return 1; }
  
    ptr = 0;

    for (j = 0; j < numterm2; j++)
      for (i = 0; i < numterm2; i++) {
	if (sscanf(s+ptr, "%lf%n", &d, &len) != 1) {
	  fprintf(stderr, "Not enough double group eigen matrix entries\n");
	  return 1; }
	dgrp->eigen->entry[j][i] = Dbl2Edbl(d);

	ptr += len;
      }
    
    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "Too many eigen matrix entries in reading double group\n");
      return 1; }

    /* Diagonal matrix */

    if ((s = GetConnect(fp)) == NULL) {
      fprintf(stderr, "Error in reading double group diag matrix\n");
      return 1; }

    ptr = 0;

    for (j = 0; j < numterm2; j++)
      for (i = 0; i < numterm2; i++) {
	if (sscanf(s+ptr, "%lf%n", &d, &len) != 1) {
	  fprintf(stderr, "Not enough double group diag matrix entries\n");
	  return 1; }
	dgrp->diag->entry[j][i] = Dbl2Edbl(d);

	ptr += len;
      }
    
    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "Too many diag matrix entries in reading double group\n");
      return 1; }

    /* Inverse eigenvector matrix */

    if ((s = GetConnect(fp)) == NULL) {
      fprintf(stderr, "Error in reading double group inveigen matrix\n");
      return 1; }

    ptr = 0;

    for (j = 0; j < numterm2; j++)
      for (i = 0; i < numterm2; i++) {
	if (sscanf(s+ptr, "%lf%n", &d, &len) != 1) {
	  fprintf(stderr, "Not enough double group inveigen matrix entries\n");
	  return 1; }
	dgrp->inveigen->entry[j][i] = Dbl2Edbl(d);

	ptr += len;
      }

    while (isspace(s[ptr]))
      ptr++;

    if (s[ptr] != '\0') {
      fprintf(stderr, "Too many inveigen matrix entries in reading double group\n");
      return 1; }

  }
  /*
  PrintMatrix(stdout, dgrp->eigen);
  PrintMatrix(stdout, dgrp->diag);
  PrintMatrix(stdout, dgrp->inveigen);
  PrintMatrix(stdout, MulMatrix(dgrp->eigen, dgrp->inveigen));
  PrintMatrix(stdout, ExpMatrix(0.01, dgrp->eigen, dgrp->diag, dgrp->inveigen));
  PrintMatrix(stdout, ExpMatrix(100, dgrp->eigen, dgrp->diag, dgrp->inveigen));
  a = MakeMatrix(16, 1);
  a->entry[0][0] = Dbl2Edbl(1);
  a->entry[1][0] = Dbl2Edbl(1);
  a->entry[2][0] = Dbl2Edbl(1);
  a->entry[3][0] = Dbl2Edbl(1);
  a->entry[4][0] = Dbl2Edbl(1);
  a->entry[5][0] = Dbl2Edbl(1);
  a->entry[6][0] = Dbl2Edbl(1);
  a->entry[7][0] = Dbl2Edbl(1);
  a->entry[8][0] = Dbl2Edbl(1);
  a->entry[9][0] = Dbl2Edbl(1);
  a->entry[10][0] = Dbl2Edbl(1);
  a->entry[11][0] = Dbl2Edbl(1);
  a->entry[12][0] = Dbl2Edbl(1);
  a->entry[13][0] = Dbl2Edbl(1);
  a->entry[14][0] = Dbl2Edbl(1);
  a->entry[15][0] = Dbl2Edbl(1);
  PrintMatrix(stdout, TransposeMatrix(a));
  PrintMatrix(stdout, TransposeMatrix(MulMatrix(ExpMatrix(0.01, dgrp->eigen, dgrp->diag, dgrp->inveigen), a)));
  PrintMatrix(stdout, TransposeMatrix(MulMatrix(ExpMatrix(0.1, dgrp->eigen, dgrp->diag, dgrp->inveigen), a)));
  PrintMatrix(stdout, TransposeMatrix(MulMatrix(ExpMatrix(1, dgrp->eigen, dgrp->diag, dgrp->inveigen), a)));
  PrintMatrix(stdout, TransposeMatrix(MulMatrix(ExpMatrix(10, dgrp->eigen, dgrp->diag, dgrp->inveigen), a)));
  PrintMatrix(stdout, TransposeMatrix(MulMatrix(ExpMatrix(100, dgrp->eigen, dgrp->diag, dgrp->inveigen), a)));
  PrintMatrix(stdout, TransposeMatrix(MulMatrix(ExpMatrix(1000, dgrp->eigen, dgrp->diag, dgrp->inveigen), a)));
*/
  PutConnect(s);

  return 0;
}

int readprod(FILE *fp, Grammar *grammar)
{
  char c[2];
  char s[4];  /* room for \0 */
  double d;
  Ruls *ruls;
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Ruldd *ruldd;
  Rulnn *rulnn;
  int l;
  int lnt, no1, no2;
  int len;
  char *t;
  int ptr;

  if ((t = GetConnect(fp)) == NULL) {
    fprintf(stderr, "Error in reading production rules\n");
    return 1; }

  ptr = 0;

  while (sscanf(t+ptr, "%1s -> %3s (%lf)%n", c, s, &d, &len) == 3) {
    ptr += len;

    if ((lnt = FindNont(*c, grammar->nont)) == -1) {
      printf("readprod: '%c' not nonterminal\n", *c);
      return 1;}

    l = strlen(s);
    if (l == 1)
      if ((no1 = FindSgrp(s[0], grammar->sgrp)) != -1) {
	/* Rule of type N -> s */
	ruls = (Ruls *)malloc(sizeof(Ruls));
	ruls->lnt = lnt;
	ruls->rsg = no1;
	ruls->prob = Dbl2Edbl(d);
	ruls->fprob = Dbl2Fdbl(d);
	Enqueue(grammar->ruls, ruls);
      }
      else if ((no1 = FindNont(s[0], grammar->nont)) != -1) {
	/* Rule of type N -> N */
	ruln = (Ruln *)malloc(sizeof(Ruln));
	ruln->lnt = lnt;
	ruln->rnt = no1;
	ruln->prob = Dbl2Edbl(d);
	ruln->fprob = Dbl2Fdbl(d);
	Enqueue(grammar->ruln, ruln);
      }
      else {
	printf("readprod: '%c' not terminal or nonterminal\n", s[0]);
	return 1;}
    else if (l == 2)
      if ((no1 = FindSgrp(s[0], grammar->sgrp)) != -1) {
	/* Rule of type N -> sN */
	if ((no2 = FindNont(s[1], grammar->nont)) == -1) {
	  printf("readprod: '%c' not nonterminal\n", s[1]);
	  return 1;}
	rulsn = (Rulsn *)malloc(sizeof(Rulsn));
	rulsn->lnt = lnt;
	rulsn->rsg = no1;
	rulsn->rnt = no2;
	rulsn->prob = Dbl2Edbl(d);
	rulsn->fprob = Dbl2Fdbl(d);
	Enqueue(grammar->rulsn, rulsn);
      }
      else if ((no1 = FindNont(s[0], grammar->nont)) != -1) {
	if ((no2 = FindNont(s[1], grammar->nont)) != -1) {
	  /* Rule of type N -> NN */
	  rulnn = (Rulnn *)malloc(sizeof(Rulnn));
	  rulnn->lnt = lnt;
	  rulnn->rnt1 = no1;
	  rulnn->rnt2 = no2;
	  rulnn->prob = Dbl2Edbl(d);
	  rulnn->fprob = Dbl2Fdbl(d);
	  Enqueue(grammar->rulnn, rulnn);
	}
	else if ((no2 = FindSgrp(s[1], grammar->sgrp)) != -1) {
	  /* Rule of type N -> Ns */
	  rulns = (Rulns *)malloc(sizeof(Rulns));
	  rulns->lnt = lnt;
	  rulns->rnt = no1;
	  rulns->rsg = no2;
	  rulns->prob = Dbl2Edbl(d);
	  rulns->fprob = Dbl2Fdbl(d);
	  Enqueue(grammar->rulns, rulns);
	}
	else {
	  printf("readprod: '%c' not terminal or nonterminal\n", s[1]);
	  return 1;}
      }
      else if ((no1 = FindDgrp(s[0], grammar->dgrp)) != -1) {
	/* Rule of type N -> dd */
	if (s[1] != s[0]) {
	  printf("readprod: double groups '%c' and '%c not identical\n",\
		 s[0], s[1]);
	  return 1;}
	ruldd = (Ruldd *)malloc(sizeof(Ruldd));
	ruldd->lnt = lnt;
	ruldd->rdb = no1;
	ruldd->prob = Dbl2Edbl(d);
	ruldd->fprob = Dbl2Fdbl(d);
	Enqueue(grammar->ruldd, ruldd);
      }
      else {
	printf("readprod: '%c' not terminal or nonterminal\n", s[0]);
	return 1;}
    else if (l == 3) {
      /* Rule of type N -> dNd */
      if ((no1 = FindDgrp(s[0], grammar->dgrp)) == -1) {
	printf("readprod: '%c' not terminal or nonterminal\n", s[0]);
	return 1;}
      if ((no2 = FindNont(s[1], grammar->nont)) == -1) {
	printf("readprod: '%c' not nonterminal\n", s[1]);
	return 1;}
      if (s[2] != s[0]) {
	printf("readprod: double groups '%c' and '%c not identical\n",\
	       s[0], s[2]);
	return 1;}
      ruldnd = (Ruldnd *)malloc(sizeof(Ruldnd));
      ruldnd->lnt = lnt;
      ruldnd->rdb = no1;
      ruldnd->rnt = no2;
      ruldnd->prob = Dbl2Edbl(d);
      ruldnd->fprob = Dbl2Fdbl(d);
      Enqueue(grammar->ruldnd, ruldnd);
    }
    else {
      printf("readprod: empty production\n");
      return 1;}
  }

  while (isspace(t[ptr]))
    ptr++;

  if (t[ptr] != '\0') {
    fprintf(stderr, "readprod: not finished '%s'\n", &t[ptr]);
    return 1;
  }

  return 0;
}

void PrintGrammar(FILE *fp, Grammar *grammar)
{
  LListCounter *nontcnt;
  LListCounter *termcnt;
  LListCounter *iprtcnt;
  LListCounter *trnscnt;
  LListCounter *sgrpcnt;
  LListCounter *dgrpcnt;
  LListCounter *lcount;
  Nont *nont;
  Term *term;
  Iprt *iprt;
  Trns *trns;
  Sgrp *sgrp;
  Dgrp *dgrp;
  Ruls *ruls;
  Ruln *ruln;
  Rulsn *rulsn;
  Rulns *rulns;
  Ruldnd *ruldnd;
  Ruldd *ruldd;
  Rulnn *rulnn;
  int nt = grammar->term->count;  /* number of terminals */
  int i, j;

  fprintf(fp, "> Nonterminal:\n");
  nontcnt = MakeCounter(grammar->nont, FIRST);
  while ((nont = Next(nontcnt)) != NULL)
    fprintf(fp, " %c", nont->sym);
  fprintf(fp, "\n\n");

  fprintf(fp, "> Terminal:\n");
  termcnt = MakeCounter(grammar->term, FIRST);
  while ((term = Next(termcnt)) != NULL)
    fprintf(fp, " %c", term->sym);
  fprintf(fp, "\n\n");
  
  fprintf(fp, "> Interpret:\n");
  iprtcnt = MakeCounter(grammar->iprt, FIRST);
  while ((iprt = Next(iprtcnt)) != NULL) {
    fprintf(fp, " %c", iprt->sym);
    for (i = 0; i < nt; i++)
      fprintf(fp, " %5.2f", Edbl2Dbl(iprt->freq->entry[0][i]));
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
  
  fprintf(fp, "> Translate:\n");
  trnscnt = MakeCounter(grammar->trns, FIRST);
  while ((trns = Next(trnscnt)) != NULL) {
    fprintf(fp, " %c", trns->newsym);
    fprintf(fp, " %c", trns->sym);
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
  
  fprintf(fp, "> Single group:\n");
  sgrpcnt = MakeCounter(grammar->sgrp, FIRST);
  while ((sgrp = Next(sgrpcnt)) != NULL) {
    fprintf(fp, " %c\n", sgrp->sym);
    for (i = 0; i < nt; i++)
      fprintf(fp, " %5.2f", Edbl2Dbl(sgrp->freq->entry[0][i]));
    fprintf(fp, "\n\n");
    for (j = 0; j < nt; j++) {
      for (i = 0; i < nt; i++)
	fprintf(fp, " %5.2f", Edbl2Dbl(sgrp->eigen->entry[j][i]));
      if (j < nt-1)
	fprintf(fp, " \\\n");
    }
    fprintf(fp, "\n\n");
    for (j = 0; j < nt; j++) {
      for (i = 0; i < nt; i++)
	fprintf(fp, " %5.2f", Edbl2Dbl(sgrp->diag->entry[j][i]));
      if (j < nt-1)
	fprintf(fp, " \\\n");
    }
    fprintf(fp, "\n\n");
    for (j = 0; j < nt; j++) {
      for (i = 0; i < nt; i++)
	fprintf(fp, " %5.2f", Edbl2Dbl(sgrp->inveigen->entry[j][i]));
      if (j < nt-1)
	fprintf(fp, " \\\n");
    }
    fprintf(fp, "\n\n");

    /*
    temp = MulMatrix(sgrp->eigen, sgrp->diag);
    temp2 = MulMatrix(temp, sgrp->inveigen);
    PrintMatrix(stdout, temp2);
    printf("\n");
    FreeMatrix(temp);
    FreeMatrix(temp2);

    temp = MulMatrix(sgrp->eigen, sgrp->inveigen);
    PrintMatrix(stdout, temp);
    printf("\n");
    FreeMatrix(temp);
    */
  }
  
  fprintf(fp, "> Double group:\n");
  dgrpcnt = MakeCounter(grammar->dgrp, FIRST);
  while ((dgrp = Next(dgrpcnt)) != NULL) {
    fprintf(fp, " %c\n", dgrp->sym);
    for (i = 0; i < nt*nt; i++)
      fprintf(fp, " %5.2f", Edbl2Dbl(dgrp->freq->entry[0][i]));
    fprintf(fp, "\n\n");
    for (j = 0; j < nt*nt; j++) {
      for (i = 0; i < nt*nt; i++)
	fprintf(fp, " %5.2f", Edbl2Dbl(dgrp->eigen->entry[j][i]));
      if (j < nt*nt-1)
	fprintf(fp, " \\\n");
    }
    fprintf(fp, "\n\n");
    for (j = 0; j < nt*nt; j++) {
      for (i = 0; i < nt*nt; i++)
	fprintf(fp, " %5.2f", Edbl2Dbl(dgrp->diag->entry[j][i]));
      if (j < nt*nt-1)
	fprintf(fp, " \\\n");
    }
    fprintf(fp, "\n\n");
    for (j = 0; j < nt*nt; j++) {
      for (i = 0; i < nt*nt; i++)
	fprintf(fp, " %5.2f", Edbl2Dbl(dgrp->inveigen->entry[j][i]));
      if (j < nt*nt-1)
	fprintf(fp, " \\\n");
    }
    fprintf(fp, "\n\n");
  }

  fprintf(fp, "> Production:\n");

  lcount = MakeCounter(grammar->ruls, FIRST);
  while ((ruls = Next(lcount)) != NULL) {
    fprintf(fp, " %c -> ", ((Nont *)SetCounter(nontcnt, ruls->lnt))->sym);
    fprintf(fp, "%c   ", ((Sgrp *)SetCounter(sgrpcnt, ruls->rsg))->sym);
    fprintf(fp, "(%8.6f)", Edbl2Dbl(ruls->prob));
  }

  InitCounter(lcount, grammar->ruln, FIRST);
  while ((ruln = Next(lcount)) != NULL) {
    fprintf(fp, " %c -> ", ((Nont *)SetCounter(nontcnt, ruln->lnt))->sym);
    fprintf(fp, "%c   ", ((Nont *)SetCounter(nontcnt, ruln->rnt))->sym);
    fprintf(fp, "(%8.6f)", Edbl2Dbl(ruln->prob));
  }

  InitCounter(lcount, grammar->rulsn, FIRST);
  while ((rulsn = Next(lcount)) != NULL) {
    fprintf(fp, " %c -> ", ((Nont *)SetCounter(nontcnt, rulsn->lnt))->sym);
    fprintf(fp, "%c", ((Sgrp *)SetCounter(sgrpcnt, rulsn->rsg))->sym);
    fprintf(fp, "%c  ", ((Nont *)SetCounter(nontcnt, rulsn->rnt))->sym);
    fprintf(fp, "(%8.6f)", Edbl2Dbl(rulsn->prob));
  }

  InitCounter(lcount, grammar->rulns, FIRST);
  while ((rulns = Next(lcount)) != NULL) {
    fprintf(fp, " %c -> ", ((Nont *)SetCounter(nontcnt, rulns->lnt))->sym);
    fprintf(fp, "%c", ((Nont *)SetCounter(nontcnt, rulns->rnt))->sym);
    fprintf(fp, "%c  ", ((Sgrp *)SetCounter(sgrpcnt, rulns->rsg))->sym);
    fprintf(fp, "(%8.6f)", Edbl2Dbl(rulns->prob));
  }

  InitCounter(lcount, grammar->ruldnd, FIRST);
  while ((ruldnd = Next(lcount)) != NULL) {
    fprintf(fp, " %c -> ", ((Nont *)SetCounter(nontcnt, ruldnd->lnt))->sym);
    fprintf(fp, "%c", ((Dgrp *)SetCounter(dgrpcnt, ruldnd->rdb))->sym);
    fprintf(fp, "%c", ((Nont *)SetCounter(nontcnt, ruldnd->rnt))->sym);
    fprintf(fp, "%c ", ((Dgrp *)SetCounter(dgrpcnt, ruldnd->rdb))->sym);
    fprintf(fp, "(%8.6f)", Edbl2Dbl(ruldnd->prob));
  }

  InitCounter(lcount, grammar->ruldd, FIRST);
  while ((ruldd = Next(lcount)) != NULL) {
    fprintf(fp, " %c -> ", ((Nont *)SetCounter(nontcnt, ruldd->lnt))->sym);
    fprintf(fp, "%c", ((Dgrp *)SetCounter(dgrpcnt, ruldd->rdb))->sym);
    fprintf(fp, "%c  ", ((Dgrp *)SetCounter(dgrpcnt, ruldd->rdb))->sym);
    fprintf(fp, "(%8.6f)", Edbl2Dbl(ruldd->prob));
  }

  InitCounter(lcount, grammar->rulnn, FIRST);
  while ((rulnn = Next(lcount)) != NULL) {
    fprintf(fp, " %c -> ", ((Nont *)SetCounter(nontcnt, rulnn->lnt))->sym);
    fprintf(fp, "%c", ((Nont *)SetCounter(nontcnt, rulnn->rnt1))->sym);
    fprintf(fp, "%c  ", ((Nont *)SetCounter(nontcnt, rulnn->rnt2))->sym);
    fprintf(fp, "(%8.6f)", Edbl2Dbl(rulnn->prob));
  }

  fprintf(fp, "\n");
}

void RobustGrammar(Grammar *grammar, double robust)
{
  int i, j;
  int numterm;
  
  numterm = grammar->term->count;

  for (i = 0; i < numterm+grammar->iprt->count; i++)
    for (j = 0; j < numterm; j++)
      if (Edbl2Dbl(grammar->quickdist[i][j]) < robust) 
	grammar->quickdist[i][j] = Dbl2Edbl(robust);
}

int FindNont(char sym, LList *nont)
{
  LListCounter *lcount = MakeCounter(nont, FIRST);
  Nont *elm;
  int num = -1;

  while ((elm = Next(lcount)) != NULL)
    if (elm->sym == sym)
      num = lcount->pos;

  free(lcount);

  return num;
}

/*------------------------------------------------------------------*/

int FindTerm(char sym, LList *term)
{
  LListCounter *lcount = MakeCounter(term, FIRST);
  Term *elm;
  int num = -1;

  while ((elm = Next(lcount)) != NULL)
    if (elm->sym == sym)
      num = lcount->pos;

  free(lcount);

  return num;
}

/*------------------------------------------------------------------*/

int FindIprt(char sym, LList *iprt)
{
  LListCounter *lcount = MakeCounter(iprt, FIRST);
  Iprt *elm;
  int num = -1;

  while ((elm = Next(lcount)) != NULL)
    if (elm->sym == sym)
      num = lcount->pos;

  free(lcount);

  return num;
}

/*------------------------------------------------------------------*/

int FindTrns(char sym, LList *trns)
{
  LListCounter *lcount = MakeCounter(trns, FIRST);
  Trns *elm;
  int num = -1;

  while ((elm = Next(lcount)) != NULL)
    if (elm->sym == sym)
      num = lcount->pos;

  free(lcount);

  return num;
}

/*------------------------------------------------------------------*/

int FindSgrp(char sym, LList *sgrp)
{
  LListCounter *lcount = MakeCounter(sgrp, FIRST);
  Sgrp *elm;
  int num = -1;

  while ((elm = Next(lcount)) != NULL)
    if (elm->sym == sym)
      num = lcount->pos;

  free(lcount);

  return num;
}

/*------------------------------------------------------------------*/

int FindDgrp(char sym, LList *dgrp)
{
  LListCounter *lcount = MakeCounter(dgrp, FIRST);
  Dgrp *elm;
  int num = -1;

  while ((elm = Next(lcount)) != NULL)
    if (elm->sym == sym)
      num = lcount->pos;

  free(lcount);

  return num;
}

/*------------------------------------------------------------------*/

int FindSym(char sym, Grammar *grammar)
{
  int num = -1;

  sym = Translate(sym, grammar->trns);

  if ((num = FindTerm(sym, grammar->term)) != -1)
    return num;

  if ((num = FindIprt(sym, grammar->iprt)) != -1)
    return grammar->term->count+num;

  /* Symbol not found, trying '-' */
  if (sym != '-') {
    /*    fprintf(stderr, "Symbol '%c' not found, using gap, '-'\n", sym);*/
    num = FindSym('-', grammar);
    return num;
  }

  fprintf(stderr, "Symbol '%c' not found\n", sym);
  exit(1);
}

/*------------------------------------------------------------------*/

char Translate(char sym, LList *trns)
{
  LListCounter *lcount = MakeCounter(trns, FIRST);
  Trns *elm;
  int num = -1;

  while ((elm = Next(lcount)) != NULL)
    if (elm->sym == sym)
      sym = elm->newsym;

  if (num != -1)
    sym = elm->newsym;

  return sym;
}
