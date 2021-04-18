#include "grammar.h"
#include "newcolprob.h"
#include "inout.h"
#include "file.h"
#include "optimize.h"

void usage(void);

int main(int argc, char **argv)
{
  FILE *gramfp, *alignfp, *treefp;         /* For input */
  FILE *pairprobfp;
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* Strng for arguments */
  Alignstr *alignstr;
  Phyl *phyl;
  Grammar *grammar;
  char *tree;
  char *ppfile, *treefile;
  int option_tree;
  int option_treeinfile;
  int option_opt;
  int option_fast;
  int option_wlen;
  int wlen;
  Colprob *cp;
  Inside *e;
  Outside *f;
  Pairprob *pp;
  Pairprob *pp2;
  PpCyk *ppcyk;
  Cyk *cyk;
  FColprob *fcp;
  FCyk *fcyk;
  double branchmul;
  int len;
  double gaplimit;
  double logp_limit;
  Header *header;
  Entry *entry;
  Entry *str_entry;
  Entry **entry_list;
  int size;
  int i, j;
  int align_bp_col, alignpos_col, seq_bp_col, seqpos_col, res_col;
  int read_error;   /* For keeping track of errors in reading entries */
  char field[MAXCOLW];
  char name[MAXCOLW];
  double a, b, x, f_val;  /* For minimizing */
  int error;
  double robust;
  int *cykstr;
  int p;
  char *structure;
  int par;
  int entry_offset;

  cmdarg = InitArgument(argc, argv);

  entry_offset = 0;
  str_entry = NULL;
  tree = NULL;
  branchmul= 1;
  gaplimit = 0.75;
  logp_limit = -3;
  option_tree = 0;
  option_treeinfile = 0;
  treefile = NULL;
  ppfile = NULL;
  option_opt = 0;
  robust = 0.01;
  option_fast = 0;
  option_wlen = 0;

  while ((s = GetArgument(cmdarg)) != NULL)
    if (strncmp(s, "-tree=", 6) == 0) {
      tree = s+6;
      option_tree = 1;
    }
    else if (strcmp(s, "-optimize") == 0) {
      option_opt = 1;
    }
    else if (strcmp(s, "-treefile") == 0) {
      if ((s = GetFilename(cmdarg)) == NULL) {
	usage();
	return 1; }
      treefile = s;
      option_tree = 1;
    }
    else if (strcmp(s, "-treeinfile") == 0) {
      option_treeinfile = 1;
      option_tree = 1;
    }
    else if (strncmp(s, "-branchmul=", 11) == 0) {
      if (sscanf(&s[11], "%lf%n",
		 &branchmul, &len) != 1 ||
	  len+11 != strlen(s)) {
	usage();
	return 1; }
    }
    else if (strncmp(s, "-robust=", 8) == 0) {
      if (sscanf(&s[8], "%lf%n",
		 &robust, &len) != 1 ||
	  len+8 != strlen(s)) {
	usage();
	return 1; }
    }
    else if (strncmp(s, "-gaplimit=", 10) == 0) {
      if (sscanf(&s[10], "%lf%n",
		 &gaplimit, &len) != 1 ||
	  len+10 != strlen(s)) {
	usage();
	return 1; }
    }
    else if (strncmp(s, "-logplimit=", 11) == 0) {
      if (sscanf(&s[11], "%lf%n",
		 &logp_limit, &len) != 1 ||
	  len+11 != strlen(s)) {
	usage();
	return 1; }
    }
    else if (strncmp(s, "-wlen=", 6) == 0) {
      if (sscanf(&s[6], "%d%n",
		 &wlen, &len) != 1 ||
	  len+6 != strlen(s)) {
	usage();
	return 1; }
      option_wlen = 1;
    }
    else if (strcmp(s, "-ppfile") == 0) {
      if ((s = GetFilename(cmdarg)) == NULL) {
	usage();
	return 1; }
      ppfile = s;
    }
    else if (strcmp(s, "-fast") == 0) {
      option_fast = 1;
    }
    else {
      usage();
      return 1; }

  if ((s = GetFilename(cmdarg)) == NULL) {
    usage();
    return 1; }
  else if ((gramfp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "scfg: Error in opening file '%s'\n", s);
    return 1; }

  if ((s = GetFilename(cmdarg)) == NULL)
    alignfp = stdin;
  else if (GetFilename(cmdarg) != NULL) {
    usage();
    return 1; }
  else if ((alignfp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "scfg: Error in opening file '%s'\n", s);
    return 1; }

  if (treefile != NULL) {
    if ((treefp = fopen(treefile, "r")) == NULL) {
      fprintf(stderr, "scfg: Error in opening file '%s'\n", s);
      return 1; }

    header = MakeHeader();
    entry = MakeEntry();

    if (ReadHeader(treefp, header) != 0)
      return 1;

    phyl = NULL;

    while ((read_error = ReadEntry(treefp, entry)) == 0) {
      if (!ReadType(entry, "TREE"))
	continue;

      phyl = ReadColPhyl(entry);

      break;
    }
    
    if (phyl == NULL) {
      fprintf(stderr, "scfg:: No tree read\n");
      return 1;
    }

    if (fclose(treefp) != 0) {
      fprintf(stderr, "scfg: Error in closing tree file\n");
      return 1; }

    if (branchmul != 1)
      fprintf(stderr, "Warning: branchmul ignored, since treefile is used\n");
  }

  if ((option_fast == 1 || option_wlen == 1) &&
      (ppfile != NULL || option_opt == 1)) {
    usage();
    return 1;
  }

  header = MakeHeader();

  if (ReadHeader(alignfp, header) != 0) {
    fprintf(stderr, "scfg: Could not read header\n");
    exit(1); }

  size = 1;

  entry_list = (Entry **)malloc(size*sizeof(Entry *));

  for (i = 0; (entry = MakeEntry()) != NULL &&
	      (read_error = ReadEntry(alignfp, entry)) == 0; i++) {
    ReadText(entry, "ENTRY", MAXCOLW, name);
    if (StrCmp(name, "structure") == 0) {
      str_entry = entry;
      i--;
      continue;
    }
    if (i >= size) {
      size *= 2;
      entry_list = (Entry **)realloc(entry_list, size*sizeof(Entry *));
    }
    entry_list[i] = entry;
  }

  if (alignfp != stdin && fclose(alignfp) != 0) {
    fprintf(stderr, "scfg: Error in closing col file\n");
    return 1; }

  if (i+1 != size)
    entry_list = (Entry **)realloc(entry_list, (i+1)*sizeof(Entry *));

  entry_list[i] = NULL;

  if (read_error == 1)
    exit(1);

  if (option_treeinfile == 1) {
    for (i = 0; entry_list[i] != NULL; i++) {
      if (!ReadType(entry_list[i], "TREE"))
	continue;

      phyl = ReadColPhyl(entry_list[i]);
      
      break;
    }

    if (phyl == NULL) {
      fprintf(stderr, "scfg: No tree read\n");
      return 1;
    }

    if (branchmul != 1)
      fprintf(stderr, "Warning: branchmul ignored, since treefile is used\n");
  }

  for (i = 0; entry_list[i] != NULL; i++) {
    if (ReadColno(entry_list[i], "residue") == 0)
      continue;
    break;
  }

  entry_offset = i;

  if (option_tree == 0) {
    if(i > 1) {
      fprintf(stderr, "scfg: a tree is need for more than one sequence\n");
      exit(1);
    }
    else {
      ReadText(entry_list[entry_offset+0], "ENTRY", MAXCOLW, name);
      tree = (char *)malloc((strlen(name)+4) * sizeof(char));
      sprintf(tree, "(%s);", name);
    }
  }

  grammar = ReadGrammar(gramfp);

  grammar->mindist = 4;

  if (fclose(gramfp) != 0) {
    fprintf(stderr, "scfg: Error in closing grammar file\n");
    return 1; }

  if (robust != 0)
    RobustGrammar(grammar, robust);

  if (treefile == NULL && option_treeinfile == 0)
    phyl = ReadPhyl(tree, branchmul);
  FixPhyl(phyl);

  alignstr = (Alignstr *)malloc(sizeof(Alignstr));

  alignstr->align = Col2Align(entry_list, (int (*)(char, void *))FindSym, (void *)grammar);

  len = EntryLength(entry_list[entry_offset+0]);
  alignstr->str = (int *)malloc(len * sizeof(int));
    
  if (str_entry == NULL) {
    align_bp_col = ReadColno(entry_list[entry_offset+0], "align_bp");
    alignpos_col = ReadColno(entry_list[entry_offset+0], "alignpos");
    seq_bp_col = ReadColno(entry_list[entry_offset+0], "seq_bp");
    seqpos_col = ReadColno(entry_list[entry_offset+0], "seqpos");

    if (align_bp_col == 0 && seq_bp_col == 0)
      for (i = 0; i < len; i++)
	alignstr->str[i] = -1;
    else
      for (i = 0; i < len; i++) {
	if (align_bp_col != 0)
	  GetField(field, entry_list[entry_offset+0], i+1, align_bp_col);
	else
	  GetField(field, entry_list[entry_offset+0], i+1, seq_bp_col);
	
	if (StrCmp(field, "s") == 0)
	  alignstr->str[i] = SINGLE;
	else if (StrCmp(field, "d") == 0)
	  alignstr->str[i] = DOUBLE;
	else {
	  alignstr->str[i] = FindPair(entry_list[entry_offset+0], i+1,
				      align_bp_col, alignpos_col,
				      seq_bp_col, seqpos_col)-1;
	}
      }
  }
  else {
    res_col = ReadColno(entry_list[entry_offset+0], "residue");
    if (res_col == 0)
      for (i = 0; i < len; i++)
	alignstr->str[i] = -1;
    else {
      structure = (char *)malloc(len * sizeof(char));
      for (i = 0; i < len; i++) {
	GetField(field, str_entry, i+1, res_col);
	structure[i] = field[0];
      }
      for (i = 0; i < len; i++) {
	alignstr->str[i] = -1;
	if (structure[i] == '(') {
	  par = 1;
	  for (j = i+1; j < len; j++) {
	    if (structure[j] == '(') par++;
	    if (structure[j] == ')') par--;
	    if (par == 0) break;
	  }
	  if (par == 0 && j >= i+grammar->mindist)
	    alignstr->str[i] = j;
	}
	else if (structure[i] == ')') {
	  par = 1;
	  for (j = i-1; j >= 0; j--) {
	    if (structure[j] == '(') par--;
	    if (structure[j] == ')') par++;
	    if (par == 0) break;
	  }
	  if (par == 0 && j+grammar->mindist <= i)
	    alignstr->str[i] = j;
	}
	else if (tolower(structure[i]) == 's')
	  alignstr->str[i] = SINGLE;
	else if (tolower(structure[i]) == 'd')
	  alignstr->str[i] = DOUBLE;
      }
    }
  }

  alignstr->map = (int *)malloc((len+1) * sizeof(int));

  for (i = 0; i < len; i++)
    alignstr->map[i] = i+1;
  alignstr->map[i] = len;

  alignstr = Removegap(grammar, alignstr, gaplimit);

  if (option_wlen == 0)
    wlen = alignstr->align->len;

  if (option_fast == 1) {
    fcp = MakeFColprob(grammar, wlen);
    fcyk = MakeFCyk(grammar, wlen);
  }
  else if (option_wlen == 1) {
    cp = MakeColprob(grammar, wlen);
    cyk = MakeCyk(grammar, wlen);
  }
  else {
    e = MakeInside(grammar, wlen);
    f = MakeOutside(grammar, wlen);
    cp = MakeColprob(grammar, wlen);
    pp = MakePairprob(wlen);
    //    pp2 = MakePairprob(wlen);
    ppcyk = MakePpCyk(wlen);
    //    cyk = MakeCyk(grammar, wlen);
  }

  if (option_opt == 1) {
    a = 1*log(branchmul);
    b = log(branchmul)+.1;
      
    InitMinimize(a, b, 0.001);
    while ((error = Minimize(&x, &f_val)) == 0) {
      branchmul = exp(x);
      phyl = ReadPhyl(tree, branchmul);
      FixPhyl(phyl);
      InitCol(phyl, grammar, wlen, alignstr, logp_limit);
      InitColprob(cp, alignstr);
      InitInside(e, cp);
      FinishCol(phyl, grammar);
      f_val = -Edbl2Dbl(LogEdouble(
        e->prob[FindNont('S', grammar->nont)][0][wlen-1]));
      fprintf(stderr, "eval %24.20f %24.20f\n", x, f_val);
    }
  
    fprintf(stderr, "evalf %24.20f %24.20f\n", branchmul, f_val);
      
    branchmul = exp(x);
    phyl = ReadPhyl(tree, branchmul);
      
    if (error == 1) {
      fprintf(stderr, "maximization problems\n");
      return 1;
    }
  }

  if (option_wlen == 1 && option_fast == 0) {
    InitCol(phyl, grammar, wlen, alignstr, logp_limit);
    cp = MakeColprob(grammar, wlen);
    InitColprob(cp, alignstr);
    cyk = MakeCyk(grammar, wlen);
    InitCyk(cyk, cp);
    cykstr = CykStr(cyk);
    for (i = 0; i < wlen; i++)
      alignstr->str[i] = cykstr[i];

    for (p = wlen; p < alignstr->align->len; p++) {
      if (p % 100 == 0)
	fprintf(stderr, ".");
      MoveColprob(cp, alignstr, p);
      MoveCyk(cyk, cp, p);
      if (cyk->numpair >= 1) {
	fprintf(stderr, " %d", p);
	cykstr = CykStr(cyk);
	for (i = 0; i < wlen; i++) {
	  if (cykstr[i] != -1) {
	    if (alignstr->str[i+p-wlen+1] == -1)
	      alignstr->str[i+p-wlen+1] = cykstr[i]+p-wlen+1;
	    else if (alignstr->str[i+p-wlen+1] != cykstr[i]+p-wlen+1)
	      alignstr->str[i+p-wlen+1] = DOUBLE;
	  }
	}
      }
    }
    fprintf(stderr, "\n");

    len = alignstr->align->len;

    for (i = 0; i < len; i++)
      if (alignstr->str[i] == SINGLE)
	alignstr->str[i] = -1;

    for (i = 0; i < len; i++) {
      if (alignstr->str[i] >= 0 && alignstr->str[alignstr->str[i]] != i) {
	alignstr->str[i] = -1;
	alignstr->str[alignstr->str[i]] = -1;
      }
    }
    
    for (i = 0; entry_list[entry_offset+i] != NULL; i++)
      AddStr(entry_list[entry_offset+i], alignstr, NULL);
  }
  else if (option_wlen == 1 && option_fast == 1) {
    InitCol(phyl, grammar, wlen, alignstr, logp_limit);
    fcp = MakeFColprob(grammar, wlen);
    InitFColprob(fcp, alignstr);
    fcyk = MakeFCyk(grammar, wlen);
    InitFCyk(fcyk, fcp);
    cykstr = FCykStr(fcyk);
    for (i = 0; i < wlen; i++)
      alignstr->str[i] = cykstr[i];

    for (p = wlen; p < alignstr->align->len; p++) {
      if (p % 100 == 0)
	fprintf(stderr, ".");
      MoveFColprob(fcp, alignstr, p);
      MoveFCyk(fcyk, fcp, p);
      if (fcyk->numpair >= 20) {
	cykstr = FCykStr(fcyk);
	for (i = 0; i < wlen; i++)
	  if (cykstr[i] != -1) {
	    if (alignstr->str[i+p-wlen+1] == -1)
	      alignstr->str[i+p-wlen+1] = cykstr[i]+p-wlen+1;
	    else if (alignstr->str[i+p-wlen+1] != cykstr[i]+p-wlen+1)
	      alignstr->str[i+p-wlen+1] = SINGLE;
	  }
      }
    }
    fprintf(stderr, "\n");

    for (i = 0; i < len; i++)
      if (alignstr->str[i] == SINGLE)
	alignstr->str[i] = -1;
    
    for (i = 0; i < len; i++)
      if (alignstr->str[i] != -1 && alignstr->str[alignstr->str[i]] != i) {
	alignstr->str[i] = -1;
	alignstr->str[alignstr->str[i]] = -1;
      }
    
    for (i = 0; entry_list[entry_offset+i] != NULL; i++)
      AddStr(entry_list[entry_offset+i], alignstr, NULL);
  }
  else if (option_fast == 1) {
    InitCol(phyl, grammar, wlen, alignstr, logp_limit);
    fprintf(stderr, "colprob\n");
    InitFColprob(fcp, alignstr);
    fprintf(stderr, "fast cyk\n");
    InitFCyk(fcyk, fcp);

    alignstr->str = FCykStr(fcyk);

    for (i = 0; entry_list[entry_offset+i] != NULL; i++)
      AddStr(entry_list[entry_offset+i], alignstr, NULL);
  }
  else {
    if (wlen > 0) {
      InitCol(phyl, grammar, wlen, alignstr, logp_limit);
      InitColprob(cp, alignstr);
      InitInside(e, cp);
      InitOutside(e, f, cp);
      InitPairprob(pp, e, f, cp);
      /*    InitPairaff(pp2, e, f, cp);*/
      InitPpCyk(ppcyk, pp);
      alignstr->str = PpCykStr(ppcyk);
    }
    for (i = 0; entry_list[entry_offset+i] != NULL; i++)
      AddStr(entry_list[entry_offset+i], alignstr, pp);

    //    printf("%.4f\n", -Edbl2Dbl(LogEdouble(
    //        e->prob[FindNont('S', grammar->nont)][0][wlen-1])));
  }

  PrintHeader(stdout, header);
  for (i = 0; entry_list[i] != NULL; i++)
    PrintEntry(stdout, entry_list[i]);

  if (ppfile != NULL) {
    if ((pairprobfp = fopen(ppfile, "w")) == NULL) {
      fprintf(stderr, "scfg: Error in opening file '%s'\n", s);
      return 1; }
    PrintPairprob(pairprobfp, pp, pp, alignstr);
    if (fclose(pairprobfp) != 0) {
      fprintf(stderr, "scfg: Error in closing pairprob file\n");
      return 1; }
  }

  return 0;
}

void usage(void)
{
  fprintf(stderr,
	  "usage: scfg [--treeinfile | --treefile <treefile> | --tree=<tree>] <grammarfile> [<colfile>]\n");
}
