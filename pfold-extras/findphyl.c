#include "file.h"
#include "phyl.h"

void usage(void);

int main(int argc, char **argv)
{
  FILE *fp, *grammarfp;
  int i;
  Phyl *phyl, *phyl2;
  Grammar *grammar;
  Phyl **phyllist;
  Dist *dist;
  char **leaf;
  Align *align;
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */
  unsigned int len;
  char *tree;
  int option_file;
  int option_tree;
  int option_newick;
  int option_fast;
  double limit;
  Entry *entry;
  int read_error;   /* For keeping track of errors in reading entries */
  int size;
  Entry **entry_list;
  Header *header;

  cmdarg = InitArgument(argc, argv);

  option_tree = 0;
  option_newick = 0;
  option_fast = 0;
  option_file = 1;
  limit = -1;

  while ((s = GetArgument(cmdarg)) != NULL)
    if (strncmp(s, "-tree=", 6) == 0 && option_tree == 0) {
      tree = s+6;
      option_tree = 1;
      option_file = 0;
    }
    else if (strcmp(s, "-newick") == 0) {
      option_newick = 1;
    }
    else if (strcmp(s, "-fast") == 0) {
      option_fast = 1;
    }
    else if (strcmp(s, "-treefile") == 0 && option_tree == 0) {
      option_tree = 1;
    }
    else if (strncmp(s, "-limit=", 7) == 0) {
      if (sscanf(&s[7], "%lf%n",
		 &limit, &len) != 1 ||
	  len+7 != strlen(s)) {
	usage();
	return 1; }
    }
    else {
      usage();
      return 1;
    }

  if (option_tree == 0) {
    if ((s = GetFilename(cmdarg)) == NULL) {
      usage();
      return 1; }
    else if ((grammarfp = fopen(s, "r")) == NULL) {
      fprintf(stderr, "findphyl: Error in opening file '%s'\n", s);
      return 1; }
    else if ((s = GetFilename(cmdarg)) == NULL)
      fp = stdin;
    else if (GetFilename(cmdarg) != NULL) {
      usage();
      return 1; }
    else if ((fp = fopen(s, "r")) == NULL) {
      fprintf(stderr, "findphyl: Error in opening file '%s'\n", s);
      return 1; }
  }
  else if (option_file == 0) {
    if (GetFilename(cmdarg) != NULL) {
      usage();
      return 1; }
  }
  else {
    if ((s = GetFilename(cmdarg)) == NULL)
      fp = stdin;
    else if (GetFilename(cmdarg) != NULL) {
      usage();
      return 1; }
    else if ((fp = fopen(s, "r")) == NULL) {
      fprintf(stderr, "findphyl: Error in opening file '%s'\n", s);
      return 1; }
    tree = GetLine(fp);
  }

  if (option_tree == 0) {
    grammar = ReadGrammar(grammarfp);

    /* Read the sequences */
    size = 1;

    entry_list = (Entry **)malloc(size*sizeof(Entry *));

    header = MakeHeader();

    if (ReadHeader(fp, header) != 0)
      return 1;

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

    align = Col2Align(entry_list, (int (*)(char, void *))FindSym, (void *)grammar);
    if (option_fast == 0)
      dist = AlignDist(align, grammar);
    else
      dist = FastAlignDist(align, grammar);
    phyl = Neighbour(dist);
  }
  else {
    phyl = ReadPhyl(tree, 1.);
    dist = PhylDist(phyl);
  }
  /*
  if (CountLeaf(phyl->root) < 5)
    PrintDist(stderr, dist);
  */
  if (limit == -1) {
    FixPhyl(phyl);
    if (option_newick == 1)
      PrintPhyl(stdout, phyl);
    else {
      entry = PhylEntry(phyl, "tree");
      PrintHeader(stdout, header);
      PrintEntry(stdout, entry);
      for (i = 0; entry_list[i] != NULL; i++)
	PrintEntry(stdout, entry_list[i]);
    }
  }
  else {
    phyllist = UPGMAlimit(dist, limit);
    for (i = 0; phyllist[i] != NULL; i++) {
      leaf = PhylLeaves(phyllist[i]);
      phyl2 = SubPhyl(phyl, leaf);
      FixPhyl(phyl2);
      if (option_newick == 1)
	PrintPhyl(stdout, phyl2);
      else {
	entry = PhylEntry(phyl2, "tree");
	PrintHeader(stdout, header);
	PrintEntry(stdout, entry);
	for (i = 0; entry_list[i] != NULL; i++)
	  PrintEntry(stdout, entry_list[i]);
      }
    }
  }
  
  return 0;
}

void usage(void)
{
  fprintf(stderr, "usage: findphyl [--limit=LIMIT] RATEFILE [COLFILE]\n");
  fprintf(stderr, "       findphyl --tree=<tree> [--limit=LIMIT]\n");
  fprintf(stderr, "       findphyl --treefile [--limit=LIMIT] [TREEFILE]\n");
}
