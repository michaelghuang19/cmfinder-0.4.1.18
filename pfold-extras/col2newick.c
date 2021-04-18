#include "file.h"
#include "phyl.h"

void usage(void);

int main(int argc, char **argv)
{
  FILE *fp;
  Header *header;
  Entry *entry;
  int read_error;   /* For keeping track of errors in reading entries */
  Phyl *phyl;
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */

  cmdarg = InitArgument(argc, argv);

  while ((s = GetArgument(cmdarg)) != NULL) {
    usage();
    return 1;
  }

  if ((s = GetFilename(cmdarg)) == NULL)
    fp = stdin;
  else if ((fp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "col2newick: Error in opening file '%s'\n", s);
    return 1; }

  header = MakeHeader();
  entry = MakeEntry();

  if (ReadHeader(fp, header) != 0)
    return 1;

  while ((read_error = ReadEntry(fp, entry)) == 0) {
    if (!ReadType(entry, "TREE"))
      continue;

    phyl = ReadColPhyl(entry);

    PrintPhyl(stdout, phyl);
  }

  if (fp != stdin && fclose(fp) != 0) {
    fprintf(stderr, "stdpair: Error in closing file\n");
    return 1; }

  if (read_error == 1)
    return 1;

  return 0;
}

void usage(void)
{
  fprintf(stderr,
	  "usage: col2newick [<file>]\n");
}
