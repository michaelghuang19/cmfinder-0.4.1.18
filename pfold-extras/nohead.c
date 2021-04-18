/*********************************************************************

  nohead.c

  usage: nohead [<file>]

  See man page for more info

  000922 Bjarne Knudsen (bk@daimi.au.dk)

*********************************************************************/

#include "col.h"
#include "file.h"

void usage(void);

int main(int argc, char **argv)
{
  FILE *fp;         /* For input */
  Header *header;
  Entry *entry;
  int read_error;   /* For keeping track of errors in reading entries */
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */

  cmdarg = InitArgument(argc, argv);

  if ((s = GetArgument(cmdarg)) != NULL) {
    usage();
    return 1; }

  if ((s = GetFilename(cmdarg)) == NULL)
    fp = stdin;
  else if (GetFilename(cmdarg) != NULL) {
    usage();
    return 1; }
  else if ((fp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "nohead: Error in opening file '%s'\n", s);
    return 1; }

  header = MakeHeader();
  entry = MakeEntry();

  if (ReadHeader(fp, header) != 0)
    return 1;

  while ((read_error = ReadEntry(fp, entry)) == 0) {
    PrintEntry(stdout, entry);
  }

  if (fp != stdin && fclose(fp) != 0) {
    fprintf(stderr, "nohead: Error in closing file\n");
    return 1; }

  if (read_error == 1)
    return 1;
  
  return 0;
}

void usage(void)
{
  fprintf(stderr, "usage: nohead [<file>]\n");
}
