/*********************************************************************

  col2fasta.c

  usage: col2fasta -m [FILE]

  This program makes a faste file from a col file. See also man page.

  000919 Bjarne Knudsen (bk@daimi.au.dk)

  Copyright (C) 2000 Bjarne Knudsen

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.

*********************************************************************/

#include "col.h"
#include "file.h"

#define MAXNAME 1000

void usage(void);

int main(int argc, char **argv)
{
  FILE *fp;         /* For input */
  Header *header;
  Entry *entry;
  int len, seqlen;  /* Sequence length */
  char field[MAXCOLW];
  int i;
  char name[MAXNAME];  /* Sequence name */
  char version[MAXNAME];  /* Sequence version */
  char version2[MAXNAME];  /* Sequence version */
  int nuc_col;      /* Column numbers */
  int read_error;   /* For keeping track of errors in reading entries */
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */
  int option_m;

  cmdarg = InitArgument(argc, argv);

  option_m = 0;

  while ((s = GetArgument(cmdarg)) != NULL)
    if (strcmp(s, "m") == 0)
      option_m = 1;
    else if (strcmp(s, "-minus") == 0)
      option_m = 1;
    else {
      usage();
      return 1; }

  if ((s = GetFilename(cmdarg)) == NULL)
    fp = stdin;
  else if (GetFilename(cmdarg) != NULL) {
    usage();
    return 1; }
  else if ((fp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "col2fasta: Error in opening file '%s'\n", s);
    return 1; }

  header = MakeHeader();
  entry = MakeEntry();

  if (ReadHeader(fp, header) != 0)
    return 1;

  while ((read_error = ReadEntry(fp, entry)) == 0) {
    if ((nuc_col = ReadColno(entry, "residue")) == 0 &&
	(nuc_col = ReadColno(entry, "nucleotide")) == 0 &&
	(nuc_col = ReadColno(entry, "aminoacid")) == 0 &&
	(nuc_col = ReadColno(entry, "code")) == 0) {
      fprintf(stderr,
	      "col2fasta: Warning: Ignoring sequence, no column info\n");
      continue; }

    if (ReadText(entry, "ENTRY", MAXNAME, name) != 0 &&
	ReadText(entry, "SEQUENCE", MAXNAME, name) != 0)
      name[0] = '\0';
    if (ReadText(entry, "VERSION", MAXNAME, version) != 0)
      version[0] = '\0';
    if (ReadText(entry, "VERSION2", MAXNAME, version2) != 0)
      version2[0] = '\0';

    for (i = 0; name[i] != '\0'; i++)
      if (name[i] == ' ')
	name[i] = '_';

    printf(">%s", name);
    if (version[0] != '\0')
      printf("_%s", version);
    if (version2[0] != '\0')
      printf("_%s", version2);
    printf("\n");

    len = EntryLength(entry);
    seqlen = 0;
    if (option_m == 1)
      for (i = len; i >= 1; i--) {
	GetField(field, entry, i, nuc_col);
	seqlen++;
	if (seqlen%60 == 1 && seqlen != 1)
	  printf("\n");
	switch (toupper(field[0])) {
	case 'A':
	  printf("T");
	  break;
	case 'C':
	  printf("G");
	  break;
	case 'G':
	  printf("C");
	  break;
	case 'T':
	case 'U':
	  printf("A");
	  break;
	case '-':
	  printf("-");
	  break;
        default:
	  printf("X");
	}
      }
    else
      for (i = 1; i <= len; i++) {
	GetField(field, entry, i, nuc_col);
	seqlen++;
	if (seqlen%60 == 1 && seqlen != 1)
	  printf("\n");
	printf("%c", toupper(field[0]));
      }
    printf("\n");
  }

  if (fp != stdin && fclose(fp) != 0) {
    fprintf(stderr, "col2fasta: Error in closing file\n");
    return 1; }

  if (read_error == 1)
    return 1;

  return 0;
}

void usage(void)
{
  fprintf(stderr,
	  "usage: col2fasta [-m | --minus] [<file]\n");
}
