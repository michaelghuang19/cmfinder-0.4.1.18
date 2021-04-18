/*********************************************************************

  col2txtalign.c

  usage: col2txtalign
              [-r<range> | -range=<range>] [--namewidth=<namewidth>]
	      [--textwidth=<textwidth>] [--space] [<file>]

  See man page for more info

  000928 Bjarne Knudsen (bk@daimi.au.dk)

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
#include <ctype.h>

#define MAXNAME 1000

void usage(void);
int inrange(int num, char *range);

int main(int argc, char **argv)
{
  Header *header;
  Entry *entry;
  FILE *fp;
  int codecol, symcol;
  int read_error;   /* for keeping track of errors in reading entries */
  int i, j;
  int width;
  char *range;
  int textlen;
  char **text;
  char **name;
  int maxlen;      /* Dynamic line allocation maximum */
  int numseq;
  int ptr;         /* Pointer in text string */
  int seqlen;
  char *writecode;
  char *symstr;
  int last_inrange;
  int curr_width, lastspace;
  int startpos;
  int startptr;
  char **sym;
  char field[MAXCOLW];
  int newblock;
  int codelen;
  int blockpage, block;
  int page;
  int nametxtwidth;
  unsigned int len; /* Used with options */
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */
  int option_space;
  int lineno;
  int colno, maxcolno;
  int txtwidth;
  int spaceptr1, spaceptr2, diff;
  char format[80];

  /* Postscript variables */
  option_space = 0;
  nametxtwidth = 10;
  txtwidth = 60;

  range = NULL;

  cmdarg = InitArgument(argc, argv);

  while ((s = GetArgument(cmdarg)) != NULL)
    if (strncmp(s, "r", 1) == 0)
      range = &s[1];
    else if (strncmp(s, "-range=", 7) == 0)
      range = &s[7];
    else if (strcmp(s, "-space") == 0)
      option_space = 1;
    else if (strncmp(s, "-namewidth=", 11) == 0) {
      if (sscanf(s, "-namewidth=%d%n", &nametxtwidth, &len) != 1 ||
	  len != strlen(s)) {
	usage();
	return 1; }
    }
    else if (strncmp(s, "-textwidth=", 11) == 0) {
      if (sscanf(s, "-textwidth=%d%n", &txtwidth, &len) != 1 ||
	  len != strlen(s)) {
	usage();
	return 1; }
    }
    else {
      usage();
      return 1; }
  
  if ((s = GetFilename(cmdarg)) == NULL)
    fp = stdin;
  else if (GetFilename(cmdarg) != NULL) {
    usage();
    return 1; }
  else if ((fp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "col2txtalign: Error in opening file '%s'\n", s);
    return 1; }

  header = MakeHeader();
  entry = MakeEntry();

  if (ReadHeader(fp, header) != 0 ||
      (read_error = ReadEntry(fp, entry)) == 1)
    return 1;
  else if (read_error == -1) {
    fprintf(stderr, "col2txtalign: Warning, no entries\n");
    return 0; }
    
  seqlen = EntryLength(entry);
  writecode = (char *)malloc((6*(seqlen-1)/2+2) * sizeof(char));

  last_inrange = 1;
  textlen = 0;

  for (ptr = i = 0; i < seqlen; i++) {
    if (inrange(i+1, range)) {
      if ((option_space == 1 && i%10 == 0 && last_inrange && i != 0) ||
	  (option_space == 0 && txtwidth != 0 && 
	   ptr%(txtwidth+1) == txtwidth))
	writecode[ptr++] = ' ';
      else if (!last_inrange) {
	strcpy(&writecode[ptr], "... ");
	ptr += 4;
      }
      writecode[ptr++] = 'c';
      textlen++;
    }
    else if (i != 0)
      if (last_inrange)
	writecode[ptr++] = ' ';
    last_inrange = inrange(i+1, range);
  }

  if (!last_inrange) {
    strcpy(&writecode[ptr], "...");
    ptr += 3;
  }

  writecode[ptr++] = '\n';
  writecode[ptr++] = '\0';

  codelen = ptr;
  writecode = (char *)realloc(writecode, ptr * sizeof(char));
  symstr = (char *)malloc(ptr * sizeof(char));

  if (txtwidth != 0)
    width = txtwidth;

  lastspace = -1;
  curr_width = 0;

  for (i = 0; writecode[i] != '\n'; i++) { 
    if (writecode[i] == ' ') {
      if (curr_width > width) {
	if (lastspace != -1) {
	  writecode[lastspace] = '\n';
	  curr_width = i - lastspace;
	  lastspace = i;
	  continue; }
	else {
	  writecode[i] = '\n';
	  curr_width = 0;
	  lastspace = -1;
	  continue; }
      }
      lastspace = i;
    }
    curr_width++;
  }  

  if (curr_width > width)
    if (lastspace != -1)
      writecode[lastspace] = '\n';

  sym = (char **)malloc(textlen * sizeof(char *));
  for (i = 0; i < textlen; i++)
    sym[i] = (char *)malloc(MAXCOLW * sizeof(char));

  if ((symcol = ReadColno(entry, "alignpos")) == 0 &&
      (symcol = ReadColno(entry, "seqpos")) == 0) {
    fprintf(stderr,
	    "Warning: no symbol info\n");
    for (i = 0; i < textlen; i++)
      sym[i][0] = '\0';
  }
  else {
    ptr = 0;
    for (i = 1; i <= seqlen; i++)
      if (inrange(i, range)) {
	GetField(field, entry, i, symcol);
	strcpy(sym[ptr++], field);
      }
  }

  for (i = 0; i < codelen; i++)
    symstr[i] = ' ';

  ptr = 0;
  newblock = 1;
  for (j = 0; j < codelen; j++) {
    if (writecode[j] == 'c' && newblock == 1) {
      if (j == 0 || symstr[j-1] == ' ')
	for (i = 0; sym[ptr][i] != '\0'; i++)
	  symstr[j+i] = sym[ptr][i];
      newblock = 0; }
    if (writecode[j] == '\n' || writecode[j] == '.')
      newblock = 1;
    if (writecode[j] == 'c')
      ptr++;
  }
  
  ptr -= 1;
  newblock = 1;
  for (j = codelen-1; j >= 0; j--) {
    if (writecode[j] == 'c' && newblock == 1) {
      if (symstr[j-strlen(sym[ptr])] == ' ')
	for (i = 0; sym[ptr][i] != '\0'; i++)
	  symstr[1+j+i-strlen(sym[ptr])] = sym[ptr][i];
      newblock = 0; }
    if (writecode[j] == '\n' || writecode[j] == '.')
      newblock = 1;
    if (writecode[j] == 'c')
      ptr--;
  }

  i = 0;

  maxlen = 30;     /* Start with 20 lines */

  text = (char **)malloc(maxlen * sizeof(char *));
  name = (char **)malloc(maxlen * sizeof(char *));

  do {
    if (i == maxlen) {
      maxlen *= 2;
      text = (char **)realloc(text, maxlen * sizeof(char *));
      name = (char **)realloc(name, maxlen * sizeof(char *));
                                                     /* Get more space */
    }

    if ((codecol = ReadColno(entry, "residue")) == 0 &&
	(codecol = ReadColno(entry, "nucleotide")) == 0 &&
	(codecol = ReadColno(entry, "aminoacid")) == 0 &&
	(codecol = ReadColno(entry, "code")) == 0) {
      fprintf(stderr,
	      "Warning: Ignoring sequence, no column info\n");
      continue; }

    text[i] = (char *)malloc((textlen + 1) * sizeof(char));
    name[i] = (char *)malloc(MAXNAME * sizeof(char));
    if (ReadText(entry, "SEQUENCE_NAME", MAXNAME, name[i]) != 0 &&
	ReadText(entry, "SEQUENCE", MAXNAME, name[i]) != 0 &&
	ReadText(entry, "ENTRY", MAXNAME, name[i]) != 0)
      name[i][0] = '\0';

    ptr = 0;
    for (j = 1; j <= seqlen; j++)
      if (inrange(j, range)) {
	GetField(field, entry, j, codecol);
	text[i][ptr++] = field[0];
      }

    text[i][textlen] = '\0';

    i++;
  } while ((read_error = ReadEntry(fp, entry)) == 0);

  if (read_error == 1)
    return 1;

  numseq = i;
  text = (char **)realloc(text, numseq * sizeof(char *));
                                             /* Free excess space */

  blockpage = 10000;
  if (blockpage == 0)
    fprintf(stderr, "col2txtalign: Warning: too many sequences on page\n");
  startpos = startptr = 0;
  block = 0;
  page = 1;
  lineno = 0;
  maxcolno = 0;
  do {
    newblock = 1;
    ptr = startptr;
    for (j = startpos; writecode[j] != '\n'; j++) {
      if (writecode[j] == 'c' && newblock == 1)
	newblock = 0;
      if (writecode[j] == 'c' &&
	  ((writecode[j+1] == ' ' && writecode[j+2] == '.') ||
	   writecode[j+1] == '\n')) {
	newblock = 0; }
      if (writecode[j] == '.')
	newblock = 1;
      if (writecode[j] == 'c')
	ptr++;
    }

    lineno++;

    for (i = 0; i < numseq; i++) {
      ptr = startptr;
      for (j = startpos; writecode[j] != '\n'; j++) {
	if (writecode[j] == 'c') {
	  ptr++;
	}
      }
    }
    for (i = 0; i < numseq; i++) {
      ptr = startptr;
      colno = 0;
      for (j = startpos; writecode[j] != '\n'; j++) {
	colno++;
	if (writecode[j] == 'c')
	  ptr++;
      }
      if (colno > maxcolno)
	maxcolno = colno;   /* for BoundingBox */
      lineno++;
    }
    startpos = j+1;
    startptr = ptr;
    lineno++;
    block++;
    if (block >= blockpage) {
      page++;
      block = 0;
    }
  } while (writecode[startpos] != '\0');

  startpos = startptr = 0;
  block = 0;
  page = 1;
  lineno = 0;
  maxcolno = 0;
  do {
    newblock = 1;
    ptr = startptr;
    for (i = 0; i <= nametxtwidth; i++)
      printf(" ");
    spaceptr1 = spaceptr2 = 0;
    for (j = startpos; writecode[j] != '\n'; j++) {
      if (writecode[j] == 'c' && newblock == 1) {
	diff = spaceptr2-spaceptr1;
	if (diff >= 1 || (spaceptr1 == 0 && spaceptr2 == 0)) {
	  for (i = 0; i < diff; i++)
	    printf(" ");
	  printf("%s", sym[ptr]);
	  spaceptr1 = strlen(sym[ptr]);
	  spaceptr2 = 0;
	}
	newblock = 0; }
      spaceptr2++;
      if (writecode[j] == 'c' &&
	  ((writecode[j+1] == ' ' && writecode[j+2] == '.') ||
	   writecode[j+1] == '\n')) {
	diff = spaceptr2-spaceptr1-strlen(sym[ptr]);
	if (diff >= 1) {
	  for (i = 0; i < diff; i++)
	    printf(" ");
	  printf("%s", sym[ptr]);
	  spaceptr1 = 0;
	  spaceptr2 = 0;
	}
	newblock = 0; }
      if (writecode[j] == '.')
	newblock = 1;
      if (writecode[j] == 'c')
	ptr++;
    }

    lineno++;

    printf("\n");
    for (i = 0; i < numseq; i++) {
      sprintf(format, "%%-%d.%ds ", nametxtwidth, nametxtwidth);
      printf(format, name[i]);
      ptr = startptr;
      colno = 0;
      for (j = startpos; writecode[j] != '\n'; j++) {
	if (writecode[j] == 'c') {
	  if (text[i][ptr] == '(' || text[i][ptr] == ')')
	    printf("%c", text[i][ptr]);
	  else
	    printf("%c", text[i][ptr]);
	  colno++;
	}
	else {
	  if (writecode[j] == '(' || writecode[j] == ')')
	    printf("%c", writecode[j]);
	  else
	    printf("%c", writecode[j]);
	  colno++;
	}
	if (writecode[j] == 'c')
	  ptr++;
      }
      if (colno > maxcolno)
	maxcolno = colno;   /* for BoundingBox */
      printf("\n");
      lineno++;
    }
    startpos = j+1;
    startptr = ptr;
    printf("\n");
    lineno++;
    block++;
    if (block >= blockpage) {
      page++;
      block = 0;
    }
  } while (writecode[startpos] != '\0');

  return 0;
}

void usage(void)
{
  fprintf(stderr,
	  "usage: col2txtalign\n"
	  "            [-r<range> | -range=<range>] [--namewidth=<namewidth>]\n"
	  "            [--textwidth=<textwidth>] [--space] [<file>]\n");
}

/* returns true if num is in range */
int inrange(int num, char *range)
{
  int ptr, ptr2;
  int low, high;

  if (range == NULL)
    return 1;

  for (ptr = 0; range[ptr] != '\0';) {
    for (ptr2 = ptr; range[ptr2] != '\0' && range[ptr2] != ','; ptr2++)
      ;
    if (sscanf(&range[ptr], " %d - %d", &low, &high) != 2) {
      if (sscanf(&range[ptr], " %d", &low) == 1)
	high = low;
      else
	continue;  /* Just continue if error in range */
    }
    if (low <= num && num <= high)
      return 1;
    ptr = ptr2;
    if (range[ptr] == ',')
      ptr++;
  }

  return 0;
}

