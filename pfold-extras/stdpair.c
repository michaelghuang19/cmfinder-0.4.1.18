/*********************************************************************

  stdpair.c

  usage: stdpair
              [-c | --color | -C<color_r>,<color_g>,<color_b> |
              --color=<color_r>,<color_g>,<color_b>] [-f |
	      --fgcolor | -F<color_r>,<color_g>,<color_b> |
              --fgcolor=<color_r>,<color_g>,<color_b>] [<file>]

  See man page for more info

  000630 Bjarne Knudsen (bk@daimi.au.dk)

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
#include "rna.h"
#include "file.h"

typedef struct tagColno {
  int align_bp, alignpos;
  int seq_bp, seqpos;
  int nuc;
  int red, green, blue;
  int red2, green2, blue2;
} Colno;

typedef struct tagColor {
  double r, g, b;
} Color;

void usage(void);
void stdpair(Entry *entry, Colno *colno);
void colorpair(Entry *entry, Colno *colno, Color *color, Color *color2);

int main(int argc, char **argv)
{
  FILE *fp;         /* For input */
  Header *header;
  Entry *entry;
  Colno *colno;     /* Column numbers */
  int read_error;   /* For keeping track of errors in reading entries */
  int option_c;     /* = 1 if option c is active */
  int option_fg;    /* = 1 if option f is active */
  unsigned int len; /* Used with -C option */
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */
  Color *color;     /* The colour of non-standard pairs */
  Color *color2;    /* The foreground colour of non-standard pairs */

  colno = (Colno *) malloc(sizeof(Colno));
  color = (Color *) malloc(sizeof(Color));
  color2 = (Color *) malloc(sizeof(Color));

  /* default options */
  option_c = 0;
  color->r = 1.0;
  color->g = 0.0;
  color->b = 0.0;

  option_fg = 0;
  color2->r = 1.0;
  color2->g = 0.0;
  color2->b = 0.0;

  cmdarg = InitArgument(argc, argv);

  while ((s = GetArgument(cmdarg)) != NULL)
    if (strcmp(s, "c") == 0)
      option_c = 1;
    else if (strcmp(s, "-color") == 0)
      option_c = 1;
    else if (strcmp(s, "f") == 0)
      option_fg = 1;
    else if (strcmp(s, "-fgcolor") == 0)
      option_fg = 1;
    else if (strncmp(s, "C", 1) == 0) {
      if (sscanf(&s[1], "%lf,%lf,%lf%n",
		 &color->r, &color->g, &color->b, &len) != 3 ||
	  len+1 != strlen(s)) {
	usage();
	return 1; }
      option_c = 1;
    }
    else if (strncmp(s, "-color=", 7) == 0) {
      if (sscanf(&s[7], "%lf,%lf,%lf%n",
		 &color->r, &color->g, &color->b, &len) != 3 ||
	  len+7 != strlen(s)) {
	usage();
	return 1; }
      option_c = 1;
    }
    else if (strncmp(s, "F", 1) == 0) {
      if (sscanf(&s[1], "%lf,%lf,%lf%n",
		 &color2->r, &color2->g, &color2->b, &len) != 3 ||
	  len+1 != strlen(s)) {
	usage();
	return 1; }
      option_fg = 1;
    }
    else if (strncmp(s, "-fgcolor=", 9) == 0) {
      if (sscanf(&s[9], "%lf,%lf,%lf%n",
		 &color2->r, &color2->g, &color2->b, &len) != 3 ||
	  len+9 != strlen(s)) {
	usage();
	return 1; }
      option_fg = 1;
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
    fprintf(stderr, "stdpair: Error in opening file '%s'\n", s);
    return 1; }

  header = MakeHeader();
  entry = MakeEntry();

  if (ReadHeader(fp, header) != 0)
    return 1;

  AddHeaderInfo(header, argc, argv);

  PrintHeader(stdout, header);

  while ((read_error = ReadEntry(fp, entry)) == 0) {
    if (!ReadType(entry, "RNA") && !ReadType(entry, "DNA")) {
      PrintEntry(stdout, entry);
      continue;
    }
    if ((colno->nuc = ReadColno(entry, "residue")) == 0)
      colno->nuc = ReadColno(entry, "nucleotide");
    colno->align_bp = ReadColno(entry, "align_bp");
    colno->alignpos = ReadColno(entry, "alignpos");
    colno->seq_bp = ReadColno(entry, "seq_bp");
    colno->seqpos = ReadColno(entry, "seqpos");
    if (colno->nuc == 0 ||
	((colno->align_bp == 0 || colno->alignpos == 0) &&
	 (colno->seq_bp == 0 || colno->seqpos == 0))) {
      fprintf(stderr, "stdpair: Warning: Ignoring sequence,"
	                " insufficient column info\n");
      PrintEntry(stdout, entry);
      continue; }

    colno->red = colno->red2 = 0;
    if (option_c == 1) {
      colno->red = EnsureCol(entry, "color_r", "1.000");
      colno->green = EnsureCol(entry, "color_g", "1.000");
      colno->blue = EnsureCol(entry, "color_b", "1.000");
    }
    if (option_fg == 1) {
      colno->red2 = EnsureCol(entry, "color2_r", "0.000");
      colno->green2 = EnsureCol(entry, "color2_g", "0.000");
      colno->blue2 = EnsureCol(entry, "color2_b", "0.000");
    }
    if (option_c == 1 || option_fg == 1)
      colorpair(entry, colno, color, color2);
    else
      stdpair(entry, colno);

    PrintEntry(stdout, entry);
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
	  "usage: stdpair\n"
	  "            [-c | --color | -C<color_r>,<color_g>,<color_b> |\n"
	  "            --color=<color_r>,<color_g>,<color_b>] [-f |\n"
	  "            --fgcolor | -F<color_r>,<color_g>,<color_b> |\n"
	  "            --fgcolor=<color_r>,<color_g>,<color_b>] [<file>]\n");
}

void stdpair(Entry *entry, Colno *colno)
{
  int i;
  char field[MAXCOLW];
  int pair;
  char nuc1, nuc2;
  int len;

  len = EntryLength(entry);

  for (i = 1; i <= len; i++) {
    pair = FindPair(entry, i, colno->align_bp, colno->alignpos,
		    colno->seq_bp, colno->seqpos);
    if (pair != 0) {  /* A pair is present */
      GetField(field, entry, i, colno->nuc);
      nuc1 = field[0];
      GetField(field, entry, pair, colno->nuc);
      nuc2 = field[0];
      if (!StdPair6(nuc1, nuc2))
	SetPair(entry, i, 0, colno->align_bp, colno->alignpos,
		colno->seq_bp, colno->seqpos);
    }
  }
}

void colorpair(Entry *entry, Colno *colno, Color *color, Color *color2)
{
  int i;
  char field[MAXCOLW];
  int pair;
  char nuc1, nuc2;
  int len;

  len = EntryLength(entry);

  for (i = 1; i <= len; i++) {
    pair = FindPair(entry, i, colno->align_bp, colno->alignpos,
		    colno->seq_bp, colno->seqpos);
    if (pair != 0) {  /* A pair is present */
      GetField(field, entry, i, colno->nuc);
      nuc1 = field[0];
      GetField(field, entry, pair, colno->nuc);
      nuc2 = field[0];
      if (!StdPair6(nuc1, nuc2)) {
	if (colno->red != 0) {
	  ChgField(entry, i, colno->red, "%5.3f", color->r);
	  ChgField(entry, i, colno->green, "%5.3f", color->g);
	  ChgField(entry, i, colno->blue, "%5.3f", color->b);
	}
	if (colno->red2 != 0) {
	  ChgField(entry, i, colno->red2, "%5.3f", color2->r);
	  ChgField(entry, i, colno->green2, "%5.3f", color2->g);
	  ChgField(entry, i, colno->blue2, "%5.3f", color2->b);
	}
      }
    }
  }
}
