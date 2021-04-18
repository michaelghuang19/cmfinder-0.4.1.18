/*********************************************************************

  extendstem.c

  usage: extendstem
            [-c   |   --color   |  -C<color_r>,<color_g>,<color_b>  |
            --color=<color_r>,<color_g>,<color_b>] [-D<color_r>,<co-
            lor_g>,<color_b> | --color2=<color_r>,<color_g>,<color_b>]
            [-f   |  --fgcolor  |  -F<color_r>,<color_g>,<color_b>  |
            --fgcolor=<color_r>,<color_g>,<color_b>] [-G<color_r>,<co-
            lor_g>,<color_b> | --fgcolor2=<color_r>,<color_g>,<color_b>]
            [-g |  --gupair]  [-u  |  --toupper]  [-s<support>  |
            --support=<support>] [<file>]

  See man page for more info

  000215 Bjarne Knudsen (bk@daimi.au.dk)

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

typedef struct tagColno {
  int align_bp, alignpos;
  int seq_bp, seqpos;
  int nuc;
  int sup;
  int red, green, blue;
  int red2, green2, blue2;
} Colno;

typedef struct tagColor {
  double r, g, b;
} Color;

void usage(void);
int prevpair(Entry *entry, int *pos, int *pair, Colno *colno);
int validpair(Entry *entry, int pos, int pair, Colno *colno);
int support(Entry *entry, int pos, int pair, int option_sup, Colno *colno);

int main(int argc, char **argv)
{
  FILE *fp;
  Header *header;
  Entry *entry;
  Colno *colno;     /* For column numbers */
  int read_error;   /* For keeping track of errors in reading entries */
  int option_c, option_g, option_u, option_fg;
                    /* = 1 if option is active */
  int i;
  int *newpair1, *newpair2, *newpair3;
                    /* For new structure. Two arrays are needed as new
                       structure may be ambiguous */
  int len;          /* For sequence length */
  char field[MAXCOLW];
  int pos, pair, oldpos, oldpair;
  char nuc1, nuc2;
  int change;
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */
  Color *color, *color2;
  Color *fgcolor, *fgcolor2;
  int option_sup;
  int pp;

  colno = (Colno *) malloc(sizeof(Colno));
  color = (Color *) malloc(sizeof(Color));
  color2 = (Color *) malloc(sizeof(Color));
  fgcolor = (Color *) malloc(sizeof(Color));
  fgcolor2 = (Color *) malloc(sizeof(Color));

  /* default options */
  option_c = 0;
  option_g = 0;
  option_u = 0;
  option_fg = 0;
  color->r = 0.5;
  color->g = 0.5;
  color->b = 1.0;

  color2->r = 0.0;
  color2->g = 1.0;
  color2->b = 1.0;

  fgcolor->r = 0.5;
  fgcolor->g = 0.5;
  fgcolor->b = 1.0;

  fgcolor2->r = 0.0;
  fgcolor2->g = 1.0;
  fgcolor2->b = 1.0;
  option_sup = -1;

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
    else if (strcmp(s, "g") == 0)
      option_g = 1;
    else if (strcmp(s, "-gupair") == 0)
      option_g = 1;
    else if (strcmp(s, "u") == 0)
      option_u = 1;
    else if (strcmp(s, "-toupper") == 0)
      option_u = 1;
    else if (strncmp(s, "C", 1) == 0) {
      if (sscanf(&s[1], "%lf,%lf,%lf%n",
		 &color->r, &color->g, &color->b, &i) != 3 ||
	  i+1 != strlen(s)) {
	usage();
	return 1; }
      option_c = 1;
    } 
    else if (strncmp(s, "-color=", 7) == 0) {
      if (sscanf(&s[7], "%lf,%lf,%lf%n",
		 &color->r, &color->g, &color->b, &i) != 3 ||
	  i+7 != strlen(s)) {
	usage();
	return 1; }
      option_c = 1;
    }
    else if (strncmp(s, "F", 1) == 0) {
      if (sscanf(&s[1], "%lf,%lf,%lf%n",
		 &fgcolor->r, &fgcolor->g, &fgcolor->b, &i) != 3 ||
	  i+1 != strlen(s)) {
	usage();
	return 1; }
      option_fg = 1;
    } 
    else if (strncmp(s, "-fgcolor=", 9) == 0) {
      if (sscanf(&s[9], "%lf,%lf,%lf%n",
		 &fgcolor->r, &fgcolor->g, &fgcolor->b, &i) != 3 ||
	  i+9 != strlen(s)) {
	usage();
	return 1; }
      option_fg = 1;
    }
    else if (strncmp(s, "s", 1) == 0) {
      if (sscanf(s, "s%d%n", &option_sup, &i) != 1 ||
	  i != strlen(s)) {
	usage();
	return 1; }
    }
    else if (strncmp(s, "-support=", 1) == 0) {
      if (sscanf(&s[9], "%d%n", &option_sup, &i) != 1 ||
	  i+9 != strlen(s)) {
	usage();
	return 1; }
    }
    else if (strncmp(s, "D", 1) == 0) {
      if (sscanf(&s[1], "%lf,%lf,%lf%n",
		 &color2->r, &color2->g, &color2->b, &i) != 3 ||
	  i+1 != strlen(s)) {
	usage();
	return 1; }
      option_c = 1;
    }
    else if (strncmp(s, "-color2=", 8) == 0) {
      if (sscanf(&s[8], "%lf,%lf,%lf%n",
		 &color2->r, &color2->g, &color2->b, &i) != 3 ||
	  i+8 != strlen(s)) {
	usage();
	return 1; }
      option_c = 1;
    }
    else if (strncmp(s, "G", 1) == 0) {
      if (sscanf(&s[1], "%lf,%lf,%lf%n",
		 &fgcolor2->r, &fgcolor2->g, &fgcolor2->b, &i) != 3 ||
	  i+1 != strlen(s)) {
	usage();
	return 1; }
      option_fg = 1;
    }
    else if (strncmp(s, "-fgcolor2=", 10) == 0) {
      if (sscanf(&s[10], "%lf,%lf,%lf%n",
		 &fgcolor2->r, &fgcolor2->g, &fgcolor2->b, &i) != 3 ||
	  i+10 != strlen(s)) {
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
    fprintf(stderr, "extendstem: Error in opening file '%s'\n", s);
    return 1; }

  header = MakeHeader();
  entry = MakeEntry();

  if (ReadHeader(fp, header) != 0)
    return 1;

  AddHeaderInfo(header, argc, argv);

  PrintHeader(stdout, header);

  newpair1 = newpair2 = newpair3 = NULL;

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
    colno->sup = ReadColno(entry, "support");
    if (colno->nuc == 0 ||
	((colno->align_bp == 0 || colno->alignpos == 0) &&
	 (colno->seq_bp == 0 || colno->seqpos == 0)) ||
	(option_sup != -1 && colno->sup == 0)) {
	  fprintf(stderr,
	      "extendstem: Warning: Ignoring sequence, no column info\n");
      PrintEntry(stdout, entry);
      continue; }

    len = EntryLength(entry);

    newpair1 = (int *)realloc(newpair1, (len+1) * sizeof(int));
    newpair2 = (int *)realloc(newpair2, (len+1) * sizeof(int));
    newpair3 = (int *)realloc(newpair3, (len+1) * sizeof(int));
                           /* len+1 since 1 to len is used */

    for (i = 1; i <= len; i++)
      newpair1[i] = newpair2[i] = newpair3[i] = 0;

    for (i = 1; i <= len; i++) {
      pos = i;
      pair = 0;
      change = 1;
      while (prevpair(entry, &pos, &pair, colno) == 0 &&
	     change == 1) {
	change = 0;
	GetField(field, entry, pos, colno->nuc);
	nuc1 = field[0];
	GetField(field, entry, pair, colno->nuc);
	nuc2 = field[0];
	if (StdPair4(nuc1, nuc2) ||
	    (option_g == 1 && StdPair6(nuc1, nuc2))) {
	  /* ...A(((...)))U... */
	  if (pair > pos &&
	      (option_sup == -1 ||
	       support(entry, pos, pair, option_sup, colno))) {
	    newpair1[pos] = pair;
	    newpair3[pair] = pos;
	    change = 1;
	  }
	  else if (validpair(entry, pair, pos, colno) &&
		   (option_sup == -1 ||
		    support(entry, pos, pair, option_sup, colno))) {
	    newpair2[pos] = pair;
	    newpair2[pair] = pos;
	    change = 1;
	  }
	}
	else if (StdPair6(nuc1, nuc2)) {
	  /* ...G(((...)))U... */
	  oldpos = pos;
	  oldpair = pair;
	  if ((pp = prevpair(entry, &pos, &pair, colno)) == 0) {
	    GetField(field, entry, pos, colno->nuc);
	    nuc1 = field[0];
	    GetField(field, entry, pair, colno->nuc);
	    nuc2 = field[0];
	    if (StdPair4(nuc1, nuc2)) {
	      /* ..CG(((...)))UG.. */
	      if (pair > pos) {
		if (option_sup == -1 ||
		    support(entry, oldpos, oldpair, option_sup, colno)) {
		  newpair1[oldpos] = oldpair;
		  newpair3[oldpair] = oldpos;
		  change = 1;
		}
		if (option_sup == -1 ||
		    support(entry, pos, pair, option_sup, colno)) {
		  newpair1[pos] = pair;
		  newpair3[pair] = pos;
		  change = 1;
		}
	      }
	      else if (validpair(entry, pair, pos, colno)) {
		if (option_sup == -1 ||
		    support(entry, oldpos, oldpair, option_sup, colno)) {
		  newpair2[oldpos] = oldpair;
		  newpair2[oldpair] = oldpos;
		  change = 1;
		}
		if (option_sup == -1 ||
		    support(entry, pos, pair, option_sup, colno)) {
		  newpair2[pos] = pair;
		  newpair2[pair] = pos;
		  change = 1;
		}
	      }
	    }
	  }
	  else if (pp == 1) {
	    /* ..(G(((...)))U).. */
	    if (oldpair > oldpos &&
		(option_sup == -1 ||
		 support(entry, oldpos, oldpair, option_sup, colno))) {
	      newpair1[oldpos] = oldpair;
	      newpair3[oldpair] = oldpos;
	      change = 1;
	    }
	    else if (validpair(entry, oldpair, oldpos, colno) &&
		     (option_sup == -1 ||
		      support(entry, oldpos, oldpair, option_sup, colno))) {
	      newpair2[oldpos] = oldpair;
	      newpair2[oldpair] = oldpos;
	      change = 1;
	    }
	  }
	}
      }
    }

    if (option_c == 1 || option_fg == 1) {
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
      for (i = 1; i <= len; i++) {
	if (newpair1[i] != 0)
	  if ((newpair2[i] != 0 && newpair1[i] != newpair2[i]) ||
	      (newpair3[i] != 0 && newpair1[i] != newpair3[i])) {
	    if (option_c == 1) {
	      ChgField(entry, i, colno->red, "%5.3f", color2->r);
	      ChgField(entry, i, colno->green, "%5.3f", color2->g);
	      ChgField(entry, i, colno->blue, "%5.3f", color2->b); 
	    }
	    if (option_fg == 1) {
	      ChgField(entry, i, colno->red2, "%5.3f", fgcolor2->r);
	      ChgField(entry, i, colno->green2, "%5.3f", fgcolor2->g);
	      ChgField(entry, i, colno->blue2, "%5.3f", fgcolor2->b); 
	    }
	  }
	  else {
	    if (option_c == 1) {
	      ChgField(entry, i, colno->red, "%5.3f", color->r);
	      ChgField(entry, i, colno->green, "%5.3f", color->g);
	      ChgField(entry, i, colno->blue, "%5.3f", color->b); 
	    }
	    if (option_fg == 1) {
	      ChgField(entry, i, colno->red2, "%5.3f", fgcolor->r);
	      ChgField(entry, i, colno->green2, "%5.3f", fgcolor->g);
	      ChgField(entry, i, colno->blue2, "%5.3f", fgcolor->b); 
	    }
	  }
	else if (newpair2[i] != 0) {
	  if (newpair3[i] != 0 && newpair2[i] != newpair3[i]) {
	    if (option_c == 1) {
	      ChgField(entry, i, colno->red, "%5.3f", color2->r);
	      ChgField(entry, i, colno->green, "%5.3f", color2->g);
	      ChgField(entry, i, colno->blue, "%5.3f", color2->b); 
	    }
	    if (option_fg == 1) {
	      ChgField(entry, i, colno->red2, "%5.3f", fgcolor2->r);
	      ChgField(entry, i, colno->green2, "%5.3f", fgcolor2->g);
	      ChgField(entry, i, colno->blue2, "%5.3f", fgcolor2->b); 
	    }
	  }
	  else {
	    if (option_c == 1) {
	      ChgField(entry, i, colno->red, "%5.3f", color->r);
	      ChgField(entry, i, colno->green, "%5.3f", color->g);
	      ChgField(entry, i, colno->blue, "%5.3f", color->b); 
	    }
	    if (option_fg == 1) {
	      ChgField(entry, i, colno->red2, "%5.3f", fgcolor->r);
	      ChgField(entry, i, colno->green2, "%5.3f", fgcolor->g);
	      ChgField(entry, i, colno->blue2, "%5.3f", fgcolor->b); 
	    }
	  }
	}
	else if (newpair3[i] != 0) {
	  if (option_c == 1) {
	    ChgField(entry, i, colno->red, "%5.3f", color->r);
	    ChgField(entry, i, colno->green, "%5.3f", color->g);
	    ChgField(entry, i, colno->blue, "%5.3f", color->b); 
	  }
	  if (option_fg == 1) {
	    ChgField(entry, i, colno->red2, "%5.3f", fgcolor->r);
	    ChgField(entry, i, colno->green2, "%5.3f", fgcolor->g);
	    ChgField(entry, i, colno->blue2, "%5.3f", fgcolor->b); 
	  }
	}
      }
    }
    else
      for (i = 1; i <= len; i++) {
	if (newpair1[i] != 0) {
	  if ((newpair2[i] != 0 && newpair1[i] != newpair2[i]) ||
	      (newpair3[i] != 0 && newpair1[i] != newpair3[i]))
	    ;
	  else if ((newpair1[newpair1[i]] == 0 ||
		    newpair1[newpair1[i]] == i) &&
		   (newpair2[newpair1[i]] == 0 ||
		    newpair2[newpair1[i]] == i) &&
		   (newpair3[newpair1[i]] == 0 ||
		    newpair3[newpair1[i]] == i)) {
	    SetPair(entry, i, newpair1[i], colno->align_bp, colno->alignpos,
		    colno->seq_bp, colno->seqpos);
	    if (option_u == 1) {
	      GetField(field, entry, i, colno->nuc);
	      field[0] = toupper(field[0]);
	      ChgField(entry, i, colno->nuc, "%s", field);
	    }
	  }
	}
	else if (newpair2[i] != 0) {
	  if ((newpair3[i] != 0 && newpair2[i] != newpair3[i]))
	    ;
	  else if ((newpair1[newpair2[i]] == 0 ||
		    newpair1[newpair2[i]] == i) &&
		   (newpair2[newpair2[i]] == 0 ||
		    newpair2[newpair2[i]] == i) &&
		   (newpair3[newpair2[i]] == 0 ||
		    newpair3[newpair2[i]] == i)) {
	    SetPair(entry, i, newpair2[i], colno->align_bp, colno->alignpos,
		    colno->seq_bp, colno->seqpos);
	    if (option_u == 1) {
	      GetField(field, entry, i, colno->nuc);
	      field[0] = toupper(field[0]);
	      ChgField(entry, i, colno->nuc, "%s", field);
	    }
	  }
	}
	else if (newpair3[i] != 0 &&
		 (newpair1[newpair3[i]] == 0 ||
		  newpair1[newpair3[i]] == i) && 
		 (newpair2[newpair3[i]] == 0 ||
		  newpair2[newpair3[i]] == i) &&
		 (newpair3[newpair3[i]] == 0 ||
		  newpair3[newpair3[i]] == i)) {
	  SetPair(entry, i, newpair3[i], colno->align_bp, colno->alignpos,
		  colno->seq_bp, colno->seqpos);
	  if (option_u == 1) {
	    GetField(field, entry, i, colno->nuc);
	    field[0] = toupper(field[0]);
	    ChgField(entry, i, colno->nuc, "%s", field);
	  }
	}
      }
    PrintEntry(stdout, entry);
  }

  if (fp != stdin && fclose(fp) != 0) {
    fprintf(stderr, "extendstem: Error in closing file\n");
    return 1; }

  if (read_error == 1)
    return 1;

  return 0;
}

void usage(void)
{
  fprintf(stderr,
	  "usage: extendstem\n"
	  "            [-c   |   --color   |  -C<color_r>,<color_g>,<color_b>  |\n"
	  "            --color=<color_r>,<color_g>,<color_b>] [-D<color_r>,<co-\n"
	  "            lor_g>,<color_b> | --color2=<color_r>,<color_g>,<color_b>]\n"
	  "            [-f   |  --fgcolor  |  -F<color_r>,<color_g>,<color_b>  |\n"
	  "            --fgcolor=<color_r>,<color_g>,<color_b>] [-G<color_r>,<co-\n"
	  "            lor_g>,<color_b> | --fgcolor2=<color_r>,<color_g>,<color_b>]\n"
	  "            [-g |  --gupair]  [-u  |  --toupper]  [-s<support>  |\n"
	  "            --support=<support>] [<file>]\n");
}

int prevpos(Entry *entry, int *pos, Colno *colno);
int nextpos(Entry *entry, int *pos, Colno *colno);

/*
   prevpair has two modes:

   *pair = 0 on call: if pos has no pair, return -1. Otherwise set
                      *pair to its pair and call prevpair.

   Otherwise: return -1 if pos >= pair, move pos to previous position
              and pair to next position (non gaps). Return 0 if none
              of them pair, 1 if they pair with eachother,-1
              otherwise.
*/
int prevpair(Entry *entry, int *pos, int *pair, Colno *colno)
{
  int p;

  if (*pair == 0) {
    p = FindPair(entry, *pos, colno->align_bp, colno->alignpos,
	                      colno->seq_bp, colno->seqpos);
    if (p != 0)    /* A pair is present */
      *pair = p;
    else
      return -1;
  }

  if (prevpos(entry, pos, colno) != 0)
    return -1;

  if (nextpos(entry, pair, colno) != 0)
    return -1;


  p = FindPair(entry, *pos, colno->align_bp, colno->alignpos,
	                    colno->seq_bp, colno->seqpos);
  if (p != 0) {    /* A pair is present */
    if (p == *pair)
      return 1;
    else
      return -1;
  }
  else {
    p = FindPair(entry, *pair, colno->align_bp, colno->alignpos,
	                       colno->seq_bp, colno->seqpos);
    if (p != 0)
      return -1;
    else
      return 0;
  }
}

/*
   Returns true if pos and pair has at least 3 nucletides between
   them. pos must be less than pair.
*/
int validpair(Entry *entry, int pos, int pair, Colno *colno)
{
  if (nextpos(entry, &pos, colno) == 1)
    return 0;

  if (nextpos(entry, &pos, colno) == 1)
    return 0;

  if (nextpos(entry, &pos, colno) == 1)
    return 0;

  if (pos >= pair)
    return 0;

  return 1;
}

int support(Entry *entry, int pos, int pair, int option_sup, Colno *colno)
{
  char field[MAXCOLW];
  int i;
  int pair2;
  int p;

  if (option_sup == -1)
    return 1;

  pair2 = pair;

  /* How far can pair go */
  for (i = 0; i < option_sup && prevpos(entry, &pair2, colno) == 0; i++)
    ;

  for (i = 0; i < option_sup; i++) {
    if (nextpos(entry, &pos, colno) != 0)
      return 0;
    GetField(field, entry, pos, colno->sup);
    if (strcmp(field, "S") == 0) { /* support is found */
      p = FindPair(entry, pos, colno->align_bp, colno->alignpos,
	                       colno->seq_bp, colno->seqpos);
      if (p < pair && p >= pair2) /* within bound */
	return 1;
    }
  }

  return 0;
}

/*
   Moves pos to previous non-gap position. Returns 0 if ok, 1 if out
   of sequence.
*/
int prevpos(Entry *entry, int *pos, Colno *colno)
{
  char field[MAXCOLW];
  int oldpos;

  oldpos = *pos;

  for (*pos -= 1; InSeq(entry, *pos); (*pos)--) {
    GetField(field, entry, *pos, colno->nuc);
    if (field[0] != '-')
      return 0;
  }

  *pos = oldpos;

  return 1;
}

/*
   Moves pos to next non-gap position. Returns 0 if ok, 1 if out of
   sequence.
*/
int nextpos(Entry *entry, int *pos, Colno *colno)
{
  char field[MAXCOLW];
  int oldpos;

  oldpos = *pos;

  for (*pos += 1; InSeq(entry, *pos); (*pos)++) {
    GetField(field, entry, *pos, colno->nuc);
    if (field[0] != '-')
      return 0;
  }

  *pos = oldpos;

  return 1;
}
