/*********************************************************************

  col2psalign.c

  usage: col2psalign
              [-ple] [-r<range> | --range=<range>] [-f<fontsize>
	      | --fit] [--space] [--width=<width>] [--height=<height>]
	      [--margin=<margin>] [--topmargin=<topmargin>]
	      [--sidemargin=<sidemargin>] [--namewidth=<namewidth>]
	      [--textwidth=<textwidth>] [--lineheight=<lineheight>]
	      [--figure] [<file>]

  See man page for more info

  000922 Bjarne Knudsen (bk@daimi.au.dk)

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
void PrintString(FILE *fp, char *s);

int main(int argc, char **argv)
{
  Header *header;
  Entry *entry;
  FILE *fp;
  int codecol, symcol, red_col, green_col, blue_col;
  int red_col2, green_col2, blue_col2;
  int read_error;   /* for keeping track of errors in reading entries */
  int i, j, k;
  int width, height;
  char *range;
  int textlen;
  char ***text;
  double **color_r;
  double **color_g;
  double **color_b;
  double **fg_color_r;
  double **fg_color_g;
  double **fg_color_b;
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
  double temp;
  double topmargin, sidemargin, pagewidth, pageheight;
  double charwidth, tlineheight, fontsize, namewidth;
  double linewidth;
  double lineheight;
  double nametxtwidth;
  unsigned int len; /* Used with options */
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */
  int option_space;
  int option_fill;
  int option_fit;
  int option_figure;
  int option_manynum;
  int landscape;
  int lineno;
  int colno, maxcolno;
  int txtwidth;
  int new_showpage;

  /* Postscript variables */
  topmargin = 20;
  sidemargin = 20;
  pagewidth = 21.0/2.54*72;   /* A4 is 21.0 cm wide */
  pageheight = 29.7/2.54*72;  /* A4 is 29.7 cm high */
  fontsize = 8;
  option_space = 0;
  option_fill = 1;
  landscape = 1;
  option_fit = 0;
  option_figure = 0;
  option_manynum = 0;
  nametxtwidth = 8;
  txtwidth = 0;
  lineheight = 1;

  range = NULL;

  cmdarg = InitArgument(argc, argv);

  while ((s = GetArgument(cmdarg)) != NULL)
    if (strncmp(s, "r", 1) == 0)
      range = &s[1];
    else if (strncmp(s, "-range=", 7) == 0)
      range = &s[7];
    else if (strcmp(s, "-empty") == 0)
      option_fill = 0;
    else if (strcmp(s, "e") == 0)
      option_fill = 0;
    else if (strcmp(s, "-portrait") == 0)
      landscape = 0;
    else if (strcmp(s, "p") == 0)
      landscape = 0;
    else if (strcmp(s, "-space") == 0)
      option_space = 1;
    else if (strcmp(s, "-fit") == 0)
      option_fit = 1;
    else if (strcmp(s, "-figure") == 0)
      option_figure = 1;
    else if (strcmp(s, "-manynum") == 0)
      option_manynum = 1;
    else if (strcmp(s, "-letter") == 0) {
      pagewidth = 8.5*72;   /* letter is 8.5 inch wide */
      pageheight = 11*72;  /* A4 is 11 inch high */
    }
    else if (strcmp(s, "l") == 0) {
      pagewidth = 8.5*72;   /* letter is 8.5 inch wide */
      pageheight = 11*72;  /* A4 is 11 inch high */
    }
    else if (strncmp(s, "f", 1) == 0) {
      if (sscanf(s, "f%lf%n", &fontsize, &len) != 1 ||
	  len != strlen(s)) {
	usage();
	return 1; }
    }  
    else if (strncmp(s, "-fontsize=", 10) == 0) {
      if (sscanf(s, "-fontsize=%lf%n", &fontsize, &len) != 1 ||
	  len != strlen(s)) {
	usage();
	return 1; }
    }  
    else if (strncmp(s, "-lineheight=", 12) == 0) {
      if (sscanf(s, "-lineheight=%lf%n", &lineheight, &len) != 1 ||
	  len != strlen(s)) {
	usage();
	return 1; }
    }
    else if (strncmp(s, "-namewidth=", 11) == 0) {
      if (sscanf(s, "-namewidth=%lf%n", &nametxtwidth, &len) != 1 ||
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
    else if (strncmp(s, "-width=", 7) == 0) {
      if (sscanf(s, "-width=%lf%n", &pagewidth, &len) != 1) {
	usage();
	return 1; }
      if (len != strlen(s)) {
	if (strcmp(&s[len], "in") == 0 || strcmp(&s[len], "inch") == 0)
	  pagewidth *= 72;    /* 72 points per inch */ 
	else if (strcmp(&s[len], "cm") == 0)
	  pagewidth *= 72 / 2.54;    /* 2.54 cm per inch */ 
      }
    }
    else if (strncmp(s, "-height=", 7) == 0) {
      if (sscanf(s, "-height=%lf%n", &pageheight, &len) != 1) {
	usage();
	return 1; }
      if (len != strlen(s)) {
	if (strcmp(&s[len], "in") == 0 || strcmp(&s[len], "inch") == 0)
	  pageheight *= 72;    /* 72 points per inch */ 
	else if (strcmp(&s[len], "cm") == 0)
	  pageheight *= 72 / 2.54;    /* 2.54 cm per inch */ 
      }
    }
    else if (strncmp(s, "-sidemargin=", 7) == 0) {
      if (sscanf(s, "-sidemargin=%lf%n", &sidemargin, &len) != 1) {
	usage();
	return 1; }
      if (len != strlen(s)) {
	if (strcmp(&s[len], "in") == 0 || strcmp(&s[len], "inch") == 0)
	  sidemargin *= 72;    /* 72 points per inch */ 
	else if (strcmp(&s[len], "cm") == 0)
	  sidemargin *= 72 / 2.54;    /* 2.54 cm per inch */ 
      }
    }
    else if (strncmp(s, "-topmargin=", 7) == 0) {
      if (sscanf(s, "-topmargin=%lf%n", &topmargin, &len) != 1) {
	usage();
	return 1; }
      if (len != strlen(s)) {
	if (strcmp(&s[len], "in") == 0 || strcmp(&s[len], "inch") == 0)
	  topmargin *= 72;    /* 72 points per inch */ 
	else if (strcmp(&s[len], "cm") == 0)
	  topmargin *= 72 / 2.54;    /* 2.54 cm per inch */ 
      }
    }
    else if (strncmp(s, "-margin=", 7) == 0) {
      if (sscanf(s, "-margin=%lf%n", &topmargin, &len) != 1) {
	usage();
	return 1; }
      if (len != strlen(s)) {
	if (strcmp(&s[len], "in") == 0 || strcmp(&s[len], "inch") == 0)
	  topmargin *= 72;    /* 72 points per inch */ 
	else if (strcmp(&s[len], "cm") == 0)
	  topmargin *= 72 / 2.54;    /* 2.54 cm per inch */ 
      }
      sidemargin = topmargin;
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
    fprintf(stderr, "col2psalign: Error in opening file '%s'\n", s);
    return 1; }

  if (option_figure == 1)
    landscape = 0;

  if (landscape == 1) {
    temp = pagewidth;
    pagewidth = pageheight;
    pageheight = temp;
    temp = sidemargin;
    sidemargin = topmargin;
    topmargin = temp;
  }

  if (option_figure == 0 && txtwidth == 0) {
    txtwidth = (pagewidth - nametxtwidth*fontsize - 2*sidemargin)/(0.75*fontsize);
    /*    printf("%f %f %f %d %d\n", pagewidth, sidemargin, fontsize, nametxtwidth, txtwidth);*/
  }  
  header = MakeHeader();
  entry = MakeEntry();

  if (ReadHeader(fp, header) != 0)
    return 1;

  while ((read_error = ReadEntry(fp, entry)) == 0) {
    if (!ReadType(entry, "TREE"))
      break;
  }

  if (read_error != 0) {
    fprintf(stderr, "col2psalign: Warning, no entries\n");
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

  if (option_fit == 1) {
    if (txtwidth == 0)
      fontsize = (pagewidth-2*sidemargin) / ((strlen(writecode)-1)*.75+nametxtwidth);
    else
      fontsize = (pagewidth-2*sidemargin) / (txtwidth*.75+nametxtwidth);
  }
  
  charwidth = .75*fontsize;
  tlineheight = fontsize*lineheight;
  namewidth = nametxtwidth*fontsize;
  linewidth = fontsize/8;

  if (txtwidth != 0)
    width = txtwidth;
  else {
    width = (pagewidth-2*sidemargin-namewidth)/charwidth;
    if (option_fit == 1)
      width += 1;
  }

  height = (pageheight-2*topmargin)/tlineheight;

  printf("%%!PS-Adobe-2.0\n");
    
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

  text = (char ***)malloc(maxlen * sizeof(char **));
  color_r = (double **)malloc(maxlen * sizeof(double *));
  color_g = (double **)malloc(maxlen * sizeof(double *));
  color_b = (double **)malloc(maxlen * sizeof(double *));
  fg_color_r = (double **)malloc(maxlen * sizeof(double *));
  fg_color_g = (double **)malloc(maxlen * sizeof(double *));
  fg_color_b = (double **)malloc(maxlen * sizeof(double *));
  name = (char **)malloc(maxlen * sizeof(char *));

  do {
    if (i == maxlen) {
      maxlen *= 2;
      text = (char ***)realloc(text, maxlen * sizeof(char **));
      color_r = (double **)realloc(color_r, maxlen * sizeof(double *));
      color_g = (double **)realloc(color_g, maxlen * sizeof(double *));
      color_b = (double **)realloc(color_b, maxlen * sizeof(double *));
      fg_color_r = (double **)realloc(fg_color_r, maxlen * sizeof(double *));
      fg_color_g = (double **)realloc(fg_color_g, maxlen * sizeof(double *));
      fg_color_b = (double **)realloc(fg_color_b, maxlen * sizeof(double *));
      name = (char **)realloc(name, maxlen * sizeof(char *));
                                                     /* Get more space */
    }
    if ((red_col = ReadColno(entry, "color_r")) == 0 ||
	(green_col = ReadColno(entry, "color_g")) == 0 ||
	(blue_col = ReadColno(entry, "color_b")) == 0)
      red_col = 0;  /* set red column to zero if color is missing */

    if ((red_col2 = ReadColno(entry, "color2_r")) == 0 ||
	(green_col2 = ReadColno(entry, "color2_g")) == 0 ||
	(blue_col2 = ReadColno(entry, "color2_b")) == 0)
      red_col2 = 0;  /* set red column to zero if color is missing */

    if ((codecol = ReadColno(entry, "residue")) == 0 &&
	(codecol = ReadColno(entry, "nucleotide")) == 0 &&
	(codecol = ReadColno(entry, "aminoacid")) == 0 &&
	(codecol = ReadColno(entry, "code")) == 0) {
      fprintf(stderr,
	      "Warning: Ignoring sequence, no column info\n");
      continue; }

    text[i] = (char **)malloc((textlen + 1) * sizeof(char *));
    color_r[i] = (double *)malloc(textlen * sizeof(double));
    color_g[i] = (double *)malloc(textlen * sizeof(double));
    color_b[i] = (double *)malloc(textlen * sizeof(double));
    fg_color_r[i] = (double *)malloc(textlen * sizeof(double));
    fg_color_g[i] = (double *)malloc(textlen * sizeof(double));
    fg_color_b[i] = (double *)malloc(textlen * sizeof(double));
    name[i] = (char *)malloc(MAXNAME * sizeof(char));
    if (ReadText(entry, "SEQUENCE_NAME", MAXNAME, name[i]) != 0 &&
	ReadText(entry, "SEQUENCE", MAXNAME, name[i]) != 0 &&
	ReadText(entry, "ENTRY", MAXNAME, name[i]) != 0)
      name[i][0] = '\0';

    ptr = 0;
    for (j = 1; j <= seqlen; j++)
      if (inrange(j, range)) {
	if (red_col != 0) {
	  GetField(field, entry, j, red_col);
	  color_r[i][ptr] = atof(field);
	  GetField(field, entry, j, green_col);
	  color_g[i][ptr] = atof(field);
	  GetField(field, entry, j, blue_col);
	  color_b[i][ptr] = atof(field);
	}
	else {
	  color_r[i][ptr] = 1;
	  color_g[i][ptr] = 1;
	  color_b[i][ptr] = 1;
	}
	if (red_col2 != 0) {
	  GetField(field, entry, j, red_col2);
	  fg_color_r[i][ptr] = atof(field);
	  GetField(field, entry, j, green_col2);
	  fg_color_g[i][ptr] = atof(field);
	  GetField(field, entry, j, blue_col2);
	  fg_color_b[i][ptr] = atof(field);
	}
	else {
	  fg_color_r[i][ptr] = 0;
	  fg_color_g[i][ptr] = 0;
	  fg_color_b[i][ptr] = 0;
	}
	GetField(field, entry, j, codecol);
	text[i][ptr] = (char *)malloc((strlen(field)+1) * sizeof(char));
	strcpy(text[i][ptr], field);
	ptr++;
      }

    text[i][textlen] = NULL;

    i++;
  } while ((read_error = ReadEntry(fp, entry)) == 0);

  if (read_error == 1)
    return 1;

  numseq = i;
  text = (char ***)realloc(text, numseq * sizeof(char **));
                                             /* Free excess space */

  blockpage = (height+1)/(numseq+2);
  if (blockpage == 0)
    fprintf(stderr, "col2psalign: Warning: too many sequences on page\n");
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
    if (option_figure == 0 && block >= blockpage) {
      page++;
      block = 0;
    }
  } while (writecode[startpos] != '\0');
  if (option_figure == 1) {
    if (txtwidth == 0) 
      printf("%%%%BoundingBox: %f %f %f %f\n",
	     0., pageheight-(lineno-1)*tlineheight-2*topmargin,
	     2*sidemargin+maxcolno*charwidth+namewidth, pageheight);
    else
      printf("%%%%BoundingBox: %f %f %f %f\n",
	     0., pageheight-(lineno-1)*tlineheight-2*topmargin,
	     2*sidemargin+txtwidth*charwidth+namewidth, pageheight);
  }
  else {
    if (block != 0) {
      page++;
    }
    printf("%%%%Pages: %d\n", page-1);
  }
  printf("%%!\n");
  printf("/b {%f 0 rmoveto gsave dup stringwidth pop 2 div neg 0 rmoveto 0.00 0.00 0.00 setrgbcolor show grestore %f 0 rmoveto} def\n", charwidth/2, charwidth/2);
  printf("/bcol {%f 0 rmoveto gsave 4 -1 roll dup 5 1 roll stringwidth pop 2 div neg 0 rmoveto setrgbcolor show grestore %f 0 rmoveto} def\n", charwidth/2, charwidth/2);
  printf("/x {gsave 0 %f rmoveto dup 0.00 0.00 0.00 setrgbcolor 0.9 0.9 scale -90 rotate show grestore %f 0 rmoveto} def\n", fontsize*0.7, charwidth);
  printf("/xcol {gsave 4 -1 roll dup 5 1 roll stringwidth pop 2 div neg 0 rmoveto setrgbcolor 0.9 0.9 scale -90 rotate show grestore %f 0 rmoveto} def\n", charwidth);
  printf("/l {gsave 0.00 0.00 0.00 setrgbcolor show grestore} def\n");
  printf("/r {gsave dup stringwidth pop neg 0 rmoveto 0.00 0.00 0.00 setrgbcolor show grestore} def\n");
  printf("/ot {gsave 0 %f rmoveto %f 0 rlineto 0 %f rlineto %f 0 rlineto 0 %f rlineto clip 0.00 0.00 0.00 setrgbcolor show grestore %f 0 rmoveto} def\n", tlineheight*1.5, namewidth-charwidth, -tlineheight*2, -(namewidth-charwidth), tlineheight*.5, namewidth);
  printf("/w {%f (...) stringwidth pop 2 div sub} def\n", (namewidth-charwidth)/2);
  printf("/c {gsave 0 %f rmoveto w 0 rlineto 0 %f rlineto w neg 0 rlineto 0 %f rlineto clip 0.00 0.00 0.00 setrgbcolor dup show grestore w 0 rmoveto (...) show gsave 0 %f rmoveto w 0 rlineto 0 %f rlineto w neg 0 rlineto 0 %f rlineto clip 0.00 0.00 0.00 setrgbcolor dup stringwidth pop w sub neg 0 rmoveto show grestore w 0 rmoveto %f 0 rmoveto} def\n", tlineheight*1.5, -tlineheight*2, tlineheight*.5, tlineheight*1.5, -tlineheight*2, tlineheight*.5, charwidth);
  printf("/t {dup stringwidth pop %f lt {ot} {c} ifelse} def\n", namewidth-charwidth);
  printf("/ret {currentpoint %f sub exch pop %f exch moveto} def\n", tlineheight, sidemargin);
  if (option_fill == 1)
    printf("/col {gsave gsave newpath grestore 0 %f rmoveto %f 0 rlineto 0 %f rlineto -%f 0 rlineto closepath setrgbcolor fill grestore} def\n", -tlineheight/6-(lineheight-1)*fontsize/2, charwidth, tlineheight, charwidth);
  else
    printf("/col {gsave gsave newpath grestore %f %f rmoveto %f 0 rlineto 0 %f rlineto -%f 0 rlineto closepath setrgbcolor %f setlinewidth stroke grestore} def\n", linewidth/2, -tlineheight/6+linewidth/2-(lineheight-1)*fontsize/2,  charwidth-linewidth, tlineheight-linewidth, charwidth-linewidth, linewidth);
  printf("/Helvetica-Bold findfont\n");
  printf("%f scalefont\n", fontsize);
  printf("setfont\n");

  startpos = startptr = 0;
  block = 0;
  page = 1;
  lineno = 0;
  maxcolno = 0;
  do {
    if (block == 0) {
      if (option_figure == 0) {
	printf("%%%%Page: %d %d\n", page, page);
	if (landscape == 1) {
	  printf("%f 0 translate\n", pageheight);
	  printf("90 rotate\n"); 
	}
      }
      printf("%f %f moveto\n\n",
	     sidemargin, pageheight-topmargin-tlineheight);
    }
    printf("() t\n");
    newblock = 1;
    ptr = startptr;
    for (j = startpos; writecode[j] != '\n'; j++) {
      if (option_manynum == 1) {
	if ((writecode[j] == 'c' && newblock == 1) ||
	    (writecode[j] == 'c' && writecode[j-1] == ' ')) {
	  printf("(");
	  PrintString(stdout, sym[ptr]);
	  printf(") l ");
	  newblock = 0; }
	printf("() b ");
	if (writecode[j] == 'c' &&
	    (writecode[j+1] == ' ' ||
	     writecode[j+1] == '\n')) {
	  printf("(");
	  PrintString(stdout, sym[ptr]);
	  printf(") r ");
	  newblock = 0; }
      }
      else {
	if (writecode[j] == 'c' && newblock == 1) {
	  printf("(");
	  PrintString(stdout, sym[ptr]);
	  printf(") l ");
	  newblock = 0; }
	printf("() b ");
	if (writecode[j] == 'c' &&
	    ((writecode[j+1] == ' ' && writecode[j+2] == '.') ||
	     writecode[j+1] == '\n')) {
	  printf("(");
	  PrintString(stdout, sym[ptr]);
	  printf(") r ");
	  newblock = 0; }
      }
      if (writecode[j] == '.')
	newblock = 1;
      if (writecode[j] == 'c')
	ptr++;
    }

    printf("ret\n");
    lineno++;

    printf("gsave\n");
    for (i = 0; i < numseq; i++) {
      printf("() t\n");
      ptr = startptr;
      for (j = startpos; writecode[j] != '\n'; j++) {
	if (writecode[j] == 'c') {
	  if (color_r[i][ptr] != 1 ||
	      color_g[i][ptr] != 1 ||
	      color_b[i][ptr] != 1)
	    printf("%.2f %.2f %.2f col ",
		   color_r[i][ptr], color_g[i][ptr], color_b[i][ptr]);
	  printf("() b ");
	  ptr++;
	}
	else
	  printf("() b ");
      }
      printf("ret\n");
    }
    printf("grestore\n");
    for (i = 0; i < numseq; i++) {
      printf("(");
      PrintString(stdout, name[i]);
      printf(") t\n");
      ptr = startptr;
      colno = 0;
      maxlen = 1;
      for (j = startpos; writecode[j] != '\n'; j++) {
	if (writecode[j] == 'c') {
	  printf("(");
	  for (k = 0; text[i][ptr][k] != '\0'; k++) {
	    if (text[i][ptr][k] == '(' || text[i][ptr][k] == ')')
	      printf("\\%c", text[i][ptr][k]);
	    else
	      printf("%c", text[i][ptr][k]);
	  }
	  printf(") ");
	  colno++;
	}
	else {
	  if (writecode[j] == '(' || writecode[j] == ')')
	    printf("(\\%c) ", writecode[j]);
	  else
	    printf("(%c) ", writecode[j]);
	  colno++;
	}
	if (writecode[j] != 'c')
	    printf("b ");
	else {
	  if (strlen(text[i][ptr]) == 1) {
	    if ((fg_color_r[i][ptr] != 0 ||
		 fg_color_g[i][ptr] != 0 ||
		 fg_color_b[i][ptr] != 0))
	      printf("%.2f %.2f %.2f bcol ",
		     fg_color_r[i][ptr], fg_color_g[i][ptr], fg_color_b[i][ptr]);
	    else
	      printf("b ");
	  }
	  else {
	    if (strlen(text[i][ptr]) > maxlen)
	      maxlen = strlen(text[i][ptr]);
	    if ((fg_color_r[i][ptr] != 0 ||
		 fg_color_g[i][ptr] != 0 ||
		 fg_color_b[i][ptr] != 0))
	      printf("%.2f %.2f %.2f xcol ",
		     fg_color_r[i][ptr], fg_color_g[i][ptr], fg_color_b[i][ptr]);
	    else
	      printf("x ");
	  }
	}
	if (writecode[j] == 'c')
	  ptr++;
      }
      if (colno > maxcolno)
	maxcolno = colno;   /* for BoundingBox */
      if (maxlen == 1)
	printf("ret\n");
      else
	for (k = 0; k < (maxlen/2); k++)
	  printf("ret\n");
      lineno++;
    }
    startpos = j+1;
    startptr = ptr;
    printf("ret\n");
    lineno++;
    block++;
    new_showpage = 1;
    if (option_figure == 0 && block >= blockpage) {
      printf("showpage\n");
      new_showpage = 0;
      page++;
      block = 0;
    }
  } while (writecode[startpos] != '\0');
  if (new_showpage == 1)
    printf("showpage\n");

  return 0;
}

void usage(void)
{
  fprintf(stderr,
	  "usage: col2psalign\n"
	  "            [-ple] [-r<range> | --range=<range>] [-f<fontsize>\n"
	  "            | --fit] [--space] [--width=<width>] [--height=<height>]\n"
	  "            [--margin=<margin>] [--topmargin=<topmargin>]\n"
	  "            [--sidemargin=<sidemargin>] [--namewidth=<namewidth>]\n"
	  "            [--textwidth=<textwidth>] [--lineheight=<lineheight>]\n"
	  "            [--figure] [<file>]\n");
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

void PrintString(FILE *fp, char *s)
{
  int i;

  for (i = 0; s[i] != '\0'; i++)
    if (s[i] == '(')
      printf("\\(");
    else if (s[i] == ')')
      printf("\\)");
    else
      printf("%c", s[i]);
}
