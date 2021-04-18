#include "file.h"

void usage(void);

int main(int argc, char **argv)
{
  FILE *fp, *fp2;         /* For input */
  CmdArg *cmdarg;   /* Command line arguments */
  char *s, *s2;          /* String for arguments */
  unsigned int len;
  int i, j;
  int size, size2;
  double scale;
  double prob, prob2;
  double width;
  int option_figure;
  int option_gray;
  int option_line;
  int grid;
  int mode;

  cmdarg = InitArgument(argc, argv);

  option_gray = 0;
  option_line = 0;
  option_figure = 0;
  width = 545.;
  grid = 10;

  while ((s = GetArgument(cmdarg)) != NULL)
    if (strncmp(s, "-width=", 7) == 0) {
      if (sscanf(&s[7], "%lf%n",
		 &width, &len) != 1 ||
	  len+7 != strlen(s)) {
	usage();
	return 1; }
      option_figure = 1;
    }
    else if (strncmp(s, "-grid=", 6) == 0) {
      if (sscanf(&s[6], "%d%n",
		 &grid, &len) != 1 ||
	  len+6 != strlen(s)) {
	usage();
	return 1; }
    }
    else if (strcmp(s, "-gray") == 0)
      option_gray = 1;
    else if (strcmp(s, "-line") == 0)
      option_line = 1;
    else {
      usage();
      return 1;
    }

  if (option_gray == 1 && option_line == 1) {
    usage();
    return 1;
  }

  mode = 1;

  if ((s = GetFilename(cmdarg)) == NULL)
    fp = stdin;
  else if ((fp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "drawdot: Error in opening file '%s'\n", s);
    return 1; }

  if ((s2 = GetFilename(cmdarg)) != NULL) {
    mode = 2;
    if (GetFilename(cmdarg) != NULL) {
      usage();
      return 1; }
    else if ((fp2 = fopen(s2, "r")) == NULL) {
      fprintf(stderr, "drawdot: Error in opening file '%s'\n", s);
      return 1; }
    fprintf(stderr, "Orange: %s\n", s);
    fprintf(stderr, "Blue:   %s\n", s2);
  }

  if (fscanf(fp, " %d", &size) != 1) {
    fprintf(stderr, "drawdot: Error in reading size\n");
    return 1;
  }

  if (mode == 2 && (fscanf(fp2, " %d", &size2) != 1 || size != size2)) {
    fprintf(stderr, "drawdot: Error in reading size\n");
    return 1;
  }

  scale = (width-10)/(size+2.1);

  printf("%%!PS-Adobe-2.0\n");
  if (option_figure == 1)
    printf("%%%%BoundingBox: %f %f %f %f\n", -2.05*scale-5, -2.05*scale-5, (size+0.05)*scale+5, (size+0.05)*scale+5);
  printf("%%!\n");
  printf("\n");
  if (option_figure != 1)
    printf("%f %f translate\n", (595-width)/2+2.05*scale+5, (842-width)/2+2.05*scale+5);
  printf("%f %f scale\n", scale, scale);
  printf("\n");
  if (option_gray == 1 && mode == 1) {
    printf("/b {\n");
    printf("  newpath\n");
    printf("  1 exch sub setgray\n");
    printf("  moveto\n");
    printf("  0.05 0.05 rmoveto\n");
    printf("  0 0.9 rlineto\n");
    printf("  0.9 0 rlineto\n");
    printf("  0 -0.9 rlineto\n");
    printf("  fill\n");
    printf("  stroke\n");
    printf("} def\n");
    printf("\n");
    printf("/s {b} def\n");
    printf("/t {b} def\n");
  }
  else if (option_line == 1) {
    printf("/b {\n");
    printf("  newpath\n");
    printf("  3 1 roll\n");
    printf("  moveto\n");
    printf("  1 exch sub sqrt 0.9 mul dup dup dup dup dup\n");
    printf("  0.05 add 0.05 rmoveto\n");
    printf("  0.9 exch sub 0 rlineto\n");
    printf("  0.9 exch sub 0 exch rlineto\n");
    printf("  dup neg exch rlineto\n");
    printf("  0.9 sub 0 rlineto\n");
    printf("  0.9 sub 0 exch rlineto\n");
    printf("  fill\n");
    printf("  stroke\n");
    printf("} def\n");
    printf("\n");
    printf("/s {\n");
    printf("  newpath\n");
    printf("  3 1 roll\n");
    printf("  moveto\n");
    printf("  0.9 mul dup dup\n");
    printf("  2 div 0.5 exch sub 0.05 exch rmoveto\n");
    printf("  0 exch rlineto\n");
    printf("  0.9 0 rlineto\n");
    printf("  0 exch neg rlineto\n");
    printf("  fill\n");
    printf("  stroke\n");
    printf("} def\n");
    printf("\n");
    printf("/t {\n");
    printf("  newpath\n");
    printf("  3 1 roll\n");
    printf("  moveto\n");
    printf("  0.9 mul dup dup\n");
    printf("  2 div 0.5 exch sub 0.05 rmoveto\n");
    printf("  0 rlineto\n");
    printf("  0 0.9 rlineto\n");
    printf("  neg 0 rlineto\n");
    printf("  fill\n");
    printf("  stroke\n");
    printf("} def\n");
  }
  else {
    printf("/b {\n");
    printf("  newpath\n");
    printf("  3 1 roll\n");
    printf("  moveto\n");
    printf("  sqrt 0.9 mul dup dup dup\n");
    printf("  2 div 0.5 exch sub dup rmoveto\n");
    printf("  0 exch rlineto\n");
    printf("  0 rlineto\n");
    printf("  0 exch neg rlineto\n");
    printf("  fill\n");
    printf("  stroke\n");
    printf("} def\n");
    printf("\n");
    printf("/s {b} def\n");
    printf("/t {b} def\n");
  }

  printf("\n");
  printf("0.1 setlinewidth\n");
  printf("\n");
  printf("newpath\n");
  printf("0 0 moveto\n");
  printf("%d 0 rlineto\n", size);
  printf("0 %d rlineto\n", size);
  printf("-%d 0 rlineto\n", size);
  printf("closepath\n");
  printf("stroke\n");
  printf("\n");
  printf("newpath\n");
  printf("0 0 moveto\n");
  printf("%d %d rlineto\n", size, size);

  if (grid != 0)
    for (i = grid; i < size; i += grid) {
      printf("0 %d moveto\n", i);
      printf("%d 0 rlineto\n", size);
      printf("%d 0 moveto\n", i);
      printf("0 %d rlineto\n", size);
    }

  printf("stroke\n");
  printf("\n");

  if (mode == 1)
    for (i = 0; i < size; i++) {
      for (j = 0; j < size; j++) {
	if (fscanf(fp, " %lf", &prob) != 1) {
	  fprintf(stderr, "drawdot: Error in reading prob\n");
	  return 1;
	}
	if (prob > 0.002)
	  printf("%d %d %8.6f b\n", i, j, prob);
      }
    }
  else
    for (i = 0; i < size; i++) {
      for (j = 0; j < size; j++) {
	if (fscanf(fp, " %lf", &prob) != 1) {
	  fprintf(stderr, "drawdot: Error in reading prob\n");
	  return 1;
	}
	if (fscanf(fp2, " %lf", &prob2) != 1) {
	  fprintf(stderr, "drawdot: Error in reading prob\n");
	  return 1;
	}
	if (option_gray == 1) {
	  if (prob > 0.02 | prob2 > 0.02) {
	    if (prob > prob2)
	      printf("%f %f %f setrgbcolor\n",
		     1-0.2*prob2, 1-0.5*prob-0.3*prob2, 1-prob+0.8*prob2);
	    else
	      printf("%f %f %f setrgbcolor\n",
		     1-prob2+0.8*prob, 1-0.5*prob2-0.3*prob, 1-0.2*prob);
	    printf("%d %d 1 b\n", i, j);
	  }
	}
	else {
	  if (prob > prob2) {
	    if (prob > 0.02) {
	      printf("1 0.5 0 setrgbcolor\n");
	      printf("%d %d %8.6f b\n", i, j, prob);
	    }
	    if (prob2 > 0.02) {
	      printf("0.8 0.2 0.8 setrgbcolor\n");
	      printf("%d %d %8.6f b\n", i, j, prob2);
	    }
	  }
	  else {
	    if (prob2 > 0.02) {
	      printf("0 0.5 1 setrgbcolor\n");
	      printf("%d %d %8.6f b\n", i, j, prob2);
	    }
	    if (prob > 0.02) {
	      printf("0.8 0.2 0.8 setrgbcolor\n");
	      printf("%d %d %8.6f b\n", i, j, prob);
	    }
	  }
	}
      }
    }

  printf("\n");

  if (mode == 1)
    for (i = 0; i < size; i++) {
      if (fscanf(fp, " %lf", &prob) != 1) {
	fprintf(stderr, "drawdot: Error in reading prob\n");
	return 1;
      }
      if (prob > 0.02) {
	printf("%d -2 %8.6f s\n", i, prob);
	printf("-2 %d %8.6f t\n", i, prob);
      }
    }
  else
    for (i = 0; i < size; i++) {
      if (fscanf(fp, " %lf", &prob) != 1) {
	fprintf(stderr, "drawdot: Error in reading prob\n");
	return 1;
      }
      if (fscanf(fp2, " %lf", &prob2) != 1) {
	fprintf(stderr, "drawdot: Error in reading prob\n");
	return 1;
      }
      if (option_gray == 1) {
	if (prob > 0.02 | prob2 > 0.02) {
	  if (prob > prob2)
	    printf("%f %f %f setrgbcolor\n",
		   1-0.2*prob2, 1-0.5*prob-0.3*prob2, 1-prob+0.8*prob2);
	  else
	    printf("%f %f %f setrgbcolor\n",
		   1-prob2+0.8*prob, 1-0.5*prob2-0.3*prob, 1-0.2*prob);
	  printf("%d -2 1 s\n", i, prob);
	  printf("-2 %d 1 t\n", i, prob);
	}
      }
      else {
	if (prob > prob2) {
	  if (prob > 0.02) {
	    printf("1 0.5 0 setrgbcolor\n");
	    printf("%d -2 %8.6f s\n", i, prob);
	    printf("-2 %d %8.6f t\n", i, prob);
	  }
	  if (prob2 > 0.02) {
	    printf("0.8 0.2 0.8 setrgbcolor\n");
	    printf("%d -2 %8.6f s\n", i, prob2);
	    printf("-2 %d %8.6f t\n", i, prob2);
	  }
	}
	else {
	  if (prob2 > 0.02) {
	    printf("0 0.5 1 setrgbcolor\n");
	    printf("%d -2 %8.6f s\n", i, prob2);
	    printf("-2 %d %8.6f t\n", i, prob2);
	  }
	  if (prob > 0.02) {
	    printf("0.8 0.2 0.8 setrgbcolor\n");
	    printf("%d -2 %8.6f s\n", i, prob);
	    printf("-2 %d %8.6f t\n", i, prob);
	  }
	}
      }
    }

  printf("\nshowpage\n");

}

void usage()
{
  fprintf(stderr,
	  "usage: drawdot [--gray | --line] [--width=<width>] [<file>]\n");
}
