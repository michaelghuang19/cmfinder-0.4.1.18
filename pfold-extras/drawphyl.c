#include "file.h"
#include "phyl.h"

typedef struct tagColno {
  int color_r;
  int color_g;
  int color_b;
  int width;
  int dash;
  int bootstrap;
} Colno;

double factor;
int option_nodenum;

void usage(void);
double drawphyl(FILE *fp, Entry *entry, PhylNode *pnode,
		double xpos, double *ypos, Colno *colno);
void PrintString(FILE *fp, char *s);

int main(int argc, char **argv)
{
  FILE *fp;
  Header *header;
  Entry *entry;
  int read_error;   /* For keeping track of errors in reading entries */
  Phyl *phyl;
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */
  unsigned int len;
  double *pos;
  Colno *colno;
  int option_figure;

  colno = (Colno *) malloc(sizeof(Colno));

  pos = (double *)malloc(sizeof(double));

  cmdarg = InitArgument(argc, argv);

  factor = 1;
  option_nodenum = 0;
  option_figure = 0;

  while ((s = GetArgument(cmdarg)) != NULL)
    if (strncmp(s, "-mult=", 6) == 0) {
      if (sscanf(&s[6], "%lf%n", &factor, &len) != 1 ||
	  len+6 != strlen(s)) {
	usage();
	return 1; }
    }
    else if (strcmp(s, "-figure") == 0)
      option_figure = 1;
    else if (strcmp(s, "-nodenum") == 0)
      option_nodenum = 1;
    else {
      usage();
      return 1; }

  if ((s = GetFilename(cmdarg)) == NULL)
    fp = stdin;
  else if (GetFilename(cmdarg) != NULL) {
    usage();
    return 1; }
  else if ((fp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "drawphyl: Error in opening file '%s'\n", s);
    return 1; }

  factor *= 72;

  header = MakeHeader();
  entry = MakeEntry();

  if (ReadHeader(fp, header) != 0)
    return 1;

  phyl = NULL;

  while ((read_error = ReadEntry(fp, entry)) == 0) {
    if (!ReadType(entry, "TREE"))
      continue;

    phyl = ReadColPhyl(entry);

    break;
  }

  if (phyl == NULL) {
    fprintf(stderr, "drawphyl: No tree read\n");
    return 1;
  }

  colno->color_r = ReadColno(entry, "color_r");
  colno->color_g = ReadColno(entry, "color_g");
  colno->color_b = ReadColno(entry, "color_b");
  colno->width = ReadColno(entry, "width");
  colno->dash = ReadColno(entry, "dash");
  colno->bootstrap = ReadColno(entry, "bootstrap");

  if (colno->bootstrap == 0)
    colno->bootstrap = ReadColno(entry, "probability");

  printf("%%!PS-Adobe-2.0\n");
  if (option_figure == 1)
    printf("%%%%BoundingBox: 15 %f 500 710\n", 670.-CountLeaf(phyl->root)*10.);
  printf("\n");
  printf("/Helvetica findfont\n");
  printf("10 scalefont\n");
  printf("setfont\n");
  printf("\n");
  printf("/halfheight\n");
  printf("  newpath 0 0 moveto (M) true charpath flattenpath pathbbox\n");
  printf("  exch pop exch sub exch pop 2 div newpath\n");
  printf("def\n");
  printf("\n");
  printf("/rshow {\n");
  printf("  0 halfheight neg rmoveto show\n");
  printf("} def\n");
  printf("\n");
  printf("/lshow {\n");
  printf("  dup stringwidth pop neg halfheight neg rmoveto show\n");
  printf("} def\n");
  printf("\n");
  printf("2 setlinecap\n");
  printf("\n");
  printf("gsave\n");

  *pos = 700;

  drawphyl(stdout, entry, phyl->root, 20, pos, colno);

  printf("grestore\n");
  printf("\n");
  printf("30 %f moveto\n", *pos-20);
  printf("0 -8 rlineto\n");
  printf("0 4 rmoveto\n");
  printf("144 0 rlineto\n");
  printf("0 4 rmoveto\n");
  printf("0 -8 rlineto\n");
  printf("-72 -4.5 rmoveto\n");
  printf("(%.2f units) dup stringwidth pop neg 2 div 0 rmoveto rshow\n", 144/factor);
  printf("stroke\n");
  printf("\n");

  printf("showpage\n");

  return 0;
}

void usage(void)
{
  fprintf(stderr,
	  "usage: drawphyl [--tree=<tree> | <file>]\n");
}

void setcolor(FILE *fp, Entry *entry, int number, Colno *colno);

/*
   drawphyl draws the subphylogeny below pnode. upon calling, *ypos is
   the y-position on the page to start from. When returning *ypos is
   the ending y-position. xpos is the x position to start from. The
   return value is the y-position of the branch to the subtree just
   drawn.
*/
double drawphyl(FILE *fp, Entry *entry, PhylNode *pnode,
		double xpos, double *ypos, Colno *colno)
{
  PhylNode *child;
  double startypos, endypos;
  int first;
  char field[MAXCOLW];
  double offset;

  if (pnode->child == NULL) { /* Leaf */
    fprintf(fp, "%f %f moveto\n", xpos, *ypos);
    fprintf(fp, "%f 0 rlineto\n", pnode->uplen*factor);
    fprintf(fp, "5 0 rmoveto\n");
    setcolor(fp, entry, pnode->number, colno);
    if (pnode->elm != NULL) {
      fprintf(fp, "(");
      PrintString(fp, (char *)pnode->elm);
      fprintf(fp, ") rshow\n");
    }
    fprintf(fp, "stroke\n\n");
    return *ypos;
  }

  /* not leaf */

  first = 1;
  startypos = *ypos;

  for (child = pnode->child; child != NULL; child = child->brother) {
    if (first == 1) {
      startypos = drawphyl(fp, entry, child,
			   xpos + pnode->uplen*factor, ypos, colno);
      endypos = startypos;
      first = 0;
    }
    else {
      *ypos -= 10;
      endypos = drawphyl(fp, entry, child,
			 xpos + pnode->uplen*factor, ypos, colno);
    }
  }

  child = pnode->child;

  fprintf(fp, "%f %f moveto\n", xpos + pnode->uplen*factor, startypos);
  fprintf(fp, "0 %f rlineto\n", (endypos - startypos)/2);
  setcolor(fp, entry, child->number, colno);
  fprintf(fp, "currentpoint\n");
  fprintf(fp, "stroke\n\n");

  /* Find last child */
  while (child->brother != NULL)
    child = child->brother;

  fprintf(fp, "moveto\n");
  fprintf(fp, "0 %f rlineto\n", (endypos - startypos)/2);
  setcolor(fp, entry, child->number, colno);
  fprintf(fp, "stroke\n\n");
  
  fprintf(fp, "%f %f moveto\n", xpos, (startypos + endypos)/2);
  fprintf(fp, "%f 0 rlineto\n", pnode->uplen*factor);
  setcolor(fp, entry, pnode->number, colno);
  fprintf(fp, "stroke\n\n");

  if (option_nodenum == 1) {
    fprintf(fp, "%f %f moveto\n", xpos + pnode->uplen*factor + 5,
	    (startypos + endypos)/2);
    setcolor(fp, entry, pnode->number, colno);
    fprintf(fp, "(%d) rshow\n", pnode->number);
    fprintf(fp, "stroke\n\n");
  }

  if (colno->bootstrap != 0) {
    GetField(field, entry, pnode->number, colno->bootstrap);
    if (strcmp(field, ".") != 0) {
      if (pnode->brother == NULL)
	offset = -5;
      else
	offset = 5;
      fprintf(fp, "%f %f moveto\n", xpos + pnode->uplen*factor-2.5,
	      (startypos + endypos)/2+offset);
      fprintf(fp, "0 0 0 setrgbcolor\n");
      fprintf(fp, "(%s) lshow\n", field);
    }
  }
    
  return (startypos + endypos)/2;
}

void setcolor(FILE *fp, Entry *entry, int number, Colno *colno)
{
  char field[MAXCOLW];

  if (colno->color_r != 0) {
    GetField(field, entry, number, colno->color_r);
    fprintf(fp, "%s ", field);
    GetField(field, entry, number, colno->color_g);
    fprintf(fp, "%s ", field);
    GetField(field, entry, number, colno->color_b);
    fprintf(fp, "%s ", field);
    fprintf(fp, "setrgbcolor\n");
  }

  if (colno->width != 0) {
    GetField(field, entry, number, colno->width);
    fprintf(fp, "%s setlinewidth\n", field);
  }

  if (colno->dash != 0) {
    GetField(field, entry, number, colno->dash);
    if (atoi(field) == 0)
      fprintf(fp, "[] 0 setdash\n");
    else if (atoi(field) == 1)
      fprintf(fp, "[1 3] 0.5 setdash\n");
    else if (atoi(field) == 2)
      fprintf(fp, "[3 5] 1.5 setdash\n");
    else
      fprintf(fp, "[3 3] 1.5 setdash\n");
  }
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
