#include "phyl.h"

typedef struct tagColno {
  int label;
  int number;
  int uplen;
  int child;
  int brother;
  int name;
} Colno;

typedef struct tagFDdata {
  Matrix **pmatrix;
  int size;
  int** count;
  double *f;
} FDdata;

FDdata *fdd;
static int MAXINT = 100;

/* Allocate and initialize Phyl */
Phyl *MakePhyl(void)
{
  Phyl *phyl;

  phyl = (Phyl *)malloc(sizeof(Phyl));
  phyl->root = (PhylNode *)malloc(sizeof(PhylNode));

  /* Initialize phyl */
  phyl->root->elm = NULL;
  phyl->root->uplen = 0;
  phyl->root->brother = NULL;
  phyl->root->child = NULL;
  phyl->root->parent = NULL;

  return phyl;
}

static void FreeNode(PhylNode *pnode, void *arg);

/* Print phylogeny in Newick format to fp */
void FreePhyl(Phyl *phyl)
{
  PostOrderTraversePhyl(phyl, FreeNode, NULL);

  free(phyl);
}

/* Free a node */
static void FreeNode(PhylNode *pnode, void *arg)
{
  free(pnode);
}

/* Insert node under pnode */
void AddChild(PhylNode *pnode, void *elm, double uplen)
{
  PhylNode *newnode;

  newnode = (PhylNode *)malloc(sizeof(PhylNode));

  /* initialize newnode */
  newnode->elm = elm;
  newnode->uplen = uplen;
  newnode->brother = pnode->child;
  newnode->parent = pnode;

  /* insert it */
  pnode->child = newnode;
}

/* Insert phyl under pnode */
void ConnectPhyl(PhylNode *pnode, Phyl *phyl)
{
  /* initialize phyl */
  phyl->root->brother = pnode->child;
  phyl->root->parent = pnode;

  /* insert it */
  pnode->child = phyl->root;

  /* Remove phyl (but not its contents) */
  free(phyl);
}

Phyl *__RPrec(char *newick, double mul);

/* Reads phylogeny in Newick format */
Phyl *ReadPhyl(char *newick, double mul)
{
  Phyl *phyl;

  /* Initialize */
  __RPrec(NULL, 0.);

  phyl = __RPrec(newick, mul);

  return phyl;
}

/* Recursive function used by ReadPhyl */
Phyl *__RPrec(char *newick, double mul)
{
  static int ptr;  /* Where we are now in the string */
  int len;       /* Name and number lengths */
  Phyl *phyl;    /* The sub phylogeny */

  if (newick == NULL) {
    ptr = 0;
    return NULL; }

  phyl = MakePhyl();

  if (isalpha(newick[ptr])) {     /* we have a leaf */
    /* Read name */
    for (len = 0; newick[ptr+len] != ':' && newick[ptr+len] != ',' &&
	          newick[ptr+len] != ')' && newick[ptr+len] != ';'; len++)
      ;
    phyl->root->elm = (void *)malloc(len+1 * sizeof(char)); 
    strncpy((char *)phyl->root->elm, newick+ptr, len);
    ((char *)phyl->root->elm)[len] = '\0';
    ptr += len;
  }
  else if (newick[ptr] == '(') {
    /* Add sub phylogenies */
    while (newick[ptr] == '(' || newick[ptr] == ',') {
      ptr++;
      ConnectPhyl(phyl->root, __RPrec(newick, mul));
      if (newick[ptr] == '(') {
	fprintf(stderr, "phyl: Expected a ',' in position %d\n", ptr);
	exit(1);
      }
    }
    if (newick[ptr] != ')') {
      fprintf(stderr, "phyl: Expected a ')' in position %d\n", ptr);
      exit(1);
    }
    ptr++;
  }
  else {
    fprintf(stderr, "phyl: Expected a '(' in position %d\n", ptr);
    exit(1);
  }

  /* Read up-length */
  if (newick[ptr] == ';' || newick[ptr] == ',' || newick[ptr] == ')')
    phyl->root->uplen = 0;
  else if (newick[ptr] == ':') {
    ptr++;
    sscanf(newick+ptr, "%lf%n", &phyl->root->uplen, &len);
    phyl->root->uplen *= mul;
    ptr += len;
  }
  else {
    fprintf(stderr, "phyl: Expected ':', ',', ')' or ';' in position %d\n", ptr);
    exit(1);
  }
    
  return phyl;
}

static void RCPrec(Entry *entry, PhylNode *pnode, int number, Colno *colno);

/* Reads phylogeny in Newick format */
Phyl *ReadColPhyl(Entry *entry)
{
  Phyl *phyl;
  Colno *colno;
  char field[MAXCOLW];

  phyl = MakePhyl();

  colno = (Colno *) malloc(sizeof(Colno));

  colno->number = ReadColno(entry, "number");
  colno->uplen = ReadColno(entry, "uplen");
  colno->child = ReadColno(entry, "child");
  colno->brother = ReadColno(entry, "brother");
  colno->name = ReadColno(entry, "name");

  if (colno->number == 0 || colno->uplen == 0 ||
      colno->child == 0 || colno->brother == 0 ||
      colno->name == 0)
    return NULL;

  if (ReadText(entry, "root", MAXCOLW, field) != 0)
    return NULL;

  RCPrec(entry, phyl->root, atoi(field), colno);

  return phyl;
}

/* Recursive function used by ReadColPhyl */
static void RCPrec(Entry *entry, PhylNode *pnode, int number, Colno *colno)
{
  char field[MAXCOLW];

  pnode->number = number;

  GetField(field, entry, number, colno->name);
  if (strcmp(field, ".") != 0) { /* Name present */
    pnode->elm = malloc((strlen(field) + 1) * sizeof(char));
    strcpy((char *)pnode->elm, field);
  }
  else
    pnode->elm = NULL;

  GetField(field, entry, number, colno->uplen);
  pnode->uplen = atof(field);

  /* Add child */

  GetField(field, entry, number, colno->child);

  if (strcmp(field, ".") != 0) {  /* We have a child */
    pnode->child = (PhylNode *)malloc(sizeof(PhylNode));
    pnode->child->parent = pnode;
    RCPrec(entry, pnode->child, atoi(field), colno);
  }
  else
    pnode->child = NULL;

  /* Add brother */

  GetField(field, entry, number, colno->brother);

  if (strcmp(field, ".") != 0) {  /* We have a brother */
    pnode->brother = (PhylNode *)malloc(sizeof(PhylNode));
    pnode->brother->parent = pnode->parent;
    RCPrec(entry, pnode->brother, atoi(field), colno);
  }
  else
    pnode->brother = NULL;
}

static int PCPrec(Entry *entry, PhylNode *pnode, int num, Colno *colno);

/* Make col format entry from phyl structure */
Entry *PhylEntry(Phyl *phyl, char *name)
{
  Entry *entry;
  Colno *colno;

  colno = (Colno *)malloc(sizeof(Colno));

  entry = NewEntry("TREE", name, CountNode(phyl->root));

  colno->label = EnsureCol(entry, "label", "N");
  colno->number = EnsureCol(entry, "number", "    .");
  colno->name = EnsureCol(entry, "name", "           .");
  colno->uplen = EnsureCol(entry, "uplen", "0.0000");
  colno->child = EnsureCol(entry, "child", "    .");
  colno->brother = EnsureCol(entry, "brother", "    .");

  AddEntryText(entry, "; root              1\n");

  PCPrec(entry, phyl->root, 1, colno);

  return entry;
}

/* Recursive function used by PhylEntry */
static int PCPrec(Entry *entry, PhylNode *pnode, int num, Colno *colno)
{
  PhylNode *child;
  int i;
  int newnum;

  ChgField(entry, num, colno->number, "%d", num);

  if (pnode->elm != NULL)
    ChgField(entry, num, colno->name, "%s", (char *)pnode->elm);

  ChgField(entry, num, colno->uplen, "%8.6f", pnode->uplen);

  newnum = num + 1;

  for (child = pnode->child, i = 0; child != NULL;
       child = child->brother, i++) {
    if (i == 0) {
      ChgField(entry, num, colno->child, "%d", newnum);
      num = newnum;
    }
    else {
      ChgField(entry, num, colno->brother, "%d", newnum);
      num = newnum;
    }
    newnum = PCPrec(entry, child, newnum, colno);
  }

  return newnum;
}

static void UPErec(PhylNode *pnode, Entry *entry, Colno *colno);

/* Print phylogeny in Newick format to fp */
void UpdatePhylEntry(Phyl *phyl, Entry *entry)
{
  Colno *colno;
  char s[MAXCOLW];
  int i;

  colno = (Colno *)malloc(sizeof(Colno));

  colno->uplen = ReadColno(entry, "uplen");
  colno->child = ReadColno(entry, "child");
  colno->brother = ReadColno(entry, "brother");
  colno->number = ReadColno(entry, "number");
  colno->name = ReadColno(entry, "name");

  for (i = EntryLength(entry); i < CountNode(phyl->root); i++) {
    AddEntryLine(entry);
    ChgField(entry, i+1, colno->number, "%d", i+1);
    ChgField(entry, i+1, colno->name, ".");
  }

  UPErec(phyl->root, entry, colno);

  sprintf(s, "%d", phyl->root->number);

  ChgEntryText(entry, "root", s);
}

/* Recursive function used by PhylEntry */
static void UPErec(PhylNode *pnode, Entry *entry, Colno *colno)
{
  PhylNode *child;
  int i;
  int newnum;

  ChgField(entry, pnode->number,
	   colno->uplen, "%8.6f", pnode->uplen);


  if (pnode->brother == NULL)
    ChgField(entry, pnode->number, colno->brother, ".");
  else
    ChgField(entry, pnode->number,
	     colno->brother, "%d", pnode->brother->number);

  if (pnode->child == NULL)
    ChgField(entry, pnode->number, colno->child, ".");
  else
    ChgField(entry, pnode->number,
	     colno->child, "%d", pnode->child->number);

  for (child = pnode->child; child != NULL; child = child->brother)
    UPErec(child, entry, colno);
}

void __PPrec(FILE *fp, PhylNode *pnode);

/* Print phylogeny in Newick format to fp */
void PrintPhyl(FILE *fp, Phyl *phyl)
{
  __PPrec(fp, phyl->root);
  fprintf(fp, ";\n");
}

/* Recursive function used by PrintPhyl */
void __PPrec(FILE *fp, PhylNode *pnode)
{
  PhylNode *pnode2;

  if (pnode->child == NULL) { /* We have a leaf */
    fprintf(fp, "%s", (char *)pnode->elm);
  }
  else {
    fprintf(fp, "(");
    if (pnode->elm != NULL) { /* We have an element */
      fprintf(fp, "%s,", (char *)pnode->elm);
    }
    for (pnode2 = pnode->child; pnode2 != NULL; pnode2 = pnode2->brother) {
      __PPrec(fp, pnode2);
      if (pnode2->brother != NULL)
	fprintf(fp, ",");
    }
    fprintf(fp, ")");
  }
  if (pnode->uplen != 0)
    fprintf(fp, ":%g", pnode->uplen);
}

void __FPrec(PhylNode *pnode);

/* this function got a seg fault when compiling with -O3 and also -O1 and -O2, but not with -O0.  It doesn't seem like a super time-critical function (I hope), so I'll just disable optimizations, and leave it.  see the 'bug' directory off of the root.  */
#pragma GCC push_options
#pragma GCC optimize ("O0")

/* Simplify phylogeny */
void FixPhyl(Phyl *phyl)
{
  PhylNode *temp;

  /* Fix everything but binary root */
  __FPrec(phyl->root);

  if (phyl->root != NULL && phyl->root->elm == NULL &&
      phyl->root->child != NULL &&
      phyl->root->child->brother != NULL &&
      phyl->root->child->brother->brother == NULL) {
    /* Root has only two children */
    if (phyl->root->child->brother->child != NULL) {
      /* Last child has children */
      /* Move last child to root, and remove root */
      temp = phyl->root->child;
      free(phyl->root);
      phyl->root = phyl->root->child->brother;
      temp->uplen += phyl->root->uplen;
      temp->brother = phyl->root->child;
      phyl->root->child = temp;
      phyl->root->parent = NULL;
    }
    else if (phyl->root->child->child != NULL) {
      /* First child has children */
      /* Exchange children */
      phyl->root->child->brother->brother = phyl->root->child;
      phyl->root->child = phyl->root->child->brother;
      phyl->root->child->brother->brother = NULL;
      /* Move last child to root, and remove root */
      temp = phyl->root->child;
      free(phyl->root);
      phyl->root = phyl->root->child->brother;
      temp->uplen += phyl->root->uplen;
      temp->brother = phyl->root->child;
      phyl->root->child = temp;
      phyl->root->parent = NULL;
    }
  }

  phyl->root->uplen = 0;
}

#pragma GCC pop_options

/* Recursive function used by FixPhyl */
void __FPrec(PhylNode *pnode)
{
  PhylNode *pnode2;

  if (pnode->child != NULL) { /* Not leaf */
    for (pnode2 = pnode->child; pnode2 != NULL; pnode2 = pnode2->brother)
      __FPrec(pnode2);
    if (pnode->elm == NULL &&
	pnode->child->brother == NULL) { /* Node has only one child */
      /* Remove pnode->child */
      pnode2 = pnode->child;
      pnode->elm = pnode->child->elm;
      pnode->uplen += pnode->child->uplen;
      pnode->child = pnode->child->child;
      free(pnode2);
    }
  }
}

void __CPrec(PhylNode *pnode, PhylNode *newnode);

/* Copy phyl and return the copy */
Phyl *CopyPhyl(Phyl *phyl)
{
  Phyl *newphyl;

  newphyl = MakePhyl();

  __CPrec(phyl->root, newphyl->root);

  return newphyl;
}

PhylNode *Child(PhylNode *pnode)
{
  return pnode->child;
}

PhylNode *Brother(PhylNode *pnode)
{
  return pnode->brother;
}

PhylNode *AddNodeAbove(Phyl *phyl, PhylNode *pnode)
{
  PhylNode *new, *temp;

  new = (PhylNode *)malloc(sizeof(PhylNode));

  if (phyl->root == pnode)
    phyl->root = new;

  if (pnode->parent->child == pnode)
    pnode->parent->child = new;
  else
    for (temp = pnode->parent->child; temp->brother != NULL;
	 temp = temp->brother)
      if (temp->brother == pnode) {
	temp->brother = new;
	break;
      }

  new->parent = pnode->parent;
  new->child = pnode;
  new->brother = pnode->brother;
  new->elm = NULL;
  pnode->parent = new;
  pnode->brother = NULL;

  pnode->uplen /= 2;
  new->uplen = pnode->uplen;

  new->number = CountNode(phyl->root);

  return new;
}


void __CPrec(PhylNode *pnode, PhylNode *newnode)
{
  PhylNode *child, *newchild;

  newnode->elm = pnode->elm; 
  newnode->uplen = pnode->uplen;
  if (pnode->brother == NULL)
    newnode->brother = NULL;
  else {
    newnode->brother = (PhylNode *)malloc(sizeof(PhylNode));
    newnode->brother->parent = newnode->parent;
  }
  if (pnode->child == NULL)
    newnode->child = NULL;
  else {
    newnode->child = (PhylNode *)malloc(sizeof(PhylNode));
    newnode->child->parent = newnode;
  }

  for (child = pnode->child, newchild = newnode->child;
       child != NULL;
       child = child->brother, newchild = newchild->brother)
    __CPrec(child, newchild);
}

void __POrec(PhylNode *pnode, void (*f)(PhylNode *, void *), void *arg);

/* Do a postorder tree traversal */
void PostOrderTraversePhyl(Phyl *phyl, void (*f)(PhylNode *, void *), void *arg)
{
  __POrec(phyl->root, f, arg);
}

/* Recursive algorithm for traversal */
void __POrec(PhylNode *pnode, void (*f)(PhylNode *, void *), void *arg)
{
  PhylNode *child;

  for (child = pnode->child; child != NULL; child = child->brother)
    __POrec(child, f, arg);

  f(pnode, arg);
}

void __PRrec(PhylNode *pnode, void (*f)(PhylNode *, void *), void *arg);

/* Do a preorder tree traversal */
void PreOrderTraversePhyl(Phyl *phyl, void (*f)(PhylNode *, void *), void *arg)
{
  __PRrec(phyl->root, f, arg);
}

/* Recursive algorithm for traversal */
void __PRrec(PhylNode *pnode, void (*f)(PhylNode *, void *), void *arg)
{
  PhylNode *child;

  f(pnode, arg);

  for (child = pnode->child; child != NULL; child = child->brother)
    __PRrec(child, f, arg);
}

static void PPrec(PhylNode *pnode, int (*f_pre)(PhylNode *, void *),
		  void *arg_pre,
		  void (*f_post)(PhylNode *, void *),
		  void *arg_post);

/*
  Do both a preorder and postorder traversal at the same time. If
  return value of f_pre is zero, do not traverse the subtree from that
  point.
*/
void PrePostOrderTraversePhyl(Phyl *phyl,
			      int (*f_pre)(PhylNode *, void *),
			      void *arg_pre,
			      void (*f_post)(PhylNode *, void *),
			      void *arg_post)
{
  PPrec(phyl->root, f_pre, arg_pre, f_post, arg_post);
}

/* Recursive algorithm for pre-post traversal */
static void PPrec(PhylNode *pnode, int (*f_pre)(PhylNode *, void *),
		  void *arg_pre,
		  void (*f_post)(PhylNode *, void *),
		  void *arg_post)
{
  PhylNode *child;

  if (f_pre(pnode, arg_pre) != 0)
    for (child = pnode->child; child != NULL; child = child->brother)
      PPrec(child, f_pre, arg_pre, f_post, arg_post);

  f_post(pnode, arg_post);
}


/* Algorithm for child traversal */
void TraversePhylChild(PhylNode *pnode,
		       void (*f)(PhylNode *, PhylNode *, void *),
		       void *arg)
{
  PhylNode *child;

  for (child = pnode->child; child != NULL; child = child->brother)
    f(child, pnode, arg);
}

static void leafcount(PhylNode *pnode, void *arg);

/* Find the number of leaves from this node and down */
int CountLeaf(PhylNode *pnode)
{
  int num;

  num = 0;

  __PRrec(pnode, leafcount, (void *)&num);

  return num;
}

static void leafcount(PhylNode *pnode, void *arg)
{
  if (pnode->child == NULL)
    (*(int *)arg)++;
}

static void nodecount(PhylNode *pnode, void *arg);

/* Find the number of nodes from this node and down */
int CountNode(PhylNode *pnode)
{
  int num;

  num = 0;

  __PRrec(pnode, nodecount, (void *)&num);

  return num;
}

double PhylLength(PhylNode *pnode)
{
  PhylNode *child;
  double len;

  len = 0;
  for (child = Child(pnode); child != NULL; child = Brother(child))
    len += PhylLength(child);

  len += pnode->uplen;

  return len;
}


static void nodecount(PhylNode *pnode, void *arg)
{
  (*(int *)arg)++;
}

Dist *MakeDist(int numseq)
{
  int i, j;
  Dist *dist;

  dist = (Dist *)malloc(sizeof(Dist));

  dist->numleaf = numseq;
  dist->name = (char **)malloc(numseq*sizeof(char *));
  dist->matrix = (double **)malloc(numseq*sizeof(double *));

  for (i = 0; i < numseq; i++)
    dist->matrix[i] = (double *)malloc(numseq*sizeof(double));

  for (i = 0; i < numseq; i++)
    for (j = 0; j < numseq; j++)
      dist->matrix[i][j] = 0;

  return dist;
}

void FreeDist(Dist *dist)
{
  int i;

  for (i = 0; i < dist->numleaf; i++)
    free(dist->matrix[i]);

  free(dist->matrix);
  free(dist->name);

  free(dist);
}

Dist *CopyDist(Dist *dist)
{
  int totleaf;
  int i, j;
  Dist *newdist;

  newdist = (Dist *)malloc(sizeof(Dist));

  totleaf = dist->numleaf;

  newdist->numleaf = totleaf;
  newdist->name = (char **)malloc(totleaf*sizeof(char *));
  newdist->matrix = (double **)malloc(totleaf*sizeof(double *));

  for (i = 0; i < totleaf; i++) {
    newdist->name[i] = (char *)malloc((strlen(dist->name[i])+1)*sizeof(char));
    strcpy(newdist->name[i], dist->name[i]);
  }

  for (i = 0; i < totleaf; i++)
    newdist->matrix[i] = (double *)malloc(totleaf*sizeof(double));

  for (i = 0; i < totleaf; i++)
    for (j = 0; j < totleaf; j++)
      newdist->matrix[i][j] = dist->matrix[i][j];

  return newdist;
}

void PrintDist(FILE *fp, Dist *dist)
{
  int i, j;

  for (i = 0; i < dist->numleaf; i++) {
    fprintf(fp, "%-10.10s", dist->name[i]);
    for (j = 0; j < dist->numleaf; j++)
      fprintf(fp, " %5.3f", dist->matrix[i][j]);
    fprintf(fp, "\n");
  }
}

typedef struct DM_arg {
  int num;
  Dist *dist;
} DM_arg;

typedef struct DM_info {
  int numleaf;
  int *leaflist;
} DM_info;

static void distcalc(PhylNode *pnode, void *dm_arg);

/* Find distance matrix from tree */
Dist *PhylDist(Phyl *original)
{
  int totleaf;
  int i, j;
  Phyl *phyl;
  DM_arg *dm_arg;
  Dist *dist;

  dm_arg = (DM_arg *)malloc(sizeof(DM_arg));
  dist = (Dist *)malloc(sizeof(Dist));

  phyl = CopyPhyl(original);

  totleaf = CountLeaf(phyl->root);

  dist->numleaf = totleaf;
  dist->name = (char **)malloc(totleaf*sizeof(char *));
  dist->matrix = (double **)malloc(totleaf*sizeof(double *));

  for (i = 0; i < totleaf; i++)
    dist->matrix[i] = (double *)malloc(totleaf*sizeof(double));

  for (i = 0; i < totleaf; i++)
    for (j = 0; j < totleaf; j++)
      dist->matrix[i][j] = 0;

  dm_arg->num = 0;
  dm_arg->dist = dist;

  PostOrderTraversePhyl(phyl, distcalc, (void *)dm_arg);

  return dist;
}

static void distcalc(PhylNode *pnode, void *arg)
{
  int totleaf;
  char **name;
  double **matrix;
  DM_info *p_inf;
  PhylNode *child;
  int inlist1, inlist2;
  int i, j, k;

  name = ((DM_arg *)arg)->dist->name;
  matrix = ((DM_arg *)arg)->dist->matrix;  

  if (pnode->child == NULL)  /* leaf */
    name[((DM_arg *)arg)->num] = pnode->elm; /* set name */

  pnode->elm = (void *)malloc(sizeof(DM_info));

  p_inf = (DM_info *)(pnode->elm);
  totleaf = ((DM_arg *)arg)->dist->numleaf;

  if (pnode->child == NULL) {  /* leaf */
    p_inf->numleaf = 1;
    p_inf->leaflist = (int *)malloc(sizeof(int));
    p_inf->leaflist[0] = (((DM_arg *)arg)->num)++;
  }
  else {
    p_inf->numleaf = 0;
    for (child = Child(pnode); child != NULL; child = Brother(child))
      p_inf->numleaf += ((DM_info *)(child->elm))->numleaf;
    p_inf->leaflist = (int *)malloc(p_inf->numleaf * sizeof(int));
    for (child = Child(pnode), i = 0; child != NULL; child = Brother(child))
      for (j = 0; j < ((DM_info *)(child->elm))->numleaf; j++, i++)
	p_inf->leaflist[i] = ((DM_info *)(child->elm))->leaflist[j];
  }

  for (i = 0; i < totleaf; i++)
    for (j = 0; j < totleaf; j++) {
      inlist1 = 0;
      inlist2 = 0;
      for (k = 0; k < p_inf->numleaf; k++) {
	if (p_inf->leaflist[k] == i)
	  inlist1 = 1;
	if (p_inf->leaflist[k] == j)
	  inlist2 = 1;
      }
      if (inlist1+inlist2 == 1)
	matrix[i][j] += pnode->uplen;
    }
}

/* Create a phylogeny from a distance mwtrix, using the UPGMA method */
Phyl *UPGMA(Dist *original)
{
  Phyl **phyl;
  Phyl *temp;
  double *height;
  double newheight;
  int i, j, min_i, min_j;
  int num;
  double min;
  Dist *dist;

  dist = CopyDist(original);

  num = dist->numleaf;

  phyl = (Phyl **)malloc(num * sizeof(Phyl *));
  height = (double *)malloc(num * sizeof(double));

  for (i = 0; i < num; i++) {
    phyl[i] = MakePhyl();
    phyl[i]->root->elm = (void *)dist->name[i];
    height[i] = 0;
    dist->matrix[i][i] = -1;
  }

  temp = phyl[0];

  for (;;) {
    min = -1;
        
    for (i = 0; i < num; i++)
      for (j = 0; j < num; j++)
	if (dist->matrix[i][j] != -1)
	  if (min == -1 || dist->matrix[i][j] < min) {
	    min = dist->matrix[i][j];
	    min_i = i;
	    min_j = j;
	  }    

    if (min == -1)
      break;
    else {
      newheight = dist->matrix[min_i][min_j]/2;
      for (i = 0; i < num; i++)
	if (dist->matrix[min_i][i] != -1) {
	  dist->matrix[min_i][i] =
	    (dist->matrix[min_i][i]*CountLeaf(phyl[min_i]->root) +
	     dist->matrix[min_j][i]*CountLeaf(phyl[min_j]->root)) /
	    (CountLeaf(phyl[min_i]->root) + CountLeaf(phyl[min_j]->root));
	  dist->matrix[i][min_i] = dist->matrix[min_i][i];
	}
      for (i = 0; i < num; i++) {
	dist->matrix[i][min_j] = -1;
	dist->matrix[min_j][i] = -1;
      }
      temp = MakePhyl();
      phyl[min_i]->root->uplen = newheight-height[min_i];
      phyl[min_j]->root->uplen = newheight-height[min_j];
      ConnectPhyl(temp->root, phyl[min_i]);
      ConnectPhyl(temp->root, phyl[min_j]);
      phyl[min_i] = temp;
      phyl[min_j] = NULL;
      height[min_i] = newheight;
    }
  }

  return temp;  /* the last phylogeny to be made */
}

/* Create phylogenies from a distance mwtrix, using the UPGMA method */
Phyl **UPGMAlimit(Dist *original, double limit)
{
  Phyl **phyl;
  Phyl *temp;
  double *height;
  double newheight;
  int i, j, min_i, min_j;
  int num;
  double min;
  Dist *dist;

  dist = CopyDist(original);

  num = dist->numleaf;

  phyl = (Phyl **)malloc(num * sizeof(Phyl *));
  height = (double *)malloc(num * sizeof(double));

  for (i = 0; i < num; i++) {
    phyl[i] = MakePhyl();
    phyl[i]->root->elm = (void *)dist->name[i];
    height[i] = 0;
    dist->matrix[i][i] = -1;
  }

  for (i = 0; i < num; i++)
    for (j = 0; j < num; j++)
      if (dist->matrix[i][j] > limit)
	dist->matrix[i][j] = -1;

  for (;;) {
    min = -1;
        
    for (i = 0; i < num; i++)
      for (j = 0; j < num; j++)
	if (dist->matrix[i][j] != -1)
	  if (min == -1 || dist->matrix[i][j] < min) {
	    min = dist->matrix[i][j];
	    min_i = i;
	    min_j = j;
	  }    

    if (min == -1)
      break;
    else {
      newheight = dist->matrix[min_i][min_j]/2;
      for (i = 0; i < num; i++) {
	if (dist->matrix[min_j][i] == -1)
	  dist->matrix[min_i][i] = -1;
	if (dist->matrix[min_i][i] != -1) {
	  dist->matrix[min_i][i] =
	    (dist->matrix[min_i][i]*CountLeaf(phyl[min_i]->root) +
	     dist->matrix[min_j][i]*CountLeaf(phyl[min_j]->root)) /
	    (CountLeaf(phyl[min_i]->root) + CountLeaf(phyl[min_j]->root));
	}
	dist->matrix[i][min_i] = dist->matrix[min_i][i];
      }
      for (i = 0; i < num; i++) {
	dist->matrix[i][min_j] = -1;
	dist->matrix[min_j][i] = -1;
      }
      temp = MakePhyl();
      phyl[min_i]->root->uplen = newheight-height[min_i];
      phyl[min_j]->root->uplen = newheight-height[min_j];
      ConnectPhyl(temp->root, phyl[min_i]);
      ConnectPhyl(temp->root, phyl[min_j]);
      phyl[min_i] = temp;
      phyl[min_j] = NULL;
      height[min_i] = newheight;
    }
  }

  for (i = 0, j = 0; i < num; i++)
    if (phyl[i] != NULL)
      phyl[j++] = phyl[i];

  phyl = (Phyl **)realloc(phyl, (j+1)*sizeof(Phyl *));
  
  phyl[j] = NULL;

  return phyl;
}

Phyl *Neighbour(Dist *original)
{
  Phyl **phyl;
  Phyl *temp;
  double *avg;
  int i, j, min_i, min_j;
  int num;
  double min;
  double d;
  Dist *dist, *dist_NJ;
  int found;
  int numclu;

  /* distance matrix */
  dist = CopyDist(original);

  /* distance matrix with Neighbour Joining distances */
  dist_NJ = CopyDist(original);

  num = dist->numleaf;

  avg = (double *)malloc(num * sizeof(double));
  
  for (i = 0; i < num; i++)
    dist->matrix[i][i] = -1;

  /* Make a Phylogeny for each sequence */

  phyl = (Phyl **)malloc(num * sizeof(Phyl *));

  for (i = 0; i < num; i++) {
    phyl[i] = MakePhyl();
    phyl[i]->root->elm = (void *)dist->name[i];
  }

  if (num == 1)
    return phyl[0];

  /* Do the Neighbour Joining */

  for (numclu = num; numclu > 1; numclu--) {
    /* New averages */
    for (i = 0; i < num; i++) {
      avg[i] = 0;
      for (j = 0; j < num; j++)
	if (dist->matrix[i][j] != -1) {
	  avg[i] += dist->matrix[i][j];
	}
      if (numclu > 2)
	avg[i] /= numclu-2;
    }

    /* Calculate Neighbour Joining distances */
    for (i = 0; i < num; i++)
      for (j = 0; j < num; j++)
	dist_NJ->matrix[i][j] = dist->matrix[i][j] - (avg[i] + avg[j]);

    /* Find minimum */

    found = 0;
        
    for (i = 0; i < num; i++)
      for (j = 0; j < num; j++)
	if (dist->matrix[i][j] != -1)
	  if (found == 0 || dist_NJ->matrix[i][j] < min) {
	    min = dist_NJ->matrix[i][j];
	    min_i = i;
	    min_j = j;
	    found = 1;
	  }

    if (found == 0)
      break;
    else {
      d = dist->matrix[min_i][min_j];
      for (i = 0; i < num; i++)
	if (dist->matrix[min_i][i] != -1) {
	  dist->matrix[min_i][i] =
	    0.5*(dist->matrix[min_i][i] + dist->matrix[min_j][i] - d);
	  dist->matrix[i][min_i] = dist->matrix[min_i][i];
	}
      for (i = 0; i < num; i++) {
	dist->matrix[i][min_j] = -1;
	dist->matrix[min_j][i] = -1;
      }
      temp = MakePhyl();
      if (d+avg[min_i]-avg[min_j] < 0)
	phyl[min_i]->root->uplen = 0;
      else
	phyl[min_i]->root->uplen = 0.5*(d+avg[min_i]-avg[min_j]);
      if (d+avg[min_j]-avg[min_i] < 0)
	phyl[min_j]->root->uplen = 0;
      else
	phyl[min_j]->root->uplen = 0.5*(d+avg[min_j]-avg[min_i]);
      ConnectPhyl(temp->root, phyl[min_i]);
      ConnectPhyl(temp->root, phyl[min_j]);
      phyl[min_i] = temp;
      phyl[min_j] = NULL;
    }
  }

  return temp;  /* the last phylogeny to be made */
}

typedef struct PL_arg {
  int num;
  char **leaf;
} PL_arg;

static void PLrec(PhylNode *pnode, void *leaf);

char **PhylLeaves(Phyl *phyl)
{
  PL_arg *pl_arg;

  pl_arg = (PL_arg *)malloc(sizeof(PL_arg));
  
  pl_arg->num = 0;
  pl_arg->leaf = (char **)malloc((CountLeaf(phyl->root)+1)*sizeof(char *));

  PostOrderTraversePhyl(phyl, PLrec, pl_arg);

  pl_arg->leaf[pl_arg->num] = NULL;

  return pl_arg->leaf;
}

static void PLrec(PhylNode *pnode, void *arg)
{
  if (pnode->child == NULL)
    ((PL_arg *)arg)->leaf[((PL_arg *)arg)->num++] = pnode->elm;
}

typedef struct CD_arg {
  int numname;
  char **name;    /* Sequence names to keep */
  PhylNode **lastpointer;
} CD_arg;

static void prune(PhylNode *pnode, void *arg);

Phyl *SubPhyl(Phyl *original, char **leaf)
{
  int i;
  CD_arg *cd_arg;
  Phyl *phyl;

  cd_arg = (CD_arg *)malloc(sizeof(CD_arg));

  phyl = CopyPhyl(original);

  for (i = 0; leaf[i] != NULL; i++)
    ;

  cd_arg->numname = i;
  cd_arg->name = leaf;

  PostOrderTraversePhyl(phyl, prune, (void *)cd_arg);

  FixPhyl(phyl);

  return phyl;
}

static void childprune(PhylNode *pnode, PhylNode *parent, void *arg);

/* sets names and makes room for NodeInfo */
static void prune(PhylNode *pnode, void *arg)
{
  ((CD_arg *)arg)->lastpointer = &pnode->child;
  TraversePhylChild(pnode, childprune, arg);
}

static void childprune(PhylNode *pnode, PhylNode *parent, void *arg)
{
  int i;
  int keep;

  if (pnode->elm != NULL) { /* a leaf */
    keep = 0;
    for (i = 0; i < ((CD_arg *)arg)->numname; i++)
      if (strcmp((char *)pnode->elm, ((CD_arg *)arg)->name[i]) == 0) {
	keep = 1;
	((CD_arg *)arg)->lastpointer = &pnode->brother;
	break;
      }
    if (keep == 0) {
      *((CD_arg *)arg)->lastpointer = pnode->brother;
      free(pnode);
    }
  }
  else {
    if (pnode->child != NULL) {  /* an inner node with children */
      ((CD_arg *)arg)->lastpointer = &pnode->brother;
    }
    else { /* an inner node with no children */
      *((CD_arg *)arg)->lastpointer = pnode->brother;
      free(pnode);
    }
  }
}

Dist *AlignDist(Align *align, Grammar *grammar)
{
  Dist *dist;
  int i, j;
  int numseq;

  numseq = align->numseq;

  dist = MakeDist(numseq);

  for (i = 0; i < numseq; i++)
    dist->name[i] = align->name[i];

  for (i = 0; i < numseq; i++)
    dist->matrix[i][i] = 0;

  for (i = 0; i < numseq; i++)
    for (j = 0; j < i; j++) {
      dist->matrix[i][j] =
	Distance(grammar, align->seq[i], align->seq[j], align->len);
      dist->matrix[j][i] = dist->matrix[i][j];
    }

  return dist;
}

Dist *FastAlignDist(Align *align, Grammar *grammar)
{
  Dist *dist;
  int i, j;
  int numseq;

  numseq = align->numseq;

  dist = MakeDist(numseq);

  for (i = 0; i < numseq; i++)
    dist->name[i] = align->name[i];

  for (i = 0; i < numseq; i++)
    dist->matrix[i][i] = 0;

  InitFastDist(grammar);

  for (i = 0; i < numseq; i++)
    for (j = 0; j < i; j++) {
      dist->matrix[i][j] =
	FastDist(grammar, align->seq[i], align->seq[j], align->len);
      dist->matrix[j][i] = dist->matrix[i][j];
    }

  return dist;
}

Dist *PartDist(Dist *dist, int *seqs)
{
  Dist *temp;
  int i, j, i2, j2;
  int size;

  size = 0;

  for (i = 0; i < dist->numleaf; i++)
    if (seqs[i] == 1)
      size++;

  temp = MakeDist(size);

  for (i = 0, i2 = 0; i < size; i++, i2++) {
    while (seqs[i2] == 0)
      i2++;
    for (j = 0, j2 = 0; j < size; j++, j2++) {
      while (seqs[j2] == 0)
	j2++;
      temp->matrix[i][j] = dist->matrix[i2][j2];
    }
  }

  for (i = 0, i2 = 0; i < size; i++, i2++) {
    while (seqs[i2] == 0)
      i2++;
    temp->name[i] = dist->name[i2];
  }

  return temp;
}


double Distance(Grammar *grammar, int *seq1, int *seq2, int len)
{
  int size;
  int i, j;
  double logdist;
  double f;
  Matrix *pmatrix, *temp;
  int **count;
  int ok, diff;
  Sgrp *sgrp;
  LListCounter *lcount = MakeCounter(grammar->sgrp, FIRST);

  sgrp = Next(lcount);

  size = sgrp->freq->cols;

  count = (int **)malloc(size * sizeof(int *));

  for (i = 0 ; i < size; i++)
    count[i] = (int *)malloc(size * sizeof(int));

  for (i = 0 ; i < size; i++)
    for (j = 0 ; j < size; j++)
      count[i][j] = 0;

  ok = 0;
  diff = 0;

  for (i = 0; i < len; i++)
    if (seq1[i] >= 0 && seq1[i] < size &&
	seq2[i] >= 0 && seq2[i] < size) {
      if (seq1[i] != seq2[i])
	diff = 1;
      ok = 1;
      count[seq1[i]][seq2[i]]++;
    }

  if (ok == 0) {
    fprintf(stderr, "Warning: sequences do not overlap\n");
    return 100;
  }

  if (diff == 0)
    return 0;

  logdist = 0;

  InitMinimize(-1., 1., 0.001);

  while (Minimize(&logdist, &f) == 0) {
    if (logdist > 10) break;
    pmatrix = 
      TransposeMatrix(temp = ExpMatrix(exp(logdist), sgrp->eigen,
				sgrp->diag, sgrp->inveigen));

    f = 0;
    for (i = 0; i < size; i++)
      for (j = 0; j < size; j++)
	if (count[i][j] != 0)
	  f -= count[i][j] * Edbl2Dbl(LogEdouble(pmatrix->entry[i][j]));
    FreeMatrix(temp);
    FreeMatrix(pmatrix);
  }

  for (i = 0 ; i < size; i++)
    free(count[i]);

  free(count);

  if (exp(logdist) > 5)
    return 5;

  return exp(logdist);
}

static double Minx(double y0, double y1, double y2);
static double Int2Dbl(double i);

double InitFastDist(Grammar *grammar)
{
  int i;
  Matrix *temp;
  Sgrp *sgrp;
  LListCounter *lcount = MakeCounter(grammar->sgrp, FIRST);

  sgrp = Next(lcount);

  fdd = (FDdata *)malloc(sizeof(FDdata));

  fdd->size = sgrp->freq->cols;

  fdd->pmatrix = (Matrix **)malloc(MAXINT * sizeof(Matrix *));

  for (i = 0; i < MAXINT; i ++) {
    fdd->pmatrix[i] = 
      TransposeMatrix(temp = ExpMatrix(Int2Dbl(i), sgrp->eigen,
				       sgrp->diag, sgrp->inveigen));
    FreeMatrix(temp);
  }

  fdd->count = (int **)malloc(fdd->size * sizeof(int *));

  fdd->f = (double *)malloc(MAXINT * sizeof(double));
}

double FastDist(Grammar *grammar, int *seq1, int *seq2, int len)
{
  int i, j, k, mink, startk, endk;
  double *f, minf;
  int ok, diff;
  int stp;
  int first;

  for (i = 0 ; i < fdd->size; i++)
    fdd->count[i] = (int *)malloc(fdd->size * sizeof(int));

  for (i = 0 ; i < fdd->size; i++)
    for (j = 0 ; j < fdd->size; j++)
      fdd->count[i][j] = 0;

  ok = 0;
  diff = 0;

  for (i = 0; i < len; i++)
    if (seq1[i] >= 0 && seq1[i] < fdd->size &&
	seq2[i] >= 0 && seq2[i] < fdd->size) {
      if (seq1[i] != seq2[i])
	diff = 1;
      ok = 1;
      fdd->count[seq1[i]][seq2[i]]++;
    }

  if (ok == 0) {
    fprintf(stderr, "Warning: sequences do not overlap\n");
    return 5.;
  }

  if (diff == 0)
    return 0;

  stp = 10;
  first = 1;

  for (k = 0; k < MAXINT; k += stp) {
    fdd->f[k] = 0;
    for (i = 0; i < fdd->size; i++)
      for (j = 0; j < fdd->size; j++)
	if (fdd->count[i][j] != 0)
	  fdd->f[k] -= fdd->count[i][j] * Edbl2Dbl(LogEdouble(fdd->pmatrix[k]->entry[i][j]));
    if (fdd->f[k] < minf || k == 0) {
      minf = fdd->f[k];
      mink = k;
    }
    else
      break;
  }

  startk = mink-stp+1;
  if (startk < 0) startk = 0;
  endk = mink+stp;
  if (endk > MAXINT) endk = MAXINT;
  
  stp = 3;

  for (k = startk; k < endk; k += stp) {
    fdd->f[k] = 0;
    for (i = 0; i < fdd->size; i++)
      for (j = 0; j < fdd->size; j++)
	if (fdd->count[i][j] != 0)
	  fdd->f[k] -= fdd->count[i][j] * Edbl2Dbl(LogEdouble(fdd->pmatrix[k]->entry[i][j]));
    if (fdd->f[k] < minf || k == startk) {
      minf = fdd->f[k];
      mink = k;
    }
    else
      break;
  }

  startk = mink-stp+1;
  if (startk < 0) startk = 0;
  endk = mink+stp;
  if (endk > MAXINT) endk = MAXINT;
  
  stp = 1;

  for (k = startk; k < endk; k += stp) {
    fdd->f[k] = 0;
    for (i = 0; i < fdd->size; i++)
      for (j = 0; j < fdd->size; j++)
	if (fdd->count[i][j] != 0)
	  fdd->f[k] -= fdd->count[i][j] * Edbl2Dbl(LogEdouble(fdd->pmatrix[k]->entry[i][j]));
    if (fdd->f[k] < minf || k == startk) {
      minf = fdd->f[k];
      mink = k;
    }
    else
      break;
  }

  if (mink == 0 || mink == MAXINT -1)
    return Int2Dbl(mink);
  else
    return Int2Dbl(mink+Minx(fdd->f[mink-1], fdd->f[mink], fdd->f[mink+1]));
}

static double Int2Dbl(double i)
{
  return exp((10.*(i - MAXINT + 1))/MAXINT)*5;
}

static double Minx(double y0, double y1, double y2)
{
  if (2*y0-4*y1+2*y2 == 0)
    return y1;

  return (y0-y2)/(2*y0-4*y1+2*y2);
}

