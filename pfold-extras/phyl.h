/*********************************************************************

  phyl.h

  Contains general functions for handling the phylogenies.

  000701 Bjarne Knudsen (bk@daimi.au.dk)

*********************************************************************/

#ifndef __phyl_h__
#define __phyl_h__

#include "matrix.h"
#include "optimize.h"
#include "llist.h"
#include "file.h"
#include "align.h"
#include "grammar.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct tagPhylNode {
  void *elm;       /* Pointer to element at node */
  double uplen;    /* Length to parent */
  struct tagPhylNode *brother;   /* Pointer to brother */
  struct tagPhylNode *child;     /* Pointer to first child */
  struct tagPhylNode *parent;    /* Pointer to parent */
  int number;      /* Number of node */
} PhylNode;

typedef struct tagPhyl {
  PhylNode *root;
} Phyl;

typedef struct Dist {
  int numleaf;
  char **name;
  double **matrix;
} Dist;

Phyl *MakePhyl(void);
void FreePhyl(Phyl *phyl);
void AddChild(PhylNode *pnode, void *elm, double uplen);
void ConnectPhyl(PhylNode *pnode, Phyl *phyl);

Phyl *ReadColPhyl(Entry *entry);
Phyl *ReadPhyl(char *newick, double mul);
Entry *PhylEntry(Phyl *phyl, char *name);
void UpdatePhylEntry(Phyl *phyl, Entry *entry);
void PrintPhyl(FILE *fp, Phyl *phyl);
void FixPhyl(Phyl *phyl);
Phyl *CopyPhyl(Phyl *phyl);

PhylNode *Child(PhylNode *pnode);
PhylNode *Brother(PhylNode *pnode);

PhylNode *AddNodeAbove(Phyl *phyl, PhylNode *pnode);

void PostOrderTraversePhyl(Phyl *phyl, void (*f)(PhylNode *, void *),
			   void *arg);
void PreOrderTraversePhyl(Phyl *phyl, void (*f)(PhylNode *, void *),
			  void *arg);
void PrePostOrderTraversePhyl(Phyl *phyl,
			      int (*f_pre)(PhylNode *, void *),
			      void *arg_pre,
			      void (*f_post)(PhylNode *, void *),
			      void *arg_post);
void TraversePhylChild(PhylNode *pnode,
		       void (*f)(PhylNode *, PhylNode *, void *),
		       void *arg);

int CountLeaf(PhylNode *pnode);
int CountNode(PhylNode *pnode);
double PhylLength(PhylNode *pnode);

Dist *MakeDist(int numseq);
void FreeDist(Dist *dist);
Dist *CopyDist(Dist *dist);
void PrintDist(FILE *fp, Dist *dist);
Dist *PhylDist(Phyl *phyl);

Phyl *UPGMA(Dist *dist);
Phyl **UPGMAlimit(Dist *dist, double limit);

Phyl *Neighbour(Dist *dist);

char **PhylLeaves(Phyl *phyl);
Phyl *SubPhyl(Phyl *phyl, char **leaf);

Dist *AlignDist(Align *align, Grammar *grammar);
Dist *FastAlignDist(Align *align, Grammar *grammar);
Dist *PartDist(Dist *dist, int *seqs);
double Distance(Grammar *grammar, int *seq1, int *seq2, int len);

double InitFastDist(Grammar *grammar);
double FastDist(Grammar *grammar, int *seq1, int *seq2, int len);

#endif
