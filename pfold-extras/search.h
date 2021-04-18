/*********************************************************************

  search.h

  Contains general functions for searching.

  000701 Bjarne Knudsen (bk@daimi.au.dk)a

*********************************************************************/

#ifndef __search_h__
#define __search_h__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct tagSuffixNode {
  void *elm;       /* Pointer to element at node */
  int depth;       /* Tree depth at node */
  struct tagSuffixNode **child;     /* Pointer to children */
} SuffixNode;

typedef struct tagSuffixTree {
  SuffixNode *root;
  int alphabet_size;
} SuffixTree;

SuffixTree *InitSuffixSearch(int alphabet_size);
void FinishSuffixSearch(SuffixTree *stree);
void AddSuffixNode(SuffixTree * stree, SuffixNode* snode,
		   int letter, void *elm);
SuffixNode *SuffixSearch(SuffixTree* stree, int *string);
void SuffixInsert(SuffixTree* stree, int *string, void *elm);

#endif
