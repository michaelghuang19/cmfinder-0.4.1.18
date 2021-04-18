#include "search.h"
#include "matrix.h"

SuffixTree *InitSuffixSearch(int alphabet_size)
{
  int i;

  SuffixTree *stree;

  stree = (SuffixTree *)malloc(sizeof(SuffixTree));
  stree->root = (SuffixNode *)malloc(sizeof(SuffixNode));

  /* Initialize stree */
  stree->alphabet_size = alphabet_size;
  stree->root->elm = NULL;
  stree->root->depth = 0;

  /* make room for children */
  stree->root->child = (SuffixNode **)malloc(stree->alphabet_size *
					     sizeof(SuffixNode *));

  for (i = 0; i < stree->alphabet_size; i++)
    stree->root->child[i] = NULL;
  
  return stree;
}

static void freesuffix(SuffixTree *stree, SuffixNode *snode);

void FinishSuffixSearch(SuffixTree *stree)
{
  freesuffix(stree, stree->root);
}

static void freesuffix(SuffixTree *stree, SuffixNode *snode)
{
  int i;

  for (i = 0; i < stree->alphabet_size; i++)
    if (snode->child[i] != NULL)
      freesuffix(stree, snode->child[i]);

  free(snode->child);
  FreeMatrix(snode->elm);
  free(snode);
}

void AddSuffixNode(SuffixTree *stree, SuffixNode *snode,
		   int letter, void *elm)
{
  int i;

  snode->child[letter] = (SuffixNode *)malloc(sizeof(SuffixNode));;

  snode->child[letter]->elm = elm;
  snode->child[letter]->depth = snode->depth+1;

  /* make room for children */
  snode->child[letter]->child = (SuffixNode **)malloc(stree->alphabet_size *
					     sizeof(SuffixNode *));

  for (i = 0; i < stree->alphabet_size; i++)
    snode->child[letter]->child[i] = NULL;
}

SuffixNode *SuffixSearch(SuffixTree* stree, int *string)
{
  SuffixNode *snode;
  int letter;

  for (snode = stree->root;
       (letter = string[snode->depth]) != -1 &&
	 snode->child[letter] != NULL;
       snode = snode->child[letter])
    ;

  return snode;
}

void SuffixInsert(SuffixTree* stree, int *string, void *elm)
{
  SuffixNode *snode;
  int letter;

  for (snode = stree->root;
       (letter = string[snode->depth]) != -1;
       snode = snode->child[letter])
     if (snode->child[letter] == NULL)
       AddSuffixNode(stree, snode, string[snode->depth], NULL);

  snode->elm = elm;
}
