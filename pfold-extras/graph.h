/*********************************************************************

  graph.h

  Contains general functions for handling graphs.

  000926 Bjarne Knudsen (bk@daimi.au.dk)

*********************************************************************/

#ifndef __graph_h__
#define __graph_h__

#include "llist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct tagGraphNode {
  void *elm;       /* Pointer to element at node */
  LList *neighbour;
  int visit;       /* Used when traversing graphs */
} GraphNode;

typedef struct tagGraph {
  LList *nodes;
} Graph;

Graph *MakeGraph(void);
void FreeGraph(Graph *graph);
GraphNode *AddGNode(Graph *graph, void *elm);
void ConnectGNode(GraphNode *gnode1, GraphNode *gnode2);
void *Element(GraphNode *gnode);
LList *Neighbour(GraphNode *gnode);

void TraverseGraph(Graph *graph, void (*f)(GraphNode *, void *),
			   void *arg);

LList *SplitGraph(Graph *graph);

#endif
