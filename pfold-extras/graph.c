/*********************************************************************

  dna.c

  Contains general functions for handling graphs.

  001016 Bjarne Knudsen (bk@daimi.au.dk)

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

#include "graph.h"

/* Allocate and initialize Graph */
Graph *MakeGraph(void)
{
  Graph *graph;

  graph = (Graph *)malloc(sizeof(Graph));

  /* Initialize graph to an empty graph */
  graph->nodes = MakeLList();

  return graph;
}

void FreeGraph(Graph *graph)
{
  LListCounter *lcount;
  GraphNode *gnode;

  /* No nodes have been visited yet */
  lcount = MakeCounter(graph->nodes, FIRST);
  while ((gnode = Next(lcount)) != NULL) {
    DestroyLList(gnode->neighbour);
    free(gnode);
  }

  DestroyLList(graph->nodes);

  free(graph);

  free(lcount);
}

GraphNode *AddGNode(Graph *graph, void *elm)
{
  GraphNode *newnode;

  newnode = (GraphNode *)malloc(sizeof(GraphNode));

  /* initialize newnode */
  newnode->elm = elm;
  newnode->neighbour = MakeLList();

  /* insert it */
  Enqueue(graph->nodes, newnode);

  return newnode;
}

void ConnectGNode(GraphNode *gnode1, GraphNode *gnode2)
{
  /* Make nodes each others neighbours */

  Enqueue(gnode1->neighbour, gnode2);
  Enqueue(gnode2->neighbour, gnode1);
}

void *Element(GraphNode *gnode)
{
  return gnode->elm;
}

LList *Neighbour(GraphNode *gnode)
{
  return gnode->neighbour;
}

static void __TGrec(GraphNode *gnode,
		   void (*f)(GraphNode *, void *), void *arg);

/*
   Do a graph traversal. The nodes visited will be neighbours of
   previously visited nodes, if possible.
*/
void TraverseGraph(Graph *graph, void (*f)(GraphNode *, void *),
			   void *arg)
{
  LListCounter *lcount;
  GraphNode *gnode;

  /* No nodes have been visited yet */
  lcount = MakeCounter(graph->nodes, FIRST);
  while ((gnode = Next(lcount)) != NULL)
    gnode->visit = 0;

  /* Visit nodes in order of connection */
  InitCounter(lcount, graph->nodes, FIRST);
  while ((gnode = Next(lcount)) != NULL)
    if (gnode->visit == 0)
      __TGrec(gnode, f, arg);
}

/* Recursive algorithm for traversal */
static void __TGrec(GraphNode *gnode,
		   void (*f)(GraphNode *, void *), void *arg)
{
  LListCounter *lcount;
  GraphNode *nextnode;

  f(gnode, arg);
  gnode->visit = 1;

  /* Visit the neighbours */
  lcount = MakeCounter(gnode->neighbour, FIRST);
  while ((nextnode = Next(lcount)) != NULL) {
    if (nextnode->visit == 0)
      __TGrec(nextnode, f, arg);
  }
}

static void __SGrec(GraphNode *gnode, Graph *newgraph);

LList *SplitGraph(Graph *graph)
{
  LListCounter *lcount;
  GraphNode *gnode;
  LList *glist;
  Graph *newgraph;

  glist = MakeLList();

  /* No nodes have been visited yet */
  lcount = MakeCounter(graph->nodes, FIRST);
  while ((gnode = Next(lcount)) != NULL)
    gnode->visit = 0;

  /* Visit nodes in order of connection */
  InitCounter(lcount, graph->nodes, FIRST);
  while ((gnode = Next(lcount)) != NULL)
    if (gnode->visit == 0) {
      newgraph = MakeGraph();
      __SGrec(gnode, newgraph);
      Enqueue(glist, newgraph);
    }

  return glist;
}

/* Recursive algorithm for traversal when splitting*/
static void __SGrec(GraphNode *gnode, Graph *newgraph)
{
  LListCounter *lcount;
  GraphNode *nextnode;

  Enqueue(newgraph->nodes, gnode);

  gnode->visit = 1;

  /* Visit the neighbours */
  lcount = MakeCounter(gnode->neighbour, FIRST);
  while ((nextnode = Next(lcount)) != NULL) {
    if (nextnode->visit == 0)
      __SGrec(nextnode, newgraph);
  }
}

