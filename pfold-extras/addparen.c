/*********************************************************************

  addparen.c

  usage: addparen [FILE]

  000927 Bjarne Knudsen (bk@daimi.au.dk)

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
#include "graph.h"

typedef struct tagColno {
  int align_bp, alignpos;
  int seq_bp, seqpos;
  int nuc;
  int red, green, blue;
  int red2, green2, blue2;
} Colno;

typedef struct tagNodeInfo {
  int pos;
  int pair;
  int color;
} NodeInfo;

typedef struct tagColor_arg {
  int error;
  int count1;
  int count2;
} Color_arg;

void usage(void);
void color(GraphNode *gnode, void *arg);
void invert(GraphNode *gnode, void *arg);

int main(int argc, char **argv)
{
  FILE *fp;         /* For input */
  Header *header;
  Entry *entry;
  int read_error;   /* For keeping track of errors in reading entries */
  CmdArg *cmdarg;   /* Command line arguments */
  char *s;          /* String for arguments */
  int i, j;
  int pair, pair2;
  int len;
  /*  char field[MAXCOLW];*/
  Colno *colno;     /* Column numbers */
  NodeInfo *nodeinfo;
  GraphNode **gnode;
  Graph *graph;
  Graph *newgraph;
  LList *glist;
  LListCounter *lcount;
  char lpar[] = {'(', '['};
  char rpar[] = {')', ']'};
  Color_arg *color_arg;
  color_arg = (Color_arg *) malloc(sizeof(Color_arg));

  colno = (Colno *) malloc(sizeof(Colno));

  cmdarg = InitArgument(argc, argv);

  if ((s = GetArgument(cmdarg)) != NULL) {
    usage();
    return 1; }

  if ((s = GetFilename(cmdarg)) == NULL)
    fp = stdin;
  else if (GetFilename(cmdarg) != NULL) {
    usage();
    return 1; }
  else if ((fp = fopen(s, "r")) == NULL) {
    fprintf(stderr, "addparen: Error in opening file '%s'\n", s);
    return 1; }

  header = MakeHeader();
  entry = MakeEntry();

  if (ReadHeader(fp, header) != 0)
    return 1;

  AddHeaderInfo(header, argc, argv);

  PrintHeader(stdout, header);

  while ((read_error = ReadEntry(fp, entry)) == 0) {
    PrintEntry(stdout, entry);
    if (!ReadType(entry, "RNA") && !ReadType(entry, "DNA"))
      continue;
    if ((colno->nuc = ReadColno(entry, "residue")) == 0)
      colno->nuc = ReadColno(entry, "nucleotide");
    colno->align_bp = ReadColno(entry, "align_bp");
    colno->alignpos = ReadColno(entry, "alignpos");
    colno->seq_bp = ReadColno(entry, "seq_bp");
    colno->seqpos = ReadColno(entry, "seqpos");
    if (colno->nuc == 0 ||
	((colno->align_bp == 0 || colno->alignpos == 0) &&
	 (colno->seq_bp == 0 || colno->seqpos == 0))) {
      fprintf(stderr, "stdpair: Warning: Ignoring sequence,"
	               " insufficient column info\n");
      continue; }
    len = EntryLength(entry);

    gnode = (GraphNode **)malloc(len * sizeof(GraphNode *));    
    graph = MakeGraph();
    for (i = 1; i <= len ; i++) {
      pair = FindPair(entry, i, colno->align_bp, colno->alignpos,
		      colno->seq_bp, colno->seqpos);
      if (pair >= i) {  /* Left side of pair */
	nodeinfo = (NodeInfo *)malloc(sizeof(NodeInfo));
	nodeinfo->pos = i;
	nodeinfo->pair = pair;
	nodeinfo->color = 0;
	gnode[i-1] = AddGNode(graph, nodeinfo);
      }
      else
	gnode[i-1] = NULL;
    }

    for (i = 1; i <= len ; i++) {
      pair = FindPair(entry, i, colno->align_bp, colno->alignpos,
		      colno->seq_bp, colno->seqpos);
      if (pair > i)   /* Left side of pair */
	for (j = i+1; j < pair; j++) {
	  pair2 = FindPair(entry, j, colno->align_bp, colno->alignpos,
			   colno->seq_bp, colno->seqpos);
	  if (pair2 == 0) /* no pair */
	    continue;
	  else if (pair2 > pair) /* a cross */
	    ConnectGNode(gnode[i-1], gnode[j-1]);
	  else if (pair2 < i) /* a cross */
	    ConnectGNode(gnode[i-1], gnode[pair2-1]);
	}
    }

    glist = SplitGraph(graph);

    color_arg->error = 0;

    lcount = MakeCounter(glist, FIRST);
    while ((newgraph = Next(lcount)) != NULL) {
      color_arg->error = 0;
      color_arg->count1 = 0;
      color_arg->count2 = 0;
      TraverseGraph(newgraph, color, color_arg);
      if (color_arg->error == -1) {
	fprintf(stderr, "addparen: structure too knotted, ignoring\n");
	break;
      }
      if (color_arg->count2 > color_arg->count1)
	TraverseGraph(newgraph, invert, NULL);
    }

    if (color_arg->error == 0) {
      for (i = 1; i <= len ; i++)
	ChgField(entry, i, colno->nuc, ".");
      for (i = 1; i <= len ; i++)
	if (gnode[i-1] != NULL) {
	  if (((NodeInfo *)gnode[i-1]->elm)->pos ==
	      ((NodeInfo *)gnode[i-1]->elm)->pair)
	    ChgField(entry, i, colno->nuc, "x");
	  else {
	    ChgField(entry, i, colno->nuc, "%c", 
		     lpar[((NodeInfo *)gnode[i-1]->elm)->color-1]);
	    ChgField(entry, ((NodeInfo *)gnode[i-1]->elm)->pair,
		     colno->nuc, "%c",
		     rpar[((NodeInfo *)gnode[i-1]->elm)->color-1]);
	  }
	}
      /*
      printf("; TYPE              parenthesis\n");
      printf("; COL 1             label\n");
      printf("; COL 2             residue\n");
      printf("; COL 3             seqpos\n");
      printf("; ENTRY             -\n");
      printf("; ----------\n");

      for (i = 1; i <= len; i++) {
	GetField(field, entry, i, colno->nuc);
	printf("P     %s %5d\n", field, i);
      }
      printf("; **********\n");
      */
      PrintEntry(stdout, entry);
    }

    free(lcount);
      
    for (i = 0; i < len ; i++)
      if (gnode[i] != NULL)
	free((NodeInfo *)Element(gnode[i]));
      
    FreeGraph(graph);
      
    free(gnode);

  }

  if (fp != stdin && fclose(fp) != 0) {
    fprintf(stderr, "stdpair: Error in closing file\n");
    return 1; }

  if (read_error == 1)
    return 1;

  return 0;
}

void usage(void)
{
  fprintf(stderr, "usage: addparen [FILE]\n");
}

void color(GraphNode *gnode, void *arg)
{
  LListCounter *lcount;
  GraphNode *newgnode;
  int newcolor;

  if (((Color_arg *)arg)->error == -1)
    return;

  newcolor = 0;

  lcount = MakeCounter(Neighbour(gnode), FIRST);
  while ((newgnode = Next(lcount)) != NULL) {
    if (((NodeInfo *)newgnode->elm)->color == 1) {
      if (newcolor == 1) {  /* impossible coloration */
	newcolor = -1;
	((Color_arg *)arg)->error = -1;
	break;
      }
      else if (newcolor == 0)  /* just fine */
	newcolor = 2;
    }
    else if (((NodeInfo *)newgnode->elm)->color == 2) {
      if (newcolor == 2) {  /* impossible coloration */
	newcolor = -1;
	((Color_arg *)arg)->error = -1;
	break;
      }
      else if (newcolor == 0)  /* just fine */
	newcolor = 1;
    }
  }

  free(lcount);

  if (newcolor == 0) /* neighbours not colored */
    newcolor = 1;

  ((NodeInfo *)gnode->elm)->color = newcolor;
  if (newcolor == 1) {
    ((Color_arg *)arg)->count1++;
  }
  else if (newcolor == 2)
    ((Color_arg *)arg)->count2++;
}

void invert(GraphNode *gnode, void *arg)
{
  if (((NodeInfo *)gnode->elm)->color == 1)
    ((NodeInfo *)gnode->elm)->color = 2;
  else if (((NodeInfo *)gnode->elm)->color == 2)
    ((NodeInfo *)gnode->elm)->color = 1;
}
