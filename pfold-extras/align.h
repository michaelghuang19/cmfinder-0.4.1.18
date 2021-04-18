#include "col.h"
#include "file.h"

#ifndef __align_h__   /* Only define the following once */
#define __align_h__

typedef struct tagAlign {
  int len;
  int numseq;
  char **name;
  int **seq;
} Align;

typedef struct tagSeqlist {
  int *len;
  int numseq;
  char **name;
  int **seq;
} Seqlist;

Align *MakeAlign(int numseq);
Align *ReadAlign(FILE *fp, int (*f)(char, void *), void *arg);
Align *Col2Align(Entry **entry_list, int (*f)(char, void *), void *arg);
int SeqNumber(Align *align, char *name);
void PrintAlign(FILE *fp, Align* align);

Seqlist *MakeSeqlist(int numseq);
Seqlist *Col2Seqlist(Entry **entry_list, int (*f)(char, void *), void *arg);

#endif
