#include "align.h"

Align *MakeAlign(int numseq)
{
  Align *align;
  int i;

  align = (Align *)malloc(sizeof(Align));

  align->name = (char **)malloc(numseq * sizeof(char *));
  align->seq = (int **)malloc(numseq * sizeof(int *));
  align->numseq = numseq;

  for (i = 0; i < numseq; i++) {
    align->name[i] = NULL;
    align->seq[i] = NULL;
  }

  return align;
}

/* Reads alignment from file */
Align *ReadAlign(FILE *fp, int (*f)(char, void *), void *arg)
{
  Align *align;
  int i;
  int size;
  Header *header;
  Entry *entry;
  Entry **entry_list;
  int read_error;

  header = MakeHeader();

  if (ReadHeader(fp, header) != 0) {
    fprintf(stderr, "align: Could not read header\n");
    exit(1); }

  size = 1;

  entry_list = (Entry **)malloc(size*sizeof(Entry *));

  for (i = 0; (entry = MakeEntry()) != NULL &&
	      (read_error = ReadEntry(fp, entry)) == 0; i++) {
    if (i >= size) {
      size *= 2;
      entry_list = (Entry **)realloc(entry_list, size*sizeof(Entry *));
    }
    entry_list[i] = entry;
  }

  if (i+1 != size)
    entry_list = (Entry **)realloc(entry_list, (i+1)*sizeof(Entry *));

  entry_list[i] = NULL;

  if (read_error == 1)
    exit(1);

  align = Col2Align(entry_list, f, arg);

  return align;
}

Align *Col2Align(Entry **entry_list, int (*f)(char, void *), void *arg)
{
  int len;
  char name[80];    /* Names are at most 80 positions */
  char field[MAXCOLW];
  int nuc_col;
  int i;
  int entnum, num;
  Align *align;

  num = 0;

  for (entnum = 0; entry_list[entnum] != NULL; entnum++) {
    nuc_col = ReadColno(entry_list[entnum], "residue");
    if (nuc_col == 0)
      continue;
    num++;
  }

  align = MakeAlign(num);

  len = -1;

  num = 0;

  for (entnum = 0; entry_list[entnum] != NULL; entnum++) {
    nuc_col = ReadColno(entry_list[entnum], "residue");
    if (nuc_col == 0) {
      fprintf(stderr, "align: Warning: Ignoring sequence,"
	                " insufficient column info\n");
      continue; }

    if (len == -1)
      len = EntryLength(entry_list[entnum]);
    else if (len != EntryLength(entry_list[entnum])) {
      fprintf(stderr, "align: Sequences not the same length\n");
      exit(1); }

    ReadText(entry_list[entnum], "ENTRY", 80, name);

    align->name[num] = (char *)malloc((strlen(name)+1) * sizeof(char));

    strcpy(align->name[num], name);

    align->seq[num] = (int *)malloc(len * sizeof(int));

    for (i = 0; i < len; i++) {
      GetField(field, entry_list[entnum], i+1, nuc_col);
      align->seq[num][i] = f(field[0], arg);
    }

    num++;
  }

  align->len = len;

  return align;
}

int SeqNumber(Align *align, char *name)
{
  int i;

  if (name == NULL)
    return -1;

  for (i = 0; i < align->numseq; i++)
    if (strcmp(align->name[i], name) == 0)
      return i;

  return -1;
}

void PrintAlign(FILE *fp, Align* align)
{
  int i, j;
  
  for (j = 0; j < align->numseq; j++) {
    fprintf(fp, "%-20s ", align->name[j]);
    for(i = 0; i < align->len; i++)
      fprintf(fp, "%d", align->seq[j][i]);
    fprintf(fp, "\n");
  }
}

Seqlist *MakeSeqlist(int numseq)
{
  Seqlist *seqlist;
  int i;

  seqlist = (Seqlist *)malloc(sizeof(Seqlist));

  seqlist->name = (char **)malloc(numseq * sizeof(char *));
  seqlist->seq = (int **)malloc(numseq * sizeof(int *));
  seqlist->len = (int *)malloc(numseq * sizeof(int));
  seqlist->numseq = numseq;

  for (i = 0; i < numseq; i++) {
    seqlist->name[i] = NULL;
    seqlist->seq[i] = NULL;
  }

  return seqlist;
}

Seqlist *Col2Seqlist(Entry **entry_list, int (*f)(char, void *), void *arg)
{
  int len;
  char name[80];    /* Names are at most 80 positions */
  char field[MAXCOLW];
  int nuc_col;
  int i;
  int num;
  Seqlist *seqlist;

  for (num = 0; entry_list[num] != NULL; num++)
    ;

  seqlist = MakeSeqlist(num);

  for (num = 0; entry_list[num] != NULL; num++) {
    nuc_col = ReadColno(entry_list[num], "residue");
    if (nuc_col == 0) {
      fprintf(stderr, "align: Warning: Ignoring sequence,"
	                " insufficient column info\n");
      continue; }

    len = EntryLength(entry_list[num]);

    seqlist->len[num] = len;

    ReadText(entry_list[num], "ENTRY", 80, name);

    seqlist->name[num] = (char *)malloc((strlen(name)+1) * sizeof(char));

    strcpy(seqlist->name[num], name);

    seqlist->seq[num] = (int *)malloc(len * sizeof(int));

    for (i = 0; i < len; i++) {
      GetField(field, entry_list[num], i+1, nuc_col);
      seqlist->seq[num][i] = f(field[0], arg);
    }
  }

  return seqlist;
}

