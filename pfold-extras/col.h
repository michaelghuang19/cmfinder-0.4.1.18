/*********************************************************************

  col.h

  Contains general functions for handling the col format.

  000119 Bjarne Knudsen (bk@daimi.au.dk)

*********************************************************************/

#ifndef __col_h__
#define __col_h__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include "file.h"

#define MAXCOLW   1000    /* Maximum width of a single column */

typedef struct tagHeader {
  char **text;         /* Text in header */
  int len;             /* Number of header lines */
} Header;

typedef struct tagEntry {
  char **text;         /* Text in entry */
  int textlen;         /* Number of text lines */
  char **line;         /* Sequence info */
  int len;             /* Sequence length */
  char *endline;       /* For '; **********' line */
} Entry;

typedef struct tagColfile {
  Header *header;
  Entry **entry;
} Colfile;

Header *MakeHeader(void);
int ReadHeader(FILE *fp, Header *header);
void AddHeaderText(Header *header, const char *format, ...);
void AddHeaderInfo(Header *header, int argc, char **argv);
void PrintHeader(FILE *fp, Header *header);

Entry *MakeEntry(void);
Entry *NewEntry(char *type, char *name, int len);
void CopyEntry(Entry *dest, Entry *source);
int ReadEntry(FILE *fp, Entry *entry);
void PrintEntry(FILE *fp, Entry *entry);
void PrintEntryLine(FILE *fp, Entry *entry, int lineno);
char *CopyEntryLine(Entry *entry, int lineno);

Colfile *ReadColfile(FILE *fp);

int ReadEntryText(FILE *fp, Entry *entry);
int ReadEntryLines(FILE *fp, Entry *entry, int num);

void PrintEntryText(FILE *fp, Entry *entry);
void PrintEntryLines(FILE *fp, Entry *entry);
void PrintEntryEnd(FILE *fp, Entry *entry);

void PrintEntryText_tex(FILE *fp, Entry *entry);

void AddEntryText(Entry *entry, const char *format, ...);
void AddEntryExt(Entry *entry, const char *format, ...);

void AddEntryLine(Entry *entry);

int EntryLength(Entry *entry);

int ReadType(Entry *entry, const char *type);
int ReadColno(Entry *entry, const char *keyword);
int ReadColnoND(Entry *entry, const char *keyword);
int ReadCollabel(Entry *entry, int colno, char *label);
int ChgCollabel(Entry *entry, int colno, const char *label);
int ReadVar(Entry *entry, const char *keyword, char *value);
int ReadText(Entry *entry, const char *keyword, int maxlen, char *value);
int ChgEntryText(Entry *entry, const char *keyword, char *newtext);
int RmvColno(Entry *entry, const char *keyword);

int CmpLine(Entry *entry1, Entry *entry2, int lineno);
void RmvLine(Entry *entry, int lineno);

char *LineAddress(Entry *entry, int lineno);

int CountField(Entry *entry, int lineno);
int GetField(char *field, Entry *entry, int lineno, int colno);
int AddField(Entry *entry, int lineno, const char *format, ...);
int ChgField(Entry *entry, int lineno, int colno, const char *format, ...);
int RmvField(Entry *entry, int lineno, int colno);

int EnsureCol(Entry *entry, const char *keyword, const char *init);
int RmvCol(Entry *entry, const char *keyword);

int InSeq(Entry *entry, int lineno);

void EntryMap(Entry *entry, void (*f)(Entry *, int, void *), void *arg);

int FindPair(Entry *entry, int lineno, int align_bp_col, int alignpos_col,
	     int seq_bp_col, int seqpos_col);
void SetPair(Entry *entry, int lineno, int pairline,
	     int align_bp_col, int alignpos_col,
	     int seq_bp_col, int seqpos_col);

#endif
