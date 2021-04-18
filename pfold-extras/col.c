/*********************************************************************

  col.c

  Contains general functions for handling the col format.

  Low memory use is preferred to fast execution.

  000111 Bjarne Knudsen (bk@daimi.au.dk)

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

/* Allocates Header */
Header *MakeHeader()
{
  Header *header;
  
  header = (Header *)malloc(sizeof(Header));
  header->text = NULL;

  return header;
}
/* 
   Reads database header from file point fp. Returns 0 on ok, 1
   otherwise.
*/
int ReadHeader(FILE *fp, Header *header)
{
  char *s;         /* String for keeping lines */
  int ptr;         /* Pointer to line number */
  int maxlen;      /* Dynamic line allocation maximum */

  if (header->text != NULL) {  /* Free old space */
    for (ptr = 0; ptr < header->len; ptr++)
      free(header->text[ptr]);
    free(header->text);
  }

  maxlen = 20;     /* Start with 20 lines */

  header->text = (char **)malloc(maxlen * sizeof(char *));

  for (ptr = 0; (s = GetLine(fp)) != NULL; ptr++) {
    if (ptr == maxlen) {
      maxlen *= 2;
      header->text = (char **)realloc(header->text, maxlen * sizeof(char *));
                                                     /* Get more space */
    }
    header->text[ptr] = s;
    if (StrnCmp(header->text[ptr], "; ==========", 12) == 0) {
      ptr++;            /* Remember last line */
      break; }
  }

  if (s == NULL) {   /* EOF before end of header */
    header->text = NULL;
    fprintf(stderr, "Error in header\n");
    return 1; }

  header->len = ptr;
  header->text = (char **)realloc(header->text, ptr * sizeof(char *));
                                             /* Free excess space */
  return 0;
}

/* Adds text before '; ==========' line in header text */
void AddHeaderText(Header *header, const char *format, ...)
{
  va_list args;
  char s[MAXCOLW*10];

  header->len++;

  header->text = (char **)realloc(header->text,
				 (header->len+1) * sizeof(char *));

  header->text[header->len-1] = header->text[header->len-2];

  va_start(args, format);
  vsprintf(s, format, args);
  va_end(args);

  header->text[header->len-2] =
              (char *)malloc((strlen(s)+1) * sizeof(char));

  strcpy(header->text[header->len-2], s);
}

/*
   Inserts info on what has been done to the file in the header text.
*/
void AddHeaderInfo(Header *header, int argc, char **argv)
{
  CmdArg *cmdarg;   /* Command line arguments */
  int ptr;          /* For making header info */
  char *t;          /* For making header info */
  char *s;

  AddHeaderText(header, ";\n");
  cmdarg = InitArgument(argc, argv);
  t = malloc((strlen(argv[0])+4000) * sizeof(char));
  sprintf(t, "; '%s", argv[0]);
  ptr = strlen(t);
  while ((s = GetArgument(cmdarg)) != NULL) {
    if (strlen(s) > 30) {
      t = realloc(t, (ptr+30+13) * sizeof(char));
      sprintf(t+ptr, " -%27.27s...", s);
    }
    else {
      t = realloc(t, (ptr+strlen(s)+13) * sizeof(char));
      sprintf(t+ptr, " -%s", s);
    }
    ptr = strlen(t);
  }
  t = realloc(t, (ptr+24) * sizeof(char));
  sprintf(t+ptr, "' was run on this file\n");
  AddHeaderText(header, "%s", t);
}

/*
   Checks the sequence type. Returns true if the type is type, false
   otherwise.
*/
int ReadType(Entry *entry, const char *type)
{
  char value[MAXCOLW];
  
  if (ReadText(entry, "TYPE", MAXCOLW, value) == 0 &&
      StrCmp(type, value) == 0)
    return 1;
  
  return 0;
}

/*
   Reads column number from entry text. The return value is the column
   number corresponding to the keyword. Returns 0 on error. No defaults.
*/
int ReadColnoND(Entry *entry, const char *keyword)
{
  int ptr;      /* For the text line number */
  int value;
  char s[MAXCOLW];

  for (ptr = 0; ptr < entry->textlen; ptr++) {
    if (sscanf(entry->text[ptr], "; COL %d %s", &value, s) == 2 &&
	StrCmp(s, keyword) == 0)
	return value;
  }

  return 0;
}

/*
   Reads column number from entry text. The return value is the column
   number corresponding to the keyword. Returns 0 on error.
*/
int ReadColno(Entry *entry, const char *keyword)
{
  int ptr;      /* For the text line number */
  int value;
  char s[MAXCOLW];
  int colinfo = 0;    /* =1 if column info is present */

  for (ptr = 0; ptr < entry->textlen; ptr++)
    if (sscanf(entry->text[ptr], "; COL %d %s", &value, s) == 2) {
      colinfo = 1;
      if (StrCmp(s, keyword) == 0)
	return value;
    }

  if (colinfo == 0) {
    if (StrCmp(keyword, "label") == 0)  /* Default value */
      return 1;
    else if (StrCmp(keyword, "nucleotide") == 0 ||
	     StrCmp(keyword, "aminoacid") == 0)  /* Default value */
      return 2;
    else if (StrCmp(keyword, "alignpos") == 0)  /* Default value */
      return 3;
    else if (StrCmp(keyword, "seqpos") == 0)  /* Default value */
      return 4;
    else if (StrCmp(keyword, "align_bp") == 0)  /* Default value */
      return 5;
  }

  return 0;
}

/*
   Reads column label from entry text. The return value is 1 on
   error. label should be able to hold the label.
*/
int ReadCollabel(Entry *entry, int colno, char *label)
{
  int ptr;      /* For the text line number */
  int value;
  char s[MAXCOLW];

  for (ptr = 0; ptr < entry->textlen; ptr++)
    if (sscanf(entry->text[ptr], "; COL %d %s", &value, s) == 2)
      if (value == colno) {
	strcpy(label, s);
	return 0; }
      
  return 1;
}

/*
   Changes column label in entry text. The return value is 1 on
   error.
*/
int ChgCollabel(Entry *entry, int colno, const char *label)
{
  int ptr;      /* For the text line number */
  int value;
  char *s;

  for (ptr = 0; ptr < entry->textlen; ptr++)
    if (sscanf(entry->text[ptr], "; COL %d", &value) == 1)
      if (value == colno) {
	s = (char *)malloc((strlen(label)+21)*sizeof(char));
	sprintf(s, "; COL %-3d           %s\n", value, label);
	entry->text[ptr] = s;
	return 0; }
      
  return 1;
}

/*
   Reads variable from entry text. Returns 0 on ok, 1 otherwise.
*/
int ReadVar(Entry *entry, const char *keyword, char *value)
{
  int ptr;      /* For the text line number */
  char s[MAXCOLW];
  
  for (ptr = 0; ptr < entry->textlen; ptr++)
    if (sscanf(entry->text[ptr], "; %s %s", s, value) == 2 &&
	StrCmp(s, keyword) == 0)
      return 0;

  return 1;
}

/*
   Reads entry text. Returns 0 on ok, 1 otherwise.
*/
int ReadText(Entry *entry, const char *keyword, int maxlen, char *value)
{
  int ptr;      /* For the text line number */
  int pos;
  int i;
  
  for (ptr = 0; ptr < entry->textlen; ptr++)
    if (sscanf(entry->text[ptr], "; %n", &pos) == 0 &&
	StrnCmp(&entry->text[ptr][pos], keyword, strlen(keyword)) == 0) {
      pos += strlen(keyword);
      while (isspace(entry->text[ptr][pos]))
	pos++;
      for (i = 0; entry->text[ptr][pos+i] != '\n' && i < maxlen-1; i++)
	value[i] = entry->text[ptr][pos+i];
      value[i] = '\0';
      return 0;
    }

  return 1;
}

/*
   Changes entry text. Returns 0 on ok, 1 otherwise.
*/
int ChgEntryText(Entry *entry, const char *keyword, char *newtext)
{
  int ptr;      /* For the text line number */
  int pos;
  int i;
  char *temp;
  
  for (ptr = 0; ptr < entry->textlen; ptr++)
    if (sscanf(entry->text[ptr], "; %n", &pos) == 0 &&
	StrnCmp(&entry->text[ptr][pos], keyword, strlen(keyword)) == 0) {
      pos += strlen(keyword);
      while (isspace(entry->text[ptr][pos]))
	pos++;
      temp = (char *)malloc((pos+strlen(newtext)+2)*sizeof(char));
      strncpy(temp, entry->text[ptr], pos);
      for (i = 0; newtext[i] != '\0'; i++)
	temp[pos+i] = newtext[i]; 
      temp[pos+i] = '\n';
      temp[pos+i+1] = '\0';
      free(entry->text[ptr]);
      entry->text[ptr] = temp;
      return 0;
    }

  return 1;
}

/*
   Removes column number info from entry text. The return value is the
   column number corresponding to the keyword. Returns 0 on error.
*/
int RmvColno(Entry *entry, const char *keyword)
{
  int ptr;      /* For the text line number */
  int value, i;
  char s[MAXCOLW];
  
  for (ptr = 0; ptr < entry->textlen; ptr++)
    if (sscanf(entry->text[ptr], "; COL %d %s", &value, s) == 2 &&
	StrCmp(s, keyword) == 0) {
      break;
    }

  if (ptr == entry->textlen)
    return 0;

  entry->text[ptr] = realloc(entry->text[ptr], 1);
  entry->text[ptr][0] = '\0';

  /* Change other column numbers */
  for (ptr = 0; ptr < entry->textlen; ptr++)
    if (sscanf(entry->text[ptr], "; COL %d %s", &i, s) == 2 &&	
	i > value)
      sprintf(entry->text[ptr], "; COL %-5d         %s\n", i-1, s);

  return value;
}


/*
   Allocates Entry using info from header. Returns NULL on error.
*/
Entry *MakeEntry()
{
  Entry *entry;

  entry = (Entry *)malloc(sizeof(Entry));
  entry->text = NULL;
  entry->line = NULL;

  return entry;
}

/*
   Makes an entry, ready for filling out.
*/
Entry *NewEntry(char *type, char *name, int len)
{
  Entry *entry;
  int i;

  entry = MakeEntry();

  entry->textlen = 3;
  entry->text = (char **)malloc(entry->textlen * sizeof(char *));

  entry->text[0] = (char *)malloc((22+strlen(type)) * sizeof(char *));
  sprintf(entry->text[0], "; TYPE              %s\n", type);
  
  entry->text[1] = (char *)malloc((22+strlen(name)) * sizeof(char *));
  sprintf(entry->text[1], "; ENTRY             %s\n", name);
  
  entry->text[2] = (char *)malloc(14 * sizeof(char *));
  sprintf(entry->text[2], "; ----------\n");

  entry->len = len;

  entry->line = (char **)malloc(len * sizeof(char *));
  for (i = 0; i < len; i++) {
    entry->line[i] = (char *)malloc(2 * sizeof(char *));
    strcpy(entry->line[i], "\n");
  }

  entry->endline = (char *)malloc(14 * sizeof(char *));
  strcpy(entry->endline, "; **********\n");

  return entry;
}

/*
   Copies source to dest entry.
*/
void CopyEntry(Entry *dest, Entry *source)
{
  int ptr;

  if (dest->text != NULL) {  /* Free old space */
    for (ptr = 0; ptr < dest->textlen; ptr++)
      free(dest->text[ptr]);
    free(dest->text);
    for (ptr = 0; ptr < dest->len; ptr++)
      free(dest->line[ptr]);
    free(dest->line);
    free(dest->endline);
  }

  dest->text = (char **)malloc(source->textlen * sizeof(char *));
  for (ptr = 0; ptr < source->textlen; ptr++) {
    dest->text[ptr] =
      (char *)malloc((strlen(source->text[ptr])+1) * sizeof(char));
    strcpy(dest->text[ptr], source->text[ptr]);
  }
  dest->textlen = source->textlen;

  dest->line = (char **)malloc(source->len * sizeof(char *));
  for (ptr = 0; ptr < source->len; ptr++) {
    dest->line[ptr] =
      (char *)malloc((strlen(source->line[ptr])+1) * sizeof(char));
    strcpy(dest->line[ptr], source->line[ptr]);
  }
  dest->len = source->len;
  
  dest->endline = (char *)malloc((strlen(source->endline)+1) * sizeof(char));
  strcpy(dest->endline, source->endline);
}

/*
   Reads entry from file point fp. Returns -1 on no more entries, 0 if
   ok, 1 otherwise.
*/
int ReadEntry(FILE *fp, Entry *entry)
{
  char *s;         /* String for keeping lines */
  int ptr;         /* Pointer to line number */
  int maxlen;      /* Dynamic line allocation maximum */

  /* Read entry text */
  
  if (entry->text != NULL) {  /* Free old space */
    for (ptr = 0; ptr < entry->textlen; ptr++)
      free(entry->text[ptr]);
    free(entry->text);
    for (ptr = 0; ptr < entry->len; ptr++)
      free(entry->line[ptr]);
    free(entry->line);
    free(entry->endline);
  }

  maxlen = 20;     /* Start with 20 lines */

  entry->text = (char **)malloc(maxlen * sizeof(char *));

  for (ptr = 0; (s = GetLine(fp)) != NULL; ptr++) {
    if (ptr == maxlen) {
      maxlen *= 2;
      entry->text = (char **)realloc(entry->text, maxlen * sizeof(char *));
                                                     /* Get more space */
    }
    entry->text[ptr] = s;
    if (StrnCmp(entry->text[ptr], "; ----------", 12) == 0) {
      ptr++;            /* Remember last line */
      break; }
  }

  if (s == NULL) {   /* EOF before end of entry text */
    if (ptr == 0) {
      entry->text = NULL;
      return -1; } /* No more entries */
    else {
      entry->text = NULL;
      fprintf(stderr, "Error in entry text\n");
      return 1; }
  }

  entry->textlen = ptr;
  entry->text = (char **)realloc(entry->text, ptr * sizeof(char *));
                                             /* Free excess space */

  /* Read sequence */

  maxlen = 100;     /* Start with 100 lines */

  entry->line = (char **)malloc(maxlen * sizeof(char *));

  for (ptr = 0; (s = GetLine(fp)) != NULL; ptr++) {
    if (StrnCmp(s, "; **********", 12) == 0) {
      entry->endline = s;
      break; }
    if (ptr == maxlen) {
      maxlen *= 2;
      entry->line = (char **)realloc(entry->line, maxlen * sizeof(char *));
                                                     /* Get more space */
    }
    entry->line[ptr] = s;
  }

  if (s == NULL) {   /* EOF before end of entry */
    entry->text = NULL;
    fprintf(stderr, "Error in entry\n");
    return 1; }

  entry->len = ptr;

  entry->line = (char **)realloc(entry->line, ptr * sizeof(char *));
                                             /* Free excess space */
  return 0;
}

/* Adds text before '; ----------' line in entry text */
void AddEntryText(Entry *entry, const char *format, ...)
{
  va_list args;
  char s[MAXCOLW];

  entry->textlen++;

  entry->text = (char **)realloc(entry->text,
				 entry->textlen * sizeof(char *));

  entry->text[entry->textlen-1] = entry->text[entry->textlen-2];

  va_start(args, format);
  vsprintf(s, format, args);
  va_end(args);

  entry->text[entry->textlen-2] =
              (char *)malloc((strlen(s)+1) * sizeof(char));

  strcpy(entry->text[entry->textlen-2], s);
}

/* 
   Adds text after '; TYPE' and '; COL' lines in entry text. The added
   text would usually be for a new column annotation.
*/
void AddEntryExt(Entry *entry, const char *format, ...)
{
  va_list args;
  char s[MAXCOLW];
  int newnum;
  int i;
  int ptr;

  newnum = 0;

  for (ptr = 0; ptr < entry->textlen; ptr++) {
    i = 0;
    sscanf(entry->text[ptr], "; TYPE %n", &i);
    sscanf(entry->text[ptr], "; COL %n", &i);
    if (i != 0)
      newnum = ptr+1;
    else
      break;
  }

  entry->textlen++;

  entry->text = (char **)realloc(entry->text,
				 entry->textlen * sizeof(char *));

  /* Move the reminder of the entry text */
  for (ptr = entry->textlen-1; ptr > newnum; ptr--)
    entry->text[ptr] = entry->text[ptr-1];

  va_start(args, format);
  vsprintf(s, format, args);
  va_end(args);

  entry->text[newnum] =
              (char *)malloc((strlen(s)+1) * sizeof(char));

  strcpy(entry->text[newnum], s);
}

/* Copies the last line in the entry */
void AddEntryLine(Entry *entry)
{
  entry->len++;

  entry->line = (char **)realloc(entry->line, entry->len * sizeof(char *));

  entry->line[entry->len-1] =
    (char *)malloc((strlen(entry->line[entry->len-2])+1) * sizeof(char));

  strcpy(entry->line[entry->len-1], entry->line[entry->len-2]);
}

/* Returns the length of the entry */
int EntryLength(Entry *entry)
{
  return entry->len;
}

/* Outputs database header to file point fp */
void PrintHeader(FILE *fp, Header *header)
{
  int i;

  for (i = 0; i < header->len; i ++)
    fprintf(fp, "%s", header->text[i]);
}

/* Outputs entry to file pointer fp */
void PrintEntry(FILE *fp, Entry *entry)
{
  int i;

  for (i = 0; i < entry->textlen; i ++)
    fprintf(fp, "%s", entry->text[i]);

  for (i = 0; i < entry->len; i ++)
    fprintf(fp, "%s", entry->line[i]);

  fprintf(fp, "%s", entry->endline);
}

/* Outputs single entry lines to file pointer fp */
void PrintEntryLine(FILE *fp, Entry *entry, int lineno)
{
  fprintf(fp, "%s", entry->line[lineno-1]);
}

char *CopyEntryLine(Entry *entry, int lineno)
{
  char *s;
  s = (char *)malloc((strlen(entry->line[lineno-1])+1)*sizeof(char));
  strcpy(s, entry->line[lineno-1]);

  return s;
}


Colfile *ReadColfile(FILE *fp)
{
  int read_error;   /* For keeping track of errors in reading entries */
  Colfile *colfile;
  Entry *entry;
  int i;

  entry = MakeEntry();

  colfile = (Colfile *)malloc(sizeof(Colfile));
  colfile->header = MakeHeader();
  colfile->entry = NULL;

  if (ReadHeader(fp, colfile->header) != 0)
    return NULL;

  for (i = 0; (read_error = ReadEntry(fp, entry)) == 0; i++) {
    colfile->entry = (Entry **)realloc(colfile->entry, (i+1)*sizeof(Entry *));
    colfile->entry[i] = MakeEntry();
    CopyEntry(colfile->entry[i], entry);
  }

  colfile->entry = (Entry **)realloc(colfile->entry, (i+1)*sizeof(Entry *));
  colfile->entry[i] = NULL;

  if (read_error == 1)
    return NULL;

  return colfile;
}

/*
   Reads entry text from file point fp. Returns -1 on no more entries,
   0 if ok, 1 otherwise.
*/
int ReadEntryText(FILE *fp, Entry *entry)
{
  char *s;         /* String for keeping lines */
  int ptr;         /* Pointer to line number */
  int maxlen;      /* Dynamic line allocation maximum */

  /* Read entry text */

  if (entry->text != NULL) {  /* Free old space */
    for (ptr = 0; ptr < entry->textlen; ptr++)
      free(entry->text[ptr]);
    free(entry->text);
    free(entry->endline);
  }

  maxlen = 20;     /* Start with 20 lines */

  entry->text = (char **)malloc(maxlen * sizeof(char *));

  for (ptr = 0; (s = GetLine(fp)) != NULL; ptr++) {
    if (ptr == maxlen) {
      maxlen *= 2;
      entry->text = (char **)realloc(entry->text, maxlen * sizeof(char *));
                                                     /* Get more space */
    }
    entry->text[ptr] = s;
    if (StrnCmp(entry->text[ptr], "; ----------", 12) == 0) {
      ptr++;            /* Remember last line */
      break; }
  }

  if (s == NULL) {   /* EOF before end of entry text */
    if (ptr == 0) {
      entry->text = NULL;
      return -1; } /* No more entries */
    else {
      entry->text = NULL;
      fprintf(stderr, "Error in entry text\n");
      return 1; }
  }

  entry->textlen = ptr;
  entry->text = (char **)realloc(entry->text, ptr * sizeof(char *));
                                             /* Free excess space */

  return 0;
}

/*
   Reads num entry lines from file point fp. Returns -1 on no more
   lines, 0 if ok, 1 otherwise. The entry will appear to have only num
   lines.
*/
int ReadEntryLines(FILE *fp, Entry *entry, int num)
{
  int ptr;         /* Pointer to line number */
  char *s;         /* String for keeping lines */
  int endofentry;

  if (entry->line != NULL) {
    for (ptr = 0; ptr < entry->len; ptr++)
      free(entry->line[ptr]);
    free(entry->line);
  }
    
  /* Read sequence */

  entry->line = (char **)malloc(num * sizeof(char *));

  endofentry = 0;

  for (ptr = 0; ptr < num && (s = GetLine(fp)) != NULL; ptr++) {
    if (StrnCmp(s, "; **********", 12) == 0) {
      entry->endline = s;
      endofentry = 1;
      break; }
    entry->line[ptr] = s;
  }

  if (s == NULL) {   /* EOF before end of entry */
    entry->line = NULL;
    fprintf(stderr, "Error in entry\n");
    return 1; }

  entry->len = ptr;

  if (ptr < num)
    entry->line = (char **)realloc(entry->line, ptr * sizeof(char *));
                                             /* Free excess space */

  if (endofentry == 1) {
    entry->line = NULL;
    return -1; }

  return 0;
}

/*
   Outputs entry text to file pointer fp.
*/
void PrintEntryText(FILE *fp, Entry *entry)
{
  int i;

  for (i = 0; i < entry->textlen; i ++) {
    fprintf(fp, "%s", entry->text[i]);
  }
}

/*
   Outputs a single entry line.
*/
void PrintEntryLines(FILE *fp, Entry *entry)
{
  int i;

  for (i = 0; i < entry->len; i ++)
    fprintf(fp, "%s", entry->line[i]);
}

/*
   Outputs entry end (; *******) to file pointer fp.
*/
void PrintEntryEnd(FILE *fp, Entry *entry)
{
  fprintf(fp, "%s", entry->endline);
}

/*
   Outputs entry text without '; %%' lines and initial ';' to file
   pointer fp.
*/
void PrintEntryText_tex(FILE *fp, Entry *entry)
{
  int i, j;
  char junks[MAXCOLW];

  for (i = 0; i < entry->textlen &&
	 sscanf(entry->text[i], "; %%%% %s", junks) != 1; i ++) {
    j = 0;
    sscanf(entry->text[i],"; %n", &j);
    fprintf(fp, "\\noindent %s\n", &entry->text[i][j]);
  }
}

/* Compares lines in different entries */
int CmpLine(Entry *entry1, Entry *entry2, int lineno)
{
  return StrCmp(entry1->line[lineno-1], entry2->line[lineno-1]);
}

/* 
  Removes a line from entry. The linenumbers of all lines are still
  unchanged.
*/
void RmvLine(Entry *entry, int lineno)
{
  entry->line[lineno-1] = realloc(entry->line[lineno-1], 1);
  entry->line[lineno-1][0] = '\0';
}

/* Returns the address of the line */
char *LineAddress(Entry *entry, int lineno)
{
  return entry->line[lineno-1];
}

/* Returns the number of columns */
int CountField(Entry *entry, int lineno)
{
  int ptr, width;
  int colno;
  char *s;

  s = (char *)malloc((strlen(entry->line[lineno-1])+1) * sizeof(char));

  ptr = 0;
  colno = 0;

  while (sscanf(&entry->line[lineno-1][ptr], " %s%n", s, &width) == 1) {
    ptr += width;
    colno++;
  }

  free(s);

  return colno;
}

/*
   The function returns a specific column from a line.

      line    the input line
      colno   the column of interest
      field   the returned field, should be allocated in advance

   Returns 0 if ok. Otherwise 1.
*/
int GetField(char *field, Entry *entry, int lineno, int colno)
{
  int ptr;
  int i, width;
  char *s;

  s = (char *)malloc((strlen(entry->line[lineno-1])+1) * sizeof(char));

  ptr = 0;

  for (i = 1; i < colno; i++) {
    if (sscanf(&entry->line[lineno-1][ptr], " %s%n", s, &width) != 1) {
      fprintf(stderr, "Column does not exist %d\n", colno);
      return 1; }
    ptr += width;
  }

  if (sscanf(&entry->line[lineno-1][ptr], " %s", field) != 1) {
    fprintf(stderr, "Column does not exist %d\n", colno);
    return 1; }

  free(s);

  return 0;
}

/*
   The function adds column(s) after a line.

      line    the input line
      format  the format the column(s) should be written in
      ...     what to put in the column(s), like printf

   The returned value is the new number of columns.
*/
int AddField(Entry *entry, int lineno, const char *format, ...)
{
  va_list args;
  char *s;
  int len;             /* Length of line */

  s = (char *)malloc((strlen(entry->line[lineno-1])+1+MAXCOLW) * sizeof(char));

  len = strlen(entry->line[lineno-1]);
  strcpy(s, entry->line[lineno-1]);

  s[len-1] = ' ';           /* Change newline to space */

  va_start(args, format);
  vsprintf(&s[len], format, args);  /* Append arguments */
  va_end(args);

  len = strlen(s);

  sprintf(&s[len], "\n");           /* Add newline */
  len++;

  entry->line[lineno-1] =
    (char *)realloc(entry->line[lineno-1], (len + 1) * sizeof(char));
                         /* Reallocate line due to greater length */

  strcpy(entry->line[lineno-1], s);

  free(s);

  return CountField(entry, lineno);
}

/*
   The function changes a column in a line.

      line    the input line
      colno   the column to be changed
      format  the format the column(s) should be written in
      value   what to put in the column(s)

   Returns 0 if ok. Otherwise 1.
*/
int ChgField(Entry *entry, int lineno, int colno, const char *format, ...)
{
  va_list args;
  int ptr;                      /* Pointer in line */
  int i;
  int width, newwidth;          /* Width of column */
  char *s, t[MAXCOLW];
  char tmpfmt[20];              /* A string used for a format */

  s = (char *)malloc((strlen(entry->line[lineno-1])+1+MAXCOLW) * sizeof(char));

  ptr = 0;

  for (i = 1; i < colno; i++) {
    if (sscanf(&entry->line[lineno-1][ptr], " %s%n", s, &width) != 1) {
      fprintf(stderr, "Column does not exist %d\n", colno);
      return 1; }
    ptr += width;
  }

  if (sscanf(&entry->line[lineno-1][ptr], " %s%n", s, &width) != 1) {
    fprintf(stderr, "Column does not exist %d\n", colno);
    return 1; }

  va_start(args, format);
  vsprintf(t, format, args);  /* Put new col(s) in t */
  va_end(args);

  newwidth = strlen(t)+1;   /* Width of new column including space */
  if (newwidth < width)
    newwidth = width;

  sprintf(tmpfmt, " %%%ds", newwidth-1);
  
  strncpy(s, entry->line[lineno-1], ptr);    /* Copy colno-1 columns */

  sprintf(&s[ptr], tmpfmt, t);   /* Insert now col(s) */
  sprintf(&s[ptr+newwidth], "%s", &(entry->line[lineno-1])[ptr+width]);
                                           /* The rest of *line */

  if (newwidth > width)
    entry->line[lineno-1] = (char *)realloc(entry->line[lineno-1],
         (strlen(s) + 1) * sizeof(char));
                         /* Reallocate line due to greater length */

  strcpy(entry->line[lineno-1], s);

  free(s);

  return 0;
}

/* This function removes a field */
int RmvField(Entry *entry, int lineno, int colno)
{
  int i;              /* Iterator */
  int ptr;            /* Pointer in line */
  int width;          /* Width of column */
  char *s;

  s = (char *)malloc((strlen(entry->line[lineno-1])+1+MAXCOLW) * sizeof(char));

  ptr = 0;

  for (i = 1; i < colno; i++) {
    if (sscanf(&entry->line[lineno-1][ptr], " %s%n", s, &width) != 1)
      return 1;
    ptr += width;
  }

  if (sscanf(&entry->line[lineno-1][ptr], " %s%n", s, &width) != 1)
    return 1;

  strncpy(s, entry->line[lineno-1], ptr);    /* Copy colno-1 columns */
  sprintf(&s[ptr], "%s", &(entry->line[lineno-1])[ptr+width]);
                                           /* The rest of *line */

  entry->line[lineno-1] = (char *)realloc(entry->line[lineno-1],
         (strlen(s) + 1) * sizeof(char));
                              /* Reallocate line due to smaller length */

  strcpy(entry->line[lineno-1], s);

  free(s);

  return 0;
}



void __ECAddField(Entry *entry, int lineno, void *init);

/*
   This function ensures the presence of a column with the given
   keyword as name. If this column is present, the function returns
   its number. Otherwise the column is added, with init as value in
   all positions.
*/
int EnsureCol(Entry *entry, const char *keyword, const char *init)
{
  int colno;

  if ((colno = ReadColnoND(entry, keyword)) != 0)
    return colno;   /* The column is present */

  EntryMap(entry, &__ECAddField, (void *)init);

  if (entry->len == 0)
    colno = 0;
  else
    colno = CountField(entry, 1);

  AddEntryExt(entry, "; COL %-3d           %s\n", colno, keyword);

  return colno;
}

/* This is used to add columns by EnsureCol */
void __ECAddField(Entry *entry, int lineno, void *init)
{
  AddField(entry, lineno, "%s", (char *)init);
}

/*
   This removes the column with the given label. Returns 0 on error,
   column number otherwise
*/
int RmvCol(Entry *entry, const char *keyword)
{
  int i;
  int colno;

  colno = RmvColno(entry, keyword);

  if (colno == 0)
    return 0;

  for (i = 1; i <= entry->len; i++)
    RmvField(entry, i, colno);

  return colno;
}

/*
   Returns true if lineno is within the sequence
*/
int InSeq(Entry *entry, int lineno)
{
  if (lineno >= 1 && lineno <= entry->len)
    return 1;

  return 0;
}

/*
   Maps the function mapfunc on all positions (in increasing order) in
   sequence in entry.

   f has three arguments, the first is the entry, the second (int) is
   the current line number. The last is a pointer to other arguments.
*/
void EntryMap(Entry *entry, void (*f)(Entry *, int, void *), void *arg)
{
  int i;

  for (i = 0; i < entry->len; i++)
    (*f)(entry, i+1, arg);
}

int __FPNextLine(Entry *entry, int *line, int search_col);
int __FPPrevLine(Entry *entry, int *line, int search_col);

/*
   Finds the linenumber of the RNA basepair corresponding to
   lineno. align_bp_col is the column number of the align_bp
   column. If this is zero, seq_bp_col is used instead.

   Return value is -1 if no basepair is found, linenumber otherwise.
*/
int FindPair(Entry *entry, int lineno, int align_bp_col, int alignpos_col,
	     int seq_bp_col, int seqpos_col)
{
  char field[MAXCOLW];
  int pair;
  int pos, pos2;
  int line, start_line, end_line;
  int pair_col;
  int search_col;
  int strand;

  if (align_bp_col != 0) {
    /* Use align_bp column */
    pair_col = align_bp_col;
    search_col = alignpos_col;
  }
  else {
    /* Use seq_bp column */
    pair_col = seq_bp_col;
    search_col = seqpos_col;
  }

  GetField(field, entry, lineno, pair_col);
  if (field[0] == '.')
    return 0;
  else
    pair = atoi(field);

  if (pair < 1)
    return 0;
  else if (pair <= entry->len) {
    /* Check if lineno pair is what we are looking for */
    GetField(field, entry, pair, search_col);
    if (field[0] != '.') {
      pos = atoi(field);
      if (pos == pair)
	return pair;
    }
  }

  /* Do a real search */

  /* establish strand */
  strand = 0;
  start_line = 1;
  if ((pos = __FPNextLine(entry, &start_line, search_col)) != 0 &&
      (start_line++, pos2 = __FPNextLine(entry, &start_line, search_col)) != 0) {
    if (pos < pos2)
      strand = 1;
    else
      strand = -1;
  }

  start_line = 1;
  end_line = entry->len;

  if ((pos = __FPNextLine(entry, &start_line, search_col)) == 0)
    /* No lines */
    return 0;

  if ((strand == -1 && pos < pair) ||
      (strand != -1 && pair < pos))
    /* pair is before start of entry */
    return 0;
  else if (pair == pos)
    /* This is it */
    return start_line;

  pos = __FPPrevLine(entry, &end_line, search_col);
  
  if ((strand == -1 && pos > pair) ||
      (strand != -1 && pair > pos))
    /* pair is after end of entry */
    return 0;
  else if (pair == pos)
    /* This is it */
    return end_line;

  for (;;) {
    if (end_line-start_line == 1)
      return 0;

    line = (start_line + end_line) / 2;
    pos = __FPPrevLine(entry, &line, search_col);
    
    if ((strand == -1 && pos < pair) ||
	(strand != -1 && pair < pos)) {
      end_line = line;
      continue;
    }
    else if (pair == pos)
      /* This is it */
      return line;

    line = (start_line + end_line) / 2;
    pos = __FPNextLine(entry, &line, search_col);
    
    if ((strand == -1 && pos > pair) ||
	(strand != -1 && pair > pos)) {
      start_line = line;
      continue;
    }
    else if (pair == pos)
      /* This is it */
      return line;

    /* Position does not exist */
    return 0;
  }
}

/*
   Returns value in search_col of next line where search_col is not
   '.', while *line is changed to this line number.  If search_col is
   not '.' in line, then the value is returned and *line is unchanged.
   If there is no more lines, where search_col is not '.', 0 is
   returned.
*/
int __FPNextLine(Entry *entry, int *line, int search_col)
{
  char field[MAXCOLW];

  for (; *line <= entry->len; (*line)++) {
    GetField(field, entry, *line, search_col);
    if (field[0] != '.')
      return atoi(field);
  }

  return 0;
}

int __FPPrevLine(Entry *entry, int *line, int search_col)
{
  char field[MAXCOLW];

  for (; *line >= 1; (*line)--) {
    GetField(field, entry, *line, search_col);
    if (field[0] != '.')
      return atoi(field);
  }

  return 0;
}

/*
   Sets the pair of lineno to the position in pairline. If pairline is
   0, the pair is set to nothing. align_bp_col and seq_bp_col are set
   as appropriate.
*/
void SetPair(Entry *entry, int lineno, int pairline,
	     int align_bp_col, int alignpos_col,
	     int seq_bp_col, int seqpos_col)
{
  char field[MAXCOLW];

  if (pairline == 0) {
    if (align_bp_col != 0)
      ChgField(entry, lineno, align_bp_col, ".");
    if (seq_bp_col != 0)
      ChgField(entry, lineno, seq_bp_col, ".");
  }
  else {
    if (align_bp_col != 0 && alignpos_col != 0) {
      GetField(field, entry, pairline, alignpos_col);
      ChgField(entry, lineno, align_bp_col, "%s", field);
    }
    if (seq_bp_col != 0 && seqpos_col != 0) {
      GetField(field, entry, pairline, seqpos_col);
      ChgField(entry, lineno, seq_bp_col, "%s", field);
    }
  }
}
