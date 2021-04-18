/*********************************************************************

  file.c

  Contains functions for handling files.

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

#include "llist.h"
#include "file.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* stack for reading and putting connected lines */
LList *_Con_stack = NULL;

/* 
   This function reads a line from fp. The line can have any length
   and a pointer to it is returned. NULL is returned on error or EOF.
*/
char *GetLine(FILE *fp)
{
  int maxlen;    /* Dynamic allocation length */
  char *s;       /* The string to hold the read line */
  int ptr;       /* Pointer to string */
  int c;         /* Read character */

  maxlen = 1024;  /* Start out with 1k */

  s = (char *)malloc((maxlen) * sizeof(char));
  ptr = 0;
  while ((c = fgetc(fp)) != EOF) {
    if (ptr == maxlen) {
      maxlen *= 2;
      s = (char *)realloc(s, maxlen * sizeof(char));  /* Get more space */
    }
    s[ptr++] = c;
    if (c == '\n')
      break;
  }

  if (ptr == 0)
    return NULL;

  s = (char *)realloc(s, (ptr+1) * sizeof(char));
                                      /* Remove excess (in most cases) */
  s[ptr] = '\0';   /* End of string */

  return s;
}

/* 
   This function reads a whole file from fp. The file can have any
   length and a pointer to it is returned. NULL is returned on error.
*/
char *GetFile(FILE *fp)
{
  int maxlen;    /* Dynamic allocation length */
  char *s;       /* The string to hold the read line */
  int ptr;       /* Pointer to string */
  int c;         /* Read character */

  maxlen = 1024;  /* Start out with 1k */

  s = (char *)malloc((maxlen) * sizeof(char));

  ptr = 0;
  while ((c = fgetc(fp)) != EOF) {
    if (ptr == maxlen) {
      maxlen *= 2;
      s = (char *)realloc(s, maxlen * sizeof(char));  /* Get more space */
    }
    s[ptr++] = c;
  }

  if (ptr == 0)
    return NULL;

  s = (char *)realloc(s, (ptr+1) * sizeof(char));
                                      /* Remove excess (in most cases) */
  s[ptr] = '\0';   /* End of string */

  return s;
}

/*
   This functions initiates setup for GetArgument and
   GetFilename. Should be used before these functions.
*/
CmdArg *InitArgument(int argc, char **argv)
{
  int i;
  CmdArg *cmdarg;

  cmdarg = (CmdArg *)malloc(sizeof (CmdArg));

  cmdarg->argc = argc;
  cmdarg->argv = argv;
  cmdarg->argc_ptr = 1;
  cmdarg->argv_ptr = 0;

  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-')
      break;
    if (strcmp(argv[i], "--") == 0) {
      i++;
      break;}
  }
  cmdarg->file_ptr = i;

  return cmdarg;
}  

/*
   This function returns a pointer to the next command line
   argument. NULL if no more arguments.

   The command line '--rec -wl45,23n34-23 -cda -- -s' gives in turn:

   -rec
   w
   l45,23
   n34-23
   c
   d
   a

   The -- stops the arguments, so '-s' is a filename.
*/
char *GetArgument(CmdArg *cmdarg)
{
  int pos;
  char *argument;
  int maxlen;        /* for dynamic memory allocation */

  if (cmdarg->argc_ptr == cmdarg->argc ||
      (cmdarg->argv_ptr == 0 &&
       (cmdarg->argv[cmdarg->argc_ptr][cmdarg->argv_ptr] != '-' ||
	strcmp(cmdarg->argv[cmdarg->argc_ptr], "--") == 0)))
    return NULL;

  if (cmdarg->argv_ptr == 0) {
    if (strncmp(cmdarg->argv[cmdarg->argc_ptr], "--", 2) == 0)
      return &cmdarg->argv[cmdarg->argc_ptr++][1];
    else
      cmdarg->argv_ptr++;                 /* move over '-' */
  }

  maxlen = 20;                  /* current maximum length of argument */

  argument = (char *)malloc(maxlen * sizeof(char));
  
  pos = 0;

  do {
    if (pos == maxlen) {
      maxlen *= 2;
      argument = (char *)realloc(argument, maxlen * sizeof(char));
                                                  /* Get more space */
    }
    argument[pos++] = cmdarg->argv[cmdarg->argc_ptr][cmdarg->argv_ptr++];
  } while (!isalpha(cmdarg->argv[cmdarg->argc_ptr][cmdarg->argv_ptr]) &&
	   cmdarg->argv[cmdarg->argc_ptr][cmdarg->argv_ptr] != '\0');

  if (cmdarg->argv[cmdarg->argc_ptr][cmdarg->argv_ptr] == '\0') {
    cmdarg->argv_ptr = 0;
    cmdarg->argc_ptr++; }

  argument = (char *)realloc(argument, (pos+1) * sizeof(char));
                                                /* adjust length */

  argument[pos] = '\0';

  return argument;
}

/* as GetArgument, gets the command line filenames */
char *GetFilename(CmdArg *cmdarg)
{
  if (cmdarg->file_ptr == cmdarg->argc)
    return NULL;

  return cmdarg->argv[cmdarg->file_ptr++];
}

/* Sets up stack for PutConnect and GetConnect */
void InitConnect(void)
{
  if (_Con_stack != NULL)
    DestroyLList(_Con_stack);

  _Con_stack = MakeLList();
}

/* Put a read line back to be read once more */
void PutConnect(char *s)
{
  Push(_Con_stack, (void *)s);
}

/*
   Reads a line from fp. A line can span several real lines that are
   ended by a '\'.
*/
char *GetConnect(FILE *fp)
{
  char *s;    /* For individual lines */
  char *line; /* For line connected with '\' */
  int ptr;    /* Pointer in line */

  if (_Con_stack->count != 0)
    /* A line has been put using PutConnect */
    return (char *)Pop(_Con_stack);

  line = NULL;
  ptr = 0;

  while ((s = GetLine(fp)) != NULL) {
    if (line == NULL)
      line = s;
    else {
      line = (char *)realloc(line, (ptr + strlen(s) + 1) * sizeof(char));
      strcpy(line+ptr, s);
    }

    for (ptr = strlen(line)-1; isspace(line[ptr]) && ptr > 0; ptr--)
      ;
    if (line[ptr] == '\\' || isspace(line[ptr]))
      ;
    else
      break;
  }

  if (s != NULL && isspace(line[ptr]))
    return NULL;

  return line;
}

int StrCmp(const char *s, const char *t)
{
  int i;

  for (i = 0; tolower(s[i]) == tolower(t[i]); i++)
    if (s[i] == '\0')
      return 0;

  return tolower(s[i])-tolower(t[i]);
}

int StrnCmp(const char *s, const char *t, int n)
{
  int i;

  for (i = 0; tolower(s[i]) == tolower(t[i]) && i < n; i++)
    if (s[i] == '\0')
      return 0;

  if (i == n)
    return 0;

  return tolower(s[i])-tolower(t[i]);
}

