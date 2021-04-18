/*********************************************************************

  file.h

  Contains functions for handling files.

*********************************************************************/

#ifndef __file_h__
#define __file_h__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "llist.h"

typedef struct tagCmdArg {
  int argc;
  char **argv;
  int argc_ptr;  /* The present argument */
  int argv_ptr;  /* The present position in this argument */
  int file_ptr;  /* Pointer to argument that is afilename */
} CmdArg;


char *GetLine(FILE *fp);
char *GetFile(FILE *fp);

CmdArg *InitArgument(int argc, char **argv);
char *GetArgument(CmdArg *cmdarg);
char *GetFilename(CmdArg *cmdarg);

void InitConnect(void);
void PutConnect(char *s);
char *GetConnect(FILE *fp);

int StrCmp(const char *s, const char *t);
int StrnCmp(const char *s, const char *t, int n);

#endif
