/*********************************************************************
                       
  FILE: grammar.h
                                             
  Grammars are defined here, and some functions for handling them.  

  990205 Bjarne Knudsen
                                                                    
**********************************************************************/

#ifndef __grammar_h__   /* Only define the following once */
#define __grammar_h__

#define prob_double

#include "llist.h"
#include "file.h"
#include "matrix.h"
#include "edouble.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* a grammar contains linked lists of nontermianls, terminals, single
   groups, double groups and various types of rules */
typedef struct tagGrammar{
  LList *nont;    /* LList of Nont */
  LList *term;    /* LList of Term */
  LList *iprt;    /* LList of interpretations */
  LList *trns;    /* LList of translations */
  LList *sgrp;    /* LList of Sgrp */
  LList *dgrp;    /* LList of Dgrp */
  LList *ruls;    /* LList of Ruls */
  LList *ruln;    /* LList of Ruln */
  LList *rulsn;   /* LList of Rulsn */
  LList *rulns;   /* LList of Rulns */
  LList *ruldnd;  /* LList of Ruldnd */
  LList *ruldd;   /* LList of Ruldd */
  LList *rulnn;   /* LList of Rulnn */
  Edouble **quickdist;
  int mindist;    /* Minimum distance for pairs */
} Grammar;

/* a nonterminal is a symbol */
typedef struct tagNont{char sym;} Nont;

/* a terminal is a symbol */
typedef struct tagTerm{char sym;} Term;

/* an interpretation is a symbol and a pointer to a list of
   probabilities */
typedef struct tagIprt{char sym; Matrix *freq;} Iprt;

/* a translation is a symbol and a its translation symbol */
typedef struct tagTrns{char sym; char newsym;} Trns;

/* an Sgrp is a symbol and a pointer to a matrix (vector) of terminal
   probabilities. Three matrices for evolution is also included. */
typedef struct tagSgrp
{
  char sym;
  Matrix *freq;
  Matrix *eigen;
  Matrix *diag;
  Matrix *inveigen;
} Sgrp;

/* a Dgrp is like Sgrp, except that numterm*numterm size matrices are
   used */
typedef struct tagDgrp
{
  char sym;
  Matrix *freq;
  Matrix *eigen;
  Matrix *diag;
  Matrix *inveigen;
} Dgrp;

/* This header file defines rules as structures. lnt is the number of
   the left nonterminal. The integers starting with 'r' are num- bers
   of produced symbols.
                                                                    
   struct    format                                                 
   ruls      N -> s                                                 
   ruln      N -> N                                                 
   rulsn     N -> sN                                                
   rulns     N -> Ns                                                
   ruldnd    N -> dNd (a)                                           
   ruldd     N -> dd (a)                                           
   rulnn     N -> NN                                                
                                                                    
  (a) the two d's are the same group, the number of which is in
  rdb. */
typedef struct tagRuls
{int lnt; int rsg; Edouble prob; Fdouble fprob;} Ruls;
typedef struct tagRuln
{int lnt; int rnt; Edouble prob; Fdouble fprob;} Ruln;
typedef struct tagRulsn
{int lnt; int rsg; int rnt; Edouble prob; Fdouble fprob;} Rulsn;
typedef struct tagRulns
{int lnt; int rsg; int rnt; Edouble prob; Fdouble fprob;} Rulns;
typedef struct tagRuldnd
{int lnt; int rdb; int rnt; Edouble prob; Fdouble fprob;} Ruldnd;
typedef struct tagRuldd
{int lnt; int rdb; Edouble prob; Fdouble fprob;} Ruldd;
typedef struct tagRulnn
{int lnt; int rnt1; int rnt2; Edouble prob; Fdouble fprob;} Rulnn;

/* prototypes */
Grammar *MakeGrammar(void);
            /* allocates and initialises Grammar */
/*Grammar *InputGrammar(char *file);*/
            /* allocates Grammar and reads it from 'file' */
Grammar *ReadGrammar(FILE *fp);
            /* reads Grammar from '*fp' and returns it */
void PrintGrammar(FILE *fp, Grammar *grammar);
            /* prints grammar to *fp */
void RobustGrammar(Grammar *grammar, double robust);
            /* makes a grammar robust */
int FindNont(char sym, LList *nont);
            /* returns number of nonterminal, -1 on error */
int FindTerm(char sym, LList *term);
            /* returns number of terminal, -1 on error */
int FindIprt(char sym, LList *iprt);
            /* returns number of interpret symbol, -1 on error */
int FindTrns(char sym, LList *trns);
            /* returns number of translate symbol, -1 on error */
int FindSgrp(char sym, LList *sgrp);
            /* returns number of single group, -1 on error */
int FindDgrp(char sym, LList *dgrp);
            /* returns number of double group, -1 on error */
int FindSym(char sym, Grammar *grammar);
            /* returns number of symbol (terminal or iprt), -1 on error */
char Translate(char sym, LList *trns);
            /* translate symbol according to trns */


/*Grammar* SingleGrammar(Grammar *grammar);*/
            /* Makes grammar without pair rules */
#endif
