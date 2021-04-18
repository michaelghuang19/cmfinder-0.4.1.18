/********************************************************************************************************
 * QRNA - Comparative analysis of biological sequences 
 *         with pair hidden Markov models, pair stochastic context-free
 *        grammars, and probabilistic evolutionary  models.
 *       
 * Version 2.0.0 (JUN 2003)
 *
 * Copyright (C) 2000-2003 Howard Hughes Medical Institute/Washington University School of Medicine
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ***********************************************************************************************************/

/* globals.h
 * Global variable definitions. 
 * This file may only be included in a main() .c file.
 *
 * ER, Fri May 28 11:33:45 CDT 1999 [St. Louis]
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "squid.h"

#define MALLOC_CHECK_ 1

/* time
 */
#define TMAX 3.0
#define TMIN 0.00005

/* accuracy
 */
#define accuracy1 0.9999999
#define accuracy  0.99
#define MARGIN    0.001
#define MARGIN2   0.01
		
/* Useful macros 
 */
#define PI           3.14159265359
#define LN2          0.69314718056
#define LN2INV       1.44269504089
#define LOG2(x)      ((x) == 0.0 ? -9999. : log(x) * LN2INV)
#define EXP2(x)      (exp((x) * LN2 )) 
#define INTSCALE     10.0
#define BIGINT       9999999	           /* prohibition in viterbi.c without overflow              */
#define BIGFLOAT     99999.999	           /* prohibition in viterbi.c without overflow              */
#define MAXRNALOOP   200                   /* maximum loop size allowed == d == (MAXRNALOOP-1)       */
#define PROHIBIT     50.0
#define GASCONST     8.31441               /* boltzman factor (J mol^{-1} K^{-1})                    */
#define GASCAL       (GASCONST*0.000239)   /* in Kcal: 1 J= 0.239 cal                                */
#define BETAINV      (GASCAL*310.15)       /* 1/beta = kT [Kcal/mol) (T=310 K)                       */
#define BETA         (1.0/BETAINV)         /* beta = 1/kT (T=310 K)                                  */
#define MAX_NUM_ENDS 50                    /* maximum number of ends reported                        */


#define CODON(x1,x2,x3)   ((x1)*16 + (x2)*4 + (x3))  /* Lookup index of a codon with bases x1x2x3 */
#define CODON5(x1,x2,x3)  ((x1)*25 + (x2)*5 + (x3))  /* Lookup index of a codon with gaps         */
#define idx(x,y)          ((x)*4+(y))                /* Lookup index of a mut prob xy             */
#define idx5(x,y)         ((x)*5+(y))                /* Lookup index of a mut5 prob xy            */

#define SYMIDX(x)    (strchr(Alphabet, (x)) - Alphabet)
#define DNAIDX(x)    (strchr(DNAAlphabet, (x)) - DNAAlphabet)

#define sreLOG2(x)   ((x) > 0 ? log(x) * 1.44269504 : -9999.)
#define sreEXP2(x)   (exp((x) * 0.69314718 )) 


#define sre_isascii(c)   (!((c) & ~0177))


extern char Alphabet[];
extern char DNAAlphabet[];
extern char SSAlphabet[];

#define RNA_ALPHABET_GAP "ACGU-"

/* additions from RJ Klein rsearch code 22 JAN 03
 */
#define RNA_ALPHABET_SIZE Alphabet_size
#define RNAPAIR_ALPHABET "AAAACCCCGGGGUUUU"
#define RNAPAIR_ALPHABET2 "ACGUACGUACGUACGU"

/* Returns true if pos. C of seq B of msa A is a gap as defined by isgap(c) 
   from squid */
#define is_rna_gap(A, B, C) (isgap(A->aseq[B][C]))

/* Returns true if pos. C of seq B of msa A is an uppercase A, C, G, T, or U */
#define is_defined_rna_nucleotide(A, B, C) (A->aseq[B][C] == 'A' || A->aseq[B][C] == 'C' || A->aseq[B][C] == 'G' || A->aseq[B][C] == 'T' || A->aseq[B][C] == 'U')
#define is_rna_nucleotide(X) (X == 'A' || X == 'C' || X == 'G' || X == 'T' || X == 'U')

/*
 * Maps to index of matrix, using binary representation of
 * nucleotides (unsorted).
 *
 * See lab book 7, p. 3-4 for details of mapping function
 */
#define matrix_index(X,Y) ((X>Y) ? X*(X+1)/2+Y: Y*(Y+1)/2+X)


extern int   Alphabet_type;
extern int   Alphabet_size;
extern int   Alphabet_iupac;











