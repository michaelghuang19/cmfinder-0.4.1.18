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


/* globals.c
 *
 *
 * E. Rivas [St. Louis]
 * 
 * 9 april 1999.
 */
                                            
#include <stdio.h>

#include "mat_globals.h"
#include "squid.h"
#include "mat_structs.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

char Alphabet[]    = AMINO_ALPHABET;
char DNAAlphabet[] = DNA_ALPHABET;
char SSAlphabet[]  = "123456789abcdefghijklmnopqrstuvwxyz";

/* additions from RJ Klein rsearch code 22 JAN 03
 */
int   Alphabet_type  = kRNA;
int   Alphabet_size  = 4;
int   Alphabet_iupac = 17;

