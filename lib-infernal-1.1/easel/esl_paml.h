/* PAML interface
 *
 *   "Phylogenetic Analysis by Maximum Likelihood"
 *   Ziheng Yang
 *   http://abacus.gene.ucl.ac.uk/software/paml.html
 *   [Yang97]
 * 
 *           incept: SRE, Tue Jul 13 13:20:08 2004 [St. Louis]
 * upgrade to Easel: SRE, Thu Mar  8 13:26:20 2007 [Janelia]
 * SVN $Id: esl_paml.h 664 2011-02-27 17:08:36Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/infernal/1.1/esl_paml.h $
 */
#ifndef eslPAML_INCLUDED
#define eslPAML_INCLUDED

#include <stdio.h>
#include <esl_dmatrix.h>

extern int esl_paml_ReadE(FILE *fp, ESL_DMATRIX *E, double *pi);


#endif /*eslPAML_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version i1.1rc2; December 2012
 * Copyright (C) 2012 HHMI Janelia Farm Research Campus
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
