/* i/o of multiple sequence alignment files in dotless UCSC A2M format
 */
#ifndef eslMSAFILE_A2M_INCLUDED
#define eslMSAFILE_A2M_INCLUDED

#include "esl_msa.h"
#include "esl_msafile.h"

extern int esl_msafile_a2m_SetInmap     (ESLX_MSAFILE *afp);
extern int esl_msafile_a2m_GuessAlphabet(ESLX_MSAFILE *afp, int *ret_type);
extern int esl_msafile_a2m_Read         (ESLX_MSAFILE *afp, ESL_MSA **ret_msa);
extern int esl_msafile_a2m_Write        (FILE *fp,    const ESL_MSA *msa);

#endif /* eslMSAFILE_A2M_INCLUDED */

/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version i1.1rc2; December 2012
 * Copyright (C) 2012 HHMI Janelia Farm Research Campus
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *
 * SVN $Id: esl_msafile_a2m.h 708 2011-07-20 12:49:10Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/infernal/1.1/esl_msafile_a2m.h $
 *****************************************************************/
