/* Generalized single linkage clustering.
 * 
 * SRE, Mon Jan  7 09:40:06 2008 [Janelia]
 * SVN $Id: esl_cluster.h 664 2011-02-27 17:08:36Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/easel/branches/infernal/1.1/esl_cluster.h $
 */
#ifndef eslCLUSTER_INCLUDED
#define eslCLUSTER_INCLUDED

extern int esl_cluster_SingleLinkage(void *base, size_t n, size_t size, 
				     int (*linkfunc)(const void *, const void *, const void *, int *), void *param,
				     int *workspace, int *assignments, int *ret_C);
#endif /*eslCLUSTER_INCLUDED*/
/*****************************************************************
 * Easel - a library of C functions for biological sequence analysis
 * Version i1.1rc2; December 2012
 * Copyright (C) 2012 HHMI Janelia Farm Research Campus
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Easel is distributed under the Janelia Farm Software License, a BSD
 * license. See the LICENSE file for more details.
 *****************************************************************/
