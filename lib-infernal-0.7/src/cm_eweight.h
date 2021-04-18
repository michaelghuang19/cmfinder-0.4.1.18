/************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 ************************************************************/

/* cm_eweight.h
 * based on:
 */
/* lsjfuncs.h
 * Declarations of external functions used in lsj_eweight.c
 * (Entropy-based sequence weighting)
 *
 * Steve Johnson
 * SVN $Id: cm_eweight.h,v 1.1 2006/03/07 07:39:51 yzizhen Exp $
 */


#include "config.h"
#include "structs.h"
#include "squid.h"
#include "msa.h"


/*
extern float CM_Eweight(struct CM_t *cm,  struct Prior_t *pri, 
		     float numb_seqs, float targetent);
*/

extern double CM_Eweight(CM_t *cm,  Prior_t *pri, 
		     float numb_seqs, float targetent);

extern void ModelContent(float *ent1, float *ent2, int M);
extern void CMRescale(CM_t *hmm, float scale);
extern double CM_Eweight_RE(CM_t *cm, Prior_t *pri, float numb_seqs, 
			    float target_relent, float *randomseq);
extern double DRelEntropy(double *p, double *f, int n);
