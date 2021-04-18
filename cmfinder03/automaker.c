/* automaker.c
* Construct a covariance model from an alignment.
* This code is adopted from COVE "fastmodelmaker.c"
* We have added structural prior in addition to mutual information to improve
* secondary structure prediction.
* In the dynamic programming algorithmg to identify a set of base pairs that optimize the
* overall probability, we add a gapcost for each helix to avoid isolated base pair.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "squid.h"
#include "structs.h"
#include "funcs.h"
#include "sre_stack.h"

#include "part_func.h"
#include "utils.h"
#include "fold_vars.h"
#include "global.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

#define MIN_HAIRPIN  3
#define PMAT         1
#define ZMAT         2
#define MINSCORE     -400
#define MAXSCORE     100
#define MUTUAL       'm'
#define ENERGY       'e'
#define COMB         'b'
#define NEGINFINITY  -9999999
#define ALPHABETSIZE 4


void mixy(char **aseq, int nseq, int alen, float* weights, int ***ret_mixy);
void merge(int **mxy, int **pxy, int alen, int ***ret_xy);



#ifdef DEBUG
static void dump_mat(int **zmat, int alen);
#endif /* DEBUG */

extern void         cm_from_guide(CM_t *cm, Parsetree_t *gtr);


/* Function: mixy()
*
* Originally from COVE package.
*
* Purpose:  given a set of N aligned sequences aseq, calculate
*           pairwise covariances (mutual information). ret_mixy
*           is allocated, filled, and returned, as a diagonal 2D
*           (NxN) matrix of values. It must be freed by
*           the caller. It is a lower diagonal matrix mxy[j][i],
*           j > i, 0..alen-1 by 0..j-1.
*
*           The values in mxy are integers. They are the average
*           secondary structure information content (i.e. weighted for
*           the number of pairs actually occurring in columns i,j)
*           in bits, to two decimal places (i.e. info*100).
*
* Returns:  mxy, which must be free'd by caller with free_mixy().
*/
void
mixy(char    **dsq,            /* array of aligned sequences, flushed right  */
	 int       nseq,		/* number of aligned sequences */
	 int       alen,		/* length of each sequence (all the same)  */
	 float*   weights,
	 int    ***ret_mxy)
{
	int    **mxy;                 /* RETURN: diagonal covariance matrix  */
	double   fx[ALPHABETSIZE];        /* singlet frequency vector        */
	double   fy[ALPHABETSIZE];	/* another singlet frequency vector    */
	double   fxy[ALPHABETSIZE][ALPHABETSIZE]; /* pairwise frequency 2D array  */
	int     idx;			/* counter for sequences               */
	int     i, j;			/* counters for columns x,y            */
	int     symi, symj;		/* counters for symbols                */
	float   pairs;		/* counter for pairs in which there are no gaps */


	mxy = IntAlloc2DArray(alen + 1);

	/* calculate mxy
	*/
	for (j = 2; j <= alen; j++){
		for (i = 1; i < j; i++)
		{
			/* zero counter array */
			for (symj = 0; symj < Alphabet_size; symj++)
			{
				fx[symj] = fy[symj] = 0.0;
				for (symi = 0; symi < Alphabet_size; symi++)
					fxy[symj][symi] = 0.0;
			}
			/* count symbols in a column */
			pairs = 0;
			for (idx = 0; idx < nseq; idx++)
			{
				/* Gaps are ignored */
				if (dsq[idx][i] == DIGITAL_GAP || dsq[idx][j] == DIGITAL_GAP )
					continue;
				symi = dsq[idx][i];
				symj = dsq[idx][j];
				if(symi >=4 || symj >=4){
					Die("Invalid symbols of seq %d at %d %d or %d %d", idx, i, symi, j, symj);
				}
				fx[symi] += weights[idx];
				fy[symj] += weights[idx];
				fxy[symi][symj] += weights[idx];
				pairs+= weights[idx];
			}

			/* convert to frequencies */
			if (pairs > 0)
				for (symi = 0; symi < Alphabet_size; symi++)
				{
					fx[symi] /=  pairs;
					fy[symi] /=  pairs;
					for (symj = 0; symj < Alphabet_size; symj++)
						fxy[symi][symj] /=  pairs;
				}
			else{
				mxy[j][i] =  NEGINFINITY;
				continue;
			}
			/* calculate mxy. 144.269504 is a conversion of ln's into
			* bits * 100: i.e. 100 * (1/log(2))
			*/
			mxy[j][i] = 0;
			for (symi = 0; symi < Alphabet_size; symi++)
				for (symj = 0; symj < Alphabet_size; symj++)
				{
					if (fxy[symi][symj] > 0.0)
						mxy[j][i] += (int) (144.269504 * fxy[symi][symj] *
						log((fxy[symi][symj] / (fx[symi] * fy[symj]))));
				}

				/* Sat Jul 17 22:17:17 1993:  We weight by pairs to get an expected score
				* over all the sequences. Fixes a problem that columns with few symbols
				* could dominate the calculation just because of noise.
				*/
				mxy[j][i] =  (int)(((float)(mxy[j][i]) * pairs) / (float)(nseq));
		}
	}

	/* dump debugging info
	*/

	*ret_mxy = mxy;
}








/*
Calculate base-pairing probabilities according to the average partition function
(which is calculated by the Avg_bppr function).
Then calculate LOD scores for favoring base pairing, based on the partition function
evidence.  Converted to the same scale as the mixy function (above).
*/
void
prxy(char    **aseq,            /* array of aligned sequences, flushed right  */
	 int       nseq,		/* number of aligned sequences */
	 int       alen,		/* length of each sequence (all the same)  */
	 double  **bp_pr,
	 float*    weights,    
	 int    ***ret_prxy)        /* RETURN: prxy array           */
{
	int    i,j;
	int    **lod;
	double *px;
	double **pxy;
	Avg_bppr(aseq,nseq, alen, bp_pr, weights, &pxy);
	px = (double*) malloc(sizeof(double) * (alen + 1));

	lod = IntAlloc2DArray(alen+1);
	for(i=0; i < alen; i++) {
		px[i+1]=0;
		for(j = 0; j < alen; j++){
			if (i != j){
				if (i < j)
					px[i+1] += pxy[j+1][i+1];
				else
					px[i+1] += pxy[i+1][j+1];
			}
		}
	}

	for(j = 1; j < alen; j++)
		for(i = 0; i < j; i++){
			lod[j+1][i+1] = (int)(0.5 * 144.269504 * log(pxy[j+1][i+1] /( (1-px[i+1]) * (1- px[j+1]))));
			if (lod[j+1][i+1] < MINSCORE) lod[j+1][i+1] = MINSCORE;
		}

		free(px);
		Free2DArray((void **)pxy, alen+1);
		*ret_prxy = lod;
}


/* combine mutual information and the prior using partition function */
void merge(int **mxy, int **pxy, int alen, int*** ret_xy)
{
	int **xy;
	int i,j;

	xy = IntAlloc2DArray(alen+1);

	for(j=2; j <= alen; j++)
		for(i=1; i < j; i++){
			xy[j][i] = mxy[j][i] + pxy[j][i] ;
			if (pxy[j][i] < 0 )
				xy[j][i] = 0;

		}
		*ret_xy = xy;
}


/* Function: zfill()
*
* Purpose:  Calculate the optimal structure for a covariance matrix
*           produced by mixy(). Uses a way-simplified form of the
*           Zuker/Nussinov dynamic programming RNA folding algorithm
*           to find the structure which a) the emitted pairs sum
*           to a maximum number of bits of covariance and b)
*           has no overlapping chords (no pseudoknots). The dynamic
*           programming matrix is allocated, filled, and returned.
*
* Returns:  ret_zmat is returned thru a passed pointer; it must be
*           free'd by the caller using free_zmat().
*/
/* Function: ztrace()
*
* Purpose:  Traceback through the dynamic programming matrix constructed
*           by zfill(). Constructs a dynamic binary tree (ztr) of z Parsetree_t
*           structures, which keep track of both topology and the
*           order in which various aligned columns are emitted.
*
*           ztr ends up being the "shell" or template upon which the
*           final model is built. ztr captures the branching structure
*           of the model tree.
*
*           Inserts are dealt with at this point. Columns with a gap frequency
*           exceeding gapthresh are excluded from ztr. The final tree ztr contains
*           MATR, MATL, MATP nodes only (with BEGIN and BIFURC of course).
*
* Data:     ztr: a traceback tree.
*                emitl = index of column emitted left (0..acol-1)
*                        or -1
*                emitr = index of column emitted right (0..acol-1)
*                        or -1
*                nodeidx = index of node in new model
*                type   = type of node (MATP_NODE, etc.)
*
* Return:   ret_ztr is allocated here and must be free'd by the caller.
*/

CM_t *
Automodelmaker(MSA *msa, char** dsq, int ** pxy, double gapthresh, int gapcost, int trim, int* constraint)
{
	int        **mxy;
	int        **xy;
	int        **pmat;
	int        **zmat;

	CM_t           *cm;		/* new covariance model                       */
	Parsetree_t    *gtr;		/* guide tree for alignment                   */
	Parsetree_t    *tr;		/* individual traces from alignment           */
	Nstack_t       *pda;		/* pushdown stack used in building gtr        */
	int            *matassign;	/* 0..alen-1 array; 0=insert col, 1=match col */
	int             idx;		/* counter over sequences in the alignment    */
	int             v;		/* index of current node                      */
	int             type;		/* type of node we're working on              */
	int  nnodes;			/* number of nodes in CM                      */
	int  nstates;			/* number of states in CM                     */

	int      i, j,li,lj;
	int      diff;
	int      mid;
	int      mat_flag = 0;
	double   sum_weight=0;
	int consensus_bp;

	/* remove existing structure, since we're going to infer the structure for the next iteration */
	if (msa->ss_cons){
		for(i=0; i < msa->alen; i++) msa->ss_cons[i] = '.';
	}

	/* build Mxy matrix, which gives mutual information for all possible pairs */
	mixy(dsq, msa->nseq, msa->alen,  msa->wgt, &mxy);   /* mutual information matrix */

	/* calculate informative priors for all possible pairs, based on Boltzmann distribution */
	if (pxy==NULL){
		prxy(msa->aseq, msa->nseq, msa->alen, NULL, msa->wgt, &pxy);
	}


	merge(mxy, pxy, msa->alen, &xy);

	matassign = MallocOrDie(sizeof(int) * (msa->alen+1));

	sum_weight = 0;
	for (i = 0; i < msa->nseq; i++) 
	{
		sum_weight += msa->wgt[i];
	}  
	for (i = 1; i <= msa->alen; i++)
	{
		double gaps = 0;
		for (idx = 0; idx < msa->nseq; idx++){
			if (dsq[idx][i] == DIGITAL_GAP)
				gaps += msa->wgt[idx];
		}
		matassign[i] = ( gaps / sum_weight > gapthresh) ? 0 : 1;
		if (constraint != NULL){
			if (constraint[i]>=0) {
				matassign[i] = 1;
				j = constraint[i];
				xy[j][i] = MAXSCORE;
			}
		}
	}

	zmat = IntAlloc2DArray(msa->alen+1);
	pmat = IntAlloc2DArray(msa->alen+1);

	/* Initialization.*/
	for(diff = 0; diff <= MIN_HAIRPIN; diff++) {
		for (i = 1; (j = i + diff) < msa->alen; i++){
			zmat[j][i] = 0;
			pmat[j][i] = MINSCORE;
		}
	}

	/* Dynamic programming stage.
	* Our recursion is:
	*    Sij = max { Si+1,j  (emit left, no covariance)
	*                Si,j-1  (emit right, no covariance)
	*                Si+1,j-1 + xy[j][i].
	*                max over mid: Si,mid + Smid+1,j (bifurcation)
	*                }
	*/
	for (diff = MIN_HAIRPIN+1; diff < msa->alen; diff++)
		for (i = 1; (j = i+diff) <= msa->alen; i++)
		{
			/* Ignore gap */
			if(!matassign[j]){
				zmat[j][i] = zmat[j-1][i];
				pmat[j][i] = pmat[j-1][i];
			}
			else if (!matassign[i]){
				zmat[j][i] = zmat[j][i+1];
				pmat[j][i] = pmat[j][i+1];
			}
			else{
				/* update pmat */
				pmat[j][i] = MINSCORE
					;
				//If it is OK to base pair i, j
				if (!constraint || constraint[j]== UNDEF_BP && constraint[i] == UNDEF_BP || constraint[j] == i){
					if (xy[j][i] > 0 || constraint && constraint[j]==i){
						pmat[j][i] =  xy[j][i];
						if (pmat[j-1][i+1] > zmat[j-1][i+1] - gapcost)
							pmat[j][i] += pmat[j-1][i+1];
						else
							pmat[j][i] += zmat[j-1][i+1] - gapcost;
					}
				}
				//If it is OK not to base pair i,j
				zmat[j][i] = MINSCORE;
				if (!constraint || constraint[j] != i){
					zmat[j][i] = zmat[j-1][i];	  
					if (zmat[j][i+1] > zmat[j][i])
						zmat[j][i] = zmat[j][i+1];
					for (mid = i+1; mid < j-1; mid++)
						if (zmat[mid][i] + zmat[j][mid+1] > zmat[j][i])
							zmat[j][i] = zmat[mid][i] + zmat[j][mid+1];
				}
				//If it is OK to base pair i,j
				if (!constraint || constraint[j]== UNDEF_BP && constraint[i] == UNDEF_BP || constraint[j] == i){
					if (constraint && constraint[j] == i || xy[j][i] > 0 && pmat[j-1][i+1] + xy[j][i] > zmat[j][i]){
						zmat[j][i] = pmat[j-1][i+1] + xy[j][i];
					}
				}
			}
		}
		//printf("zmat\n");
		//dump_mat(zmat, msa->alen+1);

		/* 4. Construct a guide tree.
		*    This code is borrowed from yarn's KHS2Trace().
		*
		*    We also keep track of how many states we'll need in the final CM,
		*    so we'll know how much to allocate -- and the number of nodes,
		*    for informational purposes.
		*/
		nstates = nnodes = 0;
		gtr = CreateParsetree();	/* the parse tree we'll grow        */
		pda = CreateNstack();		/* a pushdown stack for our indices */

		i = 1;
		j = msa->alen;
		if (trim && msa->rf){
			for(i=1; i < msa->alen && !matassign[i]; i++);
			for(j=msa->alen; j > 1 && !matassign[i]; j--);
		}

		PushNstack(pda, -1);		/* what node it's attached to */
		PushNstack(pda, i);		/* emitl */
		PushNstack(pda, j);	/* emitr */
		PushNstack(pda, ROOT_nd);	/* "state" (e.g. node type) */

		mat_flag = ZMAT;

		consensus_bp =  0;
		while (PopNstack(pda, &type))	/* pop a node type to attach */
		{
			PopNstack(pda, &j);
			PopNstack(pda, &i);	/* i..j == subseq we're responsible for */
			PopNstack(pda, &v);	/* v = index of parent node in gtr */
			for (li=i+1; li <= j; li++) if (matassign[li]) break;
			for (lj=j-1; lj >= li;lj--) if (matassign[lj]) break;
			if (i > j) {
				v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, END_nd);
				nstates += 1;		/* END_nd -> E_st */
				nnodes++;
			}
			else if (type == ROOT_nd) { /* try to push i,j; but deal with INSL and INSR */
				v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, ROOT_nd);
				for (; i <= j; i++) if (matassign[i]) break;
				for (; j >= i; j--) if (matassign[j]) break;
				PushNstack(pda, v);	/* here v==0 always. */
				PushNstack(pda, i);
				PushNstack(pda, j);
				PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
				nstates += 3;		/* ROOT_nd -> S_st, IL_st, IR_st */
				nnodes++;
			}

			else if (type == BEGL_nd) {    /* no inserts */
				v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, BEGL_nd);
				PushNstack(pda, v);
				PushNstack(pda, i);
				PushNstack(pda, j);
				PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
				nstates += 1;		/* BEGL_nd -> S_st */
				nnodes++;
			}

			else if (type == BEGR_nd)  { /* look for INSL */
				v = InsertTraceNode(gtr, v, TRACE_RIGHT_CHILD, i, j, BEGR_nd);
				for (; i <= j; i++) if (matassign[i]) break;
				PushNstack(pda, v);
				PushNstack(pda, i);
				PushNstack(pda, j);
				PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
				nstates += 2;		/* BEGR_nd -> S_st IL_st */
				nnodes++;
			}
			else if ( j - i <= MIN_HAIRPIN ||  (mat_flag == ZMAT && zmat[j][i] ==  zmat[j][li]))
			{
				/* i unpaired. This is a MATL node; allow INSL */
				v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATL_nd);
				i = li;
				PushNstack(pda, v);
				PushNstack(pda, i);
				PushNstack(pda, j);
				PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
				nstates += 3;		/* MATL_nd -> ML_st, D_st, IL_st */
				nnodes++;

				mat_flag = ZMAT;
			}
			else if ( mat_flag == ZMAT && zmat[j][i] ==  zmat[lj][i]){
				/* j unpaired. This is a MATR node; allow INSR */
				v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATR_nd);
				j = lj;
				PushNstack(pda, v);
				PushNstack(pda, i);
				PushNstack(pda, j);
				PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
				nstates += 3;		/* MATR_nd -> MR_st, D_st, IL_st */
				nnodes++;

				mat_flag = ZMAT;
			}
			else if ( mat_flag == PMAT || (mat_flag == ZMAT && zmat[j][i] == pmat[lj][li] + xy[j][i] && xy[j][i] > 0)) {
				/* i,j paired to each other. MATP. deal with INSL, INSR */
				v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, MATP_nd);
				consensus_bp ++;
				if ( mat_flag == PMAT && pmat[j][i] > pmat[lj][li] + xy[j][i]){
					mat_flag = ZMAT;
				}
				else {
					mat_flag = PMAT;
				}
				if (msa->ss_cons){
					msa->ss_cons[i-1] = '<';
					msa->ss_cons[j-1] = '>';	
				}
				i = li;
				j = lj;	
				PushNstack(pda, v);
				PushNstack(pda, i);
				PushNstack(pda, j);
				PushNstack(pda, DUMMY_nd); /* we don't know yet what the next node will be */
				nstates += 6;		/* MATP_nd -> MP_st, ML_st, MR_st, D_st, IL_st, IR_st */
				nnodes++;
			}
			else{
				v = InsertTraceNode(gtr, v, TRACE_LEFT_CHILD, i, j, BIF_nd);
				for (mid = i+1; mid < j-1; mid++)
					if (zmat[j][i] == zmat[mid][i] + zmat[j][mid+1])
					{
						/* push the right BEGIN node first */
						PushNstack(pda, v);
						PushNstack(pda, mid+1);
						PushNstack(pda, j);
						PushNstack(pda, BEGR_nd);
						/* then push the left BEGIN node */
						PushNstack(pda, v);
						PushNstack(pda, i);
						PushNstack(pda, mid);
						PushNstack(pda, BEGL_nd);
						nstates += 1;		/* BIF_nd -> B_st */
						nnodes++;
						break;
					}
					mat_flag = ZMAT;
			}
		}


		FreeNstack(pda);
		Free2DArray((void **)mxy,  msa->alen+1);
		Free2DArray((void **)pxy,  msa->alen+1);
		Free2DArray((void **)xy,  msa->alen+1);
		Free2DArray((void **)zmat,msa->alen+1);
		Free2DArray((void **)pmat,msa->alen+1);
		free(matassign);

		if (consensus_bp < 3) return NULL;

		/* OK, we've converted ct into gtr -- gtr is a tree structure telling us the
		* arrangement of consensus nodes. Now do the drill for constructing a full model
		* using this guide tree.
		*/
		cm = CreateCM(nnodes , nstates);
		cm_from_guide(cm, gtr);
		CMZero(cm);

		for (idx = 0; idx < msa->nseq; idx++)
		{
			tr = Transmogrify(cm, gtr, dsq[idx], msa->aseq[idx], msa->alen);
			ParsetreeCount(cm, tr, dsq[idx], msa->wgt[idx]);
			FreeParsetree(tr);
		}
		FreeParsetree(gtr);
		return cm;
}




static void
dump_mat(int    **mat,
		 int      alen)
{
	int i, j;
	for (j = 2; j < alen; j++)
	{
		printf("%2d", j);
		if (mat[j] != NULL){
			for (i = 1; i < j; i++)
				printf("%4d", mat[j][i]);
			puts("\n");
		}

	}
}
