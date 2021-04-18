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


#ifndef STRUCTSH_INCLUDED
#define STRUCTSH_INCLUDED

/* structs.h
 * header file for the three models (other, conding, rna).
 * Defines structures and macros.
 *
 * ER, Fri May 28 11:36:57 CDT 1999 [STL]
 *
 * 
 */ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

#include "squid.h"
#include "mat_globals.h"

/*
 * RIBOPROB Matrix type
 *
 * Contains array in one dimension (to be indexed later), matrix size,
 * H, and E.
 */
typedef struct _matrix_t {
  double *matrix;
  int     edge_size;         /* Size of one edge, e.g. 4 for 4x4 matrix */
  int     full_size;         /* Num of elements, e.g. 10 for 4x4 matirx */
  double  H;
  double  E;
} matrix_t;

typedef struct _condmatrix_t {
  double *matrix;
  int     size;         /* Size of one edge, e.g. 4 for 4x4 matrix */
} condmatrix_t;

/*
 * Full matrix definition
 */
typedef struct _fullmat_t {
  matrix_t *unpaired;
  matrix_t *paired;
  char     *name;
} fullmat_t;

typedef struct _fullcondmat_t {
  condmatrix_t *marg;
  condmatrix_t *cond;
  char         *name;
} fullcondmat_t;

/* I define a substitution process, as an evolutionary process
 * in which the probabilities P(i,j) are symmetric, and therefore they
 * can be interpred as Joint probabilities.
 * For a substitution proces, there is only one set of
 * conditional probabilties Q(i,j) = P(j|i), 
 * and one set of marginal probabilities pm(i) =sum_j P(i,j)
 * which verify the condition:
 *
 *    P(i,j) = pm(i) * Q(i,j) = pm(j) *Q(j,i)
 *
 * Because this condition is trivialy satisfied for i=j, we can use
 * the Wheelan ant Goodman process to change the marginal probabilites.
 *
 * The method is the following:
 *
 * assuming the evolutionary process: Q(t) = Qnougt * e^{ t/t^* * Rate }
 * 
 *                             where: Rate = 1/t^* log { Qnought^{-1} * Q^* }
 *
 * The new Rate matrix is given by,
 *
 * ~R(i,j) = R(i,j) * ~pm(j)/pm(j) * (~pm(j)*~pm(i) / pm(j)*pm(i) ) ^alpha(i,j)  
 *
 *                                      ,, alpha(i,j) = alpha(j,i) and in GW case they all take the value -1/2 
 *
 *
 * If Qnought is different from the identity matrix, it has to be changed. I use the Iterative  process:
 *
 *    Pnought(i,j) = p(i)*Qnought(i,j) = expS_{i,j} *  pm(i)*pm(j) --> ~Pnought(i,j) = exp{ lamda*S_{i,j} } * ~pm(i)*~pm(j) 
 *
 *                                                                ---> marginals assymtotically approximate ~pm.
 *
 *
 * Then it is easy to see that
 *
 *        ~P(i,j) = ~pm(i) * ~Q(i,j) = ~pm(j) * ~Q(j,i)  for all times, including t=0
 *
 */
struct psubs_s {
  int     L;        /* dimension      */
  double time;      /* time parameter */

  /* Joint, Conditional, and Marginal probabilities
   *
   * The relationship is,
   *
   *           P(i,j) = pm(i) * Q(i,j) = pm(j) * Q(j,i)
   *
   */
  double *P;        /* LxL Joint       probabilities -- P_{ij} = P_{ji} */
  double *Q;        /* LxL Conditional probabilities -- Q_{ij} \def P(j|i) */
  double *pm;       /* L   Marginal    probabilities -- pm_i   \def \sum_j P_{ji} */

  /* the evolutionary process is
   *
   *      Q(time) = Qnought * e^{ time * Rate }
   *
   *      Rate    = 1/time * log { Qnought^{-1} * Q(time) }
   */
  double *Qnought;  /* LxL Conditional probabilities at time zero */
  double *Rate;     /* LxL Rate Matrix */
};

/* I define a Non-substitution process, as an evolutionary process
 * in which the probabilities P(i,j) are NOT symmetric, and therefore they
 * canNOT be interpreted as Joint probabilities.
 * For a Non-substitution proces, there are TWO sets of
 * conditional probabilties 
 *                            Ql(i,j) = P(jr|il), 
 *                            Qr(i,j) = P(jl|ir), 
 * and TWO sets of marginal probabilities 
 *                            pml(i) = sum_j P(i,j)
 *                            pmr(i) = sum_j P(j,i)
 * which verify the condition:
 *
 *                            P(i,j) = pml(i) * Ql(i,j) = pmr(j) *Qr(j,i)
 *
 * Because this condition is NOT trivialy satisfied for i=j, we canNOT use
 * the Wheelan ant Goodman process to change the marginal probabilites.
 *
 *
 * To change the background basecomposition, I use the Iterative  process:
 *
 *    P(i,j) = pml(i)*Ql(i,j) = expS_{i,j} *  pml(i)*pmr(j)
 *
 *              ---> ~P(i,j) = exp{ lamda*S_{i,j} } * ~pml(i)*~pmr(j) 
 *
 *              ---> marginals assymtotically approximate ~pml and ~pmr.
 *
 *
 * Then it is easy to see that
 *
 *        ~P(i,j) = ~pml(i) * ~Ql(i,j) = ~pmr(j) * ~Qr(j,i)  for all times, including t=0
 *
 */
struct pnonsubs_s {
  int     L;        /* dimension      */
  double time;      /* time parameter */

  /* Nonsubs P's, Conditional, and Marginal probabilities
   *
   * The relationship is,
   *
   *           P(i,j) = pml(i) * Ql(i,j) = pmr(j) * Qr(j,i)
   *
   */
  double *P;        /* LxL nonsubstitution   probabilities -- P(i,j) \def P_{il,jr} neq P_{jr,il} */
  double *Ql;       /* LxL left-Conditional  probabilities -- Ql_{il,jr} \def P(jr|il) = P(i,j)/pml(i) */
  double *Qr;       /* LxL right-Conditional probabilities -- Qr_{ir,jl} \def P(jl|ir) = P(j,i)/pmr(j) */
  double *pml;      /* L   left-Marginal     probabilities -- pml_i      \def \sum_j P_{i,j} */
  double *pmr;      /* L   right-Marginal    probabilities -- pmr_i      \def \sum_j P_{j,i} */

  /* the evolutionary process is
   *
   *      Q(time) = Qnought * e^{ time * Rate }
   *
   *      Rate    = 1/time * log { Qnought^{-1} * Q(time) }
   */
  double *Qlnought;  /* LxL left-Conditional probabilities at time zero */
  double *Qrnought;  /* LxL right-Conditional probabilities at time zero */
  double *Ratel;     /* LxL left-Rate Matrix */
  double *Rater;     /* LxL rightRate Matrix */
};


/**********************************************************************
 * Time/Id correlation function ID = id_zero - a*(time-t_zero)^b for id \in [id_zero,id_cutoff]
 **********************************************************************/
struct three_times_s {
  double oth;
  double cod;
  double rna;
};

struct divergence_s {
  double id_zero;
  double id_cutoff;
  double t_zero;
  double a;
  double b;
};

struct three_divergence_s {
  struct divergence_s prxmal;  /* model id/time correlation function proximal id[100:75] */
  struct divergence_s medial;  /* model id/time correlation function medial   id[75:50]  */
  struct divergence_s distal;  /* model id/time correlation function distal   id[50:0]   */
};

extern struct three_divergence_s  othdiv;            /* OTH model id/time correlation functions */
extern struct three_divergence_s  coddiv;            /* COD model id/time correlation functions */
extern struct three_divergence_s  rnadiv_ribosum;    /* RNA model id/time correlation functions */
extern struct three_divergence_s  rnadiv_approx;     /* RNA model id/time correlation functions */
extern struct three_divergence_s  rnadiv_oldrna;     /* RNA model id/time correlation functions */

struct eigenvalue_s {
  double *real;
  double *imag;
};

struct zmatrix_s {
  double *real;
  double *imag;
};

struct complex_s {
  double real;
  double imag;
};

struct ali_s { /* arrays for an alignment of two seqs */
  char *charX;
  char *charY; /* arrays for symbols */
};

struct dos_s {                        /* list of all possible combinations of scoring systems */
  int frdiag;   /* TRUE  ==  calculates forward diagonal           */
  int frfull;   /* TRUE  ==  calculates forward full               */
  int frsemi;   /* TRUE  ==  calculates forward semi-full          */
  int fulldp;   /* TRUE  ==  calculates forward full               */
  int global;   /* TRUE  ==  calculates global alignments          */
  int nus;      /* TRUE  ==  uses NUS model for rna                */
  int semidp;   /* TRUE  ==  calculates semi-full dp scores        */
  int vidiag;   /* TRUE  ==  calculates viterbi diagonal           */
  int vifull;   /* TRUE  ==  calculates viterbi full               */
  int visemi;   /* TRUE  ==  calculates viterbi semi-full          */
  int forward;  /* TRUE  ==  calculates forward dp scores          */
  int twodiag;  /* TRUE  ==  both Forward and Viterbi diagonal dp  */
  int twofull;  /* TRUE  ==  both Forward and Viterbi full dp      */
  int twosemi;  /* TRUE  ==  both Forward and Viterbi semi-full dp */
};

struct dpd_s { /* arrays for DP matrices -  diagonal version */
  struct coddp_s  *cod;
  struct othdpd_s *oth;
  struct rnadp_s  *rna;
};

struct dpdscan_s { /* arrays for DP matrices -  scanning-diagonal version */
  struct othdpscan2_s *othscan2;
  struct coddpscan2_s *codscan2;
  struct rnadpscan2_s *rnascan2;
};

struct dpdscanfast_s { /* arrays for DP matrices -  scanning-diagonal version */
  struct othdpscanfast2_s *othscan2;
  struct coddpscanfast2_s *codscan2;
  struct rnadpscanfast2_s *rnascan2;
};

struct dpf_s { /* arrays for DP matrices -  full version */
  struct coddp_s  *cod;
  struct othdpf_s *oth;
  struct rnadp_s  *rna;
};

struct end_s { 
  int *lend;
  int *rend;
};

struct end3_s { 
  struct end_s  *cod; /* the ends for COD */
  struct end_s  *oth; /* the ends for OTH */
  struct end_s  *rna; /* the ends for RNA */
};

struct windowends_s { 
  struct  end3_s  *fwd; /* the ends for the forward strand */
  struct  end3_s  *rev; /* the ends for the reverse strand */
};

struct endscan_s { /* it will report the first local region, the last, and (MAX_NUM_ENDS-2) in between */
  int *lend;
  int *rend;
};

struct end3scan_s { 
  struct endscan_s  *cod; /* the ends for COD */
  struct endscan_s  *oth; /* the ends for OTH */
  struct endscan_s  *rna; /* the ends for RNA */
};



struct which_model_s {
  int null;
  int oth;
  int cod;
  int rna;
};

struct model_s { /* arrays for the three models + null model */
  struct codmodel_s  *cod; /* the model for COD */
  struct othmodel_s  *oth; /* the model for OTH */
  struct rnamodel_s  *rna; /* the model for RNA */

  struct nullmodel_s *null;
};

struct emodel_s { /* arrays for the three models + null model */
  struct codmodel_s  *cod; /* the model for COD */
  struct othmodel_s  *oth; /* the model for OTH */
  struct rnamodel_s  *rna; /* the model for RNA */

  struct nullmodel_s *null;

  struct which_model_s *which_model;
};

struct dpsc3_s { /* arrays for forward scores */
  double *oth;
  double *cod;
  double *rna;
};

struct sc2_s { /* to display plus and minus strand scores for a given model */
  double pl;
  double mn;
};

struct sc3_s { /* to display the three model scores */
  struct sc2_s *oth;
  struct sc2_s *cod;
  struct sc2_s *rna;
};

struct scansc_s { /* arrays for scan scores */
  double *oth;
  double *cod;
  double *rna;

  double *othrev;
  double *codrev;
  double *rnarev;
};

struct scanends_s { /* arrays for scan ends */
  struct  end3scan_s  **fwd; /* the ends for the forward strand */
  struct  end3scan_s  **rev; /* the ends for the reverse strand */
};

struct scanfast_s {  /* arrays for scan scores and corresponding local regions */
  struct scansc_s   *sc;
  struct scanends_s *ends;

};

struct scores_s { /* to display all possible scores */
  double null;                  /* null model score (given   strands)                     */
  double nullrev;               /* null model score (reverse strands)                     */
  
  struct sc3_s   *global;       /* global    dp LOD scores                                */
  struct sc3_s   *vidiag;       /* diagonal  dp LOD scores (viterbi)                      */
  struct sc3_s   *frdiag;       /* diagonal  dp LOD scores (forward)                      */
  struct sc3_s   *vifull;       /* full      dp LOD scores (viterbi)                      */
  struct sc3_s   *frfull;       /* full      dp LOD scores (forward)                      */
  struct sc3_s   *visemi;       /* semi-full dp LOD scores (viterbi)                      */
  struct sc3_s   *frsemi;       /* semi-full dp LOD scores (forward)                      */
  
  struct dpsc3_s *ardiag;       /* arrays to store scores for forward calculation         */
  struct dpsc3_s *arfull;       /* arrays to store scores for forward calculation         */
};

/**********************************************************************
 * Single-nucleotide probabilities. [To set the base composition of the different models]
 **********************************************************************/
struct singlep_s {
  double pa;                 /* p[a]  */
  double pc;                 /* p[c]  */
  double pg;                 /* p[g]  */
};
extern struct singlep_s  singlep;     /* OTH model single-emission probs */


/**********************************************************************
 * NULL-model
 **********************************************************************/
struct nullparam_s {
  double eta;                 /* S->N
                                 N->B  */
};

struct nullmodel_s {
  double  eta;          /* transition prob's.  t = eta         +*/
  double  meta;         /* transition prob's.  t = 1-eta        +*/
  double *xem;          /* RX insert emissions. xem[0..3]       +*/
  double *yem;          /* RY insert emissions. yem[0..3]       +*/
};

extern struct nullparam_s Nullparam; /* parameter for the null model */

/**********************************************************************
 * OTH-model
 **********************************************************************/
/* STATES othmodel 
 */
#define stFL      0
#define stB       1
#define stM       2
#define stX       3
#define stY       4
#define stE       5
#define stFJ      6
#define stFR      7
#define OSTATES   8

/* Indices for OTH model state transitions.
 * Used for indexing oth->t[]
 * mnemonic: Transition from Match to Match = TMM
 */
#define TFLB     0
#define TFLFR    1
#define TBM      2
#define TBX      3
#define TBY      4
#define TMM      5
#define TMX      6
#define TMY      7
#define TME      8
#define TXM      9
#define TXX     10
#define TXY     11
#define TXE     12
#define TYM     13
#define TYX     14
#define TYY     15
#define TYE     16
#define TEFJ    17
#define TEFR    18
#define TFJB    19
#define TFRFR   20
#define OTRANS  21

extern int   EmitsPerOTHState[OSTATES];    /* allowed emissions per state in the othmodel       */
extern int   IdxTransOTH[OSTATES];         /* starting point in transition indexing             */
extern int   FromStateOTH[OTRANS];         /* state from which the transition originates        */
extern int   OConnects[OSTATES][OSTATES];  /* allowed transitions in the othmodel               */
extern char *ostNAME[OSTATES];             /* ASCII states     names                            */
extern char *otrNAME[OTRANS];              /* ASCII transition names                            */
extern int   TransPerOTHState[OSTATES];    /* non-zero transitions per state in the othmodel    */
extern int   ToStateOTH[OTRANS];           /*state to which the transition goes                 */

/* Structures for othmodel 
 */
struct othparam_s {
  struct nullparam_s FLN; /* parameter for the flanking left   null model */
  struct nullparam_s FJN; /* parameter for the flanking middle null model */
  struct nullparam_s FRN; /* parameter for the flanking right  null model */

  double delta;                 /* M->X|Y         transitions   */
  double epsilon;               /* X->X   
				   Y->Y           transitions   */
  double eta;                   /* E->F_R         transitions   */
  double gamma;                 /* X->Y
				   Y->X           transitions   */
  double kappa;                 /* B->X|Y         transitions   */
  double tau;                   /* M->E 
				   X->E
				   Y->E           transitions   */
  double xi;                    /* B->E           transitions   */
};

extern struct othparam_s  OTHparam;           /* parameter for the othmodel                  */
extern struct othparam_s  OTHparam_zero;      /* parameter for the othmodel                  */
extern struct othparam_s  OTHparam_infty;     /* parameter for the othmodel                  */

struct othmodel_s {

  /*   Transition probabilities are usually accessed as a
   *   two-D array: oth->t[TMM], for instance. They are allocated
   *   such that they can also be stepped through in 1D by pointer
   *   manipulations, for efficiency in DP algorithms.
   */
  double *t;            /* transition prob's.  t[0..OTRNAS]    +*/
  double *mem;          /* match    emissions. mem[0..15]      +*/ 
  double *xem;          /* X insert emissions. xem[0..3]       +*/
  double *yem;          /* Y insert emissions. yem[0..3]       +*/

  struct nullmodel_s *FLN; /* Flanking model to the left       +*/
  struct nullmodel_s *FJN; /* Flanking model to the middle     +*/
  struct nullmodel_s *FRN; /* Flanking model to the right      +*/
};

struct othdpd_s {
  double  *flmx;
  double **fjmx;
  double  *frmx;

  double *bmx;
  double *mmx;
  double *xmx;
  double *ymx;
  double *emx;
};

struct othdpscan2_s {
  struct othdpscan_s *othscan;
  struct othdpscan_s *htoscan;

};
struct othdpscanfast2_s {
  struct othdpscanfast_s *othscan;
  struct othdpscanfast_s *htoscan;

};

struct othdpscan_s {
  double  *flmx;
  double **fjmx;
  double  *frmx;

  double *bmx;
  double *mmx;
  double *xmx;
  double *ymx;
  double *emx;
};

struct othdpscanfast_s {
  double **flmx;
  double **fjmx;
  double **frmx;

  double *bmx;
  double *mmx;
  double *xmx;
  double *ymx;
  double *emx;
};

struct othdpf_s {
  double  *flxmx;
  double **fjxmx;
  double  *frxmx;

  double  *flymx;
  double **fjymx;
  double  *frymx;

  double *bmx;
  double *mmx;
  double *xmx;
  double *ymx;
  double *emx;
};

/**********************************************************************
 * COD-model
 **********************************************************************/

/* STATES codmodel 
 */
#define stCOB     0
#define stCb      1
#define stC33     2
#define stC34     3
#define stC32     4
#define stC43     5
#define stC44     6
#define stC42     7
#define stC23     8
#define stC24     9
#define stC22    10
#define stC30    11
#define stC40    12
#define stC20    13
#define stC03    14
#define stC04    15
#define stC02    16
#define stCe     17
#define stCOJ    18
#define stCOE    19
#define CSTATES  20 

/* Indices for COD model state transitions.
 * Used for indexing cod->t[]
 * mnemonic: Transition from Match to Match = TMM
 */
#define TCOBCb    0
#define TCOBCOE   1
#define TCbC33    2
#define TCbC34    3
#define TCbC32    4
#define TCbC43    5
#define TCbC44    6
#define TCbC42    7
#define TCbC23    8
#define TCbC24    9
#define TCbC22   10
#define TCbC30   11
#define TCbC40   12
#define TCbC20   13
#define TCbC03   14
#define TCbC04   15
#define TCbC02   16
#define TCbCe    17
#define TC33Ce   18
#define TC34Ce   19
#define TC32Ce   20
#define TC43Ce   21
#define TC44Ce   22
#define TC42Ce   23
#define TC23Ce   24
#define TC24Ce   25
#define TC22Ce   26
#define TC30Ce   27
#define TC40Ce   28
#define TC20Ce   29
#define TC03Ce   30
#define TC04Ce   31
#define TC02Ce   32
#define TCeCOJ   33
#define TCeCOE   34
#define TCOJCb   35
#define TCOECOE  36
#define CTRANS   37

extern int   EmitsPerCODState[CSTATES];  /* allowed emissions per state in the codmodel     */
extern int   IdxTransCOD[CSTATES];       /* starting point in transiton indexing            */
extern char *cstNAME[CSTATES];           /* ASCII states     names                          */
extern char *ctrNAME[CTRANS];            /* ASCII transition names                          */
extern int   FromStateCOD[CTRANS];
extern int   ToStateCOD[CTRANS];
extern int   TransPerCODState[CSTATES];  /* non-zero transitions per state in the codmodel  */

/* Structures for codmodel 
 */
struct codparam_s {
  double psi;         /* COB -> Cb   */
  double eta;         /* Ce  -> COE  */
  double tau;         /* Cb  -> Ce   */
  double xi[16];      /* Cb  -> Cij  */
};

extern struct othparam_s  COBparam;        /* parameter for the othmodel O_B in codmodel  */
extern struct othparam_s  COJparam;        /* parameter for the othmodel O_J in codmodel  */
extern struct othparam_s  COEparam;        /* parameter for the othmodel O_E in codmodel  */
extern struct othparam_s  COBparam_zero;   /* parameter for the othmodel O_B in codmodel  */
extern struct othparam_s  COJparam_zero;   /* parameter for the othmodel O_J in codmodel  */
extern struct othparam_s  COEparam_zero;   /* parameter for the othmodel O_E in codmodel  */
extern struct othparam_s  COBparam_infty;  /* parameter for the othmodel O_B in codmodel  */
extern struct othparam_s  COJparam_infty;  /* parameter for the othmodel O_J in codmodel  */
extern struct othparam_s  COEparam_infty;  /* parameter for the othmodel O_E in codmodel  */

extern struct codparam_s  CODparam;       /* parameter for the codmodel                  */
extern struct codparam_s  CODparam_zero;  /* parameter for the codmodel                  */
extern struct codparam_s  CODparam_infty; /* parameter for the codmodel                  */

struct codmodel_s {
  double  *t;           /* transition prob's.  t[0..OTRNAS]            +*/
  double  *pmut;        /* mutation  prob's. pmut[16]                  +*/
  double **pcodon;      /* codon     emission  prob's. pcodon[64][64]  +*/
  double **phexa;       /* dicodon   prob's. phexa[64][64]             +*/

  struct othmodel_s *COB;
  struct othmodel_s *COJ;
  struct othmodel_s *COE;

};

struct coddp_s {

  double **COJ;

  double  *cbmx;
  double  *c33mx;
  double  *cemx;
  double  *cobmx;
  double  *cojmx;
};

struct coddpx_s {

  double **COJ;

  double  *cbmx;
  double  *c33mx;
  double  *cemx;
  double  *cobmx;
  double  *cojmx;
};

struct coddpxfast_s {

  double **COJ;

  double  *cbmx;
  double  *c33mx;
  double  *cemx;
  double  *cobmx;
  double  *cojmx;
  double  *coemx;
};

struct coddpsmart_s {

  double **COJ;

  double **cbmx;
  double **c33mx;
  double **cemx;
  double **cobmx;
  double **cojmx;
};

struct coddpscan_s {
  struct othdpscan_s *cob;
  struct othdpscan_s *coj;
  struct othdpscan_s *coe;

  struct coddpx_s  *cod;
};

struct coddpscanfast_s {
  struct othdpscanfast_s *cob;
  struct othdpscanfast_s *coj;
  struct othdpscanfast_s *coe;

  struct coddpxfast_s  *cod;
};

struct coddpscan2_s {
  struct coddpscan_s *codscan;
  struct coddpscan_s *docscan;
};

struct coddpscanfast2_s {
  struct coddpscanfast_s *codscan;
  struct coddpscanfast_s *docscan;
};

struct pamcond_s {
  double **pXY;
  double **pYX;
};

struct emitcodon_s {
  double *pcod;
  double *pcod_hexa;
  double *pcod30;
  double *pcod32;
  double *pcod34;
};

/**********************************************************************
 * RNA-model
 **********************************************************************/
 /* constraints: GAPP2GAPP >= GAPP * GAPP 
 *               INDL2INDL >= INDL * INDL 
 *
 */
#define GAPP 0.18     /* for nt-paired-to-a-gap  at t*  */
#define INDL 0.18     /* for nt-mutates-to-a-gap at t*  */

#define GAPP2GAPP 0.04    /* for gap-paired-to-a-gap  at t*  */
#define INDL2INDL 0.04     /* for gap-mutates-to-a-gap at t*  */

#define INDL_INFTY 1.0  /* probability of finding a gap at infty */
#define GAPP_INFTY 1.0  /* probability of finding a gap at infty */

/* STATES rnamodel 
 */
#define stROB    0
#define stRNA    1
#define stROJ    2
#define stROE    3
#define RSTATES  4 

/* Indices for RNA model state transitions.
 * Used for indexing rna->t[]
 * mnemonic: Transition from Match to Match = TMM
 */
#define TROBRNA  0
#define TROBROE  1
#define TRNAROJ  2
#define TRNAROE  3
#define TROJRNA  4
#define TROEROE  5
#define RTRANS   6

extern int   IdxTransRNA[RSTATES];      /* starting point in transiton indexing            */
extern char *rstNAME[RSTATES];          /* ASCII states     names                          */
extern char *rtrNAME[RTRANS];           /* ASCII transition names                          */
extern int   FromStateRNA[RTRANS];
extern int   ToStateRNA[RTRANS];
extern int   TransPerRNAState[RSTATES]; /* non-zero transitions per state in the rnamodel  */

struct rnaparam_s {
  double psi;                 /* ROB->RNA */
  double eta;                 /* RNA->ROE */
};

extern struct othparam_s  ROBparam;            /* parameter for the othmodel O_B  in rnamodel  */
extern struct othparam_s  ROJparam;            /* parameter for the othmodel O_J  in rnamodel  */
extern struct othparam_s  ROEparam;            /* parameter for the othmodel O_E  in rnamodel  */
extern struct othparam_s  Rloopparam;          /* parameter for the othmodel loop in rnamodel  */
extern struct othparam_s  ROBparam_zero;       /* parameter for the othmodel O_B  in rnamodel  */
extern struct othparam_s  ROJparam_zero;       /* parameter for the othmodel O_J  in rnamodel  */
extern struct othparam_s  ROEparam_zero;       /* parameter for the othmodel O_E  in rnamodel  */
extern struct othparam_s  Rloopparam_zero;     /* parameter for the othmodel loop in rnamodel  */
extern struct othparam_s  ROBparam_infty;      /* parameter for the othmodel O_B  in rnamodel  */
extern struct othparam_s  ROJparam_infty;      /* parameter for the othmodel O_J  in rnamodel  */
extern struct othparam_s  ROEparam_infty;      /* parameter for the othmodel O_E  in rnamodel  */
extern struct othparam_s  Rloopparam_infty;    /* parameter for the othmodel loop in rnamodel  */



extern struct rnaparam_s  RNAparam;            /* parameter for the rnamodel                  */
extern struct rnaparam_s  RNAparam_zero;       /* parameter for the rnamodel                  */
extern struct rnaparam_s  RNAparam_infty;      /* parameter for the rnamodel                  */

struct rnadp_s {

  double **ROJ;

  double  *rbmx;
  double  *rrmx;
  double  *remx;
  double  *rjmx;
};

struct rnadpx_s {

  double **ROJ;

  double  *rbmx;
  double  *rrmx;
  double  *remx;
  double  *rjmx;
};

struct rnadpxfast_s {

  double **ROJ;

  double  *rbmx;
  double  *rrmx;
  double  *remx;
  double  *rjmx;
  double  *roemx;
};

struct rnadpscan_s {
  struct othdpscan_s *rob;
  struct othdpscan_s *roj;
  struct othdpscan_s *roe;

  struct rnadpx_s  *rna;
};

struct rnadpscanfast_s {
  struct othdpscanfast_s *rob;
  struct othdpscanfast_s *roj;
  struct othdpscanfast_s *roe;

  struct rnadpxfast_s  *rna;
};

struct rnadpscanfast2_s {
  struct rnadpscanfast_s *rnascan;
  struct rnadpscanfast_s *anrscan;
};

struct rnadpscan2_s {
  struct rnadpscan_s *rnascan;
  struct rnadpscan_s *anrscan;
};

struct rnamodel_s {

  double  *t;

  struct othmodel_s *ROB;
  struct othmodel_s *ROJ;
  struct othmodel_s *ROE;

  struct nusmodel_s *nus;
  struct pi2model_s *pi2;
};

struct rnascfg_s {
  double  *sc;    /* to keep the scores */
  double  *vp;    /* to calculate internal loops in L^3 */

  struct rnamtx_s    *in;    /* Inside a lgorithm dp matrices  */
  struct rnamtx_s    *inrv;  /* Inside  algorithm dp matrices -- revese strand */

  struct rnamtx_ou_s *ou;    /* Outside algorithm dp matrices -- only one because I will do the posterior only for the winner strand */
  
  struct nrnscfg_s   *nrn;   /* matrices  for non-redundant nussinov (robin's) grammar*/

 };

struct rnascfgscan_s {
  struct rnamtxscan_s  *in;   /* Inside algorithm dp matrices  */
  struct rnamtxscan_s  *inrv; /* Inside algorithm dp matrices -- revese strand */
};

struct rnascfgscanfast_s {
  struct rnamtxscanfast_s  *in;   /* Inside algorithm dp matrices  */
  struct rnamtxscanfast_s  *inrv; /* Inside algorithm dp matrices -- revese strand */
};

struct rnamtxscan_s {
  double  *sc;    /* to keep the scores */
  double  *vp;    /* to calculate internal loops in L^3 */

  double **rnaj;          /* to calculate single stranded inside an IS with a geometrical distribution */
  struct othdpscan_s *othj;  /* weird othdp structure needed here if you want to calculate rnaj in the "exact" way with the scanning version.
			     kind of a hack*/

  /* the 3 states of the grammar 
   */
  double **vx; 
  double **wx;
  double **wbx;
};

struct rnamtxscanfast_s {
  double  *sc;    /* to keep the scores */
  double  *vp;    /* to calculate internal loops in L^3 */

  double **rnaj;          /* to calculate single stranded inside an IS with a geometrical distribution */
  struct othdpscanfast_s *othj;  /* weird othdp structure needed here if you want to calculate rnaj in the "exact" way with the scanning version.
			     kind of a hack*/

  /* the 3 states of the grammar 
   */
  double **vx; 
  double **wx;
  double **wbx;
};

struct rnamtx_s {
  double **rnaj;  /* to calculate single stranded inside an IS with a geometrical distribution */

  /* the 3 states of the grammar 
   */
  double **vx; 
  double **wx;
  double **wbx;
};

struct rnamtx_ou_s {
  double **rnaj;  /* to calculate single stranded inside an IS with a geometrical distribution */

  /* the 3 states of the grammar 
   */
  double **vx; 
  double **wx;
  double **wbx;
};

struct nrnscfg_s {
  double **sx;
  double **lx;
  double **rx;
};

/* Structures for PI2model (this is the probabilistic model for rna up to Irred surf of order 2)
 */
struct pi2model_s {

  struct othmodel_s *Rloop;

  struct rna_v5 *v;
  struct rna_w5 *w;
  struct rna_w5 *wb;

  struct rna_loop5 *is1;
  struct rna_loop5 *is2b;
  struct rna_loop5 *is2i;
};


/* Nussinov model (this is the probabilistic model for rna up to Irred surf of order 0)
 */
struct nusmodel_s {
  double tl;		    /* W->L transition           */
  double tr;		    /* W->R transition           */
  double tv;		    /* W->V transition           */
  double tw;		    /* W->WW bifurcation         */
  double te;		    /* W->e  end in a hairpin    */

  double  p[25];	    /* P(x1,y1) for left-right  singlets   */
  double  pp[25][25];       /* P(x1,y1;x2,y2) for paired emissions */
};

struct rna_w5 {
  double tl;		    /* W->L transition       */
  double tr;		    /* W->R transition       */
  double tv;		    /* W->V transition       */
  double tw;		    /* W->WW bifurcation     */

  double  pl[25];	    /* P(x1,y1) for left singlets          */
  double  pr[25];	    /* P(x1,y1) for right singlets         */
  double  pp[25][25];       /* P(x1,y1;x2,y2) for paired emissions */
};

struct rna_w {
  double tl;		    /* W->L transition       */
  double tr;		    /* W->R transition       */
  double tv;		    /* W->V transition       */
  double tw;		    /* W->WW bifurcation     */

  double  pl[16];	    /* P(x1,y1) for left singlets          */
  double  pr[16];	    /* P(x1,y1) for right singlets         */
  double  pp[16][16];       /* P(x1,y1;x2,y2) for paired emissions */
};

struct rna_v5 { 
  double t1;		          /* V->IS1 transition                 */
  double t2s;		          /* V->IS2 transition  (stems)        */
  double t2b;		          /* V->IS2 transition  (bulges)       */
  double t2i;		          /* V->IS2 transition  (intloops)     */
  double t3;		          /* V->WB WB transition               */

  double pp[25][25];              /* P(i,j) for paired emissions       */
};

struct rna_v { 
  double t1;		          /* V->IS1 transition                 */
  double t2s;		          /* V->IS2 transition  (stems)        */
  double t2b;		          /* V->IS2 transition  (bulges)       */
  double t2i;		          /* V->IS2 transition  (intloops)     */
  double t3;		          /* V->WB WB transition               */

  double pp[16][16];              /* P(i,j) for paired emissions       */
  double p2p[16][16][16][16];     /* P(i,j ; k,l) for 2 pair emissions */
};
struct rna_loop {
  double tn[MAXRNALOOP];
 
  double ps[16];	   /* P(x1,y1) for singlets in loop  */
};
struct rna_loop5 {
  double tn[MAXRNALOOP];
 
  double ps[25];	   /* P(x1,y1) for singlets in loop  */
};

struct rna_w_param_s {
  double wl;		    /* W->L transition       */
  double wr;		    /* W->R transition       */
  double wv;		    /* W->V transition       */
  double www;		    /* W->WW bifurcation     */
};

struct rna_v_param_s {
  double v1;		          /* V->IS1 transition                 */
  double v2s;		          /* V->IS2 transition  (stems)        */
  double v2b;		          /* V->IS2 transition  (bulges)       */
  double v2i;		          /* V->IS2 transition  (intloops)     */
  double v3;		          /* V->WB WB transition               */
};

struct scfg_param_s {
  struct rna_v_param_s vt;
  struct rna_w_param_s wt;
  struct rna_w_param_s wbt;
};

extern struct scfg_param_s  SCFGparam;           /* parameter for the SCFG model                */
extern struct scfg_param_s  SCFGparam_zero;      /* parameter for the SCFG model                */
extern struct scfg_param_s  SCFGparam_infty;     /* parameter for the SCFG model                */


/* non-redundant R.Dowell's nussinov grammar for posteriors
 */
struct nrn_s {
  struct nrnS_s *S;
  struct nrnL_s *L;
  struct nrnR_s *R;
};

struct nrnS_s {
  double   tss;
  double   tsl;
  double   tsr;
  double   tsb;
  double   tse;

  double **pp;
  double  *ps;
};

struct nrnL_s {
  double   tls;
  double   tll;
  double   tle;

  double **pp;
  double  *ps;
};

struct nrnR_s {
  double   trr;
  double   tre;

  double  *ps;
};

/* ungapped
 */
struct rnamodelungapped_s {

  double  *t;

  struct othmodel_s *ROB;
  struct othmodel_s *ROJ;
  struct othmodel_s *ROE;

  struct rna_v *v;
  struct rna_w *w;
  struct rna_w *wb;

  struct rna_loop *is1;
  struct rna_loop *is2b;
  struct rna_loop *is2i;
};

struct nusmodelungapped_s {
  double  *t;

  struct othmodel_s *ROB;
  struct othmodel_s *ROJ;
  struct othmodel_s *ROE;

  struct rna_w *w;
};

struct rna_s {
  double tsl;		     /* W->L transition       */
  double tsr;		     /* W->R transition       */
  double tsp;		     /* W->V transition       */
  double tss;		     /* W->WW bifurcation     */
  double tse;		     /* W->e b                */

  double  px[4][4];	     /* P(x1,y1) for singlets               */
  double  pxy[4][4][4][4];   /* P(x1,y1;x2,y2) for paired emissions */
};


/* unambiguous B.Knudsen's nussinov grammar for posteriors
 */
struct ubk_s {
  struct ubkS_s *S;
  struct ubkL_s *L;
  struct ubkF_s *F;
};

struct ubkS_s {
  double   tsl;
  double   tsb;
};

struct ubkL_s {
  double   tlf;
  double   tla;
};

struct ubkF_s {
  double   tff;
  double   tfb;
};


#endif /* STRUCTSH_INCLUDED */
















