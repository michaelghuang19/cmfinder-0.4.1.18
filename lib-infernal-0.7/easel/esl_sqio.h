/* esl_sqio.h
 * Sequence file i/o.
 * 
 * SVN $Id: esl_sqio.h 82 2005-12-13 20:21:07Z eddy $
 */
#ifndef ESL_SQIO_INCLUDED
#define ESL_SQIO_INCLUDED

#include <stdio.h>
#ifdef eslAUGMENT_MSA
#include <esl_msa.h>
#endif

/* name, accession, description, and sequence itself are of unlimited
 * length, but they are initially allocated to something sensible, as
 * set below, in the hope that any given ESL_SQ object only has to
 * make one malloc() call for them. These lengths are inclusive of the
 * \0 NUL character (so ESL_SQ_NAMELEN of 32 means we expect names <=
 * 31 chars).
 * 
 * The reallocation rule is to double the allocation every time 
 * we need to reallocate. So, with the SEQCHUNK set to 256, for example,
 * sequences may be allocated with length 256, 512, 1024, 2048... etc.
 * This results in >50% utilization of our memory, which means
 * some wastage, but we have to compromise with efficiency in
 * sequence reading speed. The esl_sq_Squeeze() function provides
 * for optimizing memory after seqs have been read, at a cost of
 * a ~20% speed hit (for memcpy()'ing within the seq to perfected space.
 */
#define eslSQ_NAMECHUNK   32	/* allocation unit for name         */
#define eslSQ_ACCCHUNK    32	/* allocation unit for accession    */
#define eslSQ_DESCCHUNK  128	/* allocation unit for description  */
#define eslSQ_SEQCHUNK   256	/* allocation unit for seqs         */

/* fread() is apparently the fastest portable way to input from disk;
 * the READBUFSIZE is the fixed size of a block to bring in at one time,
 * for character-based parsers (like the FASTA parser).
 */
#define eslREADBUFSIZE  4096

/* Unaligned file format codes
 * These codes are coordinated with the msa module.
 *   - 0 is an unknown/unassigned format (eslSQFILE_UNKNOWN, eslMSAFILE_UNKNOWN)
 *   - <=100 is reserved for sqio, for unaligned formats
 *   - >100  is reserved for msa, for aligned formats
 */
#define eslSQFILE_UNKNOWN 0
#define eslSQFILE_FASTA   1
#define eslSQFILE_EMBL    2	/* EMBL/Swissprot/TrEMBL */
#define eslSQFILE_GENBANK 3	/* Genbank */
#define eslSQFILE_DDBJ    4	/* DDBJ (currently passed to Genbank parser */
#define eslSQFILE_UNIPROT 5     /* Uniprot (passed to EMBL parser) */


typedef struct {
  FILE *fp;		      /* Open file ptr                            */
  char *filename;	      /* Name of file (for diagnostics)           */
  char *ssifile;	      /* Name of SSI index file (for diagnostics) */
  int   format;		      /* Format of this file                      */
  int   do_gzip;	      /* TRUE if we're reading from gzip -dc pipe */
  int   do_stdin;	      /* TRUE if we're reading from stdin         */
  char  errbuf[eslERRBUFSIZE];/* parse error mesg. Size must match msa.h  */

  /* Our input buffer, whether using character-based parser [fread()]
   * or line-based parser (esl_fgets()).
   */
  char *buf;		      /* buffer for fread() or fgets() input      */
  int   balloc;		      /* allocated size of buf                    */
  int   nc;		      /* #chars in buf (usually full, less at EOF)*/ 
  int   pos;		      /* current parsing position in the buffer   */
  int   linenumber;	      /* What line of the file this is (1..N)     */

  /* Format-specific information
   */
  int   addfirst;             /* TRUE to parse first line of seq record */
  int   addend;	              /* TRUE to parse last line of seq record */
  int   eof_is_ok;	      /* TRUE if record can end on EOF */
  int  (*endTest)(char *);    /* ptr to function that tests if buffer is end */
  int  *inmap;		      /* pointer to an input map, 0..255       */

  /* SSI indexing eventually goes here, including rpl,bpl counting;
   * xref squid.
   */

  /* Optional MSA augmentation: the ability to read multiple alignment
   * files as sequential seq files.
   */
#ifdef eslAUGMENT_MSA
  ESL_MSAFILE *afp;	      /* open ESL_MSAFILE for reading           */
  ESL_MSA     *msa;	      /* preloaded alignment to draw seqs from  */
  int          idx;	      /* index of next seq to return, 0..nseq-1 */
#endif /*eslAUGMENT_MSA*/
} ESL_SQFILE;




/* Allocation of name, acc, desc is all in one malloc, of length
 * nalloc+aalloc+dalloc.
 */
typedef struct {
  /*::cexcerpt::sqio_sq::begin::*/
  char *name;           /* name (mandatory)                                 */
  char *acc;            /* optional accession ("\0" if no accession)        */
  char *desc;           /* description ("\0" if no description)             */
  char *seq;            /* sequence (mandatory) [0..n-1]                    */
  char *ss;             /* secondary structure annotation [0..n-1], or NULL */
  char *dsq;            /* digitized sequence [1..n], or NULL               */
  int   n;              /* length of seq                                    */
  /*::cexcerpt::sqio_sq::end::*/

  char *optmem;         /* optimized mem storage area; see esl_sq_Squeeze() */
  int   nalloc;         /* allocated length of name */
  int   aalloc;         /* allocated length of accession */
  int   dalloc;         /* allocated length of description */
  int   salloc;         /* current allocation length for seq */
} ESL_SQ;


extern ESL_SQ *esl_sq_Create(void);
extern int     esl_sq_Reuse(ESL_SQ *sq);
extern int     esl_sq_Squeeze(ESL_SQ *sq);
extern void    esl_sq_Destroy(ESL_SQ *sq);

extern int  esl_sqfile_Open(char *seqfile, int fmt, char *env, 
			    ESL_SQFILE **ret_sqfp);
extern void esl_sqfile_Close(ESL_SQFILE *sqfp);
extern int  esl_sqfile_DetermineFormat(FILE *fp);

extern int  esl_sq_Read(ESL_SQFILE *sqfp, ESL_SQ *s);
extern int  esl_sq_Write(FILE *fp, ESL_SQ *s, int format);

#ifdef eslAUGMENT_MSA
extern int esl_sq_Dealign(char *s, char *aseq, char *gapstring, int alen);
#endif

extern int   esl_sqfile_FormatCode(char *fmtstring);
extern char *esl_sqfile_FormatString(int fmt);
extern int   esl_sqfile_IsAlignment(int fmt);

#endif /*!ESL_SQIO_INCLUDED*/

/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/
