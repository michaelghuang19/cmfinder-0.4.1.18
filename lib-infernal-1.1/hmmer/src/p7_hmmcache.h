/* A cached profile database. Used by the hmmpgmd daemon.
 */
#ifndef P7_HMMCACHE_INCLUDED
#define P7_HMMCACHE_INCLUDED

#include "esl_alphabet.h"
#include "hmmer.h"

typedef struct {
  char               *name;        /* name of the hmm database              */
  ESL_ALPHABET       *abc;         /* alphabet for database                 */

  P7_OPROFILE       **list;        /* list of profiles [0 .. n-1]           */
  uint32_t            lalloc;	   /* allocated length of <list>            */
  uint32_t            n;           /* number of entries in <list>           */
} P7_HMMCACHE;

extern int    p7_hmmcache_Open (char *hmmfile, P7_HMMCACHE **ret_cache, char *errbuf);
extern size_t p7_hmmcache_Sizeof         (P7_HMMCACHE *cache);
extern int    p7_hmmcache_SetNumericNames(P7_HMMCACHE *cache);
extern void   p7_hmmcache_Close          (P7_HMMCACHE *cache);

#endif /*P7_HMMCACHE_INCLUDED*/

/*****************************************************************
 * Infernal - inference of RNA secondary structure alignments
 * Version 1.1rc2; December 2012
 * Copyright (C) 2012 Howard Hughes Medical Institute.
 * Other copyrights also apply. See the COPYRIGHT file for a full list.
 * 
 * Infernal is distributed under the terms of the GNU General Public License
 * (GPLv3). See the LICENSE file for details.
 * 
 * SVN $Id: p7_hmmcache.h 3754 2011-11-21 14:25:31Z eddys $
 * SVN $URL: https://svn.janelia.org/eddylab/eddys/src/hmmer/branches/infernal/1.1/src/p7_hmmcache.h $
 *****************************************************************/
