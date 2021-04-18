/* keyhash.h
 * Storing keys in hash tables, similar to Perl's associative arrays.
 * 
 * SVN $Id: esl_keyhash.h 34 2005-02-24 17:25:53Z eddy $
 * SRE, Sun Jan 30 08:55:17 2005
 * from squid's gki.h, 1999.
 */

#ifndef eslKEYHASH_INCLUDED
#define eslKEYHASH_INCLUDED

#include <stdio.h>		/* for FILE */

/* esl_key_elem:
 *    key, array index pairs are kept in linked list structures.
 */
struct esl_key_elem {
  char *key;		
  int   idx;
  struct esl_key_elem *nxt;
};

/* ESL_KEYHASH:
 *    a dynamically resized hash structure; 
 *    contains a hash table and associated data
 */
typedef struct {
  struct esl_key_elem **table;
  
  int primelevel;
  int nhash;
  int nkeys;
} ESL_KEYHASH;

extern ESL_KEYHASH *esl_keyhash_Create(void);
extern void         esl_keyhash_Destroy(ESL_KEYHASH *h);
extern void         esl_keyhash_Dump(FILE *fp, ESL_KEYHASH *h);

extern int  esl_key_Store(ESL_KEYHASH *h, char *key, int *ret_index);
extern int  esl_key_Lookup(ESL_KEYHASH *h, char *key);


#endif /* eslKEYHASH_INCLUDED */

/*****************************************************************
 *    This copyrighted source code is freely distributed 
 *    under the terms of the GNU General Public License. See
 *    the files COPYRIGHT and LICENSE for details.
 *****************************************************************/
