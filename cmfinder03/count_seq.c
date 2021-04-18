#include <stdio.h>
#include <stdlib.h>

#include "squid.h"		/* general sequence analysis library    */
int main(int argc, char* argv[])
{

  char      *seqfile=NULL;        /* training sequence file                    */
  char 	    **rseqs;	          /* training sequences                        */
  SQINFO    *sqinfo;		  /* array of sqinfo structures for rseqs      */
  int  	    nseq;		  /* number of seqs */                           
  int       total_size=0;
  int       i;
  
  seqfile =argv[1];

  /* read the training seqs from file */
  if (! ReadMultipleRseqs(seqfile, SQFILE_FASTA, &rseqs, &sqinfo, &nseq))
    Die("Failed to read any sequences from file %s", seqfile);
  for(i=0; i < nseq; i++){
    total_size += sqinfo[i].len;
  }
  printf("%d %d\n", nseq, total_size / nseq);  
  for (i = 0; i < nseq; i++)
    FreeSequence(rseqs[i], &(sqinfo[i]));
  free(sqinfo);  
}
