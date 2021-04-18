#include "edit_cost.h"
#include <stdio.h>
#include <string.h>
#define GAP_SYMB '-'
#define is_gap(c) (c==GAP_SYMB || c == '_' || c== '-'|| c== '.')


int PairCodeVienna(char c1, char c2, char* b) 
{  
  /* Hairpin */
  if ( is_gap(c1)){
    *b = c1;
    return 1;
  }  
  else if ( is_gap(c2) ){
    *b = c2;
    return 1;
  }    
  else{
    *b = c1;
    *(b+1) = c2;
    return 2;
  }
}


/*---------------------------------------------------------------------------*/

void encode( int type, char* label)
{
    int   i, l;

    l = 0;
    for (i = 0; i < type; i++) {
        while (coding[l] != sep && coding[l]) l++;
        l++;
    }

    for (i = 0; coding[l+i] != sep; i++) {
        if (coding[l+i] == '\0') break;
        label[i] = coding[l+i];
    }
    label[i] = '\0';
}




int decode(char *id)
{
    int   n, quit, i;
    char  label[100], *code;

    n = 0;

    quit = 0;
    code = coding;

    while (!quit) {
        for (i = 0; code[i] != sep; i++) {
            if (code[i] == '\0') {
                quit = 1;
                break;
            }
            label[i] = code[i];
        }
        label[i] = '\0';
        if (strcmp(id, label) == 0) return (n);
        code += (i+1);
        n++;
    }

    fprintf(stderr,"Syntax error: node identifier \"%s\" not found "
		   "in coding string \"%s\"\n", id, coding);
    fprintf(stderr,"Exiting...");
}
