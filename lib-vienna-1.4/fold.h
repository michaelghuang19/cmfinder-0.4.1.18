/* function from fold.c */
extern void  comp_energy(char *sequence); 
extern void  get_structure(int start, int stop, char* seq, char* structure, char backtrack_type);
extern float  fold(char *sequence, char* structure); 
extern int *bp_energy;
extern int *ml_energy;
/* calculate mfe-structure of sequence */
extern float  energy_of_struct(char *string, char *structure);
/* calculate energy of string on structure */
extern void   free_arrays(void);           /* free arrays for mfe folding */
extern void   initialize_fold(int length); /* allocate arrays for folding */
extern void   update_fold_params(void);    /* recalculate parameters */
