/*********************************************************************
 *                                                                   *
 * This interface defines functions for probabilities. They are de-  *
 * fined in 'prob.c', and should be able to hold very small numbers. *
 *                                                                   *
 * 000712 Bjarne Knudsen                                             *
 *                                                                   *
 *********************************************************************/

#ifndef __edouble_h__
#define __edouble_h__

//#define ETYPE EDBL
#include <math.h>
#include <limits.h>
#include <stdio.h>

#ifdef ETYPE
#define FMUL 100
#define FDBL_MIN -2000000

typedef struct tagEdouble {
  double man;  /* mantissa */
  int exp;     /* exponent */
} Edouble;

typedef int Fdouble;

Edouble Dbl2Edbl(double d);            /* sets Edouble to value of double */
void AddEdouble(Edouble *e1, Edouble e2);   /* add second Edouble to first */
void SubEdouble(Edouble *e1, Edouble e2);   /* subtract second Edouble from first */
int MaxEdouble(Edouble *e1, Edouble e2);   /* set first Edouble to max */
Edouble SumEdouble(Edouble e1, Edouble e2);  /* add Edoubles */
Edouble DiffEdouble(Edouble e1, Edouble e2);  /* subtract Edoubles */
void MulEdouble(Edouble *e1, Edouble e2);  /* multiply Edoubles */
Edouble ProdEdouble(Edouble e1, Edouble e2);  /* multiply Edoubles */
void DivEdouble(Edouble *e1, Edouble e2);   /* divide Edoubles */
Edouble QoutEdouble(Edouble e1, Edouble e2);
Edouble ExpEdouble(Edouble e);              /* Exponentiate Edouble */
Edouble LogEdouble(Edouble e);              /* Returns the log of Edouble */
double Edbl2Dbl(Edouble e);              /* Convert Edouble to double */
void PrintEdouble(FILE *fp, Edouble e);
int IsZeroEdouble(Edouble e);
int PositiveEdouble(Edouble e);

Fdouble Edbl2Fdbl(Edouble e);
Fdouble Dbl2Fdbl(double d);

#define MaxFdouble(f1, f2)      ( ( ( (*(f1)) > (f2) ) ? 0 : (((*(f1))=(f2)),1) ))
#define ProdFdouble(f1, f2)     ( (f1) + (f2) )

#else

typedef double Edouble;
typedef double Fdouble;

#define Dbl2Edbl(d)            ( (d) )
#define AddEdouble(e1, e2)     ( (*e1) += (e2) )
#define SubEdouble(e1, e2)     ( (*e1) -= (e2) )
//#define MaxEdouble(e1, e2)     ( *(e1) = (*(e1) > (e2)) ? *(e1) : (e2) )
#define MaxEdouble(e1, e2)     ( (*(e1) > (e2)) ? 0 : (1, *(e1) = (e2)) )
#define SumEdouble(e1, e2)     ( (e1) + (e2) )
#define MulEdouble(e1, e2)     ( (*e1) *= (e2) )
#define ProdEdouble(e1, e2)    ( (e1) * (e2) )
#define DivEdouble(e1, e2)     ( (*e1) /= (e2) )
#define ExpEdouble(e)          (exp(e))
#define LogEdouble(e)          (log(e))
#define Edbl2Dbl(e)            (e)
#define IsZeroEdouble(e)       (e == 0)
#define PositiveEdouble(e)     (e > 0)

#define MaxFdouble(f1, f2)      ( ( ( (*(f1)) > (f2) ) ? 0 : (((*(f1))=(f2)),1) ))

#define Fdbl2Dbl(e)            (e)
#define Dbl2Fdbl(d)            ( (d) )
#define Edbl2Fdbl(e)            ( (e) )
#define ProdFdouble(e1, e2)    ( (e1) * (e2) )




#endif

#endif
