#include "edouble.h"

#ifdef ETYPE

void FixEdouble(Edouble *);

/* sets Edouble to value of double */
Edouble Dbl2Edbl(double d)
{
  Edouble e;

  e.man = d;
  e.exp = 0;

  FixEdouble(&e);

  return e;
}

/* add second Edouble to first) */
void AddEdouble(Edouble *e1, Edouble e2)
{
  int i;

  if ((e1->exp > e2.exp+60 && e1->man != 0) || e2.man == 0)
    ;
  else if (e2.exp > e1->exp+60 || e1->man == 0) {
    e1->man = e2.man;
    e1->exp = e2.exp; }
  else {
    i = e1->exp-e2.exp;
    e1->exp = e2.exp;
    
    e1->man = ldexp(e1->man, i);
    e1->man += e2.man;
    
    FixEdouble(e1);
  }
}

void SubEdouble(Edouble *e1, Edouble e2)
{
  e2.man = -e2.man;
  AddEdouble(e1, e2);
}

/* set first Edouble to max. Return true if e1 changed, false otherwise */
int MaxEdouble(Edouble *e1, Edouble e2)
{
  if ((e1->exp > e2.exp && e1->man != 0) ||
      e2.man == 0 ||
      (e1->exp == e2.exp && e1->man > e2.man))
    return 0;
  else {
    e1->man = e2.man;
    e1->exp = e2.exp;
    return 1; }
}

/* add Edoubles */
Edouble SumEdouble(Edouble e1, Edouble e2)
{
  AddEdouble(&e1, e2);

  return e1;
}

/* add Edoubles */
Edouble DiffEdouble(Edouble e1, Edouble e2)
{
  SubEdouble(&e1, e2);

  return e1;
}

/* multiply Edoubles */
void MulEdouble(Edouble *e1, Edouble e2)
{
  e1->man *= e2.man;
  e1->exp += e2.exp;
  
  FixEdouble(e1);
}

/* multiply Edoubles */
Edouble ProdEdouble(Edouble e1, Edouble e2)
{
  MulEdouble(&e1, e2);

  return e1;
}

/* divide Edoubles */
void DivEdouble(Edouble *e1, Edouble e2)
{
  e1->man /= e2.man;
  e1->exp -= e2.exp;

  FixEdouble(e1);
}

/* divide Edoubles */
Edouble QoutEdouble(Edouble e1, Edouble e2)
{
  DivEdouble(&e1, e2);

  return e1;
}

/* Logarithm of Edouble */
Edouble LogEdouble(Edouble e)
{
  Edouble a;

  a.exp = 0;
  a.man = log(e.man) + log(2.)*e.exp;

  FixEdouble(&a);

  return a;
}

/* Exponentiate Edouble */
Edouble ExpEdouble(Edouble e)
{
  Edouble a;
  double d;

  d = Edbl2Dbl(e);

  a.exp = d/log(2.);
  a.man = exp(d - log(2.)*a.exp);

  FixEdouble(&a);

  return a;
}

/* Convert Edouble to double */
double Edbl2Dbl(Edouble e)
{
  return ldexp(e.man, e.exp);
}

/* Move exponent from e->man to e->exp */
void FixEdouble(Edouble *e)
{
  int i;

  e->man = frexp(e->man, &i);

  e->exp += i;
}

/* Print Edouble to fp */
void PrintEdouble(FILE *fp, Edouble e)
{
  fprintf(fp, "%d %g\n", e.exp, e.man);
}

int IsZeroEdouble(Edouble e)
{
  return (e.man == 0);
}

int PositiveEdouble(Edouble e)
{
  return (e.man > 0);
}

Fdouble Edbl2Fdbl(Edouble e)
{
  Fdouble fd;

  fd = (Fdouble) (Edbl2Dbl(LogEdouble(e))*FMUL);

  if (fd < FDBL_MIN)
    fd = FDBL_MIN;

  return fd;
}

Fdouble Dbl2Fdbl(double d)
{
  Fdouble fd;

  fd = (d == 0) ? FDBL_MIN : (Fdouble)(log(d)*FMUL);

  if (fd < FDBL_MIN)
    fd = FDBL_MIN;

  return fd;
}

#endif

