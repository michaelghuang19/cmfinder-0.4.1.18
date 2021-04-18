/*********************************************************************

  optimize.h

  Contains general functions for optimizing.

  000919 Bjarne Knudsen (bk@daimi.au.dk)

*********************************************************************/

#ifndef __optimize_h__
#define __optimize_h__

#include <stdio.h>
#include <stdlib.h>

void InitMinimize(double ax, double bx, double tol);
int Minimize(double *x, double *f);

#endif
