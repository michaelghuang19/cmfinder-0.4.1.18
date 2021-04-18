/*********************************************************************

  matrix.h

  Contains general functions for handling the matrices.

  000701 Bjarne Knudsen (bk@daimi.au.dk)

*********************************************************************/

#ifndef __matrix_h__
#define __matrix_h__

#include "edouble.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

typedef struct tagMatrix {
  int rows;
  int cols;
  Edouble **entry;
} Matrix;

Matrix *MakeMatrix(int rows, int cols);
Matrix *ZeroMatrix(int rows, int cols);
Matrix *IdentityMatrix(int rows);
void FreeMatrix(Matrix *matrix);
Matrix *PartMatrix(Matrix *mat, int *rows);
Matrix *CopyMatrix(Matrix *mat);
void SetMatrix(Matrix *mat1, Matrix *mat2);
Matrix *MulMatrix(Matrix *mat1, Matrix *mat2);
void MulMatrix_2(Matrix *matr, Matrix *mat1, Matrix *mat2);
Matrix *ExpMatrix(double c, Matrix *eigen, Matrix *diag, Matrix *inveigen);
Matrix *SubMatrix(Matrix *mat, int row, int col);
Matrix *DiagMatrix(Matrix *mat);
Matrix *TransposeMatrix(Matrix *mat);
Edouble MatrixDet(Matrix *mat);
void PrintMatrix(FILE *fp, Matrix *mat);

#endif
