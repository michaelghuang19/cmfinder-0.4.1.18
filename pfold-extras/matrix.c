#include "matrix.h"
/* Make a matrix */
Matrix *MakeMatrix(int rows, int cols)
{
  Matrix *mat;
  int i;

  mat = (Matrix *)malloc(sizeof(Matrix));

  mat->rows = rows;
  mat->cols = cols;

  mat->entry = (Edouble **)malloc(rows * sizeof(Edouble *));
  for (i = 0; i < rows; i++)
    mat->entry[i] = (Edouble *)malloc(cols * sizeof(Edouble));

  return mat;
}

/* Return a zero matrix */
Matrix *ZeroMatrix(int rows, int cols)
{
  Matrix *mat;
  int i, j;

  mat = MakeMatrix(rows, cols);  

  for (i = 0; i < mat->rows; i++)
    for (j = 0; j < mat->cols; j++)
      mat->entry[i][j] = Dbl2Edbl(0.);

  return mat;
}

/* Return a square identity matrix */
Matrix *IdentityMatrix(int rows)
{
  Matrix *mat;
  int i;

  mat = ZeroMatrix(rows, rows);

  for (i = 0; i < mat->rows; i++)
    mat->entry[i][i] = Dbl2Edbl(1.);

  return mat;
}

/* Free a matrix from memory */
void FreeMatrix(Matrix *mat)
{
  int i;

  if (mat != NULL) {
    for (i = 0; i < mat->rows; i++)
      free(mat->entry[i]);
    free(mat->entry);

    free(mat);
  }
}

Matrix *PartMatrix(Matrix *mat, int *rows)
{
  Matrix *temp;
  int i, j, i2, j2;
  int size;

  size = 0;

  for (i = 0; i < mat->rows; i++)
    if (rows[i] == 1)
      size++;

  temp = MakeMatrix(size, size);

  for (i = 0, i2 = 0; i < size; i++, i2++) {
    while (rows[i2] == 0)
      i2++;
    for (j = 0, j2 = 0; j < size; j++, j2++) {
      while (rows[j2] == 0)
	j2++;
      temp->entry[i][j] = mat->entry[i2][j2];
    }
  }

  return temp;
}

/* Copies a matrix */
Matrix *CopyMatrix(Matrix *mat)
{
  Matrix *temp;
  int i, j;

  temp = MakeMatrix(mat->rows, mat->cols);

  for (i = 0; i < mat->rows; i++)
    for (j = 0; j < mat->cols; j++)
      temp->entry[i][j] = mat->entry[i][j];

  return temp;
}

/* Copies mat2 into mat1 */
void SetMatrix(Matrix *mat1, Matrix *mat2)
{
  int i, j;

  for (i = 0; i < mat2->rows; i++)
    for (j = 0; j < mat2->cols; j++)
      mat1->entry[i][j] = mat2->entry[i][j];
}

/* Multiplies two matrices */
Matrix *MulMatrix(Matrix *mat1, Matrix *mat2)
{
  Matrix *temp;
  int i, j, k;

  temp = MakeMatrix(mat1->rows, mat2->cols);

  for (i = 0; i < mat1->rows; i++)
    for (j = 0; j < mat2->cols; j++) {
      temp->entry[i][j] = Dbl2Edbl(0.);
      for (k = 0; k < mat1->cols; k++)
	AddEdouble(&temp->entry[i][j], ProdEdouble(mat1->entry[i][k], mat2->entry[k][j]));
      /*	temp->entry[i][j] += mat1->entry[i][k]*mat2->entry[k][j];*/
    }

  return temp;
}

/* Multiplies two matrices and puts result in matr */
void MulMatrix_2(Matrix *matr, Matrix *mat1, Matrix *mat2)
{
  int i, j, k;

  for (i = 0; i < mat1->rows; i++)
    for (j = 0; j < mat2->cols; j++) {
      matr->entry[i][j] = Dbl2Edbl(0.);
      for (k = 0; k < mat1->cols; k++)
	AddEdouble(&matr->entry[i][j], ProdEdouble(mat1->entry[i][k], mat2->entry[k][j]));
      /*	matr->entry[i][j] += mat1->entry[i][k]*mat2->entry[k][j];*/
    }
}

/*
   Returns exp(c*matrix), where eigen is the eigenvector matrix, diag
   is the diganonalized matrix and inveigen is the inverse eigenvector
   matrix.
*/
Matrix *ExpMatrix(double c, Matrix *eigen, Matrix *diag, Matrix *inveigen)
{
  Matrix *temp, *temp2;
  int i;
  int size;

  size = diag->rows;

  /* Exponentiate diag */
  temp = CopyMatrix(diag);
  for (i = 0; i < size; i++)
    temp->entry[i][i] = ExpEdouble(ProdEdouble(Dbl2Edbl(c), temp->entry[i][i]));

  temp2 = MulMatrix(eigen, temp);
  FreeMatrix(temp);
  temp = MulMatrix(temp2, inveigen);
  FreeMatrix(temp2);

  return temp;
}

/* Returns mat without row and col. */
Matrix *SubMatrix(Matrix *mat, int row, int col)
{
  Matrix *temp;
  int i, j, i2, j2;

  temp = MakeMatrix(mat->rows-1, mat->cols-1);

  for (i = i2 = 0; i < temp->rows; i++, i2++)
    for (j = j2 = 0; j < temp->cols; j++, j2++) {
      if (i2 == row)
	i2++;
      if (j2 == col)
	j2++;
      temp->entry[i][j] = mat->entry[i2][j2];
    }

  return temp;
}

/* Makes a diagonal matrix from mat (1*n or n*1) */
Matrix *DiagMatrix(Matrix *mat)
{
  Matrix *temp;
  int i, j;

  if (mat->rows == 1)
    temp = MakeMatrix(mat->cols, mat->cols);
  else
    temp = MakeMatrix(mat->rows, mat->rows);

  for (i = 0; i < mat->rows; i++)
    for (j = 0; j < mat->cols; j++)
      temp->entry[i][j] = Dbl2Edbl(0.);

  if (mat->rows == 1)
    for (i = 0; i < mat->cols; i++)
      temp->entry[i][i] = mat->entry[0][i];
  else
    for (i = 0; i < mat->rows; i++)
      temp->entry[i][i] = mat->entry[i][0];

  return temp;
}

/* Return transposed matrix */
Matrix *TransposeMatrix(Matrix *mat)
{
  Matrix *temp;
  int i, j;

  temp = MakeMatrix(mat->cols, mat->rows);

  for (i = 0; i < mat->rows; i++)
    for (j = 0; j < mat->cols; j++)
      temp->entry[j][i] = mat->entry[i][j];

  return temp;
}

/* Recursive determinant calculation */
Edouble MatrixDet(Matrix *mat)
{
  int i;
  int sign;
  Edouble det;

  det = Dbl2Edbl(0.);
  sign = 1;

  if (mat->rows == 1)
    return mat->entry[0][0];

  for (i = 0; i < mat->rows; i++) {
    AddEdouble(&det, ProdEdouble(Dbl2Edbl((double) sign), ProdEdouble(mat->entry[i][0], MatrixDet(SubMatrix(mat, i, 0)))));
    sign *= -1;
  }

  return det;
}

void PrintMatrix(FILE *fp, Matrix *mat)
{
  int i, j;

  for (i = 0; i < mat->rows; i++) {
    for (j = 0; j < mat->cols; j++)
      fprintf(fp, " %5.2g", Edbl2Dbl(mat->entry[i][j]));
    fprintf(fp, "\n");
  }
}
