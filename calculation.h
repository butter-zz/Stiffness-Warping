#ifndef __CALCULATION_H_INCLUDED__
#define __CALCULATION_H_INCLUDED__

#include <stdlib.h>

typedef struct{
	int size;
	double *v;
}Vector;

typedef struct{
	int height;
	int width;
	double *a;
}Matrix;

Vector *CreateVector(int size);
void FreeVector(Vector* vector);

Matrix *CreateMatrix(int width, int height);
void FreeMatrix(Matrix* matrix);
//ヤコビ法によって固定値をベクトルで返し、eigenvectorsに固有ベクトル
Vector *Jacobi(Matrix* matrix, Matrix* eigenvectors);

#endif /* __CALCULATION_H_INCLUDED__ */
