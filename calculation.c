#include<stdio.h>
#include<math.h>
#include"calculation.h"

#define EPS 0.0001 

//Jacobi法で使う。対角成分の上側で絶対値が一番大きな値を返し、p,qにその値の行、列をいれる。
static double GetMaxvalue(Matrix *matrix, int *p, int *q);

Vector *CreateVector(int size)
{
	Vector *temp;

	if((temp = (Vector*)malloc(sizeof(Vector))) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}

	if((temp->v = (double*)malloc(sizeof(double)*size)) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}

	temp->size = size;

	return temp;
}

void FreeVector(Vector* vector)
{
	free(vector->v);
	free(vector);
}

Matrix *CreateMatrix(int width, int height)
{
	Matrix *temp;

	if((temp = (Matrix*)malloc(sizeof(Matrix))) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}

	if((temp->a = (double*)malloc(sizeof(double)*width*height)) == NULL){
		fprintf(stderr, "Allocation Error\n");
		exit(1);
	}

	temp->width = width;
	temp->height = height;

	return temp;
}

void FreeMatrix(Matrix* matrix)
{
	free(matrix->a);
	free(matrix);
}

static double GetMaxvalue(Matrix *matrix, int *p, int *q)
{
	int i, j;
	double max;
	double temp;
	
	max = fabs(matrix->a[0*matrix->width + 1]);
	*p = 0;
	*q = 1;

	for(i=0; i<matrix->height; i++){
		for(j=i+1; j<matrix->width; j++){
			temp = fabs(matrix->a[i*matrix->width + j]);
			
			if(temp > max){
				max = temp;
				*p = i;
				*q = j;
			}
		}
	}

	return max;
}

Vector *Jacobi(Matrix* matrix, Matrix* eigenvectors)
{
	int i, j;
	int p, q;
	double max;
	double app, apq, aqq;
	double alpha, beta, gamma;
	double s, c;
	double temp;
	Vector *vec;

	for(i=0; i<eigenvectors->width; i++){
		for(j=0; j<eigenvectors->height; j++){
			if(i==j) eigenvectors->a[i*eigenvectors->width + j] = 1;
			else eigenvectors->a[i*eigenvectors->width + j] = 0;
		}
	}

	do{
		if(!(max = GetMaxvalue(matrix, &p, &q))) break;

		app = matrix->a[p*matrix->width + p];
		apq = matrix->a[p*matrix->width + q];
		aqq = matrix->a[q*matrix->width + q];

		alpha = (app - aqq)/2;
		beta = -apq;
		gamma = fabs(alpha)/sqrt(alpha*alpha + beta*beta);

		s = sqrt((1 - gamma)/2);
		c = sqrt((1 + gamma)/2);
		if(alpha*beta < 0) s = -s;
		
		for(i=0; i<matrix->height; i++){
			temp = c*matrix->a[p*matrix->width + i] - s*matrix->a[q*matrix->width + i];
			matrix->a[q*matrix->width + i] = s*matrix->a[p*matrix->width + i] + c*matrix->a[q*matrix->width + i];
			matrix->a[p*matrix->width + i] = temp;
		}

		for(i=0; i<matrix->height; i++){
			matrix->a[i*matrix->width + p] = matrix->a[p*matrix->width + i];
			matrix->a[i*matrix->width + q] = matrix->a[q*matrix->width + i];
		}

		matrix->a[p*matrix->width + p] = c*c*app + s*s*aqq - 2*s*c*apq;
		matrix->a[p*matrix->width + q] = s*c*(app-aqq) + (c*c - s*s)*apq;
		matrix->a[q*matrix->width + p] = s*c*(app-aqq) + (c*c - s*s)*apq;
		matrix->a[q*matrix->width + q] = s*s*app + c*c*aqq + 2*s*c*apq;

		for(i=0; i<matrix->height; i++){
			temp = c*eigenvectors->a[i*eigenvectors->width + p] - s*eigenvectors->a[i*eigenvectors->width + q];
			eigenvectors->a[i*eigenvectors->width + q] = s*eigenvectors->a[i*eigenvectors->width + p] + c*eigenvectors->a[i*eigenvectors->width + q];
			eigenvectors->a[i*eigenvectors->width + p] = temp;
		}
	}while(max >= EPS);

	vec = CreateVector(matrix->height);
	for(i=0; i<matrix->height; i++){
		vec->v[i] = matrix->a[i*matrix->width + i];
	}

	return vec;
}

