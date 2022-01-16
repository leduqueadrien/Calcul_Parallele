
#include "prodScalaire_dynamic_seq.h"


void showVector(double *A, int n, char* vector_name)
{
	printf("Affichage de %s\n", vector_name);
	if (n >= 100) {
		printf("Affichage de %s impossible, la vector est trop grande\n", vector_name);
		return;
	}

	int i;

	for (i=0; i<n; i++)
			printf("%6.2f%c", A[i], '\t');
	printf("\n");

}

double * allocVector(int n)
{
	double *M=NULL;

	M = (double *)malloc(n * sizeof(double));

	return M;
}

double * freeVector(double *A)
{
	free(A);
	A = NULL;

	return A;
}

int alloc2Vectors(double** A, double** B, int n) {
	*A=allocVector(n);
	if (*A != NULL) {
		*B=allocVector(n);
		if (*B == NULL) {
			*A=freeVector(*A);
			return 1;
		}
	} else {
		return 1;
	}

	return 0;
}


double vectorDotProduct(double* A, double* B, double* pd, int n)
{
	double pdsc = 0;
	int i;
	double start_time = omp_get_wtime();

	for (i=0; i<n; i++) 
		pdsc += A[i] * B[i];
	*pd = pdsc;
	return omp_get_wtime() - start_time;
}

double vectorDotProductParrallel(double* A, double* B, double* pd, int n, int nbThread)
{
	double pdsc = 0;
	int i;
	double start_time = omp_get_wtime();

	#pragma omp parallel for shared(A, B) reduction(+: pdsc) num_threads(nbThread)
	for (i=0; i<n; i++) 
		pdsc += A[i] * B[i];
	
	*pd = pdsc;
	return omp_get_wtime() - start_time;
}

