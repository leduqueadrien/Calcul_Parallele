
#include "prodMatriciel_dynamic_seq.h"


void showMatrice(double **A, int n, char* matrice_name)
{
	printf("Affichage de %s\n", matrice_name);
	if (n >= 100) {
		printf("Affichage de %s impossible, la matrice est trop grande\n", matrice_name);
		return;
	}

	int i, j;

	for (i=0; i<n; i++) {
		for (j=0; j<n; j++)
			printf("%6.2f%c", A[i][j], ((j+1)%SHOWSIZE) ? '\t' : '\n');
		printf("\n");
	}
	printf("\n");

}

double ** allocMatrice(int n)
{
	double **M=NULL;
	int i, j;

	M = (double **)malloc(n * sizeof(double *));
	if (M) {
		for (i=0; i<n; i++) {
			M[i] = (double *)malloc(n * sizeof(double));
			if (!(M[i])) {
				for (j=i-1; j>=0; j--) free(M[j]);
				free(M);
				M = NULL;
				break;
			}
		}
	}

	return M;
}

double ** freeMatrice(double **A, int n)
{
	int i;

	for (i=n-1; i>=0; i--)
		free(A[i]);
	free(A);
	A = NULL;

	return A;
}

int alloc3Matrices(double*** A, double*** B, double*** C, int n) {
	*A=allocMatrice(n);
	if (*A != NULL) {
		*B=allocMatrice(n);
		if (*B != NULL) {
			*C=allocMatrice(n);

		} else {
			*A=freeMatrice(*A, n);
			return 1;
		}
	} else {
		return 1;
	}

	if (*C == NULL) {
		*B = freeMatrice(*B, n);
		*A = freeMatrice(*A, n);

		return 1;
	}

	return 0;
}


double MatriceDotProduct(double** A, double** B, double** C, int n) {
	double start_time = omp_get_wtime();
	int i, j, k;
	for (i=0; i<n; i++) 
		for (j=0; j<n; j++) {
			C[i][j] = 0;
			for (k=0; k<n; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
	return omp_get_wtime() - start_time;
}

double MatriceDotProductParrallel(double** A, double** B, double** C, int n, int nbThread) {
	double start_time = omp_get_wtime();
	int i, j, k;
	#pragma omp parallel for private(j, k) shared(A, B, C) num_threads(nbThread)
	for (i=0; i<n; i++) 
		for (j=0; j<n; j++) {
			C[i][j] = 0;
			for (k=0; k<n; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
	return omp_get_wtime() - start_time;
}

