
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "prodMatriciel_dynamic_seq.h"

void initialisation(double** A, double** B, int n)
{
	int i, j;
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			A[i][j] = i*n+j;
			B[i][j] = i+j;
		}
	}
}

int main(int argc, char **argv)
{
	double **A=NULL, **B=NULL, **C=NULL;
	double compute_time;
	int flags;
	int size=100;
	int nb_thread = 2;

	if (argc > 1) size = atoi(argv[1]);
	if (argc > 2) nb_thread = atoi(argv[2]);


	/* Allocation de 3 matrices */
	printf("Allocation des matrices de taille %dx%d ...\n", size, size);
	flags = alloc3Matrices(&A, &B, &C, size);
	if (flags != 0) return 1;

	/* Initialisation de matrices A et B */
	printf("Initialisation des matrices ...\n");
	initialisation(A, B, size);

	showMatrice(A, size, "A");
	showMatrice(B, size, "B");

	/* Calcul du produit A*B */
	printf("Calcul de la matrice C=A*B SANS thread ...\n");
	compute_time = MatriceDotProduct(A, B, C, size);
	printf("Time : %f s\n", compute_time);

	printf("Calcul de la matrice C=A*B AVEC %d thread ...\n", nb_thread);
	compute_time = MatriceDotProductParrallel(A, B, C, size, nb_thread);
	printf("Time : %f s\n", compute_time);

	showMatrice(C, size, "C");

	// Liberation de memoire
	C = freeMatrice(C, size);
	B = freeMatrice(B, size);
	A = freeMatrice(A, size);

	return 0;
}

