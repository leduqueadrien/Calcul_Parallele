
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "prodScalaire_dynamic_seq.h"

void initialisation(double* A, double* B, int n)
{
	int i;
	for (i=0; i<n; i++) {
		A[i] = i*n;
		B[i] = i;
	}
}

int main(int argc, char **argv)
{
	double* A=NULL, *B=NULL;
	double ps;
	double compute_time;
	int flags;
	int size=100;
	int nb_thread = 2;

	if (argc > 1) size = atoi(argv[1]);
	if (argc > 2) nb_thread = atoi(argv[2]);


	/* Allocation de 3 vecteurs */
	printf("Allocation des vecteurs de taille %d ...\n", size);
	flags = alloc2Vectors(&A, &B, size);
	if (flags != 0) return 1;

	/* Initialisation des vecteurs A et B */
	printf("Initialisation des vecteurs ...\n");
	initialisation(A, B, size);

	showVector(A, size, "A");
	showVector(B, size, "B");

	/* Calcul du produit A*B */
	printf("Calcul de la vecteurs C=A*B SANS thread ...\n");
	compute_time = vectorDotProduct(A, B, &ps, size);
	printf("Time : %f s\n", compute_time);

	printf("Calcul de la vecteurs C=A*B AVEC %d thread ...\n", nb_thread);
	compute_time = vectorDotProductParrallel(A, B, &ps, size, nb_thread);
	printf("Time : %f s\n", compute_time);

	printf("Resultat du produit scalaire : %f\n", ps);

	// Liberation de memoire
	B = freeVector(B);
	A = freeVector(A);

	return 0;
}
