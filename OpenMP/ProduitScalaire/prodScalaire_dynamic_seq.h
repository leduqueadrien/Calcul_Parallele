
#ifndef PRODMATRICIEL_DYNAMIC_SEQ_H
#define PRODMATRICIEL_DYNAMIC_SEQ_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// taille maximum des matrices statiques sur etud avec 48 threads (par défaut)

#define SHOWSIZE 10

void showVector(double* A, int n, char* matrice_name);

double * allocVector(int n);

double * freeVector(double* A);

int alloc2Vectors(double** A, double** B, int size);

double vectorDotProductParrallel(double* A, double* B, double* pd, int n, int nbThread);

double vectorDotProduct(double* A, double* B, double* pd, int n);


#endif

