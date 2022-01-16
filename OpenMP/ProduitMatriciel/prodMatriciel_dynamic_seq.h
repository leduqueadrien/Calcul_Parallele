
#ifndef PRODMATRICIEL_DYNAMIC_SEQ_H
#define PRODMATRICIEL_DYNAMIC_SEQ_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// taille maximum des matrices statiques sur etud avec 48 threads (par défaut)

#define SHOWSIZE 10

void showMatrice(double** A, int n, char* matrice_name);

double ** allocMatrice(int n);

double ** freeMatrice(double** A, int n);

int alloc3Matrices(double*** A, double*** B, double*** C, int size);

double MatriceDotProductParrallel(double** A, double** B, double** C, int n, int nbThread);

double MatriceDotProduct(double** A, double** B, double** C, int n);


#endif

