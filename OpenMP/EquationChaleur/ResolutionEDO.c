
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ResolutionEDO.h"
// y : ligne
// x : colonne

int main() {
	Parametre p;
	p.a = M_PI;
	p.b = M_PI;
	p.ax = -2;
	p.bx = 2;
	p.ay = -2;
	p.by = 2;
	p.N = 10;
	p.M = 10;
	p.hx = (p.bx-p.ax)/(p.M + 1);
	p.hy = (p.by-p.ay)/(p.N + 1);

	p.Tmax = 10;
	p.T = 20;
	p.ht = p.Tmax/p.T;

	p.eps = 1e-10;
	p.maxIter = 1000;
	p.nu = 1;

	// 1 : Jacobi, 2 : Gauss-Seidel
	p.resolutionSysLin = 1;
	p.pResolution = &jacobi;
	if (p.resolutionSysLin == 2)
		p.pResolution = &gauss_Seidel;

	// 1 : euler implicite, 2 : Cranck-Nocolson
	p.schema = 2;
	p.pSchema = &eulerImplicite;
	if (p.schema == 2) {
		p.pSchema = &cranck_Nicolson;
		p.nu = p.nu/2;
	}
	
	// Calcul des constantes de l'algorithme
	p.betaX = p.nu/p.hx/p.hx;
	p.betaY = p.nu/p.hy/p.hy;
	p.lambda = 1/p.ht + 2*p.betaX + 2*p.betaY;
	p.lambdaCrankNicolson = 1/p.ht -2*p.betaX -2*p.betaY;
	p.betaXDivLambda = p.betaX/p.lambda;
	p.betaYDivLambda = p.betaY/p.lambda;
	p.invLambda = 1/p.lambda;
	p.Np2 = p.N+2;
	p.Mp2 = p.M+2;
	p.invHt = 1/p.ht;

	// Nombre de thread
	p.nbThread = 2;

	resolutionEquation(p);
	return 0;
}

void resolutionEquation(Parametre p) {
	double * u = (double *)malloc(p.Np2*p.Mp2*sizeof(double));//u courant
	double * u_nm1 = (double *)malloc(p.Np2*p.Mp2*sizeof(double));//u_n-1
	double * f = (double *)malloc(p.Np2*p.Mp2*sizeof(double));//f, second membre
	double * u_km1 = NULL; //u_k-1

	if (p.resolutionSysLin == 1)//u_k-1 non necessaire si gauss-siedel
		u_km1 = (double *)malloc(p.Np2*p.Mp2*sizeof(double));


	int nbIter = 0;

	initialisation(u, p);
	copier(u, u_nm1, p);
	if (p.resolutionSysLin == 1)
		copier(u, u_km1, p);
	
	// Affichage de l'etat initial
	afficher(u, p, nbIter, 0);

	/* Boucle de temps */
	for (int n=1; n<p.T+1; n++) {
		echangerPointeur(&u, &u_nm1);
		p.pSchema(f, u_nm1, n*p.ht, p);
		nbIter = p.pResolution(u, f, u_km1, p);
		afficher(u, p, nbIter, n*p.ht);
	}

	free(u);
	free(u_km1);
	free(u_nm1);
	free(f);
}

int jacobi(double * u, double * f, double * u_km1, Parametre p) {
	int nbIter = 0;
	double sommeDifTerme;
	double sommeTerme;
	double valeur;
	double sommeU_km1_i;
	double sommeU_km1_j;
	int i_m1;
	int i_p1;
	int j_m1;
	int j_p1;

	#pragma omp parallel private(i, i_m1, i_p1, j, j_m1, j_p1, sommeU_km1_i, sommeU_km1_j) public(p, f)
	do {
		echangerPointeur(&u, &u_km1);
		sommeDifTerme = 0;
		sommeTerme = 0;
		#pragma omp for reduction()
		for (int i=1; i<p.M+1; i++) {
			i_m1 = i-1;
			i_p1 = i+1;
			for (int j=1; j<p.N+1; j++) {
				j_m1 = j-1;
				j_p1 = j+1;
				sommeU_km1_i = u_km1[i_m1*p.Np2 + j] + u_km1[i_p1*p.Np2 + j];
				sommeU_km1_j = u_km1[i*p.Np2 + j_m1] + u_km1[i*p.Np2 + j_p1];
				valeur = p.betaYDivLambda * sommeU_km1_i
						+ p.betaXDivLambda * sommeU_km1_j
						+ p.invLambda * f[i*p.Np2 + j];
				sommeTerme += pow(valeur, 2);
				sommeDifTerme += pow(valeur - u_km1[i*p.Np2 + j], 2);
				u[i*p.Np2 + j] = valeur;
			}	
		}

		#pragma omp single
		nbIter++;
	} while (nbIter < p.maxIter && sommeDifTerme > p.eps*sommeTerme);

	return nbIter;
}

int gauss_Seidel(double * u, double * f, double * u_km1, Parametre p) {
	int nbIter = 0;
	double sommeDifTerme;
	double sommeTerme;
	double valeur;
	double sommeU_i;
	double sommeU_j;
	int i_m1;
	int i_p1;
	int j_m1;
	int j_p1;

	do {
		sommeDifTerme = 0;
		sommeTerme = 0;
		for (int i=1; i<p.M+1; i++) {
			i_m1 = i-1;
			i_p1 = i+1;
			for (int j=1; j<p.N+1; j++) {
				j_m1 = j-1;
				j_p1 = j+1;
				sommeU_i = u[i_m1*p.Np2 + j] + u[i_p1*p.Np2 + j];
				sommeU_j = u[i*p.Np2 + j_m1] + u[i*p.Np2 + j_p1];
				valeur = p.betaYDivLambda * sommeU_i
						+ p.betaXDivLambda * sommeU_j
						+ p.invLambda * f[i*p.Np2 + j];
				sommeTerme += pow(valeur, 2);
				sommeDifTerme += pow(valeur - u[i*p.Np2 + j], 2);
				u[i*p.Np2 + j] = valeur;
			}	
		}

		nbIter++;
	} while (nbIter < p.maxIter && sommeDifTerme > p.eps*sommeTerme);

	return nbIter;
}

void eulerImplicite(double * f, double * u_nm1, double t, Parametre p) {
	for (int i=1; i<p.M+1; i++) {
		for (int j=1; j<p.N+1; j++) {
			f[i*p.Np2 + j] = source(p.ax+j*p.hx, p.ay+i*p.hy, t, p) + p.invHt * u_nm1[i*p.Np2 + j];
		}
	}
}

void cranck_Nicolson(double * f, double * u_nm1, double t, Parametre p) {
	double sommeU_mn1_i;
	double sommeU_mn1_j;
	int i_m1;
	int i_p1;
	int j_m1;
	int j_p1;

	for (int i=1; i<p.M+1; i++) {
			i_m1 = i-1;
			i_p1 = i+1;
		for (int j=1; j<p.N+1; j++) {
			j_m1 = j-1;
			j_p1 = j+1;
			sommeU_mn1_i = u_nm1[i_m1*p.Np2 + j] + u_nm1[i_p1*p.Np2 + j];
			sommeU_mn1_j = u_nm1[i*p.Np2 + j_m1] + u_nm1[i*p.Np2 + j_p1];
			f[i*p.Np2 + j] = source(p.ax+j*p.hx, p.ay+i*p.hy, t, p)
				+ p.lambdaCrankNicolson * u_nm1[i*p.Np2 + j]
				+ p.betaY * sommeU_mn1_i
				+ p.betaX * sommeU_mn1_j;
		}
	}
}

double source(double x, double y, double t, Parametre p) {
	double facteur = (pow(p.a,2)+pow(p.b,2));
	return facteur * sin(p.a*x) * cos(p.b*y) + 0*t;
}

double fonctionCL(double x, double y, Parametre p) {
	return fonctionCI(x, y, p);
}

double fonctionCI(double x, double y, Parametre p) {
	return sin(p.a*x)*cos(p.b*y);
}

void initialisation(double * u, Parametre p) {
	// Condition aux limites
	for (int j=0; j<p.Np2; j++) {
		u[j] = fonctionCL(p.ax+j*p.hx, p.ay, p);
		u[(p.N+1)*p.Np2 + j] = fonctionCL(p.ax+j*p.hx, p.ay+(p.N+1)*p.hy, p);
	}
	for (int i=0; i<p.Mp2; i++) {
		u[i*p.Np2] = fonctionCL(p.ax, p.ay+i*p.hy, p);
		u[i*p.Np2 + p.N+1] = fonctionCL(p.ax+(p.N+1)*p.hx, p.ay+i*p.hy, p);
	}

	// Condition initiale
	for (int i=1; i<p.M+1; i++) {
		for (int j=1; j<p.N+1; j++) {
			u[i*p.Np2 + j] = 0;
		}
	}
}

/* Changer le palier pour modifier la frontiere de couleur */
void afficher(double * u, Parametre p, int nbIter, double t) {
	// On affiche le temps de l'iteration :
	printf("\nt = %.2fs\n", t);
	// Affichage du nombre de termes de la suite de Jacobi calcules
	printf("Nombre d'itération de %s : %d\n", p.resolutionSysLin==1?"Jacobi":"Gauss-Seidel" , nbIter);
	char blanc[6] = "\033[00m";
	// char rouge[6] = "\033[31m";
	// char vert[6] = "\033[32m";
	char jaune[6] = "\033[33m";
	// char bleuF[6] = "\033[34m";
	// char violet[6] = "\033[35m";
	// char bleuC[6] = "\033[36m";
	double valeur;
	for (int i=0; i<p.Mp2; i++) {
		for (int j=0; j<p.Np2; j++) {
			valeur = u[i*p.Np2 + j];
			if (i==0 || i==p.N+1 || j==0 || j==p.M+1) {
				printf("%s", jaune);
			} else {
				printf("%s", blanc);
			}
			// Affichage de la valeur
			printf("%.2f  ", valeur);
		}
		printf("\n");
	}
	// On remet en blanc
	printf("%s", blanc);
}

/* u contient v0 et vL, dest ne les contient pas, d'ou le dest[i-1] pour avoir la valeur du meme indice i */
void copier(double * u, double * dest, Parametre p) {
	for (int i=0; i<p.Mp2; i++) {
		for (int j=0; j<p.Np2; j++) {
			dest[i*p.Np2 + j] = u[i*p.Np2 + j];
		}
	}
}

void echangerPointeur(double ** u, double ** v) {
	double * tmp = *u;
	*u = *v;
	*v = tmp;
}
