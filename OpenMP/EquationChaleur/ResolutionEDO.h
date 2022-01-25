
#ifndef RESOLUTIONEDO_HPP
#define RESOLUTIONEDO_HPP

/* -------------------------------------------------------------------------- */
/* Parametre : Contient tous les parametres necessaires a la resolution		  */
/* -------------------------------------------------------------------------- */
struct Parametre;

/* -------------------------------------------------------------------------- */
/* ptrFonctionSchema : Pointeur sur la fonction donnant le schema utilise	  */
/* 						pour la resolution									  */
/* -------------------------------------------------------------------------- */
typedef void (*ptrFonctionSchema)(double *, double *, double, struct Parametre);

/* -------------------------------------------------------------------------- */
/* ptrFonctionSchema : Pointeur sur la fonction resolvant le systeme		  */
/* 						lineaire											  */
/* -------------------------------------------------------------------------- */
typedef int (*ptrFonctionResolution)(double *, double *, double *, struct Parametre);

typedef struct Parametre {

	/*  a et b : parametres de la fonction donant l'état initial  */
	/*  ainsi que le second membre  */
	double a;
	double b;

	/*  Espace dans lequel l'equation est resolue : [ax,bx]x[ay,by]  */
	/*  x : abscisse correspond aux colonnes  */
	/*  y : ordonnee correspond aux lignes  */
	double ax;
	double bx;
	double ay;
	double by;

	/*  Discretisation en espace : M pour x et N pour y  */
	/*  hx : pas de discretisation pour l'abscisse, hy : pas de discretisation pour l'ordonnee*/
	int N;
	int Np2;
	int M;
	int Mp2;
	double hx;
	double hy;

	/*  Discretisation en temps :  */
	/*  Tmax : temps maximal  */
	/*  T : nombre de pas de temps  */
	/*  ht : pas de la discretisation en temps*/
	double Tmax;
	int T;
	double ht;

	/*  parametre de l'equation devant le laplacien  */
	double nu;

	/*  parametres stockes pour etre calcules qu'une seule fois  */
	double lambda;
	double betaX;
	double betaY;
	double invLambda;
	double betaXDivLambda;
	double betaYDivLambda;
	double lambdaCrankNicolson;
	double invHt;

	/*  Methode de resolution de l'equation lineaire :  */
	/*  eps : precision a atteindre pour la resolution du systeme lineaire  */
	/*  maxIter : nombre d'iteration maximal pour la resolution du systeme lineaire  */
	/*  resolutionSysLin : correspond a la methode utilisee pour la resolution du systeme lineaire  */
	/*  					1 : Jacobi, 2 : Gauss-Seidel  */
	/*  pResolution : pointeur sur la fonction implementant la methode choisie  */
	double eps;
	int maxIter;
	int resolutionSysLin;
	ptrFonctionResolution pResolution;

	/*  Schema de discretisation :  */
	/*  schema : correspond au schema utilise pour la discretisation  */
	/*  			1 : euler implicite, 2 : Cranck-Nocolson  */
	/*  pSchema : pointeur sur la fonction implementant le schema choisie  */
	int schema;
	ptrFonctionSchema pSchema;

	/* Nombre de Thread a utiliser pour paralleliser la résolution du systeme lineaire */
	int nbThread;
} Parametre;

/* -------------------------------------------------------------------------- */
/* resolutionEquation : fonction qui resout et affiche l'equation  */
/* 						a chaque pas de temps   */
/* Entree : par : structure Parametre donnant tous les paramètres  */ 
/* 						necessaire a la resolution  */
/* Sortie : void  */
/* -------------------------------------------------------------------------- */
void resolutionEquation(Parametre par);

/* -------------------------------------------------------------------------- */
/* jacobi : Resout le systeme lineaire selon la methode de jacobi			  */
/* 			Calcul la limite de la sute des u_km1							  */
/* Entree : double * u : pointeur sur la solution au pas de temps precedent	  */
/* 			double * f : pointeur sur le second membre de jacobi deja calcule */
/*						 pour tout l'espace considere et au temps dont on va  */
/* 						 calculer la solution dans cette fonction			  */
/* 			double * u_km1 : pointeur sur la solution a l'iteration de		  */
/* 						jacobi precedent. Il est donne en parametre			  */
/* 						pour l'allouer qu'une seul fois						  */
/* 			par : structure Parametre donnant tous les paramètres			  */
/* 						necessaire a la resolution							  */
/* Sortie : int iter : nombre d'iteration realiser avant					  */
/*  					d'atteindre la precesion demander					  */
/* -------------------------------------------------------------------------- */	
int jacobi(double * u, double * f, double * u_km1, Parametre p);

/* -------------------------------------------------------------------------- */
/* gauss_Seidel : Resout le systeme lineaire selon la methode de gauss_Seidel */
/* 			Calcul la limite de la sute des u_km1							  */
/* Entree : double * u : pointeur sur la solution au pas de temps precedent	  */
/* 			double * f : pointeur sur le second membre de jacobi deja calcule */
/*						 pour tout l'espace considere et au temps dont on va  */
/* 						 calculer la solution dans cette fonction			  */
/* 			double * u_km1 : pointeur sur la solution a l'iteration de		  */
/* 						jacobi precedent. Il est donne en parametre			  */
/* 						pour l'allouer qu'une seul fois						  */
/* 			par : structure Parametre donnant tous les paramètres			  */
/* 						necessaire a la resolution							  */
/* Sortie : int iter : nombre d'iteration realiser avant					  */
/*  					d'atteindre la precesion demander					  */
/* -------------------------------------------------------------------------- */
int gauss_Seidel(double * u, double * f, double * u_km1, Parametre p);

/* -------------------------------------------------------------------------- */
/* eulerImplicite : Renvoit la valeur du second membre pour la resolution du  */
/* 					systeme lineaire selon le schema euler implicite		  */
/* Entree : double * f : pointeur sur le second membre a completer			  */
/* 							la fonction remplit cette variable				  */
/* 			double * u_nm1 : pointeur sur la solution au pas de temps		  */
/* 							precedent										  */
/* 			par : structure Parametre donnant tous les paramètres			  */
/* 						necessaire a la resolution							  */
/* Sortie : void															  */
/* -------------------------------------------------------------------------- */
void eulerImplicite(double * f, double * u_nm1, double t, Parametre p);

/* -------------------------------------------------------------------------- */
/* cranck_Nicolson : Renvoit la valeur du second membre pour la resolution du */
/* 					systeme lineaire selon le schema de cranck_Nicolson		  */
/* Entree : double * f : pointeur sur le second membre a completer			  */
/* 							la fonction remplit cette variable				  */
/* 			double * u_nm1 : pointeur sur la solution au pas de temps		  */
/* 							precedent										  */
/* 			par : structure Parametre donnant tous les paramètres			  */
/* 						necessaire a la resolution							  */
/* Sortie : void															  */
/* -------------------------------------------------------------------------- */
void cranck_Nicolson(double * f, double * u_nm1, double t, Parametre p);

/* -------------------------------------------------------------------------- */
/* source : renvoie la valeur de la source (second membre) de l'equation	  */
/* 			differencielle													  */
/* Entree : int x : valeur de l'abscisse ou l'on souhaite calculer sa valeur  */
/* 			int y : valeur de l'ordonnee ou l'on souhaite calculer sa valeur  */
/* 			double t : valeur du temps ou l'on souhaite calculer sa valeur	  */
/* Sortie : double : valeur de la source au point donne pour le temps donne	  */
/* -------------------------------------------------------------------------- */
double source(double x, double y, double t, Parametre p);

/* -------------------------------------------------------------------------- */
/* initialisation : initialise la solution selon la condition initiale de	  */
/* 					l'equation differentielle								  */
/* Entree : double * u : pointeur sur la solution							  */
/* 							tableau que l'on veut initialiser				  */
/* Sortie : void															  */
/* -------------------------------------------------------------------------- */
void initialisation(double * u, Parametre p);

/* -------------------------------------------------------------------------- */
/* afficher : affiche la solution dans la console							  */
/* Entree : double * u : pointeur sur la solution a afficher				  */
/* 			par : structure Parametre donnant tous les paramètres			  */
/* 						necessaire a la resolution							  */
/* 			int nbIter : nombre d'iteration de la methode de resolution		  */
/* 						du systeme lineaire									  */
/* 			double t : temps a lequelle correspond la solution a afficher	  */
/* 						necessaire a la resolution							  */
/* Sortie : void															  */
/* -------------------------------------------------------------------------- */
void afficher(double * u, Parametre p, int nbIter, double t);

/* -------------------------------------------------------------------------- */
/* copier : copie le tableau source dans le tableau destiantion				  */
/* Entree : double * u : tableau a copier, c'est la source					  */
/* 			double * dest : tableau a modifier, c'est la destination		  */
/* 			par : structure Parametre donnant la taille des tableaux		  */
/* Sortie : void															  */
/* -------------------------------------------------------------------------- */
void copier(double * u, double * dest, Parametre p);

/* -------------------------------------------------------------------------- */
/* echangerPointeur : echange le contenue de deux pointeurs					  */
/* Entree : double ** u : pointeur sur le pointeur du premier tableau		  */
/* 			double ** v : pointeur sur le pointeur du deuxieme tableau		  */
/* Sortie : void															  */
/* -------------------------------------------------------------------------- */
void echangerPointeur(double ** u, double ** v);

/* -------------------------------------------------------------------------- */
/* fonctionCL : fonction qui donne la valeur de la condition aux limites au   */
/* 				point donner												  */
/* Entree : int x : abscisse du point a evaluer								  */
/* Entree : int y : ordonnee du point a evaluer								  */
/* Sortie : double valeur : la valeur de la condition aux limites			  */
/* -------------------------------------------------------------------------- */
double fonctionCL(double x, double y, Parametre p);

/* -------------------------------------------------------------------------- */
/* fonctionCI : fonction qui donne la valeur de la condition initial au		   */
/* 				point donner												  */
/* Entree : int x : abscisse du point a evaluer								  */
/* Entree : int y : ordonnee du point a evaluer								  */
/* Sortie : double valeur : la valeur de la condition initial				  */
/* -------------------------------------------------------------------------- */
double fonctionCI(double x, double y, Parametre p);


#endif