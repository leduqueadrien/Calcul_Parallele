#ifndef __MON_FONCTIONS_H__
#define __MON_FONCTIONS_H__

/* Definition du domaine de calcul */
#define  XMIN  -5.0
#define  XMAX  5.0
#define  YMIN  -4.0
#define  YMAX  4.0

/* Erreur d'aproximation autorisee */
#define  EPSILON 0.1

/* Calcul de la vitesse de convergence de la methode pour un point donne */
/* Valeur de retour : 
 *       le nombre d'iterations necessaire pour converger a une solution
 */
unsigned int couleur1(complexe_t);

/* Calcul vers quelle solution converge la methode,  pour un point donne */
/* Valeur de retour : 
 *       le numero de la solution
 */
unsigned int couleur3(complexe_t);

#endif
