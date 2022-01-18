#include "complexe.h"
#include "fonctions.h"

/* --------------------------------------------------------- */
/* couleur1  : calcule la vitesse de convergence de la       */
/*             methode de newton pour chaque x0 initial      */
/*             se trouvant dans le domaine de definition     */
/*   entrees : l'adresse de "x=x0"                           */
/*   sortie  : couleur de "x" correspondant au nombre        */ 
/*             d'iterations necessaire pour avoir une        */
/*             solution suffisament approchee                */
/* --------------------------------------------------------- */
unsigned int couleur1(complexe_t x)
{
    complexe_t   x2, x3, xd;
    double       err;
    unsigned int iter=0;
 
     do {
      mult_comp(x, x, &x2);
      mult_comp(x, x2, &x3);
      x3.re = x3.re - 1;
      err = norme2(x3);      /* l'ecart entre x et la solution */

      xd.re = 1.0;  xd.im = 0.0;
      division(xd, x2, &x3);
      mult_scal(2./3., x, &x2);   /* x = 2x/3 + 1/3x2 */
      mult_scal(1./3., x3, &xd);
      somme(x2, xd, &x);
      iter++;
     } while (err > EPSILON);
     return (iter);
}

/* --------------------------------------------------------- */
/* couleur3  : calcule vers quelle racine converge la        */
/*             methode de newton pour chaque x0 initial      */
/*             se trouvant dans le domaine de definition     */

/*   entrees : l'adresse de "x=x0"                           */
/*   sortie  : numero de la racine approchee par la methode  */ 
/*             r0 = 1; r1 = (-1 + sqrt(3)i) / 2;             */
/*             r2 = (-1 - sqrt(3)i) / 2.                     */
/* --------------------------------------------------------- */
unsigned int couleur3(complexe_t x)
{
  complexe_t   x2, x3, xd;
  double       err;
  unsigned int iter = 0;

  do {
    mult_comp(x, x, &x2);
    mult_comp(x, x2, &x3);
    x3.re = x3.re - 1;
    err = norme2(x3); 
    
    xd.re = 1.0;  xd.im = 0.0;
    division(xd, x2, &x3);      /* x3 = 1/x^2 */
    mult_scal(2./3., x, &x2);   /* x2 = 2x/3  */
    mult_scal(1./3., x3, &xd);  /* xd = 1/3x^2 */
    somme(x2, xd, &x);          /* x = 2x/3 + 1/3x^2 */
    iter++;
  } while (err > EPSILON);

  if (x.re >0) return (0);
  else if (x.im > 0) return (1);
  else return (2);
}


