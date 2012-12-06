#include <stdio.h>
#include <complex.h>
#include "poly.h"

#define TRUE   1
#define FALSE  0

/* Calcula multiplicidade da raiz alpha com erro tao no polinomio */
/* p e armazena em Roots (fala tb se pode ser um cluster          */

void multiplicity (polynomial *p, double alpha, double tao, int ni, root **Roots) {

  int m = 0;

  do {

    lower_degree (&p);
    ++m;
    
  } while (p->order != 0 && cabs(calc_pol(*p, alpha)) > tao);

    (*Roots)->ni = ni;
    (*Roots)->mult = m;
    (*Roots)->value = alpha;
    ++(*Roots);
 
}

/* Acha anel contendo todas as raizes do polinomio */

long double find_ring (polynomial p, complex center, long double r) {

  for (r *= 1/2; has_roots (p, center, r); r *= 1/2); 
    
  return r;

}

/* UGH!!! (No comments...) */

root *lehmer (polynomial p, long double tao) {

  int jane;
  long double S, R, r, mu;
  complex alpha, beta, vi;
  root *Roots;

  Roots = (root *) malloc (p.order*sizeof (p.order));

  printf ("\nNão foi aqui!\n");

  S = r_all_zeros (p);

  printf ("\n%lf\n", S);

  if (cabsl(p.coeff[0]) <= tao) {
    multiplicity (&p, 0, tao, 0, &Roots);
  }

  printf ("\nNão foi aqui!\n");

  while (p.order != 0) {

    jane = TRUE;
    alpha = 0;
    r = S;
    printf ("\n%lf\n", r);
    
    do {
      R = find_ring (p, alpha, r);

      printf ("\nNão foi aqui e R = %lf!\n", R);

      
      if (alpha == 0) {
	S = R;
      } else if (cabsl (calc_pol(p, alpha)) > tao) {
	multiplicity (&p, alpha, tao, 1, &Roots);
	break;
      }

      vi = 5/3*R;
      mu = vi/2;

      do {
	beta = alpha + vi;

	if (alpha == beta) {
	  printf ("\nNão foi aqui!\n");
	  multiplicity (&p, alpha, tao, 1, &Roots);
	  jane = FALSE;
	  break;
	}


	if (has_roots (p, beta, mu)) {
	  alpha = beta;
	  r = mu/2;
	  break;
	} else {
	  vi = csqrtl (I*vi);
	}

      } while (1);

    } while (jane);
    
  }
  
}

int main (int argc, char **argv) {

  polynomial p;
  root *Roots;
  p.order = 2;
  p.coeff = (long double complex *) malloc ((p.order+1)*sizeof(long double complex));

  p.coeff[0] = 1;
  p.coeff[1] = 0;
  p.coeff[2] = 1;

  long double ad = 16;

  //print_poly (p);
  printf ("\n%lf\n", ad /*r_all_zeros(p)*/);

  //  Roots = lehmer (p, .00000001);

  return 0;
}
 
