#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "poly.h"

/* Calcula o valor do polinômio recursivamente p em z */

long double complex r_calc_pol (polynomial p, int order, long double complex z) {

   if (order != p.order)
      return z * r_calc_pol(p, order+1, z) + p.coeff[order];
   else
      return p.coeff[order];   
}

/* Chama funcao recursiva que resolve polinomio */

long double complex calc_pol (polynomial p, long double complex z) { 
  
  return r_calc_pol (p, 0, z); 
 
}

/*Calcula o raio do circulo centrado em zero q tem todos os zeros */

long double r_all_zeros(polynomial p){
  long double radius;
  int m;
  
  radius = cabs(p.coeff[p.order-1]);

  for (m = 2; m <= p.order; m++){
    if (radius < powl(cabs (p.coeff[p.order-m]), 1.0/(long double)(m))) {
      radius = powl (cabs (p.coeff[p.order-m]), 1.0/(long double)(m));
    }
  }

  printf ("%lf", radius);
  return 2.0*radius;
}

/* Abaixa grau do polinomio */

polynomial *lower_degree (polynomial p, long double complex root) {

  int i;

  for (i = 0; i < p.order-1; ++i) {
    p.coeff[p.order - i - 1] = root*p.coeff[p.order - i] + p.coeff[p.order - i - 1];
  }

  (p.coeff)++;
  --(p.order);

  return &p;

}

/* Imprime polinomio */

void print_poly (polynomial p) {

  int i;
  
  printf ("\n(%.20lf + i*%.20lf)", creal(p.coeff[0]), cimag(p.coeff[0]));

  for (i = 1; i <= p.order; ++i) {
    printf (" + (%.20lf + i*%.20lf)*z^%d", creal(p.coeff[i]), cimag(p.coeff[i]), i) ;
  }
  
  printf ("\n");
}

/* Copia polinomio 'a' em 'b' */

void copy_poly (polynomial *a, polynomial *b) {

  int i;

  for (i = 0; i <= a->order; ++i)
    b->coeff[i] = a->coeff[i];

  b->order = a->order;

}

/* (Acho que...) Leva em disco unitario */

polynomial *to_unit_disc (polynomial p, long double complex center, long double radius) {

  int i, j;
  polynomial *paux;
  long double aux;
  paux = (polynomial *) malloc (sizeof(polynomial));
  paux->coeff = (long double complex *) malloc ((p.order+1)*sizeof(long double complex));

  for (j = p.order; j > 0; --j) {
    for (i = 0; i < j; ++i) 
      p.coeff[p.order - 1 - i] = (center*p.coeff[p.order - i] + p.coeff[p.order - i - 1]);
  }

  p.coeff[1] *= radius;
  aux = radius;
  for (i = 2; i <= p.order; ++i){ 
    aux*=radius;
    p.coeff[i] *= aux;
  }

  copy_poly (&p, paux);

  return paux;
}

/* UMA transformada de Schur */

polynomial *one_schur_transform (polynomial p) {
  
  int i, j;
  polynomial *Tp;
  
  Tp = (polynomial *) malloc (sizeof(polynomial));

  Tp->order = p.order-1;
  Tp->coeff = (long double complex *) malloc ((p.order)*sizeof(long double complex));

  for (i = 0; i < p.order; ++i) {
    Tp->coeff[i] =(creall(p.coeff[0]) - I*cimagl(p.coeff[0])) * p.coeff[i] - p.coeff[p.order] * ( creall(p.coeff[p.order-i]) - I*cimagl(p.coeff[p.order-i]));
  } 

  return Tp;
}

/* 'k' transformada de schur */

int schur_transforms (polynomial p, long double *gamma, int k) {

  int i;
  polynomial *Tkp;

  Tkp = (polynomial *) malloc (sizeof(polynomial));
  Tkp->coeff = (long double complex *) malloc ((p.order+1)*sizeof(long double complex));

  copy_poly (&p, Tkp); 

  for (i = 0; i < k; ++i) {

    if (Tkp->order != 0)
      gamma[i] = powl (cabsl (Tkp->coeff[0]), 2) - powl (cabsl (Tkp->coeff[Tkp->order]), 2);  
    else
      gamma[i] = Tkp->coeff[0];

    if (gamma[i] < 0)
      return 1;
    
    Tkp = one_schur_transform (*Tkp);

  }

  return 0;
}

/* Verifica se p possui raizes no disco com centro em */
/* center e raio radius                               */

int has_roots (polynomial p, long double complex center, long double  radius)  {

  polynomial *g, *f;
  long double *gamma;
  gamma = (long double *) malloc (p.order*sizeof(long double));

  g = to_unit_disc (p, center, radius);
  
  if (schur_transforms (*g, gamma, g->order)) {
    return 1;
  }

  return 0;
}
