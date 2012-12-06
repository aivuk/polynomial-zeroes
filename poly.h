#include <complex.h>

typedef struct {
   long double complex *coeff;
   int order;
} polynomial;

typedef struct {
  long double complex value;
  int ni;
  int mult;
} root;
