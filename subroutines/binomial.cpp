#include "binomial.h"

double binomial(int N, int k1)
{
  double r1=lgamma_NM(N+1);
  double r2=lgamma_NM(k1+1);
  double r3=lgamma_NM(N-k1+1);

  double val=exp(r1-r2-r3);
  return val;
  
}
