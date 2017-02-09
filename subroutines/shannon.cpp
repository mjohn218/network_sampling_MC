#include "shannon.h"
#include <cmath>

double shannon(int Nbeads, int *freqhist, int Nprot, int Nsite, double *probH)
{
  int i;
  double entropy=0.0;
  double p;
  for(i=0;i<Nbeads;i++)
    {
      p=freqhist[i]*1.0/(Nprot*Nsite*1.0);
      if(p!=0)
	entropy-=p*log(p/probH[i]);
    }
  return entropy;

}
