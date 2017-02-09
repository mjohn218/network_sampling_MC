#include "pro_classes.h"
#include "calc_ratios.h"
#include "utility_calls.h"

double calc_avg_ratio(int ninterface, int **Speclist, int *numpartners)
{

  int i, j;
  int ind;
  int nbor[ninterface];
  double avgratio=0;
  double r;
  for(i=0;i<ninterface;i++){
    //now calculate neighbors
    nbor[i]=0;
    for(j=0;j<numpartners[i];j++){
      ind=Speclist[i][j];
      nbor[i]+=numpartners[ind];
    }
    //    cout <<"i: "<<i<<" partners: "<<numpartners[i]<<" nbor:  "<< nbor[i]<<endl;
    if(numpartners[i]>1)
      r=nbor[i]*1.0/(numpartners[i]*1.0);
    else
      r=1;
    avgratio+=r;
  }
  avgratio/=(1.0*ninterface);
  return avgratio;
  
}
