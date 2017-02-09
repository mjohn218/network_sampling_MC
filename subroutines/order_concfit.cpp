#include "order_concfit.h"
#include "utility_calls.h"

void order_concfit(int nstep, double *concfit, double *ordercfit)
{
  int i;
  //int nstep_opt=plist.nstep_copt;
  int nstep_opt=nstep;
  long unsigned int index[nstep_opt+1];//=new long unsigned int[nstep_opt+1];
  double value[nstep_opt+1];//=new double[nstep_opt+1];
  
  for(i=0;i<nstep_opt;i++){
    value[i+1]=concfit[i];
  
  }
  indexx(nstep_opt, &value[0], &index[0]);
  int id;
  
  for(i=0;i<nstep_opt;i++){
    id=index[i+1]-1;
    ordercfit[i]=concfit[id];
  
  }

  //delete[] index;
  //delete[] value;

}
