#include "pro_classes.h"
#include "calc_ratios.h"
#include "utility_calls.h"

void calc_ratio_mc(int ninterface, int **Speclist, int *numpartners, int *nbor)
{

  
  int i, j;
  int ind;
  for(i=0;i<ninterface;i++){
    //now calculate neighbors
    nbor[i]=0;
    for(j=0;j<numpartners[i];j++){
      ind=Speclist[i][j];
      nbor[i]+=numpartners[ind];
    }
    //    cout <<"i: "<<i<<" partners: "<<numpartners[i]<<" nbor:  "<< nbor[i]<<endl;
    
  }


  
}
