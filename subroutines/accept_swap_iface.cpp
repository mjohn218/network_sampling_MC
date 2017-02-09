#include "pro_classes.h"
#include "accept_iin_moves.h"

void accept_swap_iface(double newfit, double &oldfit, int ninterface, int *numpartners, int **Speclist, int *tmppartners, int **templist, int nwhole, Protein *wholep, Protein *wholetemp, int *p_home)
{
  //accept move
  oldfit=newfit;
  int i, j, n;
  for(i=0;i<ninterface;i++){
    n=tmppartners[i];
    numpartners[i]=n;
    for(j=0;j<n;j++)
      Speclist[i][j]=templist[i][j];
  }
  
  for(i=0;i<nwhole;i++){
    wholep[i].ninterface=wholetemp[i].ninterface;
    for(j=0;j<wholep[i].ninterface;j++){
      n=wholetemp[i].valiface[j];
      p_home[n]=i;
      wholep[i].valiface[j]=wholetemp[i].valiface[j];
    }
  }
  

} 
