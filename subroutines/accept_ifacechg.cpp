#include "pro_classes.h"
#include "accept_iin_moves.h"

void accept_ifacechg(double newfit, double &oldfit, int ntmpinterface, int &ninterface, int *numpartners, int **Speclist, int *tmppartners, int **templist, int nwhole, Protein *wholep, Protein *wholetemp)
{
  //accept move
  oldfit=newfit;
  ninterface=ntmpinterface;
 
  int i, j, n;
  for(i=0;i<ninterface;i++){
    n=tmppartners[i];
    numpartners[i]=n;
    for(j=0;j<n;j++)
      Speclist[i][j]=templist[i][j];
  }
 
  for(i=0;i<nwhole;i++){
    wholep[i].ninterface=wholetemp[i].ninterface;
    for(j=0;j<wholep[i].ninterface;j++)
      wholep[i].valiface[j]=wholetemp[i].valiface[j];
  }
  

} 
