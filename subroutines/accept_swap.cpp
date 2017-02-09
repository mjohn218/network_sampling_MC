#include "pro_classes.h"
#include "accept_iin_moves.h"

void accept_swap(double newfit, double &oldfit, int ninterface, int *numpartners, int **Speclist, int *tmppartners, int **templist)
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

} 
