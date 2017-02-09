#include "pro_classes.h"
#include "write_ppis.h"

void write_speclist(int ninterface, int *numpartners, int **Speclist )
{
  int i, j;
  for(i=0;i<ninterface;i++){
    cout <<i<<'\t'<<numpartners[i]<<'\t';
    for(j=0;j<numpartners[i];j++)
      cout <<Speclist[i][j]<<'\t';
    cout <<endl;
  }
  
}
