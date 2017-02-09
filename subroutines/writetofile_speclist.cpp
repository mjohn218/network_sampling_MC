#include "pro_classes.h"
#include "write_ppis.h"

void writetofile_speclist(ofstream &filename, int ninterface, int *numpartners, int **Speclist )
{
  int i, j;
  for(i=0;i<ninterface;i++){
    filename <<i<<'\t'<<numpartners[i]<<'\t';
    for(j=0;j<numpartners[i];j++)
      filename <<Speclist[i][j]<<'\t';
    filename <<endl;
  }
  
}
