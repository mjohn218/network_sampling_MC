#include "pro_classes.h"
#include "write_ppis.h"

void write_ppi(int nwhole, Protein *wholep)
{
  int i, j;
  for(i=0;i<nwhole;i++){
    cout <<i<<'\t'<<wholep[i].ninterface<<'\t';
    for(j=0;j<wholep[i].ninterface;j++)
      cout <<wholep[i].valiface[j]<<'\t';
    cout <<endl;
  }
}
