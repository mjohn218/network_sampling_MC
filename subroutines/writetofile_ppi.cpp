#include "pro_classes.h"
#include "write_ppis.h"

void writetofile_ppi(ofstream &filename, int nwhole, Protein *wholep)
{
  int i, j;
  for(i=0;i<nwhole;i++){
    filename <<i<<'\t'<<wholep[i].ninterface<<'\t';
    for(j=0;j<wholep[i].ninterface;j++)
      filename <<wholep[i].valiface[j]<<'\t';
    filename <<endl;
  }
}
