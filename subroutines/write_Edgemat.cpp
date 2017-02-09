#include "pro_classes.h"
#include "write_ppis.h"

void write_Edgemat(int nwhole, int **Edgemat, ppidata *ppi, ofstream &efile)
{
  int i, j, m;
  int nedge;
  for(i=0;i<nwhole;i++){
    nedge=ppi[i].nppartner;
    efile<< "p: "<<i<<endl;
    for(j=0;j<nedge;j++){
      for(m=0;m<nedge;m++){
	efile <<Edgemat[i][j*EDIM+m]<<'\t';
      }
      efile<<endl;
    }
  }


}
