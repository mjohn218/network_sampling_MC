#include "pro_classes.h"
#include "constrainParms.h"
#include "create_network.h"
#include "calc_ratios.h"
#include "utility_calls.h"

void adj_make(int nwhole, ppidata *ppi, double *Adj, ofstream &Mfile)
{
  int npart=0;
  int t, ig;
  int i, j;
  int ind;
  for(i=0;i<nwhole*nwhole;i++)
    Adj[i]=0.0;


  for(i=0;i<nwhole;i++){
    for(j=0;j<ppi[i].nppartner;j++){
      ind=ppi[i].pplist[j];
      Adj[i*nwhole+ind]=1;//i is the column, j the row, Symmetric matrix
    }
  }
  for(i=0;i<nwhole;i++){
    for(j=0;j<nwhole;j++){
      Mfile<<Adj[i*nwhole+j]<<'\t';
    }
    Mfile<<endl;
  }

}
