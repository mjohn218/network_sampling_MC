#include "net_props.h"
#include <cstdlib>

void neighbor_conn(int N, int *numpartner, int **Speclist, double  *conn, double *ctot, double &totavg, double &std, double &normavg, double &normstd)
{
  int i, j;
  totavg=0;
  double tot2=0;
  normavg=0;
  double m2=0;
  for(i=0;i<N;i++){
    conn[i]=0;
    for(j=0;j<numpartner[i];j++){
      conn[i]+=numpartner[Speclist[i][j]];
    }
    ctot[i]=conn[i];
    if(numpartner[i]>0)
      conn[i]/=(1.0*numpartner[i]);
    //cout <<"Protein: "<<i<<" numpartner: "<<numpartner[i]<<" average partner connectivity: "<<conn[i]<<" Total connectivity: "<<ctot[i]<<endl;
    totavg+=ctot[i];
    normavg+=conn[i];
    m2+=conn[i]*conn[i];
    tot2+=ctot[i]*ctot[i];
  }
  totavg/=(1.0*N);
  normavg/=(1.0*N);
  double var=tot2/(1.0*N)-totavg*totavg;
  double v2=m2/(1.0*N)-normavg*normavg;
  std=sqrt(var);
  normstd=sqrt(v2);
  //cout <<"Average total partners: "<<totavg<<" var: "<<var<<endl;
}
