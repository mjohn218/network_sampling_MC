#include "net_props.h"
#include <cstdlib>

void cluster_cof(int N, int &open, int &closed, double &avg, int  **Speclist, int *numpartner)
{
  open=0;
  closed=0;
  avg=0;
  int myopen[N];
  int myclose[N];
  double cof[N];
  int i,j, k;
  int ind;
  int t;
  double factor;
  for(i=0;i<N;i++){
    myopen[i]=0;
    myclose[i]=0;
  }
  //cout <<"In clustering coefficient "<<endl;
  //cout <<"protein "<<" coefficient "<<" nopen "<<" nclosed "<<endl;
  for(i=0;i<N;i++){
    cof[i]=0;
    if(numpartner[i]<=1){
      cof[i]=0;
    }else{
      //test for clustering
      for(j=0;j<numpartner[i];j++){
	t=0;
	ind=Speclist[i][j];
	
	while(t<numpartner[ind]){
	  
	  for(k=j+1;k<numpartner[i];k++){
	    if(Speclist[ind][t]==Speclist[i][k]){
	      cof[i]+=1;
	      closed+=1;
	      myclose[i]+=1;
	    }
	  }
	  
	  t++;
	}
      }
      factor=numpartner[i]*1.0*(numpartner[i]-1)/2.0;
      myopen[i]=factor-myclose[i];
      open+=myopen[i];
    }
			      
    if(numpartner[i]>1)
      cof[i]/=factor;
    //cout <<"protein: "<<i<<' '<<cof[i]<< ' '<<myopen[i]<<' '<<myclose[i]<<endl;
    avg+=cof[i];
  }
    //In another code which calls this function, avg is divided by N. And closed is divided by 3 to get Ntriangles. Cglobal = closed*1.0/(1.0*(open+closed). Clocal=avg/(1.0*N)
}
