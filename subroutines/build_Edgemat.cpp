#include "pro_classes.h"
#include "build_Edgemat.h"


void build_Edgemat(int Nedge, int nwhole, int **Edgemat, int *numpartner, int **Speclist, Protein *wholep, int *phome, int **ehome, ppidata *ppi)
{
  /*Create the edge matrix from the interface network
    the edges are numbered just as the protein pairs are numbered
  */
  int i, j;
  for(i=0;i<nwhole;i++){
    for(j=0;j<EDIM*EDIM;j++)
      Edgemat[i][j]=0;
  }
  int nedgep, nint;
  int nshare, p2, e1, e2;
  int s, i1, i2, t;
  //  cout <<"In Edgemat "<<endl;
  for(i=0;i<nwhole;i++){
    nint=wholep[i].ninterface;
    nedgep=ppi[i].nppartner;
    //cout <<"pro; "<<i<<" nint: "<<nint<<" nedge: "<<nedgep<<endl;
    //    if(nint<nedgep){
      //there is some edge sharing
      for(j=0;j<nint;j++){
	i1=wholep[i].valiface[j];
	nshare=numpartner[i1];
	//cout <<"interface: "<<j<<endl;
	for(t=0;t<nshare;t++){
	  i2=Speclist[i1][t];
	  p2=phome[i2];
	  e1=ehome[i][p2];
	  if(i1==i2){
	    Edgemat[i][e1*EDIM+e1]=1;//self, same interface!
	  }
	  //cout <<"in Edgemat, i1: "<<i1<<" binds i2: "<<i2<<" edge: "<<e1<<endl; 
	  for(s=t+1;s<nshare;s++){
	    i2=Speclist[i1][s];
	    p2=phome[i2];
	    e2=ehome[i][p2];
	    //cout <<"share: "<<i<<" and : "<<p2<<" edge: "<<e2<<endl; 
	    Edgemat[i][e1*EDIM+e2]=1;
	    Edgemat[i][e2*EDIM+e1]=1;
	  }
	}
	
      }
      //    }//done checking for sharing on this protein
  
  }
  //    cout <<"should be one if self-self: "<<Edgemat[1][11]<<endl;
}
