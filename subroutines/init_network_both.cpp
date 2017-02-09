#include "pro_classes.h"
#include "constrainParms.h"
#include "init_iin_net.h"

int init_network_both(int nwhole, ppidata *ppi, Protein *wholep, int *p_home, constrainParms &plist, int Nedge, int **Speclist, int *numpartners, int *e1num, int *e2num, double *abund)
{
  /*given a PPI network, create a IIN with a specific number of interfaces
    if abundance is known, give each protein 1 interface, if it's not known
    give protein new interface for each edge.
  */

  int i, j;
  int nbor[nwhole];
  int ind;
  int id;
  int icurr;
  double highratio;
  int chgprot;
  int num2, t, num, n2;
  int ncurr=nwhole;
  
  char fname[100];
  

  int ntry;
  
  cout <<"creating separate interfaces per edge: "<<endl;
  t=0;
  int nint;
  for(i=0;i<nwhole;i++){
    if(abund[i]<0){
      //give this protein a separate interface for each
      //partner
      wholep[i].ninterface=ppi[i].nppartner;
      nint=ppi[i].nppartner;
      for(j=0;j<nint;j++){
	wholep[i].valiface[j]=t;
	p_home[t]=i;
	t++;
      }
    }else{
      //give this protein only a single interface
      wholep[i].ninterface=1;
      wholep[i].valiface[0]=t;
      p_home[t]=i;
      t++;
    }
  }
  int Nif=t;
  cout <<"enumerated interfaces: "<<t<<endl; 
  //  double *A=new double[Nif*Nedge];
  int e1int[Nedge];
  int e2int[Nedge];
  /*Now need to create the speclist and numpartners array*/
  int myind[nwhole];
  for(i=0;i<nwhole;i++)
    myind[i]=0;
  
  int p1, p2;
  int i1, i2;
  for(i=0;i<Nedge;i++){
    p1=e1num[i];
    p2=e2num[i];
    i1=wholep[p1].valiface[myind[p1]];
    i2=wholep[p2].valiface[myind[p2]];
    e1int[i]=i1;
    e2int[i]=i2;
    if(p1==p2){
      if(abund[p1]<0)
	myind[p1]++;
      
    }
    else{
      if(abund[p1]<0){
	myind[p1]++;
      }
      if(abund[p2]<0){
	myind[p2]++;
      }
    }
  }

  //should be able to use e1int to create those matrices
  for(i=0;i<Nif;i++){
    numpartners[i]=0;
  }

  for(i=0;i<Nedge;i++){
    i1=e1int[i];
    i2=e2int[i];
    Speclist[i1][numpartners[i1]]=i2;
    Speclist[i2][numpartners[i2]]=i1;
    
    numpartners[i1]++;
    if(i2!=i1)
      numpartners[i2]++;
  }
 


  return Nif;
  
} 
