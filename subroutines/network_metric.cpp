#include "pro_classes.h"
#include "utility_calls.h"
#include "network_metrics2.h"
#include "matmultiply.h"

double network_metric(int nwhole, Protein *wholep, ppidata *ppi, int *numpartners, int **Speclist, Protein *refnet, int *refpartners, int **reflist, double alpha)
{
  /*Compare the current network's structure to that of the target network
    this version is simpler than above, checks if N_interface is correct
   */
  int i, j;
  int p1, p2; 
  int p1if=0;
  int p2if=0;
  int r1if=0;
  int r2if=0;
  double nd1, nd2, wd1, wd2;
  double sum1=0;
  double sum2=0;
  int nedge;
  double wref, wnet;
  //  double half=-0.5;
  //cout <<"absolute value -1/2: "<<abs(half)<<endl;
  for(i=0;i<nwhole;i++){
    /*Loop over all proteins and calculate the difference between weights for each one */
    /*Loop over the edges, since they are conserved between the two networks
     */
    p1=i;
    nedge=ppi[i].nppartner;
    //    cout <<" protein: "<<i<<" npartners: "<<ppi[i].nppartner<<endl;
    for(j=0;j<nedge;j++){
      p2=ppi[i].pplist[j];
      //cout <<j<<" partner: "<<p2<<endl;
      //so check this edge between protein p1 and p2
      wref=1.0/(1.0*refnet[p1].ninterface);
      wnet=1.0/(1.0*wholep[p1].ninterface);
      wd1=wnet-wref;
      //      cout <<"wd1: "<<wd1<<endl;
      // wref=1.0/(1.0*refnet[p2].ninterface);
//       wnet=1.0/(1.0*wholep[p2].ninterface);
//       wd2=wnet-wref;
      //now find which interface owns this edge on each protein, and how many
      //edges it controls.

      locate_interfaces(p1, p2, wholep, numpartners, Speclist, p1if, p2if);
      //cout <<"proteins: "<<p1<<"  "<<p2<<" interface in opt: "<<p1if<<' ' <<p2if<<endl;
      locate_interfaces(p1, p2, refnet, refpartners, reflist, r1if, r2if);
      //cout <<"proteins: "<<p1<<"  "<<p2<<" interface in opt: "<<r1if<<' ' <<r2if<<endl;
      nd1=numpartners[p1if]-refpartners[r1if];
      //nd2=numpartners[p2if]-refpartners[r2if];
      sum1+=abs(wd1);//+abs(wd2);//wd1*wd1+wd2*wd2;
      sum2+=nd1*nd1;//+nd2*nd2;
      //cout <<"edge; "<<j<<" 1/ref(p1).ninterface-1/opt(p1).ninterface: "<<wd1<< " same for p2: "<<wd2<<endl;
      //cout <<"edge; "<<j<<" refpartners(i1)-refpartners(i2): "<<nd1<< " same for p2: "<<nd2<<endl;
    }      	
    
  }
  double fitness=sum1+sum2*alpha;
  return fitness;
  
}
