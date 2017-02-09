#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int split_type_self(int nedge, int &flagsep)
{
  int i;
  int ptot=0;
  int b[nedge];
  double tmp;
  for(i=1;i<nedge;i++){
    tmp=binomial(nedge,i);
    b[i]=round(tmp);
    ptot+=b[i];
  } 
  int a=ptot;
  /*If there is a self interaction, also include the split self, which is the binomial
    sum with one fewer edge*/
  int b2[nedge-1];
  b2[0]=2;//there is one way to keep all and just break off the self
  ptot+=b2[0];
  for(i=1;i<nedge-1;i++){
    tmp=binomial(nedge-1,i);
    b2[i]=round(tmp);
    ptot+=b2[i];
  }
  /*just keep full distribution*/
  //ptot*=0.5;//total number of choices
  double rnum=ptot*1.0*rand_gsl();
  int w=ceil(rnum);//index between 1 and ptot, never zero 
  int sum;
  int t;
  int nswap;
  flagsep=0;
  if(w>a){
    //this is one of the self interactions
    flagsep=1;
    sum=a+b2[0];
    t=0;
    while(sum<w){
      t++;
      sum+=b2[t];
    }
    //nswap can be zero, it means splitting off the self interaction
    nswap=t;
  
    if(t>(nedge-1)/2)nswap=t-(nedge-1)/2;
  
  }else{
    //just a regular splitting
    //now figure out which split this corresponds to 
    t=1;
    sum=b[1];
    while(sum<w){
      t++;
      sum+=b[t];
    }
    nswap=t;
    if(t>nedge/2)nswap=t-nedge/2;
  
  }
  //  //cout <<"ptot: "<<ptot<<" w: "<<w<<" nedge: "<<nedge<<" n to swap: "<<t<<" what is nedge/2? "<<nedge/2<<endl;
  return nswap;
}
