#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int split_type(int nedge)
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
  /*just keep full distribution*/
  //ptot*=0.5;//total number of choices
  double rnum=ptot*1.0*rand_gsl();
  int w=ceil(rnum);//index between 1 and ptot, no zero!
  //now figure out which split this corresponds to 
  int t=1;
  int sum=b[1];
  while(sum<w){
    t++;
    sum+=b[t];
  }
  int nswap=t;
  
  if(t>nedge/2)nswap=t-nedge/2;
  return nswap;
}
