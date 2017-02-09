#include "pro_classes.h"
#include "utility_calls.h"
#include "network_metrics2.h"
#include "matmultiply.h"

void locate_interfaces(int p1, int p2, Protein *wholep, int *numpartners, int **Speclist, int &p1if, int &p2if)
{
  int s, t, k, l;
  int i1, i2;

  for(t=0;t<wholep[p1].ninterface;t++){
    i1=wholep[p1].valiface[t];
    for(s=0;s<numpartners[i1];s++){
      for(k=0;k<wholep[p2].ninterface;k++){
	i2=wholep[p2].valiface[k];
	
	if(Speclist[i1][s]==i2){
	    
	    p1if=i1;
	    p2if=i2;
	    
	    l=numpartners[i2];
	    k=wholep[p2].ninterface;
	    s=numpartners[i1];
	    t=wholep[p1].ninterface;//break from loop
	}
      }
      
    }
  } //discovered identity of interfaces on our network
}
