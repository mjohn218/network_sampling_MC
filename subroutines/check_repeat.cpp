#include "net_props.h"
#include <cstdlib>

int check_repeat(int Nmer, double **fourmer)
{
  /*check whether the most recently added fourmer is a repeat*/
  int i, j, k, l;
  int a=fourmer[Nmer-1][0];
  int b=fourmer[Nmer-1][1];
  int c=fourmer[Nmer-1][2];
  int d=fourmer[Nmer-1][3];
  int flag=0;
  for(i=0;i<Nmer-1;i++){
    if(fourmer[i][0]==a){
      if(fourmer[i][1]==b){
	if(fourmer[i][2]==c){
	  if(fourmer[i][3]==d){
	    flag=1;
	    ////cout <<"repeated fourmer number: "<<Nmer-1<<" indexes: "<<a<<' '<<b<<' '<<c<<' '<<d<<endl;
	    i=Nmer-1;
	  }
	}
      }
    }
  }
 
    
  return flag;      
  

}
