#include "enum_calls.h"
#include <fstream>
#include <iostream>

using namespace std;
void write_enum_network(int **Emat)
{
  int index=0;
  int t=0;
  int t2;
  int i, j, k;

  //loop over protein 0
  for(i=0;i<5;i++){
    //now loop over protein 1
    for(j=0;j<11;j++){
      t=j+5;
      //now loop over protein 2
      for(k=0;k<2;k++){
	t2=k+16;
	cout <<"index: "<<index<<" p0: "<<Emat[i][0]<<' '<<Emat[i][1]<<' '<<Emat[i][2]<<" p1: ";
	cout <<Emat[t][0]<<' '<<Emat[t][1]<<' '<<Emat[t][2]<<' '<<Emat[t][3]<<" p2: ";
	cout <<Emat[t2][0]<<endl;
	index++;
      }
    }
  }
    //now check protein 1
  

}
