#include "enum_calls.h"
#include <cmath>
#include <stdlib.h>

using namespace std;

int identify_network(int **Epred, int **Emat)
{
  int sep=0;
  int index=0;
  int t=0;
  int i, j;
  //protein 0
  sep=0;//0,1  0, 2 1, 2
//   cout <<"new network "<<endl;
//   for(i=0;i<3;i++){
//     for(j=0;j<3;j++){
//       cout <<Epred[1][i*EDIM+j]<<'\t';
//     }
//     cout <<endl;
//   }
//   cout <<"1, 2, 12, 11: "<<Epred[1][1]<<' '<<Epred[1][2]<<' '<<Epred[1][12]<<' '<<Epred[1][11]<<endl;
  sep+=abs(Epred[0][1]-Emat[0][0]);
  sep+=abs(Epred[0][2]-Emat[0][1]);
  sep+=abs(Epred[0][12]-Emat[0][2]);

  //check protein 0
  while(sep!=0){
    index+=22;
    t++;
    sep=0;
    sep+=abs(Epred[0][1]-Emat[t][0]);
    sep+=abs(Epred[0][2]-Emat[t][1]);
    sep+=abs(Epred[0][12]-Emat[t][2]);
  }
  //now check protein 1
  //  cout <<"current index "<<index<<" t: "<<t<<endl;
  sep=0;
  sep+=abs(Epred[1][1]-Emat[5][0]);
  sep+=abs(Epred[1][2]-Emat[5][1]);
  sep+=abs(Epred[1][12]-Emat[5][2]);
  sep+=abs(Epred[1][11]-Emat[5][3]);
  t=5;
  
  while(sep!=0){
    index+=2;
    t++;
    sep=0;
    sep+=abs(Epred[1][1]-Emat[t][0]);
    sep+=abs(Epred[1][2]-Emat[t][1]);
    sep+=abs(Epred[1][12]-Emat[t][2]);
    sep+=abs(Epred[1][11]-Emat[t][3]);
  }
  
  //now check protein 2
  sep=0;
  sep+=abs(Epred[2][1]-Emat[16][0]);
  if(sep!=0)
    index+=1;
  

  return index;
}
