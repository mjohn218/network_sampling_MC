#include "enum_calls.h"

void enumerate_network(int **Emat)
{

      //create all 5 networks for protein 0
      Emat[0][0]=0;
      Emat[0][1]=0;
      Emat[0][2]=0;
      
      Emat[1][0]=0;
      Emat[1][1]=1;
      Emat[1][2]=0;
      
      Emat[2][0]=0;
      Emat[2][1]=0;
      Emat[2][2]=1;

      Emat[3][0]=1;
      Emat[3][1]=0;
      Emat[3][2]=0;
      
      Emat[4][0]=1;
      Emat[4][1]=1;
      Emat[4][2]=1;
      
      //protein 1
      Emat[5][0]=0;
      Emat[5][1]=0;
      Emat[5][2]=0;
      Emat[5][3]=1;
      
      Emat[6][0]=0;
      Emat[6][1]=1;
      Emat[6][2]=0;
      Emat[6][3]=1;
      
      Emat[7][0]=0;
      Emat[7][1]=0;
      Emat[7][2]=1;
      Emat[7][3]=1;

      Emat[8][0]=1;
      Emat[8][1]=0;
      Emat[8][2]=0;
      Emat[8][3]=1;
      
      Emat[9][0]=1;
      Emat[9][1]=1;
      Emat[9][2]=1;
      Emat[9][3]=1;
      
      
      /*Those firt 5 all have the extra 1 in the middle, self position*/
      //1, 0, 1, 0 in the middle also possible
      //1, 1, 1, 0 in the middle possible
      //0, 0, 0, 0 in the middle possible
      //0, 1, 0, 0 2 &3 share, 1 is separate
      //1, 0, 0, 0 1&2 share, 2 not self
      //0, 0, 1, 0 2&3 share, 2 not self
      
      Emat[10][0]=0;
      Emat[10][1]=0;
      Emat[10][2]=0;
      Emat[10][3]=0;
      
      Emat[11][0]=0;
      Emat[11][1]=1;
      Emat[11][2]=0;
      Emat[11][3]=0;
      
      Emat[12][0]=0;
      Emat[12][1]=0;
      Emat[12][2]=1;
      Emat[12][3]=0;
      
      Emat[13][0]=1;
      Emat[13][1]=0;
      Emat[13][2]=0;
      Emat[13][3]=0;
      
      Emat[14][0]=1;
      Emat[14][1]=1;
      Emat[14][2]=1;
      Emat[14][3]=0;
      
      Emat[15][0]=1;
      Emat[15][1]=0;
      Emat[15][2]=1;
      Emat[15][3]=0;
      
      //protein 2
      Emat[16][0]=0;
      Emat[17][0]=1;
      //now have to give each network a number

    
}
