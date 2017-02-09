#include "net_props.h"
#include <cstdlib>

void find_type(int i, int j, int k, int l, int*Adj, int N, double **hist, int tak)
{
    
    int a,b,c,d,e,f;
    
    //tak++;
    
    a=Adj[i*N+j];
    b=Adj[i*N+k];
    c=Adj[i*N+l];
    d=Adj[j*N+k];
    e=Adj[j*N+l];
    f=Adj[k*N+l];
    
    ////cout << "The number of edges is: " << a+b+c+d+e+f << endl;
    
    //Now check to see if there are three edges
    if(a+b+c+d+e+f==3) {
        //See if it is a star
        if((a+b+c==3)||(a+d+e==3)||(b+d+f==3)||(c+e+f==3)){
            hist[3][0]+=1;
                    }
        else{ //it must be an open square
            hist[3][1]+=1;
        }
    }
    //Now see if there are four edges
    else if(a+b+c+d+e+f==4){
        //See if it is a square
        if((a+b+c==2)&&(a+d+e==2)&&(b+d+f==2)&&(c+e+f==2)){
            hist[4][1]+=1;
            ////cout << "A square has been found. The four nodes are: " << i << " " << j << " " << k << " " << l << endl;
        }
        //Otherwise it must be a star
        else{
            hist[4][0]+=1;
            
        }
    }
    else if(a+b+c+d+e+f==5){
        hist[5][0]+=1;
        

    }
    else if(a+b+c+d+e+f==6){
        hist[6][0]+=1;
    }

    
}
