#include "net_props.h"
#include <cstdlib>

double cluster_cof_mod(int i1, int N, int *numpartners, int **Speclist)
{
    
    //This code find the local clustering coefficient for i1, but excludes triangles where two or three of the nodes have a self-loop.
    
    int i,j;
    int self=0;
    double clocal;
    //Find kfirst
    int kfirst;
    
    for(i=0;i<numpartners[i1];i++){
        if(Speclist[i1][i]==i1)
            self=1;
    }
    
    kfirst=numpartners[i1];

    //if(self==1){
    //clocal=0;//Ignore triangles where there is a self-loop, since this does not represent a specificity penalty
      //kfirst--;
    //}
    if(kfirst==1){
        clocal=0;
    }
    else{
        
        //Make Adj from Speclist
        int Adj[N*N];
        //Speclist is the most up-to-date here, so remake Adj from Speclist
        for(i=0;i<N*N;i++)
            Adj[i]=0;
        for(i=0;i<N;i++){
            for(j=0;j<numpartners[i];j++){
                Adj[i*N+Speclist[i][j]]=1;
            }
        }
        
        //Must find the number of triangles
        
        int w1,w2;
        int ntri=0;
        for(i=0;i<kfirst;i++){
            for(j=i+1;j<kfirst;j++){ //Only check each pair once
                w1=Speclist[i1][i];
                w2=Speclist[i1][j];
                if(Adj[w1*N+w2]==1){
		  if(Adj[w1*N+w1]+Adj[w2*N+w2]+self<2){//Ignore if at least two interfaces have a self-loop
                    ntri++;
		  }
                }
            }
        }
        
        clocal=ntri*2.0 / (kfirst * (kfirst - 1));
        
        if(clocal>1){
            cout << "Error: cluster_cof_mod. Exiting" << endl;
            exit(1);
        }
    
    }

    return clocal;
    
}
