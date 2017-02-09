#include "net_props.h"
#include <cstdlib>

void grid_cof(double *gcs, int N, int *numpartners, int **Speclist)
{

  //July 21st, 2015 by David Holland
    //Code to find local (modified) grid coefficient for all interfaces

    //The equation used here is (Q + 1/(k+1)) / (Z + 1)
    //Where Q is the number of closed squares (4 nodes, 4 edges), k is the degree minus self-interactions (ie the number of interfaces one step away from i1), and Z is the number of "possible" closed squares, equal to k2nd * k * (k-1) / 2, where k2nd is the number of interfaces two steps away from i1.
  
  int i,j,k,l;
    int tmp,tmp2;
    
  int Adj[N*N];
  //Speclist is the most up-to-date here, so remake Adj from Speclist
  for(i=0;i<N*N;i++)
    Adj[i]=0;
  for(i=0;i<N;i++){
    for(j=0;j<numpartners[i];j++){
      Adj[i*N+Speclist[i][j]]=1;
    }
  }

  //Make sure that Adj is symmetric
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(Adj[i*N+j]!=Adj[j*N+i]){
          cout << "Adj is not symmetric (grid cof). Exiting.";
          exit(1);
      }
    }
  }

    
  
  
//Code for finding squares
    int nsquares[N];
    for(i=0;i<N;i++){
      nsquares[i]=0;
    }
    
  for(i=0;i<N;i++){
    for(j=i+1;j<N;j++){
      if(Adj[i*N+j]==1){
          for(k=i+1;k<N;k++){
              if(Adj[j*N+k]==1 && j!=k){
                  for(l=j+1;l<N;l++){
                      if(Adj[k*N+l]==1 && Adj[i*N+l]==1 && k!=l){
                          if(Adj[i*N+k]==0 && Adj[j*N+l]==0){
			    nsquares[i]++;
			    nsquares[j]++;
			    nsquares[k]++;
			    nsquares[l]++;
			    //cout<<"Square found: "<< i << ", " << j << ", " << k << ", " << l << endl;
                          }
                      }
                  }
              }
          }
      }
    }
  }

//Find interfaces one step away from i1

  int kfirst,ksec;
  int conn1[N];
  int conn2[N];
  int Q,Z;  
  for(i=0;i<N;i++){  
  
    kfirst=numpartners[i];
    for(j=0;j<N;j++){
        conn1[j]=0;
    }
    for(j=0;j<numpartners[i];j++){
        conn1[Speclist[i][j]]=1;
        if(i==Speclist[i][j]){ //There exists a self-interaction
            kfirst-=1;
            conn1[i]=0;
        }
    }
    
    
  //Next, need to find number of unique second partners for each node
    ksec=0;
    for(j=0;j<N;j++)
        conn2[j]=0;
    for(j=0;j<numpartners[i];j++){
        if(Speclist[i][j]!=i){ //Ignore self-edges
            tmp=Speclist[i][j];
            for(k=0;k<numpartners[tmp];k++){
                tmp2=Speclist[tmp][k];
                if(tmp2!=tmp) //ie don't count tmp's self interactions
                    conn2[tmp2]=1;
            }
        }
    }
    conn2[i]=0;//don't count interactions back with i
    for(j=0;j<N;j++){ //Add all unique second partners
        if(conn1[j]==0 && conn2[j]>0){ //Ignore first partners
          ksec+=1;
        }
    }
    //Now calculate the grid coefficient
    
    Z = ksec * kfirst * (kfirst-1) / 2;
    Q= nsquares[i];
    double dummy=1.0;
    gcs[i] = (Q*1.0 + dummy) / (Z*1.0 + dummy);
  }
    
}
