#include "pro_classes.h"
#include "utility_calls.h"
#include "network_metrics2.h"
#include "matmultiply.h"

int find_islands(int N, int *Adj0, int **modlist, int *modsize)
{
    //Find islands in a network, any set of nodes where it is possible to reach each node from any other node.
    
    int i,j;
    double *Adj = new double[N*N];
    for(i=0;i<N*N;i++)
        Adj[i]=double(Adj0[i]);
    
    //Make sure all diagonal elements are set to 1
    for(i=0;i<N;i++)
      Adj[i*N+i]=1.0;

    //Print out Adj
    /*cout << "Adjacency Matrix:" << endl;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            cout<<Adj[i*N+j]<<'\t';
        }
        cout<<endl;
    }*/
    
    //Now perform the algorithm. Until Adj_i+1 = Adj_i, Adj_i+1 = Adj_0 * Adj_i, set all values >1 to 1.
    double *Adjold = new double[N*N];
    double *Adjnew = new double[N*N];
    for(i=0;i<N*N;i++){
        Adjnew[i]=0;
        Adjold[i]=Adj[i]; //Remember when I do this for real to make sure all the diagonal elements equal 1!
    }
    
    bool check=false;
    bool tmpcheck=false;
    
    char transA='t';
    char transB='n';
    double alph=1;
    double beta=0;
    int lda = N;
    int ldb = N;
    int ldc = N;
    
    int stp = 1;
    while(check==false){
      //cout << "Step " << stp << endl;
        dgemm(&transA,&transB,&N,&N,&N,&alph,Adj,&lda,Adjold,&ldb,&beta,Adjnew,&ldc);
        for(i=0;i<N*N;i++){
            if(Adjnew[i]>1)
                Adjnew[i]=1;
        }
        
        
        //Print out Adjnew
        
        /*for(i=0;i<N;i++){
         for(j=0;j<N;j++){
         cout<<Adjnew[i*N+j]<<'\t';
         }
         cout<<endl;
         }*/
        
        
        tmpcheck=true;
        for(i=0;i<N*N;i++){
            if(Adjold[i]!=Adjnew[i])
                tmpcheck=false;
        }
        
        if(tmpcheck==false){
            for(i=0;i<N*N;i++)
                Adjold[i]=Adjnew[i];
        }
        stp++;
        check=tmpcheck;//will break the loop if tmpcheck is true, ie Adjold and Adjnew were equal
    }
    
    //Now print out final matrix
    /*    cout << "Final matrix:" << endl;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            cout<<Adjnew[i*N+j]<<' ';
        }
        cout<<endl;
	}*/
    
    //Record islands
    int modloc[N];
    for(i=0;i<N;i++)
        modloc[i]=-1;
    
    //modlist input must be NxN, modsize should be N
    /*int **modlist = new int*[N];
    for(i=0;i<N;i++)
        modlist[i]=new int[N];
    
    int *modsize = new int[N];*/
    for(i=0;i<N;i++)
        modsize[i]=0;
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
	modlist[i][j]=0;
      }
    }

    bool newmod=false;
    bool truenew=false;
    int modcurr=0;
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            if(modloc[j]<0 && Adjnew[i*N+j]==1){
                truenew=true;
                modloc[j]=modcurr;
                modlist[modcurr][j]=1;
            }
            if(i!=(N-1) && Adjnew[i*N+j]==1){
                if(Adjnew[(i+1)*N+j]!=Adjnew[i*N+j]){
                    newmod=true;
                    
                }
            }
            else if(modloc[j]<0 && Adjnew[i*N+j]==0)
                modlist[modcurr][j]=0;
        }
        
        if(!truenew){
            modcurr--;
            truenew=true;
        }
        
        if(newmod){
            modcurr++;
            newmod=false;
            truenew=false;
        }
    }
    
    //Create and print out modlist
    modcurr++;

    for(i=0;i<modcurr+1;i++){
        //cout<<"In island " << i << ":\t";
        for(j=0;j<N;j++){
            if(modlist[i][j]==1){
                //cout<< j << "\t";
                modsize[i]++;
            }
        }
        //cout<< "Size: " << modsize[i] << endl;
    }

    delete [] Adj;
    delete [] Adjold;
    delete [] Adjnew;
    return modcurr;
}
