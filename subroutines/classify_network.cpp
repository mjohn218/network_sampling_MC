
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <sys/time.h>
#include "pro_classes.h"
#include "classify_network.h"

using namespace std;

int classify_network(int nwhole, int ninterface, int *Adj, int **Speclist, int *p_home, Protein *wholep, int *numpartners)
{
    //Use a decision tree to determine what IIN state the following PPI network is in:
    //0 2   0   1
    //1 1   0
    //It is possible to add one extra edge, and the maximum number of interfaces is 6.
    //29 possible states. Determine which state the IIN is in via a decision tree.

  cout<<"Inside classify_network" << endl;
    int state=-1;
    
    int i,j;
    int N=ninterface;
    int self_pos=0;
    int transedgealone=0;
    //Count edges
    int nedge=0;
    int nconnect[2]; //0: 0->0 connections, 1: 0->1 connections
    nconnect[0]=0;
    nconnect[1]=0;
    for(i=0;i<N;i++){
        for(j=0;j<numpartners[i];j++){
            if(Speclist[i][j]==i){
                self_pos++;
                nedge+=10;
            }
            else
                nedge+=5;
                
            if(p_home[i]==p_home[Speclist[i][j]] && i!=Speclist[i][j])//0->0 connection
                nconnect[0]++;
	    else if(p_home[i]==p_home[Speclist[i][j]] && i==Speclist[i][j])//0->0 connection
	      nconnect[0]+=2;
            else if(p_home[i]==0)
                nconnect[1]++;
            
            if(p_home[i]==0 && p_home[Speclist[i][j]]==1 && numpartners[i]==1)
                transedgealone++;
        }
    }
    nedge/=10;
    nconnect[0]/=2;
    if(nedge!=(nconnect[0] + nconnect[1])){
        cout<<"Problem with nconnect. Neconnect equals: "<<endl;
	cout<<"O->0: " << nconnect[0]<<endl;
	cout<<"0->1: " << nconnect[1]<<endl;
	cout<<"But nedge=" << nedge<<endl;
        exit(1);
        }
    
    if(nedge==2){
        if(self_pos==1){//Is the self-edge not split?
            if(N==2)
                state=0;
            else if(N==3)
                state=1;
        }
        else{
            if(N==3)
                state=2;
            else if(N==4)
                state=3;
        }
        
    }
    else if(nedge==3){
        if(nconnect[0]==2){
            if(transedgealone==0){
                if(self_pos==0){
                    if(N==5)
                        state=4;
                    else if(N==4){
                        state=5;
                        for(i=0;i<N;i++){
                            if(numpartners[i]==3)
                                state=6;
                        }
                    }
                }
                else if(self_pos==2)
                    state=7;
                else if(self_pos==1){
		  /*bool transhasself=false;
                    for(i=0;i<N;i++){
                        for(j=0;j<numpartners[i];j++){
                            if(p_home[i]==0 && p_home[Speclist[i][j]]==1 && Adj[i*N+i]==1)
                                transhasself=true;
                        }
			}*/
                    if(N==3){
		      state=10;
		      for(i=0;i<N;i++){
			if(numpartners[i]==3)
                            state=8;
		      }
		    }
                    else if(N==4){
		      state=11;
		      for(i=0;i<N;i++){
			if(numpartners[i]==2){
			  for(j=0;j<2;j++){
			    if(Speclist[i][j]==i)
			      state=9;
			  }
			}
		      }
		      
                    }
                }
            }
            else if(transedgealone==1){
                if(self_pos==1){
                    if(N==5)
                        state=12;
                    else if(N==4)
                        state=13;
                }
                else if(self_pos==2){
                    state=14;
                }
                else if(self_pos==0){
                    if(N==5)
                        state=15;
                    else if(N==6)
                        state=16;
                }
            }
        }
        else if(nconnect[1]==2){
            if(wholep[1].ninterface==1){
	      if(self_pos==1){
		if(N==3)
		  state=17;
		else if(N==4)
		  state=29;
	      }
	      else if(self_pos==0){
		if(N==3)
		  state=18;
		else if(N==4)
		  state=19;
		else if(N==5)
		  state=30;
	      }
            }
            else if(wholep[1].ninterface==2){
                if(self_pos==1){
                    if(N==3)
                        state=20;
                    else if(N==5)
                        state=21;
                    else if(N==4){
                        if(transedgealone==0)
                            state=22;
                        else if(transedgealone==1)
                            state=23;
                    }
                }
                else if(self_pos==0){
                    bool sharedpart=false;
                    int w1=wholep[1].valiface[0];
                    int w2=wholep[1].valiface[1];
                    if(Speclist[w1][0]==Speclist[w2][0])
                        sharedpart=true;
                    if(sharedpart==true){
                        if(N==4)
                            state=24;
                        else if(N==5)
                            state=25;
                    }
                    else{
                        if(N==6)
                            state=26;
                        else if(N==4)
                            state=27;
                        else if(N==5)
                            state=28;
                    }
                }
            }
        }
    }
    else{//Error. Nedge should be 2 or 3.
        cout<<"Error. nedge is " << nedge << endl;
        exit(1);
    }
        
        if(state==-1){
            cout<< "Error. Could not find state. Exiting" << endl;
            exit(1);
        }

	        
        return state;
        
        
        
}
