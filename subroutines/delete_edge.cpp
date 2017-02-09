#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int delete_edge(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, int maxni, int *selfflag, int PPIedge)
{
    
    int chg=0;
    int N=ninterfaces;
    
    int i,j;
    double rnum1,rnum2;
    int w1, w2;
    int ind, tmp;
    int i1,i2,p1,p2;
    
    //Create Adj
    
    for(i=0;i<N*N;i++)
        Adj[i]=0;
    for(i=0;i<N;i++){
      tmp=numpartners[i];
      ////cout << i << "\t" << tmp;
        for(j=0;j<tmp;j++){
            ind=Speclist[i][j];
	    ////cout << "\t" << ind;
            Adj[i*N+ind]=1;//i is the column, j the row, Symmetric matrix
        }
	////cout << endl;
    }
    

    //First, select a protein
    //rnum1=rand()*1.0*nwhole/((double)RAND_MAX+1);
    //p1=int(rnum1);//index of first protein
    //Find a PPI partner. Self-interactions are allowed
    
    //rnum2=rand()*1.0*ppi[p1].nppartner/((double)RAND_MAX+1);
    //tmp=int(rnum2);
    //p2=ppi[p1].pplist[tmp];

    //    //cout<<"p1,p2 are " << p1 << ", "<< p2 << endl;

    //Find a PPI connection
    int NPedge=0;
    for(i=0;i<nwhole;i++)
        NPedge+=ppi[i].nppartner;
     
    //rnum1=trand()*1.0*NPedge/RAND_MAX;
    rnum1=NPedge*1.0*rand_gsl();
    tmp=int(rnum1);
    p1=-1;
    for(i=0;i<nwhole;i++){
      for(j=0;j<ppi[i].nppartner;j++){
	  if(j==tmp){
	      p1=i;
	      p2=ppi[p1].pplist[j];
            }
      }
      if(p1<0)
	  tmp-=ppi[i].nppartner;
        else
	  break;
    }

    
    
    //For testing purposes:
    // p1=0;
    //rnum1=trand()*2.0/RAND_MAX;
    ////cout<"rnum is " << int(rnum1) << endl;
    //if(int(rnum1)==0)
    //p2=1;
    //else 
    //p2=0;
         
    //cout<<"p1 is " << p1 << ", p2 is " << p2 << endl;

    //Now we need to see if nconnect>1. Also, the interfaces for each edge need to have >1 partner, or else they will become orphan edges
    int nconnect = 0;
    int ndel=0;
    if(p1!=p2){
    for(i=0;i<wholep[p1].ninterface;i++){
        for(j=0;j<wholep[p2].ninterface;j++){
            w1=wholep[p1].valiface[i];
            w2=wholep[p2].valiface[j];
            if(Adj[w1*N+w2]==1){// && numpartners[w1]>1 && numpartners[w2]>1){
                //cout << "Connection found b/t interfaces " << w1 << " and " << w2 << endl;
                nconnect++;
                if(numpartners[w1]>1 && numpartners[w2]>1)
                    ndel++;
	    }
        }
    }
    }
    else if(p1==p2){
      for(i=0;i<wholep[p1].ninterface;i++){
	for(j=i;j<wholep[p1].ninterface;j++){
	  w1=wholep[p1].valiface[i];
	  w2=wholep[p1].valiface[j];
        if(Adj[w1*N+w2]==1){// && numpartners[w1]>1 && numpartners[w2]>1)
            nconnect++;
	    //cout<<"Connection found b/t interfaces " << w1 << " and " << w2 << endl;
            if(numpartners[w1]>1 && numpartners[w2]>1)
                ndel++;
        }
      }
      }
    }

    if(nconnect>1 && ndel>0)
        chg=1;
    
    //Want to avoid accessing the same move in multiple ways to balance the probabilities
    if(p2<p1)
      chg=0;
    
    if(chg!=0){//Perform the move.
        //Now that we know there is more than one edge between these proteins that may be deleted, pick one at random
        //cout<<"There are " <<nconnect<< " connections between the proteins, of which " << ndel << " may be deleted."<<endl;
        int edgepairs[ndel][2];
        nconnect=0;
	if(p1!=p2){
        for(i=0;i<wholep[p1].ninterface;i++){
            for(j=0;j<wholep[p2].ninterface;j++){
                w1=wholep[p1].valiface[i];
                w2=wholep[p2].valiface[j];
                if(Adj[w1*N+w2]==1 && numpartners[w1]>1 && numpartners[w2]>1){
                    edgepairs[nconnect][0]=w1;
                    edgepairs[nconnect][1]=w2;
                    nconnect++;
                }
                
            }
        }
	}
	else if(p1==p2){
	  for(i=0;i<wholep[p1].ninterface;i++){
	    for(j=i;j<wholep[p1].ninterface;j++){
	      w1=wholep[p1].valiface[i];
	      w2=wholep[p1].valiface[j];
	      if(Adj[w1*N+w2]==1 && numpartners[w1]>1 && numpartners[w2]>1){
		edgepairs[nconnect][0]=w1;
		edgepairs[nconnect][1]=w2;
		nconnect++;
	      }
	    }
	  }
	}

        rnum1=rand_gsl()*1.0*nconnect;
        ind=int(rnum1);
        i1=edgepairs[ind][0];
        i2=edgepairs[ind][1];
        
        //Now find pgen_f
        int fact1,fact2,fact;
        int factf=1;
        int copy=0;
        int t1,t2,t3,t4;
        
        double pchooseedge_f = 1.0/(nconnect);
        //Need to see if edges are perfect "copies"
        for(i=0;i<nconnect;i++){
            if(i!=ind){
                t1=edgepairs[i][0];
                t2=edgepairs[ind][0];
                //cout<<"t1, t2 = " << t1 << ", " << t2 << endl;
                if(t1==t2)
                    copy=1;
                else
                    fact1=find_copies(nwhole,N,numpartners, wholep,Speclist,Adj, t1,p_home,t2,copy);
                if(copy==1){
                    copy=0;
                    t3=edgepairs[i][1];
                    t4=edgepairs[ind][1];
                    //cout << "t3, t4 = " << t3 << ", " << t4 << endl;
                    fact2=find_copies(nwhole,N,numpartners, wholep,Speclist,Adj, t3,p_home,t4,copy);
                    if(copy==1)
                        factf++;
                }
            }
        }
        double pgen_f=double(factf)*pchooseedge_f;
        //cout<< "In move 6: pchooseedge_f is " << pchooseedge_f << " and fact_f is " << factf << endl;
        
        //Now delete the edge. Need to be careful in case it is a self-edge
        //cout << "In move 6: Deleting edge connecting iface " << i1 << " and iface " << i2 << endl;

	if(i1==i2){
            numpartners[i1]--;
            Adj[i1*N+i1]=0;
            
        }
        else{
            numpartners[i1]--;
            numpartners[i2]--;
            Adj[i1*N+i2]=0;
            Adj[i2*N+i1]=0; //The Speclist will be updated below
            
        }
        
        /*update the Speclist */
	int t;
        for(i=0;i<N;i++){
            t=0;
            for(j=0;j<N;j++){
                if(Adj[i*N+j]==1){
                    Speclist[i][t]=j;
                    t++;
                }
            }
        }

	        //Now find pgen_b
        double pchoosepair_b = 1.0/(wholep[p1].ninterface * wholep[p2].ninterface);
        
       
        //pgen will be affected if either i1 or i2 has a "copy"
        
	  fact1=find_copies(nwhole,N, numpartners, wholep, Speclist, Adj, i1, p_home, i2, copy);
	  fact2=find_copies(nwhole,N, numpartners, wholep, Speclist, Adj, i2, p_home, i1, copy);
            if(copy==0)
                fact=fact1*fact2;
            else if(copy==1 && i1!=i2)//i1 and i2 are copies of each other. But see if they are identical.
                fact=fact1*(fact1-1)/2;
            else if(i1==i2)
                fact=fact1;//if x copies on protein p1, then x ways to make a new self-edge

	    if(p_home[i1]==p_home[i2] && i1!=i2)
	    fact*=2;
        
        double pgen_b=pchoosepair_b*1.0*fact;
        
	//cout<< "pchoosepair_b="<<pchoosepair_b<< " and factb=" << fact << ". pgen_b is " << pgen_b<<endl;;

	pgen_ratio=pgen_b/pgen_f;
        
            
    }
    return chg;
    
}
