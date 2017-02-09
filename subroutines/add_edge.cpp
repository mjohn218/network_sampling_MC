#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"


int add_edge(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, int maxni, int maxne, int *selfflag, double &pf, double &pb, int PPIedge)
{
    //Move created by DH
    //The purpose of this move is to add an edge between a protein and one of its PPI partners, since when you look at an IIN, there may be several edges between proteins
    //The reverse move will be delete_edge. A combine_interface move that deletes an edge will not be allowed since it is irreversible.

    int ind;
    int i, j;
    int w1,w2;
    int chg=1;
    double rnum1,rnum2;
    int sumf=0;
    int p1, p2, i1, i2;
    int tmp;
    int copy=0;
    
    double nedge=0;
    for(i=0;i<ninterfaces;i++){
      for(j=0;j<numpartners[i];j++){
	if(Speclist[i][j]!=i)
	  nedge+=0.5;
	else
	  nedge+=1.0;
      }
    }

    if(int(nedge)==maxne)
      chg=0;

    //Create Adj
    int N=ninterfaces;
    for(i=0;i<N*N;i++)
        Adj[i]=0;
    for(i=0;i<N;i++){
        for(j=0;j<numpartners[i];j++){
            ind=Speclist[i][j];
            Adj[i*N+ind]=1;//i is the column, j the row, Symmetric matrix
        }
    }
    
    //First, select a protein
    //rnum1=rand()*1.0*nwhole/((double)RAND_MAX+1);
    //p1=int(rnum1);//index of first protein
    //Find a PPI partner. Self-interactions are allowed
    ////cout<< "rnum1 was " << rnum1 << " and p1 is " << p1 << endl;
    //rnum2=rand()*1.0*ppi[p1].nppartner/((double)RAND_MAX+1);
    //tmp=int(rnum2);
    //p2=ppi[p1].pplist[tmp];
    
    //Select a PPI connection
    ////cout<<"Number of PPI edges is " << PPIedge << endl;
    
    int NPedge=0;
    for(i=0;i<nwhole;i++)
        NPedge+=ppi[i].nppartner;
    
    
    rnum1=1.0*NPedge*rand_gsl();
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

    //Want to avoid having multiple ways to access the same move to balance the probabilities
    if(p2<p1)
      chg=0;
    //For testing purposes:
    //p1=0;
    //rnum1=trand()*2.0/RAND_MAX;
    ////cout<<"rnum is " << int(rnum1) << endl;
    //if(int(rnum1)==0)
    //p2=1;
    //else
    //p2=0;
    //p2=0;


    //cout<<"p1 is " << p1 << ", its selfflag is " << selfflag[p1] << ", and p2 is " << p2 << endl;
    
    //Find nconnections between p1 and p2.
    //cout << "Looking for connections between proteins " << p1 << " and " << p2 << endl;
    int nconnect = 0;
    if(p1!=p2){
    for(i=0;i<wholep[p1].ninterface;i++){
        for(j=0;j<wholep[p2].ninterface;j++){
            w1=wholep[p1].valiface[i];
            w2=wholep[p2].valiface[j];
	    if(Adj[w1*N+w2]==1)
                nconnect++;
        }
    }
    }
    else if(p1==p2){
      for(i=0;i<wholep[p1].ninterface;i++){
	for(j=i;j<wholep[p1].ninterface;j++){
	  w1=wholep[p1].valiface[i];
	  w2=wholep[p1].valiface[j];
	if(Adj[w1*N+w2]==1)
	  nconnect++;
	}
      }
    }

    //cout << "There are " << nconnect << " connections.\n";
    if(nconnect==0){
      cout << "No connections between PPI partners\n";
      cout << "P1: " << p1 << " P2: " << p2 << endl;
      exit(1);
    }
    
    //Now we must randomly select an interface on p1 and on p2. Again, self-interactions are allowed, so if p1=p2, then it's possible to add a self-interaction, or just an interaction between two interfaces on the same protein.
    rnum1=rand_gsl()*1.0*wholep[p1].ninterface;
    w1=int(rnum1);
    rnum2=rand_gsl()*1.0*wholep[p2].ninterface;
    w2=int(rnum2);
    i1=wholep[p1].valiface[w1];
    i2=wholep[p2].valiface[w2];
    
    //If there is already a connection, don't do the move
    if(Adj[i1*N+i2]==1){
      //cout<<"Already an edge between " << i1 << " and " << i2 << endl;
      chg=0;
    }
    
    if(chg!=0){
      int fact1, fact2, fact;
      double p_selectifaces=1.0/(wholep[p1].ninterface * wholep[p2].ninterface);

    
    //pgen will be affected if either i1 or i2 has a "copy"
	
      fact1=find_copies(nwhole,N, numpartners, wholep, Speclist, Adj, i1, p_home, i2, copy);
      fact2=find_copies(nwhole,N, numpartners, wholep, Speclist, Adj, i2, p_home, i1, copy);
      if(copy==0)
          fact=fact1*fact2;
      else if(copy==1 && i1!=i2)//i1 and i2 are copies of each other. But see if they are identical.
          fact=fact1*(fact1-1)/2;
      else if(i1==i2)
          fact=fact1;//if x copies on protein p1, then x ways to make a new self-edge
    
      //if i1 and i2 are on the same protein but are different, two ways to get same result
      if(p1==p2 && i1!=i2)
      fact*=2;

      //cout << "In Move 5: Adding edge connecting " << i1 << " and " << i2 << endl;
      //cout << "p_selectifaces = " << p_selectifaces << " and fact_f = " << fact << endl;
      double pgen_f=fact*1.0*p_selectifaces; //p_selectp2 will cancel out.
    
    //Now that i1 and i2 are selected, we can add an edge
      if(i1==i2){
        int n1=numpartners[i1];
        numpartners[i1]++;
        Speclist[i1][n1]=i1;
        Adj[i1*N+i1]=1;
        nconnect++;
      }
      else{
        int n1=numpartners[i1];
        int n2=numpartners[i2];
        numpartners[i1]++;
        numpartners[i2]++;
        Speclist[i1][n1]=i2;
        Speclist[i2][n2]=i1;
        Adj[i1*N+i2]=1;
        Adj[i2*N+i1]=1;
        nconnect++;
    }
    
    //Viola! There is now an extra edge in the system. Now to calculate pgen_ratio.
    
    int edgepairs[nconnect][2];
    nconnect = 0;//Calculate again because we want to skip interfaces with one edge
    //ind=0;
    if(p1!=p2){
    for(i=0;i<wholep[p1].ninterface;i++){
        for(j=0;j<wholep[p2].ninterface;j++){
            w1=wholep[p1].valiface[i];
            w2=wholep[p2].valiface[j];
	    if(w1!=w2 && Adj[w1*N+w2]==1 && numpartners[w1]>1 && numpartners[w2]>1){
	      edgepairs[nconnect][0]=w1;
	      edgepairs[nconnect][1]=w2;
	      nconnect++;
	    }
	    if(w1==i1 && w2==i2){
	      ind=nconnect-1;
////cout<<"ind is now " << ind << endl;
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
	if(w1==i1 && w2==i2){
	  ind=nconnect-1;
////cout<<"ind is now " << ind << endl;
	}
	else if(w1==i2 && w2==i1){
	  ind=nconnect-1;
	//cout<<"ind is now " << ind << endl;
	}	
	}
      }
    }

    //cout<< "In move 5: ind is " << ind << " and nconnect is " << nconnect << endl;
    
    double pchooseedge_b = 1.0/nconnect; //p_selectp2 will cancel out. Need to select an edge between p1 and p2 to delete
        int factb=1;
	int t1,t2,t3,t4;
        //Need to see if edges are perfect "copies"
        for(i=0;i<nconnect;i++){
            if(i!=ind){
                t1=edgepairs[i][0];
                t2=edgepairs[ind][0];
		//cout<<"t1, t2 = " << t1 << ", " << t2 << endl;
                fact1=find_copies(nwhole,N,numpartners,wholep,Speclist,Adj,t1,p_home,t2,copy);
                if(copy==1){
                    copy=0;
                    t3=edgepairs[i][1];
                    t4=edgepairs[ind][1];
		    //cout << "t3, t4 = " << t3 << ", " << t4 << endl;
                    fact2=find_copies(nwhole,N,numpartners,wholep,Speclist,Adj,t3,p_home,t4,copy);
                    if(copy==1)
                        factb++;
                }
            }
        }
        
	//cout << "pchooseedge_b is " << pchooseedge_b << " and factb is " << factb << endl;
    double pgen_b=factb*1.0*pchooseedge_b;

    pgen_ratio=pgen_b/pgen_f;
    
    }
    
    return chg;
    
}
