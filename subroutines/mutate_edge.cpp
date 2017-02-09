#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int mutate_edge(int nwhole, int *numpartners, int **Speclist, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *selfflag, int *p_home)
{
  // //cout << "Inside mutate_edge" << endl;
  /*In this case first choose a protein that has more than one interface*/
  /*ADDED SELF*/
  pgen_ratio=1;
  double rnum1=rand_gsl()*1.0*nwhole;
  int p1=int(rnum1);//index of first protein
    
  int chg=1;
  int i, t, j;
  int loc;
  int i1, ind, ip1;
  double pf, pb;
  int nedge=0;
    int copy=0;
    
    //Make Adj matrix
    int N=0;
    for(i=0;i<nwhole;i++)
        N+=wholep[i].ninterface;
    
    int * Adj2 = new int[N*N];
    for(i=0;i<N*N;i++)
        Adj2[i]=0;
    for(i=0;i<N;i++){
        for(j=0;j<numpartners[i];j++){
            ind=Speclist[i][j];
            Adj2[i*N+ind]=1;//i is the column, j the row, Symmetric matrix
        }
    }
    
    
  ////cout <<"in mutate edge, ninterfaces : "<<wholep[p1].ninterface<<endl;
  chg=0;
  /*if(wholep[p1].ninterface<ppi[p1].nppartner){ //meaning: # interfaces < # of protein partners, so at least one interface has >1 partners (DH)
    chg=1; //Allow an edge movement (DH)
    
  /*If they are equal, you can still move an edge if there is a self interaction on the protein*/
  /*if(wholep[p1].ninterface==ppi[p1].nppartner &&selfflag[p1]==1){
    
    chg=1;//binds to self, so its possible to have an extra interface
  }*/
  
    //New criteria by DH
    if(wholep[p1].ninterface<ppi[p1].nppartner)
        chg=1;
    else{
    for(i=0;i<wholep[p1].ninterface;i++){
        ind=wholep[p1].valiface[i];
        if(numpartners[ind]>1)
            chg=1;
        for(j=0;j<numpartners[ind];j++){
            if(Speclist[ind][j]==ind)
                chg=1;
        }
    }
    }
    

  if(wholep[p1].ninterface==1)
    chg=0; //No edge movement allowed.


  int oself=0;
  int w, spf;
  double denom;
  int tmp;
  if(chg!=0){
    /*Pick one edge off this protein*/
      /*Maggie's code. Obsolete since there may be multiple edges connecting proteins
    if(selfflag[p1]==1){
      //find out if it is to itself, or split self
      
      for(i=0;i<wholep[p1].ninterface;i++){
	ind=wholep[p1].valiface[i]; 
          for(t=0;t<numpartners[ind];t++){
              if(Speclist[ind][t]==ind)
                  oself=1;
          } 
      }
      if(oself==0){
	
          oself=-1;//the self interaction is split
          nedge=ppi[p1].nppartner+1;
      }
      else{ //self interaction to the same interface
          nedge=ppi[p1].nppartner;
      }
      denom=1.0*nedge;
    }
    else{
      //no self
      nedge=ppi[p1].nppartner;
      denom=1.0*nedge;
    }*/
      
      //New code by DH. Count nedge by looking at each interface.
      nedge=0;
      for(i=0;i<wholep[p1].ninterface;i++)
          nedge+=numpartners[wholep[p1].valiface[i]];
    
      //cout<<"nedge is " << nedge << endl;
    
    rnum1=rand_gsl()*1.0*nedge;
    w=int(rnum1);
    i=0;
    loc=w;
    
    i1=-1;
    
    while(i<wholep[p1].ninterface){
      ind=wholep[p1].valiface[i];
    
      if(loc<numpartners[ind]){
          ip1=Speclist[ind][loc];
          i1=ind;
	
          i=wholep[p1].ninterface;//break out
      }
      else{
          loc-=numpartners[ind];
	
	  i++;
      }
    }
    
    
    /*Now pick another interface on this protein to move the edge to*/
    double rnum2=rand_gsl()*1.0*wholep[p1].ninterface;
    
    int i2=wholep[p1].valiface[int(rnum2)];
      
    while(i1==i2){
      rnum2=rand_gsl()*1.0*wholep[p1].ninterface;
      i2=wholep[p1].valiface[int(rnum2)];
    }
    //cout <<"in mutate edge: p1 "<<p1 <<" ninterface "<<wholep[p1].ninterface<<" rnum2: "<<rnum2<<" i2: "<<i2<<" i1: "<<i1<<" ip1: "<<ip1<<endl;
    //Check for shared partners, in which case stop the move (added by DH)
    for(i=0;i<numpartners[i2];i++){
      if(Speclist[i2][i]==ip1)
	chg=0;
    }
      

  
    if(chg!=0){
        
        //Want to see if any edges on other interfaces are "copies"
        //To be a copy, both interfaces the edges are connected to must be copies.
        int edgepairs[nedge][2];
        
        int w2=0;
        for(i=0;i<wholep[p1].ninterface;i++){
            tmp=wholep[p1].valiface[i];
            for(j=0;j<numpartners[tmp];j++){
                edgepairs[w2+j][0]=tmp;
                edgepairs[w2+j][1]=Speclist[tmp][j];
                if(tmp==i1 && Speclist[tmp][j]==ip1)
                    w=w2+j;
            }
            w2+=numpartners[tmp];
        }
            
        int fact1,fact2,factf;
	int t1,t2,t3,t4;
        factf=1;//include itself
        for(i=0;i<nedge;i++){
            if(i!=w){
                t1=edgepairs[i][0];
                t2=edgepairs[w][0];
                if(t1==t2)
                    copy=1;
                else
                fact1=find_copies(nwhole, N,numpartners, wholep,Speclist,Adj2, t1,p_home,t2,copy);
                if(copy==1){
                    copy=0;
                    t3=edgepairs[i][1];
                    t4=edgepairs[w][1];
                    if(t3==t4)
                        copy=1;
                    else
                    fact2=find_copies(nwhole, N,numpartners, wholep,Speclist,Adj2, t3,p_home,t4,copy);
                    if(copy==1){
                        factf++;
                        //cout<<"Edgecopy found. Edge one connects: " << t2 << " to " << t4 << ". Edge two connects " << t1 << " to " << t3 << endl;
                    }
                    
                }
            }
        }

    pb=1.0;
    pf=1.0;
    pgen_ratio=pb/pf;
    
        //cout<<"factf is now " << factf << endl;
        
        //Find copies of i2 for pgen
        int ndups,ndups2;
        if(i2==ip1)//Un-splitting a self-edge
            ndups=1;
        else{
            ndups = find_copies(nwhole, N, numpartners, wholep, Speclist, Adj2, i2, p_home, i1, copy);
            if(copy==1)
                ndups--;
            ndups2 = find_copies(nwhole, N, numpartners, wholep, Speclist, Adj2, i2, p_home, ip1, copy);
            if(copy==1 && i1!=ip1)
                ndups--;
        }
        
        
        factf*=ndups;
        //cout<<"factf is "<<factf <<". ndups was " << ndups << ", taking into account that copy="<<copy<<endl;
        
    if(ip1==i1){
      /*Self binding splitting to another edge 
	Move i1's partner ip1 to i2, keep the rest of his partners*/
    
      //pf=1.0/(1.0*ppi[p1].nppartner);
      //pb=1.0/(1.0*ppi[p1].nppartner+1.0);
        pf=1.0/(1.0*nedge);
        pb=1.0/(1.0*nedge+1.0);
      pgen_ratio=pb/pf;
      t=0;
      /*Now swap these*/
      for(i=0;i<numpartners[i1];i++){
          if(Speclist[i1][i]==ip1){
              Speclist[i1][i]=i2;
              //t++;
          }
      }
      //i1 and ip1 are the same, so updating one also updates the other
      //i1 still has a partner! it is i2!

      Speclist[i2][numpartners[i2]]=ip1;
      numpartners[i2]++;//added one
    }else if(numpartners[i1]==1){
      chg=0;
      
    }else{
      /*Move i1's partner ip1 to i2, keep the rest of his partners*/
      
      /*First check if this edge is a split self edge*/
      spf=0;
      if(p_home[ip1]==p1){
          pf=1.0/(1.0*nedge);
          pb=1.0/(1.0*nedge);
          spf=1;
      }
      t=0;
      /*Now swap these*/
      for(i=0;i<numpartners[i1];i++){
          if(Speclist[i1][i]!=ip1){
              Speclist[i1][t]=Speclist[i1][i];
              t++;
          }
      }
      numpartners[i1]--;//got rid of one partner
      if(ip1==i2){
          nedge--;
          pb=1.0/(1.0*nedge);
          //forming a self interaction from a within protein self interaction
          for(i=0;i<numpartners[i2];i++){
              if(Speclist[i2][i]==i1){
                  Speclist[i2][i]=i2;
              }
          }
	
      }else{
	/*Check about a stand alone split edge, results in 2 ways to get to same state
	  if you are moving to one or creating one, extra factor of 2*/
          /*if(spf==0 && oself==-1){
              tmp=Speclist[i2][0];
              if(numpartners[i2]==1 &&p_home[tmp]==p1){
                  if(numpartners[tmp]==1){
                      pf*=2.0;
                  }
              }
              tmp=Speclist[i1][0];
              if(numpartners[i1]==1 &&p_home[tmp]==p1){
                  if(numpartners[tmp]==1){
                      pb*=2.0;
                  }
              }
          }*/
          //Above code obsolete (DH). Will find copies of i2 above instead.
          for(i=0;i<numpartners[ip1];i++){
              if(Speclist[ip1][i]==i1)
                  Speclist[ip1][i]=i2;
          }
          Speclist[i2][numpartners[i2]]=ip1;
          numpartners[i2]++;//added one
      }
        pgen_ratio=pb/pf;
    }
   
        //Now need to get factb for pgen
        
        nedge=0;
        for(i=0;i<wholep[p1].ninterface;i++)
            nedge+=numpartners[wholep[p1].valiface[i]];
        
        //Find location of moved edge
        w2=0;
        t=0;
        for(i=0;i<wholep[p1].ninterface;i++){
            tmp=wholep[p1].valiface[i];
            for(j=0;j<numpartners[tmp];j++){
                edgepairs[w2+j][0]=tmp;
                edgepairs[w2+j][1]=Speclist[tmp][j];
                if(tmp==i2 && j==numpartners[tmp]-1)
                    w=t;
                t++;
            }
            w2+=numpartners[tmp];
        }
        
        //update Adj2
        
        for(i=0;i<N*N;i++)
            Adj2[i]=0;
        for(i=0;i<N;i++){
            for(j=0;j<numpartners[i];j++){
                ind=Speclist[i][j];
                Adj2[i*N+ind]=1;//i is the column, j the row, Symmetric matrix
            }
        }
        //Find edge copies
        copy=0;
        int factb=1;//include itself
        for(i=0;i<nedge;i++){
            if(i!=w){
                t1=edgepairs[i][0];
                t2=edgepairs[w][0];
                if(t1==t2)
                    copy=1;
                else
                    fact1=find_copies(nwhole, N,numpartners, wholep,Speclist,Adj2, t1,p_home,t2,copy);
                if(copy==1){
                    copy=0;
                    t3=edgepairs[i][1];
                    t4=edgepairs[w][1];
                    if(t3==t4)
                        copy=1;
                    else
                        fact2=find_copies(nwhole, N,numpartners, wholep,Speclist,Adj2, t3,p_home,t4,copy);
                    if(copy==1){
                        factb++;
                        //cout<<"Edgecopy found. Edge one connects: " << t2 << " to " << t4 << ". Edge two connects " << t1 << " to " << t3 << endl;
                    }
                }
            }
        }

        //cout<<"factb is now " << factb<<endl;
        
        //Find copies of i1 for pgen, since the reverse move is moving an edge back.
        if(i1==ip1)//A self-edge was split
            ndups=1;
        else{
            ndups=find_copies(nwhole,N,numpartners,wholep,Speclist,Adj2,i1,p_home,i2,copy);
            if(copy==1)
                ndups--;
            ndups2=find_copies(nwhole,N,numpartners,wholep,Speclist,Adj2,i1,p_home,ip1,copy);
            if(copy==1 && i2!=ip1)
                ndups--;
        }
        factb*=ndups;
        //cout<<"ncopies of i1: " << ndups << " taking into account that copy=" << copy<< ". factb="<<factb<<endl;
        
        pgen_ratio=pgen_ratio*1.0*factb/(1.0*factf);
    }    
    }//end if move possible

  delete [] Adj2;
  return chg;

}
