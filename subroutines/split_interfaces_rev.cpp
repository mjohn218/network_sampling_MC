#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int split_interfaces_rev(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, int maxni, int *selfflag, double &pf, double &pb)
{
  ////cout << "Inside split_interfaces" << endl;
  /*Split one interface into two (add an interface)
  */
  /*ADDED SELF*/
  int N=ninterfaces+1;//for the adjacency matrix, Make room for new interface
  /*create the adjacency matrix*/
  int ind;
  int i, j; 
  int t;
  int chg=1;
  double rnum1;
  int p1;
  
   
    numpartners[N-1]=0;//this is as yet new and undefined interface
    //Create Adj matrix
    for(i=0;i<N*N;i++)
      Adj[i]=0;
    for(i=0;i<N;i++){
      for(j=0;j<numpartners[i];j++){
          ind=Speclist[i][j];
          Adj[i*N+ind]=1;//i is the column, j the row, Symmetric matrix
      }
    }
      
      
    
    rnum1=rand_gsl()*1.0*nwhole;
    p1=int(rnum1);//index of first protein
    
  chg=0; //Check to see if an interface either has >1 partner or a self-edge
  for(i=0;i<wholep[p1].ninterface;i++){
      ind=wholep[p1].valiface[i];
    if(numpartners[ind]>1 || Adj[ind*N+ind]==1)chg=1;
  }

  if(ninterfaces==maxni)chg=0;//If maximum interfaces has been reached, don't make a change
  
    int edgeused[N];
  int flagsep;
  int sflag;
  double rnum2;
  int nswap;
  int w, i1, ni1;
  int nedge;
  double sum;
  int sp1, tmp, tmp2;
  int i2, inew;
  double psplit_f, pgen_f;
  double pcombine_b, pgen_b;
  double fact;
    int copyf, copyb, cp;
  if(chg!=0){
    //pick 1 interfaces off this protein p1, add another
    rnum1=rand_gsl()*1.0*wholep[p1].ninterface;
     w=int(rnum1);
     i1=wholep[p1].valiface[w];
    while(numpartners[i1]==1 && Speclist[i1][0]!=i1){
      rnum1=rand_gsl()*1.0*wholep[p1].ninterface;
      w=int(rnum1);
      i1=wholep[p1].valiface[w];
      ////cout << "numpartners[i1] is " << numpartners[i1] << endl;
    }
   
    i2=ninterfaces;//create new interface at edge of matrix

    //cout <<"Split interface: "<<i1<<" on protein: "<<p1<<" numpartner: "<<numpartners[i1]<<" create interface: "<<N-1<<" on protein: "<<p1<<endl;
    
    /*If we decide to perform this move
      then we need to calculate the probability of generating each move
    */
    ni1=wholep[p1].ninterface;
    nedge=numpartners[i1];
    /*keep ni1 fixed, add i2 to this protein*/
    wholep[p1].ninterface++;
    wholep[p1].valiface[ni1]=i2;
    p_home[i2]=p1;
      
      sp1=0; //for pgen_f
      for(i=0;i<ni1;i++){
          tmp=wholep[p1].valiface[i];
          if(numpartners[tmp]>1 || Adj[tmp*N+tmp]==1)
              sp1++;
      }
      
    if(nedge==1){
      //this is a self splitting
      ////cout <<"self splitting a single edge "<<endl;
        sum=1.0;
        psplit_f=1.0/(sp1*1.0*sum);
        //Also need to find copies. Any other ifaces with one self-edge and nothing else.
        copyf=0;
        for(i=0;i<ni1;i++){
            tmp=wholep[p1].valiface[i];
            if(numpartners[tmp]==1 && Adj[tmp*N+tmp]==1)
                copyf++;
        }
	//cout << "In move4: copy_f is " << copyf << " and psplit_f is " << psplit_f << endl;
        pgen_f=1.0*copyf*psplit_f;//*prob[split_move]
        pcombine_b=2.0/(1.0*(ni1+1)*ni1);
        //Need to see if there are any copies, e.g. any other pairs of split interfaces
        copyb=0;
                

       
       ////cout <<"protein: "<<p1<<" number of interfaces on i1: "<<ni1<<" and ways to choose 1 to split, afterwards sp1: "<<sp1<<endl;
    
       
       //pf=pgen_f;
       //pb=pgen_b;
       // //cout <<"in split: forward: "<<pgen_f<<" backward: "<<pgen_b<<" ratio: "<<pgen_ratio<<" protein: "<<p1<<endl;
      
       ////cout <<"nedge: "<<nedge<<" original num: "<<numpartners[i1]<<endl;
        
        numpartners[i1]=1;
        numpartners[i2]=1;
      
      /*update i2 and i1's partners: delete from i1, add to i2 */
      
        Adj[i1*N+i1]=0;//i1 no longer bound to self
            
        Adj[i2*N+i1]=1;//i2 is now bound to it.
        Adj[i1*N+i2]=1;
      //      //cout <<"self splitting, now i1: "<<i1<<" binds to i2: "<<i2<<endl;
        for(i=0;i<ni1+1;i++){
            for(j=i+1;j<ni1+1;j++){
                tmp=wholep[p1].valiface[i];
                tmp2=wholep[p1].valiface[j];
                if(Adj[tmp*N+tmp2]==1 && numpartners[tmp]==1 && numpartners[tmp2]==1)
                    copyb++;
            }
        }
        ////cout<< "In move 4: copy_b is " << copyb << endl;
        pgen_b=1.0*copyb*pcombine_b;//*prob[combine_move]
        pgen_ratio=pgen_b/pgen_f;
        
        
    /*update the speclist */
        for(i=0;i<N;i++){
            t=0;
            for(j=0;j<N;j++){
                if(Adj[i*N+j]==1){
                    Speclist[i][t]=j;
                    t++;
                }
            }
        }

    }else {
      //move involves more than just splitting a self interaction
      
        /*determine if this interface has a self interaction*/
        int pos_self=0;
        if(Adj[i1*N+i1]==1)pos_self=1;
      
        sum=0;        
        for(i=1;i<nedge;i++)
          sum+=binomial(nedge, i);
        if(pos_self==1){
          sum+=2;
	/*we have to also include possiblility of splitting self across 2 interfaces*/
            for(i=1;i<nedge-1;i++)
                sum+=binomial(nedge-1, i);
        }
        
        //cout <<"nedge: "<<nedge<<" binomial sum: "<<sum<<endl;
       //sp1=0;//because we can only split on interfaces with more than 1 partner
      //or self
        
      //Already calculated sp1 above, so commenting this part out (DH)
      //for(i=0;i<ni1;i++){
	//tmp=wholep[p1].valiface[i];
	//if(numpartners[tmp]>1 || Speclist[tmp][0]==tmp)
	  //sp1++;
      //}
      
      /*If one interaction of i1 is a self, there is 2 ways to split that, 
       if we choose that move, divide psplit_f by 2*/
      
      ////cout <<"protein: "<<p1<<" number of interfaces on i1: "<<ni1<<" and ways to choose 1 to split, afterwards sp1: "<<sp1<<endl;
      ////cout <<"total number of edges: "<<nedge<<" pos_self: "<<pos_self<<" sum: "<<sum<<endl;
      copyf=find_copies(nwhole, N, numpartners, wholep, Speclist, Adj, i1, p_home, i1, cp);//See if i1 has any copies on i1.
      
      psplit_f=1.0/(sp1*0.5*sum);
      pgen_f=1.0*copyf*psplit_f;//*prob[split_move]
      
      //cout << "In move 4. copy_f is " << copyf << " and psplit_f is " << psplit_f << " because sum was " << sum << endl;
        //Will have to alter pgen_f and calculate pgen_b later, since it gets more complicated. (DH)
      pcombine_b=2.0/(1.0*(ni1+1)*ni1);
        
        
        
        //For pgen_f: need to see if there are any edge copies. Afterwards, will need to see how they are distributed across the split interfaces.

        int edgecopy[N*N];
	for(i=0;i<N*N;i++)
	  edgecopy[i]=0;
        bool edge_copies=false;
        for(i=0;i<nedge;i++){
            tmp=Speclist[i1][i];
            if(tmp==i1)continue;
            for(j=i+1;j<nedge;j++){
                tmp2=Speclist[i1][j];
                if(tmp2==i1)continue;
                find_copies(nwhole, N, numpartners, wholep, Speclist, Adj, tmp, p_home, tmp2, cp);
                if(cp==1){
                    edgecopy[tmp*N+tmp2]=1;
                    edgecopy[tmp2*N+tmp]=1;
                    edge_copies=true;
                    //cout<<"Interface " << tmp << " is a copy of interface " << tmp2 << " on protein " << p_home[tmp]<<endl;
                }
            }
            t++;
        }
        
        
        
        
      
      ////cout <<"nedge: "<<nedge<<" original num: "<<numpartners[i1]<<endl;
    
      if(pos_self==1){
          flagsep=0;
          ////cout <<"interface : "<<i1<<" binds to itself "<<endl;
          nswap=split_type_self(nedge, flagsep);
          ////cout <<"numswap: "<<nswap<<" nedge : "<<nedge<<" flagsep: "<<flagsep<<endl;
	/*Need to choose how to split the selected interface*/
          if(flagsep==1){
	  //Csplititng self across
              Adj[i1*N+i1]=0;
              Adj[i2*N+i1]=1;
              Adj[i1*N+i2]=1;
        
        
              if(nswap==0){
	    /*in this case, we have only the self interaction being moved to the other interface*/
                  numpartners[i2]=1;
	    
              }else{
	    /*splitting self across, and also nswap edges*/
                  numpartners[i1]=nedge-nswap; //losing nswap edges, maintaining split self, which is not included in nswap
                  numpartners[i2]=nswap+1; //earning nswap, also earning split self +1
                  ////cout  <<" numpartners i1: "<<numpartners[i1]<<" numpartnesr[i2]: "<<numpartners[i2]<<endl;
                  ////cout <<"performing self split across, num to swap: "<<nswap<<endl;
	    
                  int iswap[nswap];
                  i=0;
	    /*select which edges to move off of i1 and onto i2, but don't select self, moving that automatically */
                  while(i<nswap){
                      rnum2=rand_gsl()*1.0*nedge;
                      w=int(rnum2);
                      iswap[i]=Speclist[i1][w];
                      while(iswap[i]==i1){
		//don't choose self
                          rnum2=rand_gsl()*1.0*nedge;
                          w=int(rnum2);
                          iswap[i]=Speclist[i1][w];
                      }
	      
                      if(i>0){
		//make sure you don't pick the same interface twice
                          sflag=0;
                          for(j=0;j<i;j++){
                              if(iswap[j]==iswap[i])
                                  sflag=1;
                          }
                          if(sflag==1)
                              i-=1;//reselect this interface
                      }
                      i++;
                  }
	    
	  		
	    /*update i2 and i1's partners: delete from i1, add to i2 */
                  for(t=0;t<nswap;t++){
                      inew=iswap[t];
                      Adj[i1*N+inew]=0;//i1 no longer bound to inew
                      Adj[inew*N+i1]=0;
	      
                      Adj[i2*N+inew]=1;//i2 is now bound to it.
                      Adj[inew*N+i2]=1;
                      ////cout <<"i1 no longer binds to partner  "<<inew<<" i2: "<<i2<< "binds to : "<<inew<<endl;
                  }//end looping over all the moved edges
          
          
              }//done checking for the nswap=0 case
          }else{
	  /*keeping the self interaction to itself*/
	  
              numpartners[i1]=nedge-nswap; //i1 has fewer partners
              numpartners[i2]=nswap;
              ////cout  <<" numpartners i1: "<<numpartners[i1]<<" numpartnesr[i2]: "<<numpartners[i2]<<endl;
              ////cout <<"num to swap: "<<nswap<<endl;
	  
              int iswap[nswap];
              i=0;
	  
	  /*select which edges to move off of i1 and onto i2*/
              while(i<nswap){
                  rnum2=rand_gsl()*1.0*nedge;
                  w=int(rnum2);
                  iswap[i]=Speclist[i1][w];
	    
                  if(i>0){
	      //make sure you don't pick the same interface twice
                      sflag=0;
                      for(j=0;j<i;j++){
                          if(iswap[j]==iswap[i])
                              sflag=1;
                      }
                      if(sflag==1)
                          i-=1;//reselect this interface
                  }
                  i++;
              }
	  
	  /*update i2 and i1's partners: delete from i1, add to i2 */
              for(t=0;t<nswap;t++){
                  inew=iswap[t];
                  if(inew==i1){
                      Adj[i1*N+inew]=0;//i1 no longer bound to self
                      Adj[i2*N+i2]=1;//i2 is now bound to self!!!
                      ////cout <<"in split, i1 no longer binds to self: "<<i1<<" i2 binds to self  : "<<i2<<endl;
                  }else{
                      Adj[i1*N+inew]=0;//i1 no longer bound to inew
                      Adj[inew*N+i1]=0;
	      
                      Adj[i2*N+inew]=1;//i2 is now bound to it.
                      Adj[inew*N+i2]=1;
                      ////cout <<"i1 no longer binds to partner  "<<inew<<" i2: "<<i2<< "binds to : "<<inew<<endl;
	      
                  }
	    
              }//end looping over all the moved edges
        
        
          }
	
          //Need to calculate pgen_b in here.
          
      }else{
	/*This protein does not have a self interaction*/ //pos_self=0
          nswap=split_type(nedge);
          ////cout<<"Done with split_type"<<endl;
	/*Need to choose how to split the selected interface*/
          numpartners[i1]=nedge-nswap; //i1 has fewer partners
          numpartners[i2]=nswap;
          ////cout  <<"Normal split,  numpartners i1: "<<numpartners[i1]<<" numpartnesr[i2]: "<<numpartners[i2]<<endl;
          ////cout <<"num to swap: "<<nswap<<endl;
          int iswap[nswap];
          i=0;
          int sflag;
          double rnum2;
	/*select which edges to move off of i1 and onto i2*/
          while(i<nswap){
              rnum2=rand_gsl()*1.0*nedge;
              w=int(rnum2);
              iswap[i]=Speclist[i1][w];
	  
              if(i>0){
	    //make sure you don't pick the same interface twice
                  sflag=0;
                  for(j=0;j<i;j++){
                      if(iswap[j]==iswap[i])
                          sflag=1;
                  }
                  if(sflag==1)
                      i-=1;//reselect this interface
              }
              i++;
          }
          ////cout<<"Done selecting edges off of i1 and i2"<<endl;
      
          int inew;
	/*update i2 and i1's partners: delete from i1, add to i2 */
          for(t=0;t<nswap;t++){
              inew=iswap[t];
              Adj[i1*N+inew]=0;//i1 no longer bound to inew
              Adj[inew*N+i1]=0;
	  
              Adj[i2*N+inew]=1;//i2 is now bound to it.
              Adj[inew*N+i2]=1;
	  // //cout <<"i1 no longer binds to partner  "<<inew<<" i2: "<<i2<< "binds to : "<<inew<<endl;
          }//end looping over all the moved edges
	/*test to see what the p_combine will be*/
	
         
    //Maggie's code. Outdated: must find duplicate ifaces instead.
    /*fact=1.0;
	if(numpartners[inew]==1){
	  if(numpartners[i2]==1 && p_home[inew]==p1)
	    fact=2.0;
	}
	if(numpartners[i1]==1){
	  for(i=0;i<N;i++){
	    if(Adj[i1*N+i]==1 &&p_home[i]==p1){
	      if(numpartners[i]==1)
		fact=2.0;
	    }
	  }
	}*/
	
	          
	
      }//done checking if possibility for a self move
        ////cout<<"Done with move. Now calculating pgen."<<endl;
        //Calculate pgen_b and pgen_ratio here.
        //Update the Speclist
        for(i=0;i<N;i++){
            t=0;
            for(j=0;j<N;j++){
                if(Adj[i*N+j]==1){
                    Speclist[i][t]=j;
                    t++;
                }
            }
        }

        ////cout<<"Speclist updated."<<endl;
        int fact1,fact2;
        cp=0;
        fact1=find_copies(nwhole, N, numpartners, wholep, Speclist, Adj, i1, p_home, i2, cp);//See if i1 has any copies
        fact2=find_copies(nwhole, N, numpartners, wholep, Speclist, Adj, i2, p_home, i1, cp);
        if(cp==1)
            fact=1.0*fact1*(fact1-1)/2.0;
        else
            fact=1.0*fact1*fact2;
        
        //Complications if there are >1 split self-edges, since not all combinations produce same result.
        //Check for split self-edges
        int nsplits=0;
        int splittesti1=0;
        int splittesti2=0;
        //int splitslist[wholep[p1].ninterface/2][2];
        for(i=0;i<wholep[p1].ninterface;i++){
            tmp=wholep[p1].valiface[i];
            for(j=i+1;j<wholep[p1].ninterface;j++){
                tmp2=wholep[p1].valiface[j];
                if(Adj[tmp*N+tmp2]==1){
                    //splitslist[nsplits][0]=tmp;
                    //splitslist[nsplits][1]=tmp2;
                    nsplits++;
                    if(tmp==i1 || tmp2==i1)
                        splittesti1=1;
                    if(tmp==i2 || tmp2==i2)
                        splittesti2=1;
                }
            }
        }
        if(nsplits>1){
            if(splittesti1==1 && splittesti2==1){
                //Are these two connected to each other?
                if(Adj[i1*N+i2]==1 && fact1==fact2){
                    //Then factb is just the number of i1 copies
                    fact=fact1;
                    if(cp==1)
                        fact=fact1/2;
                }
                else if(cp==1)
                    //If the two are not connected, then they are in separate splits. This will only be a problem if i1 and i2 are copies, in which case we need to subtract the cases where a self-edge is recreated
                    fact-=(fact1/2);
            }
        }
        
        //cout << "In move 4: (for b) fact1 is " << fact1 << " and fact2 is " << fact2 << " and copy is " << cp << endl;
        //cout << "In move 4: fact (for b) is " << fact << " and pcombine_b is " << pcombine_b << endl;
        pcombine_b*=fact;
        pgen_b=pcombine_b;//*prob[combine_move]
        
        //Now see if any edge copies were distributed. This affects pgen_f
        
        double fact_f2=1.0;
        //cout<<"edge_copies is " << edge_copies << " nedge is " << nedge << endl;
        
        if(edge_copies){

	  for(i=0;i<N;i++){
	    for(j=i;j<N;j++){
	      if(edgecopy[i*N+j]!=edgecopy[j*N+i]){
		//cout<<"Error with edgecopy. Not symmetric" <<endl;
		exit(1);
	      }
	    }
	  }

            for(i=0;i<N;i++)
                edgeused[i]=0;
            int i1cop[nedge];
            int i2cop[nedge];
            for(i=0;i<nedge;i++){
                i1cop[i]=0;
                i2cop[i]=0;
            }
            t=0;
            ////cout<< "i1 is " << i1 << " and i2 is " << i2 << endl;
            ////cout<<"Partners of i1:"<<endl;
            //for(i=0;i<numpartners[i1];i++)
	    //  //cout<<Speclist[i1][i]<<"\t";
            ////cout<<"\nPartners of i2:"<<endl;
            //for(i=0;i<numpartners[i2];i++)
	    //  //cout<<Speclist[i2][i]<<"\t";
            ////cout<<endl;
            for(i=0;i<numpartners[i1]+numpartners[i2];i++){
                if(i<numpartners[i1]){
                    tmp=Speclist[i1][i];
                    if(tmp!=i1 && tmp!=i2 && edgeused[tmp]==0)
                        i1cop[t]++;
                }
                else{
                    tmp=Speclist[i2][i-numpartners[i1]];
                    if(tmp!=i1 && tmp!=i2 && edgeused[tmp]==0)
                        i2cop[t]++;
                }
                if(tmp==i1 || tmp==i2 || edgeused[tmp]==1)continue;
                for(j=i+1;j<numpartners[i1]+numpartners[i2];j++){
                    if(j<numpartners[i1])
                        tmp2=Speclist[i1][j];
                    else
                        tmp2=Speclist[i2][j-numpartners[i1]];
                    if(tmp2==i1 || tmp2==i1)continue;
                    //cout<<"Was " << tmp << " a copy of " << tmp2 << "?"<<endl;
                    //cout<<edgecopy[tmp*N+tmp2] << " " << edgecopy[tmp2*N+tmp]<<endl;
                    if(edgecopy[tmp*N+tmp2]==1 || edgecopy[tmp2*N+tmp]==1){
                        //cout<<"Yes it was." << endl;
                        edgeused[tmp2]=1;
                        if(j<numpartners[i1])
                            i1cop[t]++;
                        else
                            i2cop[t]++;
                    }
                    
                }
                t++;
            }
            
            //For testing: print out i1cop and i2cop
            ////cout<<"i1cop:";
            //for(i=0;i<nedge;i++)
	    //  //cout<<"\t"<<i1cop[i];
            ////cout<<"\n"<<"i2cop:";
            //for(i=0;i<nedge;i++)
	    //  //cout<<"\t"<<i2cop[i];
            ////cout<<endl;
            
            for(i=0;i<nedge;i++){
                if(i1cop[i]>0 && i2cop[i]>0){
                    tmp=i1cop[i]+i2cop[i];
                    fact_f2*=binomial(tmp,i1cop[i]);
                    //cout<<"Increasing fact_f2 to " << fact_f2 << endl;
                }
            }
            
            if(cp==1)
                fact_f2*=0.5;
            
        }
        
        pgen_f*=fact_f2;
        //cout<<"fact_f2 is " << fact_f2 << " making the final pgen_f " << pgen_f << endl;
        pgen_ratio=pgen_b/pgen_f;
        pf=pgen_f;
        pb=pgen_b;
        
        
    }//done checking for a single edge move
  
    
    //Print out Adj for testing
    //for(i=0;i<N;i++){
    ////cout << i << "\t" << numpartners[i];
    //for(j=0;j<numpartners[i];j++)
    //    //cout << "\t" << Speclist[i][j];
    ////cout<<endl;
    //}
      
  
    
  } //done completing the move because chg==1
  ////cout <<"ninterfaces in move 4: "<<ninterfaces<<endl;
    
   
  return chg;
    
}
