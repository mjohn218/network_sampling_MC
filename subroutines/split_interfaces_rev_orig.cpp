#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int split_interfaces_rev_orig(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, int maxni, int *selfflag, double &pf, double &pb)
{

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
  if(ninterfaces==maxni)chg=0;
  if(chg!=0){
    numpartners[N-1]=0;//this is as yet new and undefined interface
    for(i=0;i<N*N;i++)
      Adj[i]=0;
    for(i=0;i<N;i++){
      for(j=0;j<numpartners[i];j++){
	ind=Speclist[i][j];
	Adj[i*N+ind]=1;//i is the column, j the row, Symmetric matrix
      }
    }
    
    rnum1=trand()*1.0*nwhole/RAND_MAX;
    p1=int(rnum1);//index of first protein
    
    chg=0;
    if(wholep[p1].ninterface<ppi[p1].nppartner)
      chg=1;//some are shared, because of self, it is possible for ninterface to be >npartner
    if(wholep[p1].ninterface==ppi[p1].nppartner &&selfflag[p1]==1){
      chg=1;
    }
  }
  int flagsep;
  int sflag;
  double rnum2;
  int nswap;
  int w, i1, ni1;
  int nedge;
  double sum;
  int sp1, tmp;
  int i2, inew;
  double psplit_f, pgen_f;
  double pcombine_b, pgen_b;
  double fact;
  if(chg!=0){
    //pick 1 interfaces off this protein p1, add another
    rnum1=trand()*1.0*wholep[p1].ninterface/RAND_MAX;
     w=int(rnum1);
     i1=wholep[p1].valiface[w];
    while(numpartners[i1]==1 && Speclist[i1][0]!=i1){
      rnum1=trand()*1.0*wholep[p1].ninterface/RAND_MAX;
      w=int(rnum1);
      i1=wholep[p1].valiface[w];
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
      
    if(nedge==1){
      //this is a self splitting
      //cout <<"self splitting a single edge "<<endl;
      sum=1.0;
      sp1=1;//because we can only split on interfaces with more than 1 partner
      tmp;
      for(i=0;i<ni1;i++){
	tmp=wholep[p1].valiface[i];
	if(numpartners[tmp]>1)
	  sp1++;
      }
      
       psplit_f=1.0/(sp1*1.0*sum);
       pgen_f=psplit_f;//*prob[split_move]
       pcombine_b=2.0/(1.0*(ni1+1)*ni1);
       pgen_b=pcombine_b;//*prob[combine_move]
       //cout <<"protein: "<<p1<<" number of interfaces on i1: "<<ni1<<" and ways to choose 1 to split, afterwards sp1: "<<sp1<<endl;

    
       pgen_ratio=pgen_b/pgen_f;
       pf=pgen_f;
       pb=pgen_b;
      //    cout <<"in split: forward: "<<pgen_f<<" backward: "<<pgen_b<<" ratio: "<<pgen_ratio<<" protein: "<<p1<<endl;
      
      //cout <<"nedge: "<<nedge<<" original num: "<<numpartners[i1]<<endl;
            
      /*Need to choose how to split the selected interface*/
      numpartners[i1]=1;//nedge-nswap; //i1 has fewer partners
      numpartners[i2]=1;//nswap;
      //    cout  <<" numpartners i1: "<<numpartners[i1]<<" numpartnesr[i2]: "<<numpartners[i2]<<endl;  
      
      
      
      /*update i2 and i1's partners: delete from i1, add to i2 */
      
      Adj[i1*N+i1]=0;//i1 no longer bound to self
            
      Adj[i2*N+i1]=1;//i2 is now bound to it.
      Adj[i1*N+i2]=1;
      //      cout <<"self splitting, now i1: "<<i1<<" binds to i2: "<<i2<<endl;
    
    }else {
      //move involves more than just splitting a self interaction
      
      sum=0;
      for(i=1;i<nedge;i++)
	sum+=binomial(nedge, i);
      
      //cout <<"nedge: "<<nedge<<" binomial sum: "<<sum<<endl;
      
      /*determine if this interface has a self interaction*/
      int pos_self=0;
      for(i=0;i<numpartners[i1];i++){
	if(Speclist[i1][i]==i1){
	  pos_self=1;
	}
      }
      if(pos_self==1){
	sum+=2;
	/*we have to also include possiblility of splitting self across 2 interfaces*/
	for(i=1;i<nedge-1;i++)
	  sum+=binomial(nedge-1, i);
      }
       sp1=0;//because we can only split on interfaces with more than 1 partner
      //or self
      
      for(i=0;i<ni1;i++){
	tmp=wholep[p1].valiface[i];
	if(numpartners[tmp]>1 || Speclist[tmp][0]==tmp)
	  sp1++;
      }
      
      /*If one interaction of i1 is a self, there is 2 ways to split that, 
       if we choose that move, divide psplit_f by 2*/
      
      //cout <<"protein: "<<p1<<" number of interfaces on i1: "<<ni1<<" and ways to choose 1 to split, afterwards sp1: "<<sp1<<endl;
      // cout <<"total number of edges: "<<nedge<<" pos_self: "<<pos_self<<" sum: "<<sum<<endl;
          
      psplit_f=1.0/(sp1*0.5*sum);
      pgen_f=psplit_f;//*prob[split_move]
      
      pcombine_b=2.0/(1.0*(ni1+1)*ni1);
      pgen_b=pcombine_b;//*prob[combine_move]
      
      pgen_ratio=pgen_b/pgen_f;
      pf=pgen_f;
      pb=pgen_b;
      // cout <<"in split: forward: "<<pgen_f<<" backward: "<<pgen_b<<" ratio: "<<pgen_ratio<<" protein: "<<p1<<endl;
      
      //cout <<"nedge: "<<nedge<<" original num: "<<numpartners[i1]<<endl;
      
      if(pos_self==1){
	
	flagsep=0;
	//	cout <<"interface : "<<i1<<" binds to itself "<<endl;
	nswap=split_type_self(nedge, flagsep);
	//cout <<"numswap: "<<nswap<<" nedge : "<<nedge<<" flagsep: "<<flagsep<<endl;
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
	    // cout  <<" numpartners i1: "<<numpartners[i1]<<" numpartnesr[i2]: "<<numpartners[i2]<<endl;  
	    //cout <<"performing self split across, num to swap: "<<nswap<<endl;
	    
	    int iswap[nswap];
	    i=0;
	    /*select which edges to move off of i1 and onto i2, but don't select self, moving that automatically */
	    while(i<nswap){
	      rnum2=trand()*1.0*nedge/RAND_MAX;
	      w=int(rnum2);
	      iswap[i]=Speclist[i1][w];
	      while(iswap[i]==i1){
		//don't choose self
		rnum2=trand()*1.0*nedge/RAND_MAX;
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
	      //cout <<"i1 no longer binds to partner  "<<inew<<" i2: "<<i2<< "binds to : "<<inew<<endl;
	    }//end looping over all the moved edges
	  
	  }//done checking for the nswap=0 case
	}else{
	  /*keeping the self interaction to itself*/
	  
	  numpartners[i1]=nedge-nswap; //i1 has fewer partners
	  numpartners[i2]=nswap;
	  //cout  <<" numpartners i1: "<<numpartners[i1]<<" numpartnesr[i2]: "<<numpartners[i2]<<endl;  
	  //cout <<"num to swap: "<<nswap<<endl;
	  
	  int iswap[nswap];
	  i=0;
	  
	  /*select which edges to move off of i1 and onto i2*/
	  while(i<nswap){
	    rnum2=trand()*1.0*nedge/RAND_MAX;
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
	      //	  cout <<"in split, i1 no longer binds to self: "<<i1<<" i2 binds to self  : "<<i2<<endl;
	    }else{
	      Adj[i1*N+inew]=0;//i1 no longer bound to inew
	      Adj[inew*N+i1]=0;
	      
	      Adj[i2*N+inew]=1;//i2 is now bound to it.
	      Adj[inew*N+i2]=1;
	      //cout <<"i1 no longer binds to partner  "<<inew<<" i2: "<<i2<< "binds to : "<<inew<<endl;
	      
	    }
	    
	  }//end looping over all the moved edges
	}
	
      }else{
	/*This protein does not have a self interaction*/
	nswap=split_type(nedge);
	/*Need to choose how to split the selected interface*/
	numpartners[i1]=nedge-nswap; //i1 has fewer partners
	numpartners[i2]=nswap;
	//cout  <<"Normal split,  numpartners i1: "<<numpartners[i1]<<" numpartnesr[i2]: "<<numpartners[i2]<<endl;  
	//      cout <<"num to swap: "<<nswap<<endl;
	int iswap[nswap];
	i=0;
	int sflag;
	double rnum2;
	/*select which edges to move off of i1 and onto i2*/
	while(i<nswap){
	  rnum2=trand()*1.0*nedge/RAND_MAX;
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
	
      
	int inew;
	/*update i2 and i1's partners: delete from i1, add to i2 */
	for(t=0;t<nswap;t++){
	  inew=iswap[t];
	  Adj[i1*N+inew]=0;//i1 no longer bound to inew
	  Adj[inew*N+i1]=0;
	  
	  Adj[i2*N+inew]=1;//i2 is now bound to it.
	  Adj[inew*N+i2]=1;
	  //cout <<"i1 no longer binds to partner  "<<inew<<" i2: "<<i2<< "binds to : "<<inew<<endl;
	}//end looping over all the moved edges
	/*test to see what the p_combine will be*/
	/*update the speclist */
	fact=1.0;
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
	}
	
	pcombine_b*=fact;
	pgen_b=pcombine_b;//*prob[combine_move]
	
	pgen_ratio=pgen_b/pgen_f;
	pf=pgen_f;
	pb=pgen_b;
	
      }//done checking if possibility for a self move
    }//done checking for a single edge move
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
    
  } //done completing the move because chg==1
  //  cout <<"ninterfaces in move 4: "<<ninterfaces<<endl;
  return chg;
}
