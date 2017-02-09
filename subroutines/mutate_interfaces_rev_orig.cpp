#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int mutate_interfaces_rev_orig(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *selfflag)
{

  /*In this case, choose first a protein with multiple interfaces, 
    remove one (and therefore add edges to the original interface
    Then choose another protein with one interface with multiple partners,
    and add it there.
  */
  int N=ninterfaces;//for the adjacency matrix
  /*create the adjacency matrix*/
  int ind;
  int i, j; 
  int t;
  for(i=0;i<N*N;i++)
    Adj[i]=0;
  for(i=0;i<N;i++){
    for(j=0;j<numpartners[i];j++){
      ind=Speclist[i][j];
      Adj[i*N+ind]=1;//i is the column, j the row, Symmetric matrix
    }
  }
  
  double rnum1=trand()*1.0*nwhole/RAND_MAX;
  int p1=int(rnum1);//index of first protein
  rnum1=trand()*1.0*nwhole/RAND_MAX;
  int p2=int(rnum1);
  while(p1==p2){
     rnum1=trand()*1.0*nwhole/RAND_MAX;
     p2=int(rnum1);
  }
  int type=2;//type is either interface combine/donate=1, or interface split/accept=2
  int chg=1;
  double ph1=1;
  double ph2=1;

  if(wholep[p1].ninterface==1){
    //p1 would have to be acceptor of interface
    ph1=1;
    if(wholep[p1].ninterface<ppi[p1].nppartner){
      chg=1;
      type=2;
      if(wholep[p2].ninterface==1)
	chg=0;//p2 can't donate
    }else if(wholep[p1].ninterface==ppi[p1].nppartner && selfflag[p1]==1){
      //then we CAN be the second protein either;
      chg=1;
      type=2;
      if(wholep[p2].ninterface==1)
	chg=0;//p2 can't donate
    }else{
      chg=0;//p1 cannot donate, and p1 cannot accept
    }
  }else{
    type=1;//p1 is the donator
    if(wholep[p2].ninterface==ppi[p2].nppartner){
      if(selfflag[p2]==0){
	//p2 cannot accept
	ph1=1;
	if(wholep[p2].ninterface==1){
	  chg=0; //p2 can't accept or donate
	}else{
	  //p2 could be the donator
	  if(wholep[p1].ninterface==ppi[p1].nppartner){
	    if(selfflag[p1]==0){
	      //p1 can't accept
	      chg=0;
	    }else{
	      type=2;
	    }
	  }else if(wholep[p1].ninterface<ppi[p1].nppartner){
	    type=2;//p1 accepts, p2 is the donor
	  }else
	    chg=0;//p1 has more interfaces than partners, no acceptance allowed!
	}
      }
      //else, p2 can accept
    }else if(wholep[p2].ninterface<ppi[p2].nppartner){
      type=1;
      //p2 can be the acceptor
      //check if we can perform the other way, just for pgen
      if(wholep[p2].ninterface!=1){
	//p2 can also donate, can p1 accept?
	if(wholep[p1].ninterface!=ppi[p1].nppartner)
	  ph1=0.5;
      }
      
    }else{
      //p2 cannot be the acceptor, more interfaces than partners
      //can it donate?
      if(wholep[p2].ninterface==1){
	chg=0;
      }else{
	//p2 can donate
	//check if p1 can accept
	ph1=1;
	if(wholep[p1].ninterface<ppi[p1].nppartner){
	  chg=1;
	  type=2;
	
	}else if(wholep[p1].ninterface==ppi[p1].nppartner && selfflag[p1]==1){
	  //then we CAN be the second protein either;
	  chg=1;
	  type=2;
	
	}else{
	  chg=0;//p1 cannot accept
	}
      }
    } 
    
    
  }
  int tmp;
  if(type==2){
    tmp=p1;//we want p1 to be an interface donor (combines 2 into 1)
    p1=p2;
    p2=tmp;
  }
  
  if(chg!=0){
    //pick 2 interfaces off this protein p1, one will be deleted
    rnum1=trand()*1.0*wholep[p1].ninterface/RAND_MAX;
    double rnum2=trand()*1.0*wholep[p1].ninterface/RAND_MAX;
    int w=int(rnum1);
    int w2=int(rnum2);
    
    while(w==w2){
      rnum2=trand()*1.0*wholep[p1].ninterface/RAND_MAX;
      w2=int(rnum2);
    }
    int i1=wholep[p1].valiface[w];
    int i2=wholep[p1].valiface[w2];
    
    
    /*Now choose one of the interfaces on this protein 2 that has multiple partners*/
    int ni_p2=wholep[p2].ninterface;
    int i3;
    rnum2=trand()*1.0*ni_p2/RAND_MAX;
    w=int(rnum2);
    i3=wholep[p2].valiface[w];
    while(numpartners[i3]==1 && Speclist[i3][0]!=i3){
      rnum2=trand()*1.0*ni_p2/RAND_MAX;
      w=int(rnum2);
      i3=wholep[p2].valiface[w];
    }  
    //        cout <<"Mutate Interface : "<<i1<<" from protein: "<<p1<< " to protein: "<<p2<<endl;
    //cout <<"and merge onto  interface: "<<i2<<" and split from interface: "<<i3<<" on: "<<p2<<endl;
    /*If we decide to perform this move
      then we need to calculate the probability of generating each move
    */
    
    if(numpartners[i3]==1){
      //in this case we are splitting a self interaction across 2 interfaces
      int ni1=wholep[p1].ninterface;
      double pcombine_f=2.0/(1.0*ni1*(ni1-1));
      int ni2=wholep[p2].ninterface;
      int sp2=0;//because we can only split on interfaces with more than 1 partner
      int tmp;
      for(i=0;i<ni2;i++){
	tmp=wholep[p2].valiface[i];
	if(numpartners[tmp]>1)
	  sp2++;
      }
      //      int nedge=numpartners[i3];
      double sum=1.0;
      //for(i=1;i<nedge;i++)
      //	sum+=binomial(nedge, i);
      
      double psplit_f=1.0/(sp2*1.0*sum);//only one way to split the self interaction
      double pgen_f=pcombine_f*psplit_f*ph1;
      //    cout <<"in swap interface move, psplit_For: "<<psplit_f<<" pcomb_f: "<<pcombine_f<<endl;
      /*Now calculate reverse move probability */
      int nedge2=numpartners[i1]+numpartners[i2];
      sum=0;
      for(i=1;i<nedge2;i++)
	sum+=binomial(nedge2, i);
      int sp1=1;//one split is possible due to i1 and i2 combining
      for(i=0;i<ni1;i++){
	tmp=wholep[p1].valiface[i];
	if(tmp!=i1 && tmp!=i2){
	  if(numpartners[tmp]>1)
	    sp1++;
	}
      }
      /*on reverse, p1 accepts, and p2 donates, check whether the 
	other direction is possible*/
      if(ni1-1!=1){
	//p1 can also donate, check if p2 can accept (or split)
	if(ppi[p2].nppartner!=ni2+1)
	  ph2=0.5;
      }
      double psplit_b=1.0/(sp1*0.5*sum);
      double pcombine_b=2.0/(1.0*(ni2+1)*ni2);
      double pgen_b=pcombine_b*psplit_b*ph2;
      //cout <<"in swap interface move, psplit_back: "<<psplit_b<<" pcomb_b: "<<pcombine_b<<" sp1: "<<sp1<<endl;
      pgen_ratio=pgen_b/pgen_f;
      //cout <<"pgen_ratio: "<<pgen_ratio<<endl;
      /*Now move interface i1 to protein p2, 
	update wholep[p2]
	speclist[i1], i3 has fewer partners
	update speclist[i3]*/
      wholep[p2].ninterface++;
      wholep[p2].valiface[ni_p2]=i1;
      
      
      int ni_3=numpartners[i3];//1
      int ni_one=numpartners[i1];
      numpartners[i2]+=numpartners[i1];
      
      //      int nswap=split_type(ni_3);
      /*Need to choose how to split the selected interface*/
      numpartners[i1]=1;//nswap; //i1 has different partners now, it is inew..
      numpartners[i3]=1;//ni_3-nswap;
      
      /*protein p1 has fewer interfaces (got rid of i1), 
	interface i2 has more partners
	update wholep[p1]
	Speclist[i2]*/
      int ni=wholep[p1].ninterface;
      wholep[p1].ninterface--;
      t=0;
      for(i=0;i<ni;i++){
	if(wholep[p1].valiface[i]!=i1){
	  wholep[p1].valiface[t]=wholep[p1].valiface[i];
	  t++;
	}
      }
      
      //update i2's partners to include i1's, careful if i1 binds to itself
      int flag=1;
      for(i=0;i<ni_one;i++){
	ind=Speclist[i1][i];
	flag=1;
	if(ind==i1){
	  Adj[i1*N+i1]=0;//i1 no longer bound to its self
	  Adj[i2*N+i2]=1;//now i2 binds to itself
	}else{
	  Adj[i1*N+ind]=0;//i1 no longer bound to its original partners
	  Adj[ind*N+i1]=0;
	  Adj[i2*N+ind]=1;//add i1's original partners to i2
	  Adj[ind*N+i2]=1;

	}
	
      }
      int inew;
      /*update i3 and i1's partners: delete from i3, add to i1 */
      Adj[i3*N+i3]=0;//i3 no longer bound to self
      Adj[i1*N+i3]=1;//i1 is now bound to i3.
      Adj[i3*N+i1]=1;
      //            cout <<"in mutate, self split, now i3: "<<i3<< " binds to i1: "<<i1<< " on p2: "<<p2<<endl;
      
      
    }else {
      
      //splitting separate edges, however, it is still possible that one is a self
      int ni1=wholep[p1].ninterface;
      double pcombine_f=2.0/(1.0*ni1*(ni1-1));
      int ni2=wholep[p2].ninterface;
      int sp2=0;//because we can only split on interfaces with more than 1 partner
      int tmp;
      for(i=0;i<ni2;i++){
	tmp=wholep[p2].valiface[i];
	if(numpartners[tmp]>1)
	  sp2++;
      }
      int nedge=numpartners[i3];
      double sum=0;
      for(i=1;i<nedge;i++)
	sum+=binomial(nedge, i);
      
      double psplit_f=1.0/(sp2*0.5*sum);
      double pgen_f=pcombine_f*psplit_f*ph1;
      //    cout <<"in swap interface move, psplit_For: "<<psplit_f<<" pcomb_f: "<<pcombine_f<<endl;
      /*Now calculate reverse move probability */
      int nedge2=numpartners[i1]+numpartners[i2];
      sum=0;
      for(i=1;i<nedge2;i++)
	sum+=binomial(nedge2, i);
      int sp1=1;//one split is possible due to i1 and i2 combining
      for(i=0;i<ni1;i++){
	tmp=wholep[p1].valiface[i];
	if(tmp!=i1 && tmp!=i2){
	  if(numpartners[tmp]>1)
	    sp1++;
	}
      }
      /*on reverse, p1 accepts, and p2 donates, check whether the 
	other direction is possible*/
      if(ni1-1!=1){
	//p1 can also donate, check if p2 can accept (or split)
	if(ppi[p2].nppartner!=ni2+1)
	  ph2=0.5;
      }
      double psplit_b=1.0/(sp1*0.5*sum);
      double pcombine_b=2.0/(1.0*(ni2+1)*ni2);
      double pgen_b=pcombine_b*psplit_b*ph2;
      //cout <<"in swap interface move, psplit_back: "<<psplit_b<<" pcomb_b: "<<pcombine_b<<" sp1: "<<sp1<<endl;
      pgen_ratio=pgen_b/pgen_f;
      //cout <<"pgen_ratio: "<<pgen_ratio<<endl;
      /*Now move interface i1 to protein p2, 
	update wholep[p2]
	speclist[i1], i3 has fewer partners
	update speclist[i3]*/
      wholep[p2].ninterface++;
      wholep[p2].valiface[ni_p2]=i1;
      
      
      int ni_3=numpartners[i3];
      int ni_one=numpartners[i1];
      numpartners[i2]+=numpartners[i1];
      
      int nswap=split_type(ni_3);
      /*Need to choose how to split the selected interface*/
      numpartners[i1]=nswap; //i1 has different partners now, it is inew..
      numpartners[i3]=ni_3-nswap;
      
      /*protein p1 has fewer interfaces (got rid of i1), 
	interface i2 has more partners
	update wholep[p1]
	Speclist[i2]*/
      int ni=wholep[p1].ninterface;
      wholep[p1].ninterface--;
      t=0;
      for(i=0;i<ni;i++){
	if(wholep[p1].valiface[i]!=i1){
	  wholep[p1].valiface[t]=wholep[p1].valiface[i];
	  t++;
	}
      }
      int iswap[nswap];
      i=0;
      int sflag;
      /*select which edges to move off of i3 and onto i1*/
      while(i<nswap){
	rnum2=trand()*1.0*ni_3/RAND_MAX;
	w=int(rnum2);
	iswap[i]=Speclist[i3][w];
	
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
      
      //update i2's partners to include i1's
      int flag=1;
      for(i=0;i<ni_one;i++){
	ind=Speclist[i1][i];
	flag=1;
	if(ind!=i3){
	  if(ind==i1){
	    Adj[i1*N+i1]=0;//i1 no longer bound to its self
	    Adj[i2*N+i2]=1;//now i2 binds to itself

	  }else{
	    if(ind==i2)numpartners[i2]--;// i2 binds to itself, from binding across protein self
	    Adj[i1*N+ind]=0;//i1 no longer bound to its original partners
	    Adj[ind*N+i1]=0;
	    Adj[i2*N+ind]=1;//add i1's original partners to i2
	    Adj[ind*N+i2]=1;
	  }
	}else{
	  //check and make sure i1 is not equal to any of the swapped edges
	  for(t=0;t<nswap;t++){
	    if(i1==iswap[t])
	      flag=0;
	  }
	  if(flag==1){
	    /*i1 would still be attached to i3, so add i3 to i2  */
	    Adj[i1*N+ind]=0;//i1 no longer bound to its original partners
	    Adj[ind*N+i1]=0;
	    Adj[i2*N+ind]=1;//add i1's original partners to i2
	    Adj[ind*N+i2]=1;
	  }
	}//end checking about i1=inew
      }
      int inew;
      /*update i3 and i1's partners: delete from i3, add to i1 */
      for(t=0;t<nswap;t++){
	inew=iswap[t];
	if(inew==i3){
	  Adj[i3*N+inew]=0;//i3 no longer bound to self
	  Adj[i1*N+i1]=1;//i1 is now bound to self!!!
	  
	}else{
	  
	  Adj[i3*N+inew]=0;//i3 no longer bound to inew
	  Adj[inew*N+i3]=0;
	  if(i1==inew){
	    /*special case where i1 is bound first to i3*/
	    Adj[i1*N+i2]=1;//because i1 has been replaced on p1 by i2
	    Adj[i2*N+i1]=1;
	    
	  }else{
	    Adj[i1*N+inew]=1;//i1 is now bound to it.
	    Adj[inew*N+i1]=1;
	    
	  }
	}
      }//end looping over all the moved edges
      
      
      
    }
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
  return chg;
}
