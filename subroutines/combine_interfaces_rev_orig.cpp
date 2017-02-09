#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int combine_interfaces_rev_orig(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, double &pf, double &pb)
{

  /*In this case we are deleting an interface
   */
  /*ADDED SELF*/
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
  
  int chg=1;
  if(wholep[p1].ninterface==1)
    chg=0;
  
  
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
    
    //    cout <<"in combine, number of interfaces: "<<wholep[p1].ninterface<<" w: "<<w<<" w2: "<<w2<<endl;
    //cout <<"Delete Interface : "<<i1<<" from protein: "<<p1<< " combine with interface i2: "<<i2<<endl;

    /*If we decide to perform this move
      then we need to calculate the probability of generating each move
    */
    int ni1=wholep[p1].ninterface;
    double pcombine_f=2.0/(1.0*ni1*(ni1-1));
    /*if we have split self, there are in fact two ways to end up in the exact same final state*/
    int tmp;
    double fact=1.0;
    if(numpartners[i1]==1){
      for(i=0;i<ni1;i++){
	tmp=wholep[p1].valiface[i];
	if(Speclist[i1][0]==tmp && tmp!=i1){
	  //you bind to another interface on this protein
	  if(numpartners[tmp]==1)
	    fact=2.0;
	}
      }
    }
    if(numpartners[i2]==1){
      for(i=0;i<ni1;i++){
	tmp=wholep[p1].valiface[i];
	if(Speclist[i2][0]==tmp && tmp!=i2){
	  //you bind to another interface on this protein
	  if(numpartners[tmp]==1)
	    fact=2.0;
	}
      }
    }
    pcombine_f*=fact;
    double pgen_f=pcombine_f;//*prob[combine_move]
    /*now do reverse move probabilit*/
    int sp1=1;//one split is possible due to i1 and i2 combining

    for(i=0;i<ni1;i++){
      tmp=wholep[p1].valiface[i];
      if(tmp!=i1 && tmp!=i2){
	if(numpartners[tmp]>1 || Speclist[tmp][0]==tmp)
	  sp1++;
      }
    }
    //    cout <<"protein: "<<p1<<" number of interfaces on i1: "<<ni1<<" and ways to choose 1 to split, afterwards sp1: "<<sp1<<endl;
    /*If any of these interactions are self, have to add to the sum, more ways to split,
      also, if i1 is bound to i2, you'll overcount the number of total edges by 1 */
    int pos_sep=0;
    int flag_self=0;
    for(i=0;i<numpartners[i1];i++){
      if(Speclist[i1][i]==i1 )
	flag_self=1;
      if(Speclist[i1][i]==i2)
	pos_sep=1;//
    }
    for(i=0;i<numpartners[i2];i++){
      if(Speclist[i2][i]==i2 )
	flag_self=1;
    }
    
    int nedge2=numpartners[i1]+numpartners[i2];
    if(pos_sep==1)nedge2-=1;
    double sum=0;
    if(nedge2==1){
      sum=2.0;//because we want 1 and its multiplied by 0.5
    }else{
      for(i=1;i<nedge2;i++)
	sum+=binomial(nedge2, i);
      /*if this new combined interface has a self interaction, add in more ways to split*/
      if(pos_sep==1 || flag_self==1){
	sum+=2;//keep it the same configuration except move off the self
	for(i=1;i<nedge2-1;i++)
	  sum+=binomial(nedge2-1,i);
      }
    }

    //    cout <<"total number of edges after combination: "<<nedge2<<"  pos_sep: "<<pos_sep<<" flag_self: "<<flag_self<<" sum: "<<sum<<endl;
    double psplit_b=1.0/(sp1*0.5*sum);//*prob[split move]
    double pgen_b=psplit_b;
    pgen_ratio=pgen_b/pgen_f;
    pf=pgen_f;
    pb=pgen_b;
    //cout <<"in combine: forward: "<<pgen_f<<" backward: "<<pgen_b<<" ratio: "<<pgen_ratio<<" protein: "<<p1<<endl;
    
    int ni_one=numpartners[i1];
    
    
    
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
    
    int self_ind=-1;
    /*test to see if i1 binds to itself*/
    for(i=0;i<numpartners[i1];i++){
      if(Speclist[i1][i]==i1)
	self_ind=i;//self binding
    }
    numpartners[i2]+=numpartners[i1];    
    //update i2's partners to include i1's
    if(self_ind<0){
      //normal, no self 
      //careful you are not createing a self interaction
      for(i=0;i<ni_one;i++){
	ind=Speclist[i1][i];
	Adj[i1*N+ind]=0;//i1 no longer bound to its original partners
	Adj[ind*N+i1]=0;
	if(ind==i2)
	  numpartners[i2]--;//creating a self interaction from cross protein self
	
	Adj[i2*N+ind]=1;//add i1's original partners to i2
	Adj[ind*N+i2]=1;
	//	cout <<"move "<<ind<<" from i1: "<<i1<<" to i2: "<<i2<<endl;
      }
      
    }else{
      //careful that i1 binds to itself
      for(i=0;i<ni_one;i++){
	ind=Speclist[i1][i];
	if(ind==i1){
	  Adj[i1*N+ind]=0;//i1 no longer bound to itself
	  Adj[i2*N+i2]=1;//now i2 binds to itself!
	  //cout <<"move self from i1; "<<i1<<" to i2: "<<i2<<endl; 
	}else{
	  Adj[i1*N+ind]=0;//i1 no longer bound to its original partners
	  Adj[ind*N+i1]=0;
	  Adj[i2*N+ind]=1;//add i1's original partners to i2
	  Adj[ind*N+i2]=1;
	  //cout <<" move i1's former interface "<<ind<<" to i2: "<<i2<<endl;
	}
      }
    }
    /*Now that we have deleted an interface, need to copy the last interface into the slot for 
      i1's position. That is, relabel the N-1 interface as the i1 interface, to keep the numbering continuous
    */
    numpartners[i1]=numpartners[N-1];
    /*If the N-1 protein binds to itself, be careful */
    //if N-1 binds to itself, now i1 binds to itself

    for(i=0;i<N;i++){
      Adj[i1*N+i]=Adj[(N-1)*N+i];
      Adj[i*N+i1]=Adj[i*N+N-1];
    }
    if(Adj[(N-1)*N+N-1]==1)Adj[i1*N+i1]=1;
    int p2=p_home[N-1];
    
    p_home[i1]=p2;
    //    cout<<"i1: "<<i1<<" p_home: "<<p2<<" N-1: "<<N-1<<endl;
    for(i=0;i<wholep[p2].ninterface;i++){
      if(wholep[p2].valiface[i]==(N-1))
	wholep[p2].valiface[i]=i1;
    }

    /*update the speclist */
    for(i=0;i<N-1;i++){
      t=0;
      for(j=0;j<N-1;j++){
	if(Adj[i*N+j]==1){
	  Speclist[i][t]=j;
	  t++;
	}
      }
    }
  } //done completing the move because chg==1
  return chg;
}
