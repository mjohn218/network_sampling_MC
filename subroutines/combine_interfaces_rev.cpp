#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int combine_interfaces_rev(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, double &pf, double &pb)
{

  /*In this case we are deleting an interface
   */
  /*ADDED SELF*/

  int N=ninterfaces;//for the adjacency matrix
  /*create the adjacency matrix*/
  int ind;
  int i, j, k; 
  int t;
  for(i=0;i<N*N;i++)
    Adj[i]=0;
  for(i=0;i<N;i++){
    for(j=0;j<numpartners[i];j++){
      ind=Speclist[i][j];
      Adj[i*N+ind]=1;//i is the column, j the row, Symmetric matrix
    }
  }
    
    
  
  double rnum1=rand_gsl()*1.0*nwhole;
  int p1=int(rnum1);//index of first protein
  
  int chg=1;
  if(wholep[p1].ninterface==1){
    //cout << "In move 3: chg=0 because only one interface\n";
    chg=0;
  }
  
  //  cout << "Move 3: Checkpoint 1" << endl;
  if(chg!=0){
    //pick 2 interfaces off this protein p1, one will be deleted
    rnum1=rand_gsl()*1.0*wholep[p1].ninterface;
    double rnum2=rand_gsl()*1.0*wholep[p1].ninterface;
    int w=int(rnum1);
    int w2=int(rnum2);
    
    while(w==w2){
      rnum2=rand_gsl()*1.0*wholep[p1].ninterface;
      w2=int(rnum2);
    }
    int i1=wholep[p1].valiface[w];
    int i2=wholep[p1].valiface[w2];
    
    //cout <<"in combine, number of interfaces: "<<wholep[p1].ninterface<<" w: "<<w<<" w2: "<<w2<<endl;
    //cout <<"Delete Interface : "<<i1<<" from protein: "<<p1<< " combine with interface i2: "<<i2<<endl;

    /*If we decide to perform this move
      then we need to calculate the probability of generating each move
    */
    int ni1=wholep[p1].ninterface;
    double pcombine_f=2.0/(1.0*ni1*(ni1-1));
    /*if we have split self, there are in fact two ways to end up in the exact same final state*/
    int tmp,tmp2;
    //double fact=1.0;
    //Maggie's code
    /*if(numpartners[i1]==1){
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
      }*/
    //DH code. The point is that if i1 or i2 is identical with a third interface on the protein, fact=2. But shared partners complicate things. Basically, if 2 or more ifaces
    //are identical, then there are multiple ways of getting the result
    //first check to see if there are "duplicates" of i1
    int ndups1, ndups2; //Since i1 itself counts
    int fact;
    int copy=0;
    int copy_orig=0;
    ndups1=find_copies(nwhole, N, numpartners, wholep, Speclist, Adj, i1, p_home, i2, copy);
    ndups2=find_copies(nwhole, N, numpartners, wholep, Speclist, Adj, i2, p_home, i1, copy);
    //cout << "ndups1 is " << ndups1 << ", ndups2 is " << ndups2 << ", and copy is " << copy << endl;
      if(copy==1){//i1 and i2 are copies
        fact=ndups1*(ndups1-1)/2;
          copy_orig=1;
      }
    else
        fact=ndups1*ndups2;
    
//More complications if there are >1 split self-edges.
      int nsplits=0;
      int splittesti1=0;
      int splittesti2=0;
      int splitslist[wholep[p1].ninterface/2][2];
      for(i=0;i<wholep[p1].ninterface;i++){
          tmp=wholep[p1].valiface[i];
          for(j=i+1;j<wholep[p1].ninterface;j++){
              tmp2=wholep[p1].valiface[j];
              if(Adj[tmp*N+tmp2]==1){
                  splitslist[nsplits][0]=tmp;
                  splitslist[nsplits][1]=tmp2;
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
              if(Adj[i1*N+i2]==1 && ndups1==ndups2){
                  //Then factf is just the number of i1 copies
                  fact=ndups1;
                  if(copy==1)
                      fact=ndups1/2;
              }
              else if(copy==1)
                  //If the two are not connected, then they are in separate splits. This will only be a problem if i1 and i2 are copies, in which case we need to subtract the cases where a self-edge is recreated
                  fact-=(ndups1/2);
          }
      }
      
      
    //cout << "In move 3: fact_f is " << fact << " and pcombine_f is " << pcombine_f << endl;
    double pgen_f=double(fact)*pcombine_f;//*prob[combine_move]
    //cout << "Move 3: checkpoint 2" << endl;
      
    int sp1=1;//one split is possible due to i1 and i2 combining

    for(i=0;i<ni1;i++){
      tmp=wholep[p1].valiface[i];
      if(tmp!=i1 && tmp!=i2){
          if(numpartners[tmp]>1 || Speclist[tmp][0]==tmp)
              sp1++;
      }
    }
    ////cout <<"protein: "<<p1<<" number of interfaces on i1: "<<ni1<<" and ways to choose 1 to split, afterwards sp1: "<<sp1<<endl;
    /*If any of these interactions are self, have to add to the sum, more ways to split,
      also, if i1 is bound to i2, you'll overcount the number of total edges by 1 */
    int pos_sep=0;
      int fs1,fs2;
      fs1=0;
      fs2=0;
    int flag_self=0;
    for(i=0;i<numpartners[i1];i++){
      if(Speclist[i1][i]==i1 )
          fs1=1;//Self-edge present
      if(Speclist[i1][i]==i2)
          pos_sep=1;//i1 binds to i2
    }
    for(i=0;i<numpartners[i2];i++){
      if(Speclist[i2][i]==i2 )
          fs2=1;
    }
      
      if(fs1==1 || fs2==1)
          flag_self=1;
      
      if(fs1==1 && fs2==1){
	//cout<< "In Move 3: making chg=0 because combining two self-edges\n";
          chg=0; //Since this will result in two self-edges
      }
    //Look for shared partners. If i1 and i2 have a shared partner, DO NOT DO THE MOVE since it is irreversible. Added by DH
    int nshared=0;
    //int tmpshare;
    //int patch=0;
    for(i=0;i<numpartners[i1];i++){
      for(j=0;j<numpartners[i2];j++){
	if(Speclist[i1][i]==Speclist[i2][j]){
	  //cout << "Share found: iface " << Speclist[i1][i] << " connected to both iface "<<i1<<" and "<< i2<< endl;
	    nshared++;
    }
      }
    }
    if(nshared>0)
      chg=0;
    
    //    cout << "Move 3: checkpoint 3" << endl;
  if(chg!=0){//no shared partners. May continue move.
  
  
    int nedge2=numpartners[i1]+numpartners[i2];
      //cout<<"In move 3: nedge2 is " << nedge2 << endl;

        //Need to record which edge was originally on each interface. May be needed for pgen_b
      int edgeloc[N];
      for(i=0;i<nedge2;i++){
          if(i<numpartners[i1])
              edgeloc[Speclist[i1][i]]=i1;
          else
              edgeloc[Speclist[i2][i-numpartners[i1]]]=i2;
      }
      
      int i1orig=i1;
      int i2orig=i2;
      
      
    if(pos_sep==1)nedge2-=1;//Meaning, i1 was a partner of i2. This will decrease the total partners by 1 when they bind. (Note by DH)
    double sum=0;
      //pgen_b will require a little more effort below.
    if(nedge2==1){
      sum=2.0;//because we want 1 and it's multiplied by 0.5
    }else{
      for(i=1;i<nedge2;i++)
        sum+=binomial(nedge2, i);
        
      //if this new combined interface has a self interaction, add in more ways to split/
      if(pos_sep==1 || flag_self==1){
          sum+=2;//keep it the same configuration except move off the self
          for(i=1;i<nedge2-1;i++)
              sum+=binomial(nedge2-1,i);
      }
    }
      
      
    // //cout <<"total number of edges after combination: "<<nedge2<<"  pos_sep: "<<pos_sep<<" flag_self: "<<flag_self<<" sum: "<<sum<<endl;
    //double psplit_b=1.0/(sp1*0.5*sum);//*prob[split move]
    
    // //cout <<"in combine: forward: "<<pgen_f<<" backward: "<<pgen_b<<" ratio: "<<pgen_ratio<<" protein: "<<p1<<endl;
    
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
      //careful you are not creating a self interaction
      for(i=0;i<ni_one;i++){
          ind=Speclist[i1][i];
          Adj[i1*N+ind]=0;//i1 no longer bound to its original partners
          Adj[ind*N+i1]=0;
          if(ind==i2)
              numpartners[i2]--;//creating a self interaction from cross protein self
	
          Adj[i2*N+ind]=1;//add i1's original partners to i2
          Adj[ind*N+i2]=1;
          ////cout <<"move "<<ind<<" from i1: "<<i1<<" to i2: "<<i2<<endl;
      }
      
                
    }else{
      //careful that i1 binds to itself
      for(i=0;i<ni_one;i++){
          ind=Speclist[i1][i];
          if(ind==i1){
            Adj[i1*N+ind]=0;//i1 no longer bound to itself
              Adj[i2*N+i2]=1;//now i2 binds to itself!
              ////cout <<"move self from i1; "<<i1<<" to i2: "<<i2<<endl; 
          }else{
              Adj[i1*N+ind]=0;//i1 no longer bound to its original partners
              Adj[ind*N+i1]=0;
              Adj[i2*N+ind]=1;//add i1's original partners to i2
              Adj[ind*N+i2]=1;
              ////cout <<" move i1's former interface "<<ind<<" to i2: "<<i2<<endl;
          }
      }
    }
      
    //    cout << "Move 3: Checkpoint 4" << endl;
      
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
    if(Adj[(N-1)*N+N-1]==1)
        Adj[i1*N+i1]=1;
      
      
      
      
    int p2=p_home[N-1];
    
    p_home[i1]=p2;
    ////cout<<"i1: "<<i1<<" p_home: "<<p2<<" N-1: "<<N-1<<endl;
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
      //Now get reverse probability
      
      //New code by DH to find edge copies when splitting
      
      
      //If i2 happened to be N-1, then it is now the number i1.
      if(i2==N-1)
          i2=i1;
      
      t=0;
      edgeloc[i1]=edgeloc[N-1];//Since what was previously N-1 is now i1.
      int edgeused[N];
      int edge_copies=0;
      for(i=0;i<N;i++)
          edgeused[N]=0;
      int i1cop[numpartners[i2]],i2cop[numpartners[i2]];
      for(i=0;i<numpartners[i2];i++){
          i1cop[i]=0;
          i2cop[i]=0;
      }
      for(i=0;i<numpartners[i2];i++){
         tmp=Speclist[i2][i];
          if(tmp==i2 || edgeused[tmp]==1)continue;
          //cout<<"The edge connecting to " << tmp << " is being checked against ";
          if(edgeloc[tmp]==i1orig)
              i1cop[t]++;
          if(edgeloc[tmp]==i2orig)
              i2cop[t]++;
          for(j=i+1;j<numpartners[i2];j++){
              tmp2=Speclist[i2][j];
              //cout << tmp2 << " to see if they are copies." << endl;
              if(tmp2==i2)continue;
              find_copies(nwhole, N, numpartners, wholep, Speclist, Adj, tmp, p_home, tmp2, copy);
              if(copy==1 && edgeloc[tmp2]==i1orig){
                  edgeused[tmp2]=1;
                  i1cop[t]++;
              }
              if(copy==1 && edgeloc[tmp2]==i2orig){
                  edgeused[tmp2]=1;
                  i2cop[t]++;
              }
              if(copy==1){
                  edge_copies=1;
                  //cout<<"They are copies."<<endl;
              }
          }
          t++;
      }
      
      //cout<<"i1cop:";
      //for(i=0;i<t;i++)
          //cout<<"\t"<<i1cop[i];
      //cout<<"\n" << "i2cop:";
      //for(i=0;i<t;i++)
          //cout<<"\t"<<i2cop[i];
      //cout<<endl;
      
      double fact_b2=1.0;
      if(edge_copies==1){
      for(i=0;i<numpartners[i2];i++){
          if(i1cop[i]>0 && i2cop[i]>0){
              tmp=i1cop[i]+i2cop[i];
              fact_b2*=binomial(tmp,i1cop[i]);
              //cout<<"Increasing fact_b2 to " << fact_b2 << endl;
          }
      }
      
      if(copy_orig==1)
          fact_b2*=0.5;
      }
          
      
      double psplit_b=1.0/(sp1*0.5*sum);
      //cout <<"total number of edges after combination: "<<nedge2<<"  pos_sep: "<<pos_sep<<" flag_self: "<<flag_self<<" sum: "<<sum<<endl;
      
      int fact_b;
      
      fact_b=find_copies(nwhole, N, numpartners, wholep, Speclist, Adj, i2, p_home, i2, copy);
      //cout<<"In move 3: psplit_b is " << psplit_b << "and fact_b is " << fact_b <<" and fact_b2 is " << fact_b2 << endl;
      
      double pgen_b=fact_b*1.0*fact_b2*psplit_b;
      //double pgen_b=fact_b*1.0*psplit_b;
      pgen_ratio=pgen_b/pgen_f;
      pf=pgen_f;
      pb=pgen_b;
  } //done completing the move because chg==1
  }
    
  //Code added by DH to test for errors
  //In this case, Speclist is altered by Adj at the end. Adj should be symmetric
  for(i=0;i<N-1;i++){
    for(j=0;j<N-1;j++){
      if(Adj[i*N+j]!=Adj[j*N+i]){
	//cout << "Inside combine_interfaces: Adj not symmetric" << endl;
	//cout << "i equals: " << i << ", j equals: " << j << endl;
	//cout << "numpartners[i]: " << numpartners[i] << ", numpartners[j]: " << numpartners[j] << endl;
      }
    }
  }
  
  //  cout << "End of move 3 reached" << endl;
  return chg;
    
    
}
