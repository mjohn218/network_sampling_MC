#include "net_props.h"
#include <cstdlib>

void three_ways(int i, int j, int k, int l, int *Adj, int N, int &t, double **fourmer, double**hist, int *origlist, int *typelist)
{
  int flag;
  int nsub=4;
  int *ordered=new int[nsub];
  int ne;
  origlist[3]=l;
  int type;
  int sign;

    
  if(Adj[j*N+l]==1){
    //l is now a partner of j and it is the fourth component
    //sort these guys to make sure it is not a repeat
    order_nodes(nsub, origlist, ordered);
		
      //if(t>=N)
        //  throw invalid_argument("t is greater than N");
      
    fourmer[t][0]=ordered[0];
    fourmer[t][1]=ordered[1];
    fourmer[t][2]=ordered[2];
    fourmer[t][3]=ordered[3];
    t++;
    /*Should establish whether this fourmer is a repeat*/
    flag=check_repeat(t, fourmer);
    if(flag==1){
      //then this fourmer is a repeat, so skip it
      t--;
    }else{
      //evaluate the number of edges between them
      //so far we know that i connects to j, to k and to l
      ne=3;
      if(Adj[j*N+k]==1)
	ne++;
      if(Adj[i*N+l]==1)
	ne++;
      if(Adj[k*N+l]==1)
	ne++;
      //because i is the hub of this 4-mer, there is only 1 unique solution for each number of edges
      type=0;
      if(ne==3) type=1;
      if(ne==4){
	//type could be 0 or 1
          if(Adj[k*N+l]==1)
	  type=1;
              
      }
      ////cout <<"t: "<<t<<" ne: "<<ne<<" type: "<<type<<" hist: "<<hist[ne][type]<<" address: " <<endl;
      hist[ne][type]+=1;
              sign=1;
        
        //if(t-1>=N){
        //    //cout << "N equals: " << N << ". But t-1 equals: " << t-1 << endl;
        //    throw invalid_argument("Error with typelist[t-1]");
        //}
        
      //if(type==1)sign=-1;
      //typelist[t-1]=ne*sign;
      
    }
  }
  else if(Adj[k*N+l]==1){
    //l is now a partner of j and it is the fourth component
    //sort these guys to make sure it is not a repeat
    order_nodes(nsub, origlist, ordered);
    
    fourmer[t][0]=ordered[0];
    fourmer[t][1]=ordered[1];
    fourmer[t][2]=ordered[2];
    fourmer[t][3]=ordered[3];
    t++;
    /*Should establish whether this fourmer is a repeat*/
    flag=check_repeat(t, fourmer);
    if(flag==1){
      //then this fourmer is a repeat, so skip it
      t--;
    }else{
      //evaluate the number of edges between them
      //so far we know that i connects to j, to k and to l
      ne=3;
      if(Adj[j*N+k]==1)
	ne++;
      if(Adj[i*N+l]==1)
	ne++;
      if(Adj[j*N+l]==1)
	ne++;
      //because i is the hub of this 4-mer, there is only 1 unique solution for each number of edges
      type=0;
      if(ne==3) type=1;
      if(ne==4){
	//type could be 0 or 1
	if(Adj[j*N+l]==1)
	  type=1;
      }
      //       //cout <<"t: "<<t<<" ne: "<<ne<<" type: "<<type<<" hist: "<<hist[ne][type]<<endl;
      hist[ne][type]+=1;
      sign=1;
      //if(type==1)sign=-1;
      //typelist[t-1]=ne*sign;
      
    }
  }
  else if(Adj[i*N+l]==1){
    //l is a partner of i and it is the fourth component
    //sort these guys to make sure it is not a repeat
    order_nodes(nsub, origlist, ordered);
    
    fourmer[t][0]=ordered[0];
    fourmer[t][1]=ordered[1];
    fourmer[t][2]=ordered[2];
    fourmer[t][3]=ordered[3];
    t++;
    /*Should establish whether this fourmer is a repeat*/
    flag=check_repeat(t, fourmer);
    if(flag==1){
      //then this fourmer is a repeat, so skip it
      t--;
    }else{
      //evaluate the number of edges between them
      //so far we know that i connects to j, to k and to l
      ne=3;
      if(Adj[j*N+k]==1)
	ne++;
      if(Adj[j*N+l]==1)
	ne++;
      if(Adj[k*N+l]==1)
	ne++;
        
    if(ne>6)
        throw invalid_argument("ne is greater than 6");
      //because i is the hub of this 4-mer, there is only 1 unique solution for each number of edges
      type=0;
      ////cout <<"t: "<<t<<" ne: "<<ne<<" type: "<<type<<" hist: "<<hist[ne][type]<<endl;
      hist[ne][type]+=1;
      sign=1;
        
      //if(type==1)sign=-1;
      //typelist[t-1]=ne*sign;
      
    }
  }
    
  delete[] ordered;
}
