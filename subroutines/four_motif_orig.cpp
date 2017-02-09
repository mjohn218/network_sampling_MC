#include "net_props.h"

void four_motif_orig(int N, int *Adj, double **fourmer, double *histfour, double**hist, int *typelist)
{
  /*There are 6 subgraphs that have 4 nodes*/
  /*two have 3 edges*/
  /*two have 4 edges*/
  /*one has 5 edges*/
  /*one has 6 edges*/
  /*All self loops are ignored.*/
  int i,j, k, l;
  int ind;
  double factor;
  int Ntype=6;
  int Nrepeats=2;//max number of graphs with given edge number

  //never less than 3 edges
  hist[3][0]=0;
  hist[3][1]=0;
  hist[4][0]=0;
  hist[4][1]=0;
  hist[5][0]=0;
  hist[6][0]=0;
  cout <<"in four subgraph motifs! "<<endl; 
  /*Try using the adjacency matrix*/
  int nsub=4;
  int *origlist=new int[nsub];
  
  int t=0;//this index keeps track of which fourmer we're on
  int flag, ne; 
  int type;
//   cout <<"Current adjacency in four motif "<<endl;
//   for(i=0;i<N;i++){
//     for(j=0;j<N;j++){
//       //cout <<Adj[i*N+j]<<'\t';
//     }
//     //cout <<endl;
//  }
  for(i=0;i<N;i++){
    origlist[0]=i;
    for(j=i+1;j<N;j++){
      origlist[1]=j;
      if(Adj[i*N+j]==1){
	//these two nodes are connected, find a third partner
	//first loop over i's partners.
	for(k=0;k<N;k++){
	  //make sure k is not the same as j
	  origlist[2]=k;
	  if(k!=j && k!=i){
	    if(Adj[i*N+k]==1){
	      //now k is the third component and is bound to i
	      /*find a fourth component by looping over partners of i, and j, and k. 
		each of those subnetworks could be distinct*/
	      for(l=0;l<N;l++){
		if(l!=k &&l!=i &&l!=j){
		  three_ways(i, j, k, l, Adj, N, t, fourmer, hist, origlist, typelist);
		}
	      }//end looping over l
	    }
	    if(Adj[j*N+k]==1){
	      //now k is the third component and is bound to J
	      /*find a fourth component by looping over partners of i, and j, and k. 
		each of those subnetworks could be distinct*/
	      for(l=0;l<N;l++){
		if(l!=k &&l!=i &&l!=j){
		  three_ways(j, i, k, l, Adj, N, t, fourmer, hist, origlist, typelist);//i and j swap positions
		}
	      }//end looping over l
	    }
	    

	  }
	}//end looping over k
      }
    }//end looping over j
  }//end looping over i
  /*Now count total number of subgraphs*/
  double total=hist[3][0]+hist[3][1]+hist[4][0]+hist[4][1]+hist[5][0]+hist[6][0];
  cout <<"total 4-mers: "<<total<<' '<<t<<endl;
  
  cout <<"Three edges, star shape: "<< hist[3][0]*1.0/(1.0*total)<<endl;
  cout <<"Three edges, square: "<<hist[3][1]*1.0/(1.0*total)<<endl;
  cout <<"Four edges, star: "<<hist[4][0]*1.0/(1.0*total)<<endl;
  cout <<"Four edges, square: "<<hist[4][1]*1.0/(1.0*total)<<endl;
  cout <<"Five edges: "<<hist[5][0]*1.0/(1.0*total)<<endl;
  cout <<"Six edges: "<<hist[6][0]*1.0/(1.0*total)<<endl;
  cout <<"List fourmers: "<<total<<endl;
  histfour[0]=hist[3][0]*1.0/(1.0*total);
  histfour[1]= hist[3][1]*1.0/(1.0*total);
  histfour[2]=hist[4][0]*1.0/(1.0*total);
  histfour[3]=hist[4][1]*1.0/(1.0*total);
  histfour[4]=hist[5][0]*1.0/(1.0*total);
  histfour[5]=hist[6][0]*1.0/(1.0*total);

  //  for(i=0;i<total;i++){
    //  cout <<fourmer[i][0]<<' '<<fourmer[i][1]<<' '<<fourmer[i][2]<<' '<<fourmer[i][3]<<' '<<typelist[i]<<endl;
  //}
  delete[] origlist;

}
