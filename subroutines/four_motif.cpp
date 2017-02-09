#include "net_props.h"
#include <cstdlib>

void four_motif(int N, int *Adj, double **fourmer, double *histfour, double**hist, int *typelist)
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
    
  /*double **hist2 = new double*[7];
    for (int jj=0;jj<7;jj++)
        hist2[jj]=new double[2];
    
    hist2[3][0]=0;
    hist2[3][1]=0;
    hist2[4][0]=0;
    hist2[4][1]=0;
    hist2[5][0]=0;
    hist2[6][0]=0;*/
    
  ////cout <<"in four subgraph motifs! "<<endl; 
  /*Try using the adjacency matrix*/
  int nsub=4;
  int *origlist=new int[nsub];
  
  int t=0;//this index keeps track of which fourmer we're on
  int flag, ne;
  int type;
//   //cout <<"Current adjacency in four motif "<<endl;
//   for(i=0;i<N;i++){
//     for(j=0;j<N;j++){
//       ////cout <<Adj[i*N+j]<<'\t';
//     }
//     ////cout <<endl;
//  }
    
  //Make adjacency matrix out of Speclist (DH)

  //Print out Adj matrix (test code)
  //for(i=0;i<N;i++){
    // //cout<<i << "\t" << numpartners[i] << "\t";
    //for(j=0;j<N;j++){
  ////cout<<Adj[i*N+j] << "\t";
  //}
  ////cout<<endl;
  //}


  //Check for symmetry
    for(i=0;i<N;i++){
        for(j=0;j<N;j++){
            if(Adj[i*N+j]!=Adj[j*N+i]){
                //cout << "Error: mismatch in Adj (four motif)" << endl;
		//cout << "I1 is: " << i << ". I2 is: " << j << endl;
                exit(1);
		//i=N;
                //break;
            }
        }
    }

   
   //Maggie's code
    /*
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
	      //find a fourth component by looping over partners of i, and j, and k. each of those subnetworks could be distinct
	      for(l=0;l<N;l++){
		if(l!=k &&l!=i &&l!=j){
		  three_ways(i, j, k, l, Adj, N, t, fourmer, hist2, origlist, typelist);
		}
	      }//end looping over l
	    }
	    if(Adj[j*N+k]==1){
	      //now k is the third component and is bound to J
	      //find a fourth component by looping over partners of i, and j, and k. each of those subnetworks could be distinct
	      for(l=0;l<N;l++){
		if(l!=k &&l!=i &&l!=j){
		  three_ways(j, i, k, l, Adj, N, t, fourmer, hist2, origlist, typelist);//i and j swap positions
		}
	      }//end looping over l
	    }
	    

	  }
        
	}//end looping over k

      }
        
    }//end looping over j
      
  }//end looping over i*/
    
    //End Maggie's code.
    
    
   
    
    int a,b,c,d,e,f;
    
    //My code for finding motifs w/o having to bother with former but avoiding repeats. I have tested this code against Maggie's 10 times and each time the numbers found for all 6 motifs were equal. It's possible there's a very rare three edge star it misses
   
    int tak;
    
    for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){
            //See if there is a connection
            if (Adj[i*N+j]==1){
                for(k=i+1;k<N;k++){
                    //Find a connection. See that this is not a connection already looked at.
                    if ((Adj[j*N+k]==1) && (k>j)){
                        for(l=i+1;l<N;l++){
                            if(Adj[k*N+l]==1){
                                ////cout << "Good so far1" << endl;
                                //tak++;
                                    //Check if a repeated connection. Else, find type.
                                
                                if(l>k){
                                    find_type(i,j,k,l,Adj,N,hist,1);//1
                                }
                                
                                else if((Adj[j*N+l]==0) && (l>j)){
                                    find_type(i,j,k,l,Adj,N,hist,2);//2
                                }
                                else if((l<j) && (Adj[i*N+l]==0) && (Adj[j*N+l]==0)){
                                    find_type(i,j,k,l,Adj,N,hist,3);//3
                                }
                            }
                            else{
                                if((Adj[i*N+l]==1) && (l>j) && (k!=l)){
                                    if(Adj[j*N+l]==0)
                                        find_type(i,j,k,l,Adj,N,hist,11);
                                    else if (l>k)
                                        find_type(i,j,k,l,Adj,N,hist,15);
                                }
                                else if((Adj[j*N+l]==1) && (l>k))
                                    find_type(i,j,k,l,Adj,N,hist,13);
                            }
                            
                        }
                    }
            
                    else if((Adj[j*N+k]==1) && (k<j) && (Adj[i*N+k]==0)){
                        for(l=i+1;l<N;l++){
                            if(Adj[k*N+l]==1){
                                ////cout << "Good so far2" << endl;
                                //tak++;
                                if(l>j){
                                    find_type(i,j,k,l,Adj,N,hist,4);//4
                                }
                                else if((l<k) && (Adj[j*N+l]==0) && (Adj[i*N+l]==0)){
                                    find_type(i,j,k,l,Adj,N,hist,5);//5
                                }
                                else if((l>k) && (l<j) && (Adj[i*N+l]==0)){
                                    find_type(i,j,k,l,Adj,N,hist,6);//6
                                }
                            }
                            else{
                                if((Adj[i*N+l]==1) && (l>k) && (l>j))
                                    find_type(i,j,k,l,Adj,N,hist,12);
                                else if((Adj[j*N+l]==1) && (l>k)){
                                    if( (l>j) && (Adj[i*N+l]==1) )
                                       find_type(i,j,k,l,Adj,N,hist,14);
                                    else if (Adj[i*N+l]==0)
                                        hist[3][0]++; //16
                                        
                                       
                                }
                            }
                        }
                    }
                    else if((Adj[j*N+k]==0) && (Adj[i*N+k]==1) && (k>j)){
                        for(l=i+1;l<N;l++){
                            if (Adj[j*N+l]==0){
                                if ((l>k) && (Adj[k*N+l]==1))
                                    find_type(i,j,k,l,Adj,N,hist,7);//7
                                else if ((l>j) && (l<k) && (Adj[i*N+l]==0) && (Adj[k*N+l]==1))
                                    find_type(i,j,k,l,Adj,N,hist,8);//8
                                else if((l<j) && (Adj[i*N+l]==0) && (Adj[k*N+l]==1))
                                    find_type(i,j,k,l,Adj,N,hist,9);//9 Should be a three edge star
                                else if((Adj[k*N+l]==0) && (Adj[i*N+l]==1) && (l>k)){
                                    hist[3][0]++;//10
                                }
                            }
                        }
                    }
                }
            }
        }
    }

                          
                
        

    //End new code by DH
    
  /*Now count total number of subgraphs*/
  double total=hist[3][0]+hist[3][1]+hist[4][0]+hist[4][1]+hist[5][0]+hist[6][0];
    //double totalm = hist2[3][0]+hist2[3][1]+hist2[4][0]+hist2[4][1]+hist2[5][0]+hist2[6][0];
    ////cout <<"total 4-mers: "<<total<< " " << totalm << endl;
  
  
  histfour[0]=hist[3][0]*1.0/(1.0*total+1E-8); //hub
  histfour[1]=hist[3][1]*1.0/(1.0*total+1E-8); //open square (chain)
  histfour[2]=hist[4][0]*1.0/(1.0*total+1E-8); //flag
  histfour[3]=hist[4][1]*1.0/(1.0*total+1E-8); //square
  histfour[4]=hist[5][0]*1.0/(1.0*total+1E-8);
  histfour[5]=hist[6][0]*1.0/(1.0*total+1E-8);
  /*
    //cout <<"Three edges, star shape: "<< hist[3][0]<<endl;//histfour[0]<<endl;
     //cout <<"Three edges, square: "<<hist[3][1] <<endl;//histfour[1]<<endl;
    //cout <<"Four edges, star: "<<hist[4][0] <<endl;//histfour[2]<<endl;
     //cout <<"Four edges, square: "<<hist[4][1] <<endl;//histfour[3]<<endl;
     //cout <<"Five edges: "<<hist[5][0] << endl;
     //cout <<"Six edges: "<<hist[6][0]<<endl;
     //cout <<"List fourmers: "<<total<<endl;
  */
    
   /* //cout <<"Do my code and maggie's code produce the same results?" << endl;
    if (hist[3][0]==hist2[3][0])
        //cout << "Yes ";
    else
        //cout << "No ";
    if (hist[3][1]==hist2[3][1])
        //cout << "Yes ";
    else
        //cout << "No ";
    if (hist[4][0]==hist2[4][0])
        //cout << "Yes ";
    else
        //cout << "No ";
    if (hist[4][1]==hist2[4][1])
        //cout << "Yes ";
    else
        //cout << "No ";
    if (hist[5][0]==hist2[5][0])
        //cout << "Yes ";
    else
        //cout << "No ";
    if (hist[6][0]==hist2[6][0])
        //cout << "Yes " << endl;
    else
        //cout << "No " << endl;*/


  //  for(i=0;i<total;i++){
    //  //cout <<fourmer[i][0]<<' '<<fourmer[i][1]<<' '<<fourmer[i][2]<<' '<<fourmer[i][3]<<' '<<typelist[i]<<endl;
  //}
  delete[] origlist;

}
