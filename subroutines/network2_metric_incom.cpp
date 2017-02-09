#include "pro_classes.h"
#include "utility_calls.h"
#include "network_metrics2.h"
#include "matmultiply.h"

double network2_metric_incom(int nwhole,  ppidata *ppi, int **Eref, int **Epred)
{
  /*Compare the current network's structure to that of the target network
   */
  int i, j;
  int nedge;

  int t=0;
  int flag;
  
  int nstep=0;
  int e1=0;
  int s, m;
  int df; 
  int nwrong[nwhole];
  long unsigned int index[nwhole+1];
  double value[nwhole+1];
  int brk;

  int ind;
  int tms;
  /*Have an 'group' matrix, if you have a 1, it means you are in that pair are in the same group  */
  //  cout <<"Network metric . "<<endl;
  for(i=0;i<nwhole;i++){
    // cout << "Testing. We are at protein number: " << i << endl;
    /*Loop over all proteins and calculate the difference between weights for each one */
    /*Loop over the edges, since they are conserved between the two networks
     */
    nedge=ppi[i].nppartner;
    flag=1;
    if(nedge==1)flag=0;
        cout <<"Protein: "<<i<<"  nedge: "<<nedge<<endl; 
    //tms=0;
    while(flag==1){
      cout <<"flag: "<<tms<<"  protein: "<<i<<endl;
      //tms++;
      for(j=0;j<nedge;j++)
	nwrong[j]=0;
      for(j=0;j<nedge;j++){
	/*In this version, we will evaluate the number of interface combines
	  or splits needed to move from the predicted network to the target
	  network
	*/
	/*loop over the elements in his row*/
	for(s=j;s<nedge;s++){
	  
	  df=Epred[i][j*EDIM+s]-Eref[i][j*EDIM+s];
	  if(abs(df)==1){
	    nwrong[j]++;
	    nwrong[s]++;
	  }
	  
	}
	cout <<"edge: "<<j<<" wrong: "<<nwrong[j]<<endl;
      }//end looping over the edges in this protein	      
      
      /*Now start fixing the network to count the number of moves needed*/
      
      for(j=0;j<nedge;j++)
	value[j+1]=nwrong[j];
      indexx(nedge, &value[0], &index[0]);//indexing is like fortran, from 1:N!!
      ind=index[1]-1;
      j=1;
      while(nwrong[ind]==0 && j<nedge){
	j++;
	ind=index[j]-1;
      }
      cout <<"index of incorrect edge: "<<ind<<" nwrong: "<<nwrong[ind]<<" nedges: "<<nedge<<endl;
      if(nwrong[ind]==0)flag=0;//all is correct
      //ind contains the index of which edge to correct
      
      if(flag==1){
	cout <<"flag 1: "<<"first fix: "<<ind<<" nwrong: "<<nwrong[ind]<<endl;
	/*Perform corrections on this edge by looping over all other edges*/
	for(s=0;s<nedge;s++){
	  //  if(s!=ind){
	    //test whether this connection is correct
	    df=Epred[i][ind*EDIM+s]-Eref[i][ind*EDIM+s];
	    cout <<"diff: "<<df<<" ind: "<<ind<<" s: "<<s<<endl;
	    if(df==1){
	      //then disconectn this edge
	      Epred[i][ind*EDIM+s]=0;
	      Epred[i][s*EDIM+ind]=0;
	      for(m=0;m<nedge;m++){
		if(Epred[i][ind*EDIM+m]==1){
		  //these are grouped, disconnect s from them
		  Epred[i][s*EDIM+m]=0;
		  Epred[i][m*EDIM+s]=0;
		  
		}
	      }
	      nstep++;
	      
	    }else if(df==-1){
	      //need to add an edge, so disconnect that edge from wherever it is
	      //s is a missing edge
	      //if s has 0 partners, then just add it, otherwise break it, then add it
	      brk=0;
	      for(m=0;m<nedge;m++){
		if(Epred[i][s*EDIM+m]==1){
		  Epred[i][s*EDIM+m]=0;
		  Epred[i][m*EDIM+s]=0;
		  brk=1;
		}
	      }
	      Epred[i][ind*EDIM+s]=1;
	      Epred[i][s*EDIM+ind]=1;
	      /*You also need to add this edge to all of ind's partners*/
	      for(m=0;m<nedge;m++){
		if(m!=ind){
		  if(Epred[i][ind*EDIM+m]==1){
		    //ind is connected to this edge, so add it also to s, unless it is s.
		    if(m!=s){
		      Epred[i][s*EDIM+m]=1;
		      Epred[i][m*EDIM+s]=1;
		    }
		  }
		}
	      }
	      nstep++;
	      if(brk==1)nstep++;//we also broke s from someone else
	      
	    }
	    
	    //}//don't check against self
	}//end looping over shared edges
      }
    
    }//done correcting this protein
    e1+=nedge;
  }
  
  double fitness=nstep;
  return fitness;
  
}
