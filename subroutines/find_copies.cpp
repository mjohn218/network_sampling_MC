#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int find_copies(int nwhole, int N, int *numpartners, Protein *wholep, int **Speclist, int *Adj, int i1, int *p_home, int bait, int &copy)
{
    /*Trying out a new, less messy algorithm
     Come up with an "edge profile" for each interface. Want to record:
     1) If there is a self-interaction
     2) If it is a split-self
     3) How many connections it has to every other protein in the system
     4) What protein is the interface on
     
     For interfaces to be identical, they and their partners must have identical profiles
     */
    ////cout << "Inside find_copies. i1 is " << i1 << endl;

    int i,j,k,l,m;
    int i2, w1, w2;
    copy=0;
    
    int p1=p_home[i1];
    //    int **profile=new int*[N];
    int profile[N][nwhole+2];
    for(i=0;i<N;i++){
      ////cout<<"First for loop. i="<<i<<endl;
      //profile[i]=new int[nwhole+2];
        for(j=0;j<nwhole+2;j++){
	  ////cout<<"Second for loop.j="<<j<<"\t";
            profile[i][j]=0;
	    ////cout<<"profile[i][j]="<<profile[i][j]<<endl;
	}
    }
    ////cout << "Checkpoint0. Okay setting up profiles\n";
    for(j=0;j<N;j++){
      ////cout<<"In first for loop." << endl;
        profile[j][nwhole+1]=p_home[j];//Record home protein
        for(i=0;i<numpartners[j];i++){
	  ////cout << "In second for loop. i="<<i;
            w1=Speclist[j][i];
            if(w1==j)
                profile[j][nwhole]=1; //Records a self-interaction as true
            else
                profile[j][p_home[w1]]++;//Records number of times it connects to a protein
        }
    }
    
    ////cout << "Checkpoint1. Okay so far.\n";
    bool match;
    int isMatch[N*N]; //Find interfaces with the same profiles. N*N matrix
    for(i=0;i<N*N;i++)
        isMatch[i]=0;
    for(i=0;i<N;i++){
        for(j=i;j<N;j++){
            if(p_home[i]!=p_home[j] || numpartners[i]!=numpartners[j])
                continue;
            match=true;
            for(k=0;k<nwhole+1;k++){
                if(profile[i][k]!=profile[j][k]){
                    match=false;
                    break;
                }
            }
            if(match==true){
                isMatch[i*N+j]=1;
                isMatch[j*N+i]=1;
            }
        }
    }
    ////cout << "Checkpoint2. Okay so far.\n";
    
    ////cout<< "N is " << N << ", i1 is " << i1 << ", bait is " << bait << endl;

    ////cout<<"Do i1 and bait have a matching profile? 1 if yes. " << isMatch[i1*N+bait] << endl;;

    //Now that we have created profiles for every interface in the system, see if there are copies on a protein
    int iface_checked[N];
    for(i=0;i<N;i++)
        iface_checked[i]=0;
    
    //int iface_checked2[N];
    //for(i=0;i<N;i++)
    //iface_checked2=0;
    bool checking, isdup;
    int nmatches[N];
    int edgeused[N];
    int toCheck[N][2];//Records pairs of ifaces to see if they have identical profiles
    int trueMatch[N];//Records if partners have same profiles
    for(i=0;i<N;i++)
        trueMatch[N]=0;
    int t,itr,iter1,iter2;
    
    int ndups=1; //Count i1 itself as a copy
    ////cout<<"Checkpoint3. Beginning for loop.\n";
    ////cout<<"In find_copies: ninterfaces is " << wholep[p1].ninterface << endl;
    for(i=0;i<wholep[p1].ninterface;i++){
        i2=wholep[p1].valiface[i];
        ////cout<<"iface being checked is " << i2 << endl;
        isdup=false;
        if(i2==i1)
            continue;
        if(isMatch[i2*N+i1]==1){//Need to check partners on other proteins recursively
            //repeat=true;
            for(j=0;j<N;j++){
                toCheck[j][0]=-1;
                toCheck[j][1]=-1;
            }
            t=0;
            itr=0;
            for(j=0;j<N;j++)
                iface_checked[j]=0;
            toCheck[0][0]=i1;
	    toCheck[0][1]=i2;
            //while(repeat=true){
                //repeat=false;
                //iter1=i1;
                //iter2=i2;
            iface_checked[i1]=1;//Meaning this is already in the toCheck array
            iface_checked[i2]=1;
                checking=true;
                while(checking==true){
		  ////cout << "In while loop" << endl;
                    checking=false;
                    nmatches[itr]=0;
                    iter1=toCheck[itr][0];
                    iter2=toCheck[itr][1];
                    if(iter1==-1){//No more sets of partners to check.
                        isdup=true;//Since we have not broken out of the loop already, all sets of partners must have the same profiles
                        break;          
                    }
                    for(j=0;j<N;j++)
                        edgeused[j]=0;
                    //iface_checked[iter1]=1;
                    //iface_checked[iter2]=1;
                    
                    if(Adj[iter1*N+iter1]==1 && Adj[iter2*N+iter2]==1){//Check to see if both have self-edges
                        nmatches[itr]++;
                        //nmatches[iter2*N+iter1]++;
                    }
                    for(j=0;j<numpartners[iter1];j++){//Choose a partner of iter1
                        w1=Speclist[iter1][j];
                        if(w1==iter1)
                            continue;
                        for(k=0;k<numpartners[iter2];k++){//See if there is a partner of iter2 with the same profile
                            w2=Speclist[iter2][k];
                            if(w2==iter2 || edgeused[w2]==1)
                                continue;
			    if(w1==w2){//Case: shared partner
			      nmatches[itr]++;
			      edgeused[w2]=1;//Now that we have used this partner of iter2, shouldn't use it more than once
			    }
			    //Shared partner case wasn't being caught below
                            else if(isMatch[w1*N+w2]==1){
                                //nmatches[iter1*N+iter2]++;
                                //nmatched[iter2*N+iter1]++;
                                nmatches[itr]++;
                                edgeused[w2]=1;//Now that we have used this partner of iter2, shouldn't use it more than once
                                if(iface_checked[w2]==0){//Add pair to array to be checked later
                                    t++;
                                    toCheck[t][0]=w1;
                                    toCheck[t][1]=w2;
                                    iface_checked[w1]=1;//Now placed in the toCheck array
                                    iface_checked[w2]=1;
				    ////cout<<"to check: pair of interfaces " << w1 << " and " << w2 << endl;
                                }
                                //iter1=w1;//Now check to see
                                //iter2=w2;
                                //nmatches[iter1*N+iter2]++;
                                //nmatched[iter2*N+iter1]++;//Already know that backwards is a match
                                //checking=true;//Keep iterating
                                //j=numpartners[iter1];//break out of for loops
                                //break;
                            }
                        }
                    }
                    if(nmatches[itr]==numpartners[iter1] && numpartners[iter1]==numpartners[iter2]){
                        trueMatch[itr]=1;//Partners of iter1 and iter2 have matching profiles
                        itr++;
			////cout << iter1 << " and " << iter2 << " are a match.\n";
                        checking=true;//If there is no match, then the while loop will break since checking=false
                    }
                
                }
            //}
            
            
            
        }
        if(isdup==true){
            ndups++;
            if(i2==bait)
                copy=1;
        }
    }
    
    if(bait==i1)
        copy=1;

    if(ndups>1){
      //cout<<"Copies of " << i1 << " found. ndups=" << ndups << endl;
    }
      //for(i=0;i<N;i++)
      //delete profile[i];

      //delete [] profile;
    //delete [] isMatch;

    return ndups;
    
}
