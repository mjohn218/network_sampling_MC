#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int mutate_edge_orig(int nwhole, int *numpartners, int **Speclist, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *selfflag, int *p_home)
{

  /*In this case first choose a protein that has more than one interface*/
  /*ADDED SELF*/
  pgen_ratio=1;
  double rnum1=trand()*1.0*nwhole/RAND_MAX;
  int p1=int(rnum1);//index of first protein
    
  int chg=1;
  int i, t;
  int loc;
  int i1, ind, ip1;
  double pf, pb;
  int nedge;
  //  cout <<"in mutate edge, ninterfaces : "<<wholep[p1].ninterface<<endl;
  chg=0;
  if(wholep[p1].ninterface<ppi[p1].nppartner)
    chg=1;
  /*If they are equal, you can still move an edge if there is a self interaction on the protein*/
  if(wholep[p1].ninterface==ppi[p1].nppartner &&selfflag[p1]==1){
    chg=1;//binds to self, so its possible to have an extra interface
  }
  if(wholep[p1].ninterface==1)
    chg=0;
  int oself=0;
  int w, spf;
  double denom;
  int tmp;
  if(chg!=0){
    /*Pick one edge off this protein*/
    if(selfflag[p1]==1){
      /*find out if it is to itself, or split self*/
      for(i=0;i<wholep[p1].ninterface;i++){
	ind=wholep[p1].valiface[i];
	for(t=0;t<numpartners[ind];t++){
	  if(Speclist[ind][t]==ind)
	    oself=1;
	}
      }
      if(oself==0){
	oself=-1;//the self interaction is split
	nedge=ppi[p1].nppartner+1;
      }else{
	nedge=ppi[p1].nppartner;
      }
      denom=1.0*nedge;
    }else{
      //no self
      nedge=ppi[p1].nppartner;
      denom=1.0*nedge;
    }
    rnum1=trand()*1.0*nedge/RAND_MAX;
    w=int(rnum1);
    i=0;
    loc=w;
    while(i<wholep[p1].ninterface){
      ind=wholep[p1].valiface[i];
      if(loc<numpartners[ind]){
	ip1=Speclist[ind][loc];
	i1=ind;
	i=wholep[p1].ninterface;//break out
      }else{
	loc-=numpartners[ind];
	i++;
      }
    }
    
    /*Now pick another interface on this protein to move the edge to*/
    double rnum2=trand()*1.0*wholep[p1].ninterface/RAND_MAX;

    int i2=wholep[p1].valiface[int(rnum2)];
    //cout <<"in mutate edge: p1 "<<p1 <<" ninterface "<<wholep[p1].ninterface<<" rnum2: "<<rnum2<<" i2: "<<i2<<" i1: "<<i1<<endl;
    while(i1==i2){
      rnum2=trand()*1.0*wholep[p1].ninterface/RAND_MAX;
      i2=wholep[p1].valiface[int(rnum2)];
    }
    
    pb=1.0;
    pf=1.0;
    pgen_ratio=pb/pf;
    if(ip1==i1){
      /*Self binding splitting to another edge 
	Move i1's partner ip1 to i2, keep the rest of his partners*/
      //      cout <<"in move edge, moving i1 self interaction: "<<i1<<" to interaction between i1 and i2: "<<i2<<endl;
      pf=1.0/(1.0*ppi[p1].nppartner);
      pb=1.0/(1.0*ppi[p1].nppartner+1.0);
      pgen_ratio=pb/pf;
      t=0;
      /*Now swap these*/
      for(i=0;i<numpartners[i1];i++){
	if(Speclist[i1][i]==ip1){
	  Speclist[i1][i]=i2;
	  //t++;
	}
      }
      //i1 and ip1 are the same, so updating one also updates the other
      //i1 still has a partner! it is i2!

      Speclist[i2][numpartners[i2]]=ip1;
      numpartners[i2]++;//added one
    }else if(numpartners[i1]==1){
      
      chg=0;
      
    }else{
      /*Move i1's partner ip1 to i2, keep the rest of his partners*/
      //      cout <<"in move edge, move i1 interaction: "<<i1<<" with his partner ip1: "<<ip1<<" to i2: "<<i2<<endl;
      /*First check if this edge is a split self edge*/
      spf=0;
      if(p_home[ip1]==p1){
	pf=1.0/(1.0*ppi[p1].nppartner+1.0);
	pb=1.0/(1.0*ppi[p1].nppartner+1.0);
	spf=1;
      }
      t=0;
      /*Now swap these*/
      for(i=0;i<numpartners[i1];i++){
	if(Speclist[i1][i]!=ip1){
	  Speclist[i1][t]=Speclist[i1][i];
	  t++;
	}
      }
      numpartners[i1]--;//got rid of one partner
      if(ip1==i2){
	pb=1.0/(1.0*ppi[p1].nppartner);
	//forming a self interaction from a within protein self interaction
	for(i=0;i<numpartners[i2];i++){
	  if(Speclist[i2][i]==i1){
	    Speclist[i2][i]=i2;
	  }
	}
	
      }else{
	/*Check about a stand alone split edge, results in 2 ways to get to same state
	  if you are moving to one or creating one, extra factor of 2*/
	if(spf==0 && oself==-1){
	  tmp=Speclist[i2][0];
	  if(numpartners[i2]==1 &&p_home[tmp]==p1){
	    if(numpartners[tmp]==1){
	      pf*=2.0;
	    }
	  }
	  tmp=Speclist[i1][0];
	  if(numpartners[i1]==1 &&p_home[tmp]==p1){
	    if(numpartners[tmp]==1){
	      pb*=2.0;
	    }
	  }
	}
	for(i=0;i<numpartners[ip1];i++){
	  if(Speclist[ip1][i]==i1)
	    Speclist[ip1][i]=i2;
	}
	Speclist[i2][numpartners[i2]]=ip1;
	numpartners[i2]++;//added one
      }
      pgen_ratio=pb/pf;
    }
    //cout <<" move! "<<i1 <<' '<<i2<<' '<<ip1<<endl;
  }//end if move possible
  
  return chg;

}
