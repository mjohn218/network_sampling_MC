#include "pro_classes.h"
#include "add_iface_to_ppi.h"
#include "utility_calls.h"
#include "binomial.h"

int mutate_connections(int nwhole, int *numpartners, int **Speclist, Protein *wholep)
{

  /*In this case first choose a protein that has more than one interface*/
  /*ADDED  SELF*/
  double rnum1=trand()*1.0*nwhole/RAND_MAX;
  int p1=int(rnum1);//index of first protein
    
//   while(wholep[p1].ninterface==1){
//     rnum1=trand()*1.0*nwhole/RAND_MAX;
//     p1=int(rnum1);
//   }
  int chg=1;
  
  if(wholep[p1].ninterface==1)
    chg=0;
  if(chg!=0){
    //now pick 2 interfaces off this protein p1, and swap their partners
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
    
    //now choose a partner of i1
    rnum2=trand()*1.0*numpartners[i1]/RAND_MAX;
    w=int(rnum2);
    int ip1=Speclist[i1][w];
    rnum2=trand()*1.0*numpartners[i2]/RAND_MAX;
    w=int(rnum2);
    int ip2=Speclist[i2][w];
    int i;
    int t;
    ////cout <<"in mutate connections: "<<i1<<" i2: "<<i2<<" ip1: "<<ip1<<" ip2: "<<ip2<<endl;
    /*Now swap these*/
    if(i1==ip1){
      //self
      ////cout <<"self reconnect: from i1 to i1: "<<i1<< " to i2 to i2: "<<i2<<endl;
      for(i=0;i<numpartners[i1];i++){
          if(Speclist[i1][i]==ip1)
              Speclist[i1][i]=ip2;
      }
      
      for(i=0;i<numpartners[i2];i++){
          if(Speclist[i2][i]==ip2)
              Speclist[i2][i]=i2;//move self to this edge
      }
      for(i=0;i<numpartners[ip2];i++){
          if(Speclist[ip2][i]==i2)
              Speclist[ip2][i]=i1;
      }
    }else if(i2==ip2){
      //self
      ////cout <<"self reconnect: from i2 to i2: "<<i2<< " to i1 to i1: "<<i1<<endl;
      for(i=0;i<numpartners[i2];i++){
	if(Speclist[i2][i]==ip2)
	  Speclist[i2][i]=ip1;
      }
      
      for(i=0;i<numpartners[i1];i++){
	if(Speclist[i1][i]==ip1)
	  Speclist[i1][i]=i1;//move self to this edge
      }
      for(i=0;i<numpartners[ip1];i++){
	if(Speclist[ip1][i]==i1)
	  Speclist[ip1][i]=i2;
      }

    }else{
      
      /*unless it is the two interfaces that bind to one another*/
      if(i1==ip2 && i2==ip1){
          chg=0;
      }else if(i1==ip2){
	//for i1, can't just replace, need to delete both ip1 and i2 
	t=0;
	//first delete ip1
	for(i=0;i<numpartners[i1];i++){
	  if(Speclist[i1][i]!=ip1){
	    Speclist[i1][t]=Speclist[i1][i];
	    t++;
	  }
	}
	//now replace i2 with ip2
	numpartners[i1]--;
	for(i=0;i<numpartners[i1];i++){
	  if(Speclist[i1][i]==i2)
	    Speclist[i1][i]=ip2;//now binds to self
	}
	for(i=0;i<numpartners[ip1];i++){
	  if(Speclist[ip1][i]==i1)
	    Speclist[ip1][i]=i2;
	}
	for(i=0;i<numpartners[i2];i++){
	  if(Speclist[i2][i]==ip2)
	    Speclist[i2][i]=ip1;
	}
	
      }else if(i2==ip1){
	//for i1, can't just replace, need to delete both ip1 and i2 
	t=0;
	//first delete ip2
	for(i=0;i<numpartners[i2];i++){
	  if(Speclist[i2][i]!=ip2){
	    Speclist[i2][t]=Speclist[i2][i];
	    t++;
	  }
	}
	//now replace i1 with ip1
	numpartners[i2]--;
	for(i=0;i<numpartners[i2];i++){
	  if(Speclist[i2][i]==i1)
	    Speclist[i2][i]=ip1;//now binds to self
	}
	for(i=0;i<numpartners[ip2];i++){
	  if(Speclist[ip2][i]==i2)
	    Speclist[ip2][i]=i1;
	}
	for(i=0;i<numpartners[i1];i++){
	  if(Speclist[i1][i]==ip1)
	    Speclist[i1][i]=ip2;
	}
	
      }else{
	//normal 
	for(i=0;i<numpartners[i1];i++){
	  if(Speclist[i1][i]==ip1)
	    Speclist[i1][i]=ip2;
	}
	for(i=0;i<numpartners[ip1];i++){
	  if(Speclist[ip1][i]==i1)
	    Speclist[ip1][i]=i2;
	}
	for(i=0;i<numpartners[i2];i++){
	  if(Speclist[i2][i]==ip2)
	    Speclist[i2][i]=ip1;
	}
	for(i=0;i<numpartners[ip2];i++){
	  if(Speclist[ip2][i]==i2)
	    Speclist[ip2][i]=i1;
	}
      }

    }
    
    ////cout <<" swap! "<<i1 <<' '<<i2<<' '<<ip1<<' '<<ip2<<endl;
  }
  return chg;
}
