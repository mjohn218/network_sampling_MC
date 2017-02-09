#include "pro_classes.h"
#include "constrainParms.h"
#include "create_network.h"
#include "calc_ratios.h"
#include "utility_calls.h"

int create_network(int nwhole, ppidata *ppi, Protein *wholep, int *p_home, constrainParms &plist, int maxinterface, int Nedge, int **Speclist, int *numpartners )
{

  /*
    given a PPI, define an initial IIN.
    One way would be to start out with same number of interfaces as proteins and
    start splitting
   A reasonable metric, in keeping with the highest specificity IINs, would be to
   count the number of partners your neighbors have, divided by the number of neighbors.
   This number should be close to 1.
  */
  int i, j;
  int nbor[nwhole];
  int ind;
  double *ratio=(double *)malloc(sizeof(double) *nwhole);
  long unsigned int *index=(long unsigned int *) malloc(sizeof(long unsigned int) *(nwhole+1));
  ppidata * iin=(ppidata *) malloc(sizeof(ppidata) *nwhole);
  ppidata * tempiin=(ppidata *) malloc(sizeof(ppidata) *nwhole);
  int id;
  int icurr;
  double highratio;
  int chgprot;
  int num2, t, num, n2;
  int ncurr=nwhole;
  double *Adj=(double *)malloc (sizeof(double) *ncurr*ncurr);
  char fname[100];
  ofstream Mfile;
  double avg_ratio[MAXP];
  int ntry;
  
  for(i=0;i<nwhole;i++){
    iin[i]=ppi[i];
    wholep[i].ninterface=1;//currently all proteins have 1 interfaces
    wholep[i].valiface[0]=i;//your interface is the same as your protein index
  }
  

  /*Now the proteins are ordered*/
  /*At this stage can add an interface on this protein, randomly, or loop through all possible
    additions, and take the best one*/
  
  int nit=0;

  avg_ratio[0]=2;
  
  while(ncurr<maxinterface && avg_ratio[0] >1){
    /*Continue looping to add more interfaces to 
      the ppi network*/
    cout <<"Niteration: "<<nit<<endl;
    calc_ratio(ncurr, iin, ratio);
    order_ratio(ncurr, ratio, index);
    
    for(i=0;i<ncurr;i++){
      id=index[i+1]-1;
      cout <<"index: "<<id<<" ratio: " <<ratio[id]<<endl;;
    }
    
    ncurr++;//current number of interfaces
    icurr=ncurr-1;//index of added interface
    chgprot=index[ncurr-1]-1;//index of final protein (highest ratio)
    cout <<"Protein we are adding the interface:  "<<chgprot<<endl;
    tempiin=(ppidata *)realloc(tempiin, sizeof(ppidata) * ncurr);
    iin=(ppidata *)realloc(iin, sizeof(ppidata) * ncurr);
    for(i=0;i<ncurr-1;i++)
      tempiin[i]=iin[i];
    index= (long unsigned int *)realloc (index, sizeof(long unsigned int) * (ncurr+1));
    ratio=(double *)realloc(ratio, sizeof(double) *ncurr);
    num=wholep[chgprot].ninterface;//current number of interfaces on this protein
    wholep[chgprot].ninterface++;
    wholep[chgprot].valiface[num]=icurr; //new interface
    n2=ncurr*ncurr;
    Adj=(double *)realloc (Adj, sizeof(double) * (n2)); 
  
    cout <<"alter network: "<<endl;
    ntry=iin[chgprot].nppartner;
    for(i=0;i<ntry;i++){
      //reset temp ppi network 
      for(j=0;j<ncurr-1;j++)
	tempiin[j]=iin[j];
      
      id=iin[chgprot].pplist[i];
      cout <<"first protein with new partner: "<<id<<" i: "<<i<<endl;
    //now split off each partner, reevaluate the ratio, and compare
      tempiin[icurr].nppartner=1;
      tempiin[icurr].pplist[0]=id;
      //adjust the changed protein, and also adjust id's partners 
      tempiin[chgprot].nppartner=iin[chgprot].nppartner-1;
      
      num2=iin[chgprot].nppartner;
      t=0;
      for(j=0;j<num2;j++){
	if(iin[chgprot].pplist[j]!=id){
	  tempiin[chgprot].pplist[t]=iin[chgprot].pplist[j];
	t++;
	}
      }
      tempiin[id].nppartner=iin[id].nppartner; //id is no longer partnered with chgprot, but it is still partnered with icurr
      num2=iin[id].nppartner;
      cout <<"protein: "<<id<<" npartners: "<<num2<<endl;
      t=0;
      for(j=0;j<num2;j++){
	if(iin[id].pplist[j]!=chgprot){
	  tempiin[id].pplist[t]=iin[id].pplist[j];
	  t++;
	}
      }
      tempiin[id].pplist[t]=icurr;
      sprintf(fname, "Adj_matrix.alt%d.out", id);
      Mfile.open(fname);
      adj_make(ncurr, tempiin, Adj, Mfile);  
      Mfile.close();
      /*now we should have a new network, save in tempiin
	so test its ratios
      */
      calc_ratio(ncurr, tempiin, ratio);
      order_ratio(ncurr, ratio, index);
      cout <<"adding interface num connecting to: "<<id<<endl;
      avg_ratio[i+1]=0;
      for(j=0;j<ncurr;j++){
	ind=index[j+1]-1;
	cout <<"index: "<<ind<<" ratio: " <<ratio[ind]<<endl;;
	avg_ratio[i+1]+=ratio[ind];
      }
      avg_ratio[i+1]/=(1.0*ncurr);
      cout <<"Avg ratio: "<<avg_ratio[i+1]<<" protein: "<<id<<endl; 
    }
    /*Tried all different partners, now choose the solution with the best answer*/
    indexx(ntry, avg_ratio, index);
    ind=index[1]-1;
    avg_ratio[0]=avg_ratio[index[1]];
    id=iin[chgprot].pplist[ind];
    cout <<"protein with lowest ratio: "<<id<<endl;
    cout <<"Min ratio: "<<avg_ratio[0]<<endl;
    /*So now we can reset the list for this protein*/

    iin[icurr].nppartner=1;
    iin[icurr].pplist[0]=id;
    
    //adjust the changed protein, and also adjust id's partners 
    iin[chgprot].nppartner=iin[chgprot].nppartner-1;
    num2=iin[chgprot].nppartner+1;
    t=0;
    for(j=0;j<num2;j++){
      if(iin[chgprot].pplist[j]!=id){
	tempiin[chgprot].pplist[t]=iin[chgprot].pplist[j];
	t++;
      }
    }
    for(j=0;j<num2-1;j++)
      iin[chgprot].pplist[j]=tempiin[chgprot].pplist[j];//copy into iin
    
    
    num2=iin[id].nppartner;
    cout <<"protein: "<<id<<" npartners: "<<num2<<endl;
    t=0;
    for(j=0;j<num2;j++){
      if(iin[id].pplist[j]!=chgprot){
	tempiin[id].pplist[t]=iin[id].pplist[j];
	t++;
      }
    }
    iin[id].pplist[t]=icurr;
    for(j=0;j<num2-1;j++)
      iin[id].pplist[j]=tempiin[id].pplist[j];//copy into iin
    
    
    sprintf(fname, "NewAdj_matrix.%d.out", nit);
    Mfile.open(fname);
    adj_make(ncurr, iin, Adj, Mfile);  
    Mfile.close();
    nit++;
  }//end while loop over changin the network
  cout <<"Final number of interfaces after initial network creation: "<<ncurr<<endl;
  cout <<"Final ratio: "<<avg_ratio[0]<<endl;
  cout <<"Final degree: "<<2.0*Nedge/(1.0*ncurr)<<" Starting degree: "<<2.0*Nedge/(1.0*nwhole)<<endl;
  int nadd=ncurr-nwhole;
  cout <<"Number interfaces added: "<< nadd<<" Percent added: "<<nadd*1.0/(1.0*nwhole)<<endl;
  /*Now we should be able to write out which interfaces belong to which proteins, and how the interfaces connect to one another*/
  sprintf(fname, "Interface_NetworkStart.%d.out",ncurr);
  ofstream intfile(fname);
  for(i=0;i<ncurr;i++){
    intfile <<i<<'\t'<<iin[i].nppartner<<'\t';
    for(j=0;j<iin[i].nppartner;j++)
      intfile<<iin[i].pplist[j]<<'\t';
    intfile<<endl;
  }
  sprintf(fname, "Protein_InterfacesStart.%d.out", ncurr);
  ofstream pfile(fname);
  for(i=0;i<nwhole;i++){
    pfile <<i<<'\t'<<wholep[i].ninterface<<'\t';
    for(j=0;j<wholep[i].ninterface;j++)
      pfile<<wholep[i].valiface[j]<<'\t';
    pfile<<endl;
  }
  for(i=0;i<ncurr;i++){
    numpartners[i]=iin[i].nppartner;
    //cout <<"i: "<<i<<" iin[i]: "<<iin[i].nppartner<<endl;
    for(j=0;j<iin[i].nppartner;j++)
      Speclist[i][j]=iin[i].pplist[j];
  }
  for(i=0;i<nwhole;i++){
    for(j=0;j<wholep[i].ninterface;j++){
      id=wholep[i].valiface[j];
      p_home[id]=i;
      //cout <<"protein: "<<i<<" interface: "<<id<<endl;
    }
  }
  
  free(ratio);
  free(index);
  free(iin);
  free(tempiin);
  free(Adj);
  return ncurr;
  
} 
