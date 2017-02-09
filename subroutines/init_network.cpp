#include "pro_classes.h"
#include "constrainParms.h"
#include "init_iin_net.h"

int init_network(int nwhole, ppidata *ppi, Protein *wholep, int *p_home, constrainParms &plist, int Nedge, int **Speclist, int *numpartners )
{
  //given a PPI network, create a IIN with a specific number of interfaces
  /*If the protein has a single partner, it only has a single interface
    if the protein has multiple partners, it can have from 0:n interfaces
  */


  /*One way would be to start out with same number of interfaces as proteins and
   start splitting
   A reasonable metric, in keeping with the highest specificity IINs, would be to
   count the number of partners your neighbors have, divided by the number of neighbors.
   This number should be close to 1.
  */
  int i, j;
  int nbor[nwhole];
  int ind;
//   double *ratio=new double[nwhole];
//   long unsigned int *index=new long unsigned int[nwhole+1];
//   ppidata * iin=new ppidata [nwhole];
//   ppidata * tempiin=new ppidata[nwhole];
  int id;
  int icurr;
  double highratio;
  int chgprot;
  int num2, t, num, n2;
  int ncurr=nwhole;

  char fname[100];
  
  double avg_ratio[MAXP];
  int ntry;
  
  for(i=0;i<nwhole;i++){
    //iin[i]=ppi[i];
    wholep[i].ninterface=1;//currently all proteins have 1 interfaces
    wholep[i].valiface[0]=i;//your interface is the same as your protein index
  }
  

  /*Now the proteins are ordered*/
  /*At this stage can add an interface on this protein, randomly, or loop through all possible
    additions, and take the best one*/
  
  int nit=0;

  avg_ratio[0]=2;
    
  /*Code added by DH to add interfaces to proteins (7-25-14). Will try to have an average of about 3.5 interfaces/protein*/
    
   /* double rnum;
    
    for(i=0;i<nwhole;i++){
        //Cannot have more interfaces than edges, so first check to see how many partners
        j=1;
        while (j<ppi[i].npartner){
            
        }
            
            
    }*/
    
  
  cout <<"Final number of interfaces after initial network creation: "<<ncurr<<endl;
  cout <<"Final ratio: "<<avg_ratio[0]<<endl;
  cout <<"Final degree: "<<2.0*Nedge/(1.0*ncurr)<<" Starting degree: "<<2.0*Nedge/(1.0*nwhole)<<endl;
  int nadd=ncurr-nwhole;
  cout <<"Number interfaces added: "<< nadd<<" Percent added: "<<nadd*1.0/(1.0*nwhole)<<endl;
  /*Now we should be able to write out which interfaces belong to which proteins, and how the interfaces connect to one another*/
  // sprintf(fname, "Interface_NetworkStart.%d.out",ncurr);
//   ofstream intfile(fname);
//   for(i=0;i<ncurr;i++){
//     intfile <<i<<'\t'<<iin[i].nppartner<<'\t';
//     for(j=0;j<iin[i].nppartner;j++)
//       intfile<<iin[i].pplist[j]<<'\t';
//     intfile<<endl;
//   }
//   sprintf(fname, "Protein_InterfacesStart.%d.out", ncurr);
//   ofstream pfile(fname);
//   for(i=0;i<nwhole;i++){
//     pfile <<i<<'\t'<<wholep[i].ninterface<<'\t';
//     for(j=0;j<wholep[i].ninterface;j++)
//       pfile<<wholep[i].valiface[j]<<'\t';
//     pfile<<endl;
//   }
  for(i=0;i<ncurr;i++){
    numpartners[i]=ppi[i].nppartner;
    //cout <<"i: "<<i<<" iin[i]: "<<iin[i].nppartner<<endl;
    for(j=0;j<ppi[i].nppartner;j++)
      Speclist[i][j]=ppi[i].pplist[j];
  }
  for(i=0;i<nwhole;i++){
    for(j=0;j<wholep[i].ninterface;j++){
      id=wholep[i].valiface[j];
      p_home[id]=i;
      //cout <<"protein: "<<i<<" interface: "<<id<<endl;
    }
  }
  

  return ncurr;
  
} 
