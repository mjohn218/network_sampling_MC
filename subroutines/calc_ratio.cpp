#include "pro_classes.h"
#include "calc_ratios.h"
#include "utility_calls.h"

void calc_ratio(int nwhole, ppidata *ppi, double *ratio)
{
  int nbor[nwhole];
  int i, j;
  int ind;
  for(i=0;i<nwhole;i++){
    if(ppi[i].nppartner==1){
      //this protein only has one partner
      ratio[i]=1;
      nbor[i]=1;
    }else{
      //this protein has multiple partners
      //count up the partners of his neighbors to get the ratio 
      //now calculate neighbors
      
      nbor[i]=0;
      for(j=0;j<ppi[i].nppartner;j++){
	ind=ppi[i].pplist[j];
	nbor[i]+=ppi[ind].nppartner;//[ind];
      }
      
      /*Consider not dividing this by the number of partners, to weight 
	nodes of high connectivity, since they are the most troublesome
	for creating a large gap
      */
      ratio[i]=nbor[i]*1.0/(ppi[i].nppartner*1.0);
      
    }
    cout <<"i: "<<i<<" partners: "<<ppi[i].nppartner<<" nbor:  "<< nbor[i]<<" ratio: "<<ratio[i]<<endl;
  }
  
}
