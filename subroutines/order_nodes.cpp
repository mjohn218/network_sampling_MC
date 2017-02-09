#include "net_props.h"
#include <cstdlib>

void order_nodes(int Nnode, int *origlist, int *ordered)
{
  int i;
  int num=Nnode+1;
  long unsigned int *index=(long unsigned int *) malloc(sizeof(long unsigned int) *(num));
  double *value=new double[num];
  for(i=0;i<Nnode;i++)
    value[i+1]=origlist[i];
  indexx(Nnode, value, index);
  int id;
  for(i=0;i<Nnode;i++){
    
    id=index[i+1]-1;//indexing is from 1:N, skips zero
    ordered[i]=origlist[id];
  }

  free(index);
  
  delete[] value;

}
