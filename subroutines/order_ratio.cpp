#include "pro_classes.h"
#include "calc_ratios.h"
#include "utility_calls.h"

void order_ratio(int nwhole, double *ratio, long unsigned int *index)
{
  int i;
  double value[nwhole+1];
  for(i=0;i<nwhole;i++)
    value[i+1]=ratio[i];
  indexx(nwhole, &value[0], index);//indexing is like fortran, from 1:N!!
  //in ascending order


}
