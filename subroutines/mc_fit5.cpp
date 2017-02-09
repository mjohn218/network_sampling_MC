#include "pro_classes.h"
#include "mc_fits.h"

double mc_fit5(int ninterface, double *grid_cofs, double *clocals, double ksi, double beta, int nwhole, Protein *wholep, double mu, int *numpartners, double edgediff)
{
    
    double fitness=0;
    int i,j;
    
    for(i=0;i<ninterface;i++){
      if(numpartners[i]>1){
        fitness+=exp(ksi*(1-grid_cofs[i]))-1; //A high grid_cof means lower value (more fit system)
      }
    }


    double tmp;
    cout<<"ksi component of fitness is: " << fitness << endl;
    tmp = fitness;

    for(i=0;i<ninterface;i++){
      if(numpartners[i]>1){
	fitness+=exp(beta*clocals[i])-1; //clocal>0 means a triangle present. Penalize this.
      }
    }
    cout << "The beta component is: " << fitness-tmp << endl;
    tmp=fitness;
    //This version penalizes the number of extra interfaces, not interfaces per protein.
    
    fitness+=exp(mu*(ninterface-nwhole));

    cout << "The mu component is: " << fitness-tmp << endl;
    tmp=fitness;
    fitness+=exp(edgediff); //Penalizes having too many edges
    cout<<"The omega component is: " << fitness-tmp << endl;
    return fitness;
    
}
