//
//  main.cpp
//  AlphaAnalyze
//
//  Created by Benjamin Shapiro on 2/20/15.
//  Copyright (c) 2015 Johns Hopkins. All rights reserved.
//
//  Edited version by David Holland on 6/20/2015

#include "findAlpha.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////
int index(int row , int column, int width)
{
    /* this function allows storage of 2d array in 1d array using index conversion */
    int indexVal;
    
    indexVal= row * width + column; // convert index using standard formula
    
    return indexVal;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////



void printarray( int *array, int length)
{
    /* prints an integer array */
    for (int i = 0; i < length; i++)
    {
        cout << array[i] << endl;
    }
    return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void printarraydoub( double *array, int length)
{
    /* prints an double array */
    for (int i = 0; i < length; i++)
    {
        cout << array[i] << endl;
    }
    return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////


void print2darray (int *TwoDarray, int rows, int columns )
{
    // prints a 2d array using Prof. Johnson conventions
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < columns; ++j)
        {
            cout << TwoDarray[index(i, j, columns)]<< ' ';
        }
        cout << endl;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<double> divideArray(vector<double> array, double divisor, int length)
{  /* This function divides a double array with a double */
  //double*  dividedArray = new double[length];  // allocate space for scale free array
  vector<double> dividedArray(length);
    for (int count = 0 ; count < length ; count++ ) // itter through every element
    {
        //cout << array[count] << endl;
        dividedArray[count] =  array[count]/divisor; // divide each element of the array
        //cout << dividedArray[count];
    }
    return dividedArray;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////


double sum_array(vector<double> a, int num_elements)
{
    /* This function returns the sum of an array */
    
    int count; // counter for looping through array
    double sum; // hold sum of array
    sum = 0; // initialize sum
    for (count = 0 ; count < num_elements ; count++ )
    {
        sum = sum + a[count];
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

int sum_arrayint(int *a, int num_elements)
{
    /* This function returns the sum of an INTEGER array */
    
    int count; // counter for looping through array
    int sum; // hold sum of array
    sum = 0; // initialize sum
    for (count = 0 ; count < num_elements ; count++ )
    {
        sum = sum + a[count];
    }
    return sum;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<double> Cum_Sum(vector<double> array, int length)
{
    /* this function returns the cumulative sum of an array and creates a pointer to a NEW variable which it returns */
    //double* cumsum = new double[length];  // allocate space for scale free array
  vector<double> cumsum(length);
  //  double currentsum=0;
    for (int count = 0 ; count < length ; count++ ) // itter through every element
    {
      //  currentsum+=array[count];
      cumsum[count] = sum_array(array, count + 1); // call sum to sum up count elements before current
        
    }
    return cumsum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
double *ScaleFree(int Nodes, double alpha)
{
    /* This function returns an array which represents a scale free network with parameter alpha */
    
    double* Pn = new double[Nodes];  // allocate space for scale free array
    
    for (int x = 0; x <= Nodes - 1 ; x++)        // fill scale free array
    {
        Pn[x] = pow((x+1),-alpha);
        // cout <<Pn[x] <<endl; // Debugging line
    }
    
    return Pn;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
int CountOrphans(int *array, int length)
{
    /* This function takes  1 x n array and outputs the corresponding number of 0 elements */
    int sum = 0; // this will keep track of the number of 0s we have found
    for (int i = 0 ; i < length; i++)
    {
        if (array[i] == 0) // do we have a hit?
        {
            sum++;
        }
    }
    return sum;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
int * ZeroArray( int * array, int length )
{
    /* This function will take all the elements of an integer array and turn them into 0s */
    for (int i = 0 ; i < length; i++)
    {
        array[i] = 0;
    }
    
    return array;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> ZeroArrayDouble( vector<double> array, int length )
{
    /* This function will take all the elements of a double array and turn them into 0s */
    for (int i = 0 ; i < length; i++)
    {
        array[i] = 0;
    }
    
    return array;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<int>  locateIndexOfZeros(vector<int> Vect, int * array, int length)
{
    /* this function takes a 1 x length array and return a VECTOR with the indexed positions of the 0s. */
    // vector<int> indexContainer; // define a vectore to hold the indices for 0s ALLOCATED MEMORY
    Vect.clear(); // clear vector
    for (int i  = 0; i < length; i++)
    {
        if (array[i] == 0 )
        {
            //indexContainer.push_back(i); // this will insert array[i] to last position
            Vect.push_back(i); // insert array[i] to last position
        }
        
    }
    
    //return indexContainer;
    return Vect;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////

void printVect(vector<int> vecArray)
{
    /* This function prints a vector of unknown length */
    
    for (int i = 0; i < vecArray.size(); i++)
    {
        cout << vecArray[i] << endl;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////

int *NeighborArray( int *TwoDArray, int Nodes)
{
    
    /* this function will sum a Two dimentional array represented as a 1 dimentional array creating a single array with the sum of each row in the previous inputed matrix */
    int * SumVector = new int[Nodes]; // allocate space for the neighborhood vector to store # neighbors per node
    SumVector  = ZeroArray(SumVector, Nodes);
    //cout << endl << "SUM";
    //printarray(SumVector, Nodes);
    //cout << endl;
    for (int i = 0 ; i < Nodes; i++) // loop through the rows
    {
        for (int j = 0; j < Nodes; j++ ) // loop through the columns
        {
            SumVector[i] = SumVector[i] + TwoDArray[index(i, j, Nodes)];
            
        }
    }
    
    return SumVector;
    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Pick(double *ProbVect, int Nodes)
{
    /* This function takes an input of two arrays of the same length. It will chose an element from Nodvect with corresponding ProbVect Probability */
    
    double randomNumber; // this will be the number we generate
    double normalizingValue; // sum of probability vector for normalizatio
    vector<double> normalizedProbabilities; //vector to hold normalized probabilities
    vector<double> cumulativeProbabilities; // holds the cumulative probabilities
    randomNumber = rand() % 10000 ;
    randomNumber = randomNumber/10000;  // uniform with discrete probabilities of 0,.0001,...,.9999
    vector<double> ProbVect2(ProbVect, ProbVect + Nodes);
    normalizingValue = sum_array(ProbVect2, Nodes); // sum the vector
    normalizedProbabilities = divideArray(ProbVect2, normalizingValue, Nodes); // normalize probability vect
    cumulativeProbabilities= Cum_Sum(normalizedProbabilities, Nodes); // array containing cumulative probs
    //delete [] normalizedProbabilities; // FREE MEMORY
    // printarray(cumulativeProbabilities, Nodes);
    for (int i = 0 ; i < Nodes ; i++)        // fill scale free array
    {
        if (randomNumber <= cumulativeProbabilities[i])
        {
	  //delete [] cumulativeProbabilities; // FREE MEMORY
            return i;
        }
        
    }
cout << "ERRRRROOOOOORRRRRRRR See Pick function";
return 10000; // THIS WILL YIELD AN ERROR. PROGRAM SHOULD NEVER REACH THIS STATEMENT
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
int leaveout(double *NormProbs, int Nodes, int Edges, double Alpha)
{

  vector<double> NormProbs2(NormProbs, NormProbs+Nodes);
    double sumPn = sum_array(NormProbs2, Nodes);
    double EOrphans = 0; // this holds the expected number of orphans
    double EOrphansPrev = 1; // stores prev. number of orphans calculated
    vector<double> Summingvector(Nodes);
    Summingvector=ZeroArrayDouble(Summingvector, Nodes);
    double temp = 0;
    
    while(abs(EOrphans - EOrphansPrev) > .01)
    {
        EOrphansPrev = EOrphans;
        
        for (int i = 0; i < Nodes; i++)
        {
            
            temp = pow(i+1, -Alpha);
            Summingvector[i] = pow((1-temp/sumPn), 2*(Edges-EOrphans));
        }
        EOrphans = sum_array(Summingvector, Nodes);
        
    }
    //    delete [] Summingvector;
    
    return (ceil(EOrphans));
    
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////


int * GenerateAlphaNetwork(int Nodes, int Edges, double alpha)
{
  //cout<<"Inside GenerateAlphaNetwork" << endl;  
  double *Pn; // define the array which will hold the scale free network
    int *Network = new int[Nodes * Nodes]; // this holds the network storage bin
    int Eorphans; // the number of orphans we will expect
    int count; // this is a counter for while loop
    int Orphans; // this will hold the current count for the number of orphans we have in the network
    int * Degrees; // this will hold an array with the degree of each node
    int a; // hold the place of chosen nodes to be connected
    int b; // hold the next connection of node
    vector<int> orphanedIndexes; // this will hold the indexes of all the orphaned nodes
    int Mcurrent = Edges; // in case we have too many orphans
    Pn = ScaleFree(Nodes,alpha);  // generate a pointer to a scale free network
    //printarray(Pn, Nodes); //Debugging line
   // Eorphans = 100; // This is a rough estimate on the number of orphans we will expect
    Eorphans = leaveout(Pn, Nodes, Edges, alpha);
    
    //cout << " LEAVING OUT: " << Eorphans << endl;
    
    count = 0;
    //print2darray(Network, Nodes, Nodes);
    
    
    Network = ZeroArray(Network, Nodes*Nodes); // Initialize Network to be 0s NEED TO FREE
    
    if (Eorphans < Edges && Eorphans >= 0 )
    {
    while (count < Edges - Eorphans ) //
    {
        a = Pick(Pn, Nodes);  // pick first node with corresponding probability
        b = Pick(Pn, Nodes); // pick second node with corresponding probability
        // cout << Network[index(a,b,Nodes)] <<endl;
        if (Network[index(a,b,Nodes)] == 0)
        {
            Network[index(a,b,Nodes)] = 1;  // set matrix value to one
            Network[index(b,a,Nodes)] = 1; // set opposing matrix value to one
            count = count + 1 ;
        }
        Degrees = NeighborArray(Network, Nodes); // returns array with varying degrees. ALLOCATED MEMORY
        Orphans = CountOrphans(Degrees, Nodes); // we now need to find the number of orphans in the degree vector!
        if(count == (Edges - Eorphans) && (Orphans > Eorphans)) // if we have more orphans then edges left to add *THIS WAS CHANGED ON MARCH 29th*
        {
            count = 0;    // trash it
            //print2darray(Network, Nodes, Nodes);
            Network = ZeroArray(Network, Nodes*Nodes); // Reset Network to be 0s NEED TO FREE
            // print2darray(Network, Nodes, Nodes);
        }
        //printarray(Degrees,Nodes);
        delete [] Degrees; // FREE SPACE allocated to degrees ****
    }
    
    Degrees = NeighborArray(Network, Nodes); // returns array with varying degrees. ALLOCATED MEMORY
    Orphans = CountOrphans(Degrees, Nodes); // we now need to find the number of orphans in the degree vector!
    //printarray(Degrees,Nodes);
    while (Orphans > 0 )  // continue until no orphans left
    {
        orphanedIndexes = locateIndexOfZeros(orphanedIndexes, Degrees, Nodes); // find indexes of orphans ALLOCATED MEMORY
        a = orphanedIndexes[0]; // chose orphaned index
        b = Pick(Pn, Nodes); // pick independent prob of another node
        if (Network[index(a,b,Nodes)] == 0) //check to make sure the nodes do not share an edge
        {
            Network[index(a,b,Nodes)] = 1;  // assign node matrix place holder
            Network[index(b,a,Nodes)] = 1;  //  assign node matrix place holder
            Eorphans  = Eorphans - 1 ; //  update edges left
        }
        delete [] Degrees; // FREE MEMORY of the old pointer *QUESTION - WHY DOES THIS NOT WORK?*
        //cout << Orphans << endl<< endl ;
        Degrees = NeighborArray(Network, Nodes); // returns array with varying updated degrees. ALLOCATED MEMORY
        //printarray(Degrees, Nodes);
        Orphans = CountOrphans(Degrees, Nodes); // we now need to find the number of orphans in the degree vector
    }
    
    // BUT WAIT! What if we have some extra edges left over... Let's add them back in.
    
    while(Eorphans > 0)
    {
        a = Pick(Pn, Nodes);  // pick first node with corresponding probability
        b = Pick(Pn, Nodes); // pick second node with corresponding probability
        if (Network[index(a,b,Nodes)] == 0) //check to make sure the nodes do not share an edge
        {
            Network[index(a,b,Nodes)] = 1;  // assign node matrix place holder
            Network[index(b,a,Nodes)] = 1;  //  assign node matrix place holder
            Eorphans = Eorphans - 1;
            
        }
    }
    delete [] Degrees;
    delete [] Pn;
    //print2darray(Network, Nodes, Nodes);  // these print the network & Degrees for debugging purposes
    //cout << endl;
    // cout << endl;
    // printarray(Degrees, Nodes);
    // cout << sum_arrayint(Degrees, Nodes);
    return Network;
        
        
        /* WHAT IF THERE ARE TOO MANY INTERFACES!?!?! WE DEAL WITH THAT HERE: */
    }
    else  // the expected value of orphans is more than the number of edges available. We attach one by one
    {
        Degrees = NeighborArray(Network, Nodes); // returns array with varying degrees. ALLOCATED MEMORY
        Orphans = CountOrphans(Degrees, Nodes); // we now need to find the number of orphans in the degree vector!
        delete [] Degrees;
        while (Orphans > 0)
        {
            
            if ( Orphans  <=  Mcurrent*2 - 2 )
            {
                a = Pick(Pn, Nodes);  // pick first node with corresponding probability
                b = Pick(Pn, Nodes); // pick second node with corresponding probability
                if (Network[index(a,b,Nodes)] == 0) //check to make sure the nodes do not share an edge
                {
                    Network[index(a,b,Nodes)] = 1;  // assign node matrix place holder
                    Network[index(b,a,Nodes)] = 1;  //  assign node matrix place holder
                    Mcurrent = Mcurrent - 1;
                }
            } else
            {
                Degrees = NeighborArray(Network, Nodes); // returns array with varying degrees. ALLOCATED MEMORY
                orphanedIndexes = locateIndexOfZeros(orphanedIndexes, Degrees, Nodes); // find indexes of orphans ALLOCATED MEMORY
                delete [] Degrees;
                a = orphanedIndexes[0]; // chose orphaned index
                b = orphanedIndexes[1]; // chose other orphaned index
                if (Network[index(a,b,Nodes)] == 0) //check to make sure the nodes do not share an edge
                {
                    Network[index(a,b,Nodes)] = 1;  // assign node matrix place holder
                    Network[index(b,a,Nodes)] = 1;  //  assign node matrix place holder
                    Mcurrent = Mcurrent - 1;
               }

            }
            Degrees = NeighborArray(Network, Nodes); // returns array with varying degrees. ALLOCATED MEMORY
            Orphans = CountOrphans(Degrees, Nodes); // we now need to find the number of orphans in the degree vector!
            delete [] Degrees;
        }
        //cout << Mcurrent << endl;

    }
                        
    
    return Network;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
double findvals(int * array, int value, int length)
{ /* this function finds the number of times a value is in an array */
    double count = 0; // keep track of the number of values we have found
    for (int i = 0; i < length; i++)
    {
        if (array[i] == value )
        {
            count = count + 1;
        }
    }

    return count;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> PDFcalculate (int * Network, int dimentions)
{
    int * PDFHolder; // define our PDF in integer form
    vector<double> PDF(dimentions); // more permanent double form
    PDFHolder = NeighborArray(Network, dimentions); // we now have a vector with all the neighbors
    //printarray(PDFHolder, dimentions);
    vector<double> tempholderfordeletion;
    for (int i = 0; i < dimentions; i++ )
    {
        PDF[i] = findvals(PDFHolder, i+1, dimentions);
    }
    
    //printarraydoub(PDF, dimentions);
    //tempholderfordeletion = PDF;
    PDF = divideArray(PDF,(double)dimentions, dimentions); // **HOW DO I FREE MEMORY FOR THIS CASE?
    //delete [] tempholderfordeletion; //deletes old address of pdf Changed on april 1
    //printarraydoub(PDF, dimentions);
    
    delete [] PDFHolder; // return memory
    
    return PDF;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

vector<double> addTwoArrays(vector<double> array1, vector<double> array2, int length)

/* this function takes two arrays and adds the elements of each together. It stores in first array */
{
    //double*  addedArray = new double[length];  // allocate space for scale free array

    for (int i = 0 ; i < length; i++ )
    {
        array1[i] = array1[i] + array2[i];
    }

    return array1;
    
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////
vector<double> AvgPDF (int Nodes, int Edges, double alpha, int numIters)
{
    /* This function takes an input of alpha, the edges of a network, and the number of nodes. It outputs an average pmf of the alpha network. numIters is the number of networks to average over */
    int *AlphaNetwork;  // this defines our network
    vector<double> pdf;
    vector<double> rollingsum(Nodes);
    //double * avgpmf = new double [Nodes];
    vector<double> avgpmf; // replaced on April 1st
    rollingsum = ZeroArrayDouble(rollingsum, Nodes);
    for (int i = 0; i < numIters; i++)
    {
        AlphaNetwork = GenerateAlphaNetwork(Nodes, Edges, alpha ); // create network
        pdf = PDFcalculate(AlphaNetwork, Nodes); // get pdf of network
        rollingsum = addTwoArrays(rollingsum,pdf, Nodes); // this function adds arrays and stores in rollingsum
        //delete [] pdf;
        delete [] AlphaNetwork;
    }
    avgpmf = divideArray(rollingsum, numIters, Nodes);
    //printarraydoub(avgpmf, Nodes);
    //delete rollingsum;
    //printarraydoub(avgpmf, Nodes); // ** TESTING **
    return avgpmf;
    
    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
double chierror(vector<double> pmf1 , vector<double> pmf2, int length)
{
  //cout << "Inside chierror\n";  
  double  rollingtracker; // keeps track of error
  //double * cdf1;
  //double * cdf2;
  vector<double> cdf1(length);
  vector<double> cdf2(length);
  cdf1 = Cum_Sum(pmf1, length); // find cdfs
  //cout <<"cdf1 found. ";
  cdf2 = Cum_Sum(pmf2, length);
  //cout << "cdf2 found.\n";
  rollingtracker = 0 ; // initialize value
  for (int i = 0; i < length; i++ ) {

    rollingtracker = rollingtracker + pow((cdf1[i]-cdf2[i]),2);

  }
  //cout << "What is rollingtracker? " << rollingtracker << endl;
  return rollingtracker;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////

int indexofSmallestElement(vector<double> array, int size)
{
    int index = 0;
    
    for(int i = 1; i < size; i++)
    {
        if(array[i] < array[index])
            index = i;
    }
    
    return index;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////
double findAlpha(double *pmfCompare , int Edges, int Nodes, int nwhole, int maxni, int PPIedge, int maxne, double **sampPMFs)
{
    /* this network takes a pmf of known network with an inputed number of edges and nodes. It then finds the alpha which corresponds to that network with an accuracy of .1 */


  double alpha;
  //double * avgpdf; // pdf to compare to
  //Using a vector instead to avoid memory errors. Edited by DH
  //vector<double> avgpdf(Nodes);
  vector<double> vec_pmf(pmfCompare, pmfCompare+Nodes);
  vector<double> ChiValues(11);
  ChiValues = ZeroArrayDouble(ChiValues, 11);
  double truealpha; // this is the return value
  //New variables by DH
  int line, j;
  int dim1 = maxne-PPIedge+1;
  for (int i = 0; i <= 10; i++)
    {
      //cout << endl;
      //cout << i ;
      alpha = .1*(i);
      //      avgpdf = AvgPDF(Nodes, Edges, alpha, 30); // call the control function which net generator will call
     
      /*Edited code by David Holland: drawing avgpdf from a file*/
      //Find the right PMF to compare to
      line = (Nodes-nwhole)*11*dim1 + (Edges-PPIedge)*11 + i;
      //cout << "line equals: " << line;
      //avgpdf = new double[Nodes];
      //avgpdf = sampPMFs[line];
      vector <double> avgpdf(sampPMFs[line], sampPMFs[line]+Nodes);

      /*End edited code*/
      
      ChiValues[i] = chierror(avgpdf,vec_pmf,Nodes);

    }

  // we now need ot find the minimum chierror
  //printarraydoub(ChiValues, 11);
  truealpha = (indexofSmallestElement(ChiValues, 11))*.1; // find smallest index and convert to alpha

  //delete []ChiValues; //DELETE CHI VALUES ** CHANGE MARCH 30th
  //No longer necessary to delete since pointer arrays have been changed into vectors (DH)
  return truealpha;
  
}

double findAlpha2(double *pmfCompare , int Edges, int Nodes)
{
    /* this network takes a pmf of known network with an inputed number of edges and nodes. It then finds the alpha which corresponds to that network with an accuracy of .1 */


  double alpha;
  //double * avgpdf; // pdf to compare to
  //Using a vector instead to avoid memory errors. Edited by DH
  vector<double> avgpdf(Nodes);
  vector<double> vec_pmf(pmfCompare, pmfCompare+Nodes);
  vector<double> ChiValues(11);
  ChiValues = ZeroArrayDouble(ChiValues, 11);
  double truealpha; // this is the return value
  
  for (int i = 0; i <= 10; i++)
    {
      //cout << endl;
      //cout << i ;
      alpha = .1*(i);
      avgpdf = AvgPDF(Nodes, Edges, alpha, 30); // call the control function which net generator will call
      ChiValues[i] = chierror(avgpdf,vec_pmf,Nodes);
      
    }

  // we now need ot find the minimum chierror
  //printarraydoub(ChiValues, 11);
  truealpha = (indexofSmallestElement(ChiValues, 11))*.1; // find smallest index and convert to alpha

  //delete []ChiValues; //DELETE CHI VALUES ** CHANGE MARCH 30th
  //No longer necessary to delete since pointer arrays have been changed into vectors (DH)
  return truealpha;
  
}

void readPMFfile(ifstream &pmf_file, double **PMFs, int lines)
{

  int nodecurr, edgecurr;
  double alphacurr;
  string line;
  int i,j;
  getline(pmf_file,line);//Skip header line
  for(i=0;i<lines;i++){
    getline(pmf_file,line);
    stringstream stream(line);
    stream >> nodecurr;
    stream >> edgecurr;
    stream >> alphacurr;
    PMFs[i]=new double[nodecurr];
    for(j=0;j<nodecurr;j++){
      stream >> PMFs[i][j];
    }
  }
  pmf_file.close();

}
