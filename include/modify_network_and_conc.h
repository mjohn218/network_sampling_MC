
void modify_network_and_conc(int nwhole, int ninterface, int *numpartners, int **Speclist, Protein *wholep, ppidata* ppi, int *p_home, constrainParms &plist,  ofstream &globalfile, int **ehome, int **Epred, int **templist, int *tmppartners, Protein *wholetemp, int * nbor, int *Adj, double *Amat, double *indivconc,double *complexconc, double *cumul, double *yeast_abund, int *constrain, string *genid, double *abund, ofstream &idfile, ofstream &matchfile);
void modify_network_and_conc_ErbB(int nwhole, int ninterface, int *numpartners, int **Speclist, Protein *wholep, ppidata* ppi, int *p_home, constrainParms &plist,  ofstream &globalfile, int **ehome, int **Epred, int **templist, int *tmppartners, Protein *wholetemp, int * nbor, int *Adj, double *Amat, double *indivconc,double *complexconc, double *cumul, double *yeast_abund, int *constrain, string *genid, double *abund, ofstream &idfile, ofstream &matchfile);