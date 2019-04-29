#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

#include "qEpitopeEval.h"   // Evaluate epitopes
#include "monteCarlo.h"     // Monte Carlo
#include "tools.h"          // Numerical tools
#include "io.h"             // Data I/O
#include "mtrnd.h"          // Random number generator


// Runs the main program

void runEpitopeEval(RunParameters &r, std::vector<std::vector<int> > &eStates, std::vector<int> &eStartSite, int nout) {

    // Define variables
    
    Vector p;
    Vector J, expJ;
    
    // Retrieve couplings from file
    
    FILE *dataIn=fopen(r.getInfile().c_str(),"r");
        
    if (dataIn!=NULL) getCouplings(dataIn,J);
    else { printf("Error reading input from file %s",r.getInfile().c_str()); exit(1); }
    
    // Resize p and expJ
    
    for (int i=0;i<J.size();i++) { p.push_back(std::vector<double>(J[i].size(),0)); expJ.push_back(std::vector<double>(J[i].size(),0)); }
    for (int i=0;i<J.size();i++) { for (int j=0;j<J[i].size();j++) expJ[i][j] = exp(J[i][j]); }
	
    // Prepare Monte Carlo
    
    int N = sizetolength(J.size());
    
    if (r.useVerbose) printf("Got N=%d, len(h[0])=%d\n",N,(int)J[0].size());
    
    // Get default starting configuration, if nontrivial
    
    std::vector<int> lattice(N);
    
    if (r.useStart) {
    
        FILE *startIn=fopen(r.getStartInfile().c_str(),"r");
        for (int i=0;i<N;i++) fscanf(startIn,"%d",&lattice[i]);
    
    }
    
    else { for (int i=0;i<N;i++) lattice[i]=(int) p[i].size(); }
    
    // Prepare to simulate
    
    srand((unsigned)time(0));
    
    // Run MC and output Potts configurations
    
    if (r.saveStates) {
        
        FILE *dataOut=fopen(r.getStatesOutfile(nout).c_str(),"w");
        getEpitopeStates(dataOut, J, expJ, eStates, eStartSite, r.b, r.runs, p, lattice);
        fclose(dataOut);
        
    }
    
    // Run MC and get error
    
    getEpitopeCorrelations(J, expJ, eStates, eStartSite, r.b, r.runs, p, lattice);
    
    // Write out correlations
    
    FILE *dataOut=NULL;
    
    if (nout>=0) dataOut=fopen(r.getCorrelationsOutfile(nout).c_str(),"w");
    else         dataOut=fopen(r.getCorrelationsOutfile().c_str(),    "w");
    
    for (int i=0;i<N;i++) {
    
        fprintf(dataOut,"%.6e",p[i][0]);
        for (int j=1;j<p[i].size();j++) fprintf(dataOut,"\t%.6e",p[i][j]);
        fprintf(dataOut,"\n");
        
    }
    
    if (r.saveP2) { for (int i=N;i<p.size();i++) {
    
        fprintf(dataOut,"%.6e",p[i][0]);
        for (int j=1;j<p[i].size();j++) fprintf(dataOut,"\t%.6e",p[i][j]);
        fprintf(dataOut,"\n");
        
    } }
    
    fclose(dataOut);
    
}

