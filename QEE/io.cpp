#include <vector>

#include "io.h"     // Data input and output
#include "tools.h"  // Numerical tools


// Read correlations and number of sites in from a file

void getCorrelations(FILE *input, Vector &p) {

    double c;
    char o;
    
    while (fscanf(input,"%le",&c)==1) {
    
        p.push_back(std::vector<double>());
        (p.back()).push_back(c);
        
        while (fscanf(input,"%c",&o)==1) {
    
            if (o=='\n' || o=='\r') break;
            
            fscanf(input,"%le",&c);
            (p.back()).push_back(c);
            
        }
        
    }
	
}


// Read couplings in from a file

void getCouplings(FILE *input, Vector &p) {

    getCorrelations(input, p);
	
}


// Read secondary structure from file 

void getSS(FILE *input, IntVector &p) {

    int c;
    char o;
    
    while (fscanf(input,"%d",&c)==1) {
    
        p.push_back(std::vector<int>());
        (p.back()).push_back(c);
        
        while (fscanf(input,"%c",&o)==1) {
    
            if (o=='\n' || o=='\r') break;
            
            fscanf(input,"%d",&c);
            (p.back()).push_back(c);
            
        }
    
    }

}


// Print final couplings out to file in a form understood by ising.cpp

void printCouplings(FILE *output, const Vector &J) {
    
    for (int i=0;i<J.size();i++) {
    
        fprintf(output,"%.6e",J[i][0]);
        
        for (int j=1;j<J[i].size();j++) fprintf(output,"\t%.6e",J[i][j]);
        
        fprintf(output,"\n");
        
    }
    
    fflush(output);
	
}


// Print supplementary information to file

void printSupplementaryOutput(FILE *output, double theta, const std::vector<double> &error, double finalS, int maxClusterSize, unsigned long numClusters, unsigned long numSignificantClusters) {
    
    fprintf(output,"%le",theta);
    for (int i=0;i<error.size();i++) fprintf(output,"\t%le",error[i]);
    fprintf(output,"\t%le\t%d\t%lu\t%lu\n",finalS,maxClusterSize,numClusters,numSignificantClusters);
    fflush(output);
	
}

