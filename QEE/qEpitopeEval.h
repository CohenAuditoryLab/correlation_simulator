#ifndef QLEARNSPARSE_H
#define QLEARNSPARSE_H


#include <iostream>
#include <string>

#include "tools.h"  // Numerical tools


// This class holds the parameters needed for running the ising program

class RunParameters {
    
public:
    
    std::string directory;          // Path to the directory where the inut file is located
                                    // Output will also be sent to this directory
    std::string infile;             // Input file from which couplings are to be read
    std::string outfile;            // Output file (prefix) where data is to be written
    
    double b;                       // Number of data points to take in a single run
    double runs;                    // Number of times to run dynamics
    
    std::vector<std::vector<std::vector<int> > > epitopeStates; // Configuration of sites in the targeted epitope
    std::vector<std::vector<int> >            epitopeStartSite; // Index of the first site in the epitope
    
    int on;                         // Start number for keeping track of all epitopes and epitope combinations
    
    bool useStart;                  // Use different starting sequence than all wild-type
    std::string startfile;          // File containing the starting sequence
    
    bool useVerbose;                // If true, print extra information while program is running
    bool saveP2;                    // If true, record pair correlations in additions to one-point correlations
    bool saveStates;                // If true, record MC sampled configurations
    
    RunParameters() {
        
        directory=".";
        infile="input";
        outfile="output";
        b=1000000;
        runs=1;
        on=0;
        
        useStart=false;
        useVerbose=false;
        saveP2=false;
        saveStates=true;
        
    }
    std::string getInfile()                   { return (directory+"/"+infile+".j");                    }
    std::string getCorrelationsOutfile()      { return (directory+"/"+outfile+".p");                   }
    std::string getCorrelationsOutfile(int i) { return (directory+"/"+outfile+"-"+tostring(i)+".p");   }
    std::string getStatesOutfile()            { return (directory+"/"+outfile+".dat");                 }
    std::string getStatesOutfile(int i)       { return (directory+"/"+outfile+"-"+tostring(i)+".dat"); }
    std::string getStartInfile()              { return (directory+"/"+startfile+".dat");               }
    ~RunParameters() {}
    
};


//void runEpitopeEval(RunParameters &);
void runEpitopeEval(RunParameters &, std::vector<std::vector<int> > &, std::vector<int> &, int);


#endif
