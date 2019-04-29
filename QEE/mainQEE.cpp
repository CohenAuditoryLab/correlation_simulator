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


/*********************************************************************
 
                   COMMAND LINE INPUT FORMAT
 
 Command line instructions tell the program where to look for input
 files and where to send output, as well as the setting of various
 parameters (gamma, theta, etc) and flags (useSparse, etc).
 Note that numerical parameters may be entered in either scientific
 (recommended) or standard decimal notation. True/false switches are
 off by default - if entered in the command line, the corresponding
 option is set to true. Defaults can be found in dataStructures.h.
 The conventions are given below:
 
 -(flag name): (type of input expected)
 
 -d: string
    Default: "." (current directory)
    Path to the directory where the data file is located, and where
    output will be written.
 
 -i: string
    Default: "input"
    The location of the file containing a set of couplings, the
    starting values for the Monte Carlo learning algorithm.
 
 -o: string
    Default: "output"
    The location of the file where output is to be sent. Each
    different type of output file will have a different file type, 
    e.g. .j for couplings.
    
 -c: string
    Default: "input"
    The location of the file containing the set of correlations to
    reproduce (i.e. the correlations obtained from the data).
    
 -s: string
    Default: none
    Starting configuration for MC simulations.
    
 -g2: real number
    Default: 0.0
    The L2 regularization strength. L2 regularization is enabled by setting the 
    regularization strength to a nonzero value using this flag, or by using the
    -ag flag below.
 
 -ag: none
    Attempt to set the L2 regularization strengths to its optimal value,
    based on the number of samples (input) in the data.
    
 -b: real number
    Default: 1.0e+4
    Number of samples used to compute the correlations. Used to 
    determine the inference error.
    
 -mcb: integer
    Default: 8.0e+5
    Number of Monte Carlo samples to take to check inference error.
    
 -mcr: integer
    Default: 1
    Number of independent Monte Carlo runs to perform.
    
 -e: real number
    Default: 1.0
    Maximum tolerable error threshold. The MC learning algorithm will 
    continue to run until the error on the one- and two-point correlations 
    falls below this level.
    
 -v: none
    Enable verbose output.
 
 *********************************************************************/


// MAIN PROGRAM

int main(int argc, char *argv[]) {
    
	RunParameters r;
    
    // Process command line input
    
    for (int i=1;i<argc;i++) {
        
        // Location of input/output files
        
        if      (strcmp(argv[i],"-d")==0)   { if (++i==argc) break; else r.directory=argv[i];                               }
        else if (strcmp(argv[i],"-i")==0)   { if (++i==argc) break; else r.infile=argv[i];                                  }
        else if (strcmp(argv[i],"-o")==0)   { if (++i==argc) break; else r.outfile=argv[i];                                 }
        else if (strcmp(argv[i],"-s")==0)   { if (++i==argc) break; else { r.useStart=true; r.startfile=argv[i]; }          }
        else if (strcmp(argv[i],"-on")==0)  { if (++i==argc) break; else r.on=strtoint(argv[i]);                            }
        
        // Epitope input
        
        else if (strcmp(argv[i],"-n")==0)   { r.epitopeStates.push_back(std::vector<std::vector<int> >());
                                              r.epitopeStartSite.push_back(std::vector<int>());                             }
        else if (strcmp(argv[i],"-e")==0)   { if (++i==argc) break;
                                              else if (strcmp(argv[i],"[")==0) {
                                                std::vector<int> tempvec;
                                                while (strcmp(argv[++i],"]")!=0) tempvec.push_back(strtoint(argv[i]));
                                                (r.epitopeStates.back()).push_back(tempvec);
                                              } }
        else if (strcmp(argv[i],"-es")==0)  { if (++i==argc) break; else (r.epitopeStartSite.back()).push_back(strtoint(argv[i])); }
        
        // Monte Carlo settings
        
        else if (strcmp(argv[i],"-mcb")==0) { if (++i==argc) break; else r.b=strtoint(argv[i]);                             }
        else if (strcmp(argv[i],"-mcr")==0) { if (++i==argc) break; else r.runs=strtoint(argv[i]);                          }
        
        // Optional output
        
        else if (strcmp(argv[i],"-p2")==0)  { r.saveP2=true;                                                                }
        else if (strcmp(argv[i],"-v")==0)   { r.useVerbose=true;                                                            }
        
        else printf("Unrecognized command! '%s'\n",argv[i]);
        
    }
    
    if (r.epitopeStates.size()==0) {
        
        std::vector<std::vector<int> > eStates;
        std::vector<int> eStartSite;
        
        runEpitopeEval(r, eStates, eStartSite, -1);
        
    }
    
    else if (r.epitopeStates.size()==1) {
    
        std::vector<std::vector<int> > eStates(r.epitopeStates[0]);
        std::vector<int> eStartSite(r.epitopeStartSite[0]);
        
        runEpitopeEval(r, eStates, eStartSite, -1);
    
    }
    
    else { for (int i=0;i<r.epitopeStates.size();i++) {
    
        std::vector<std::vector<int> > eStates(r.epitopeStates[i]);
        std::vector<int> eStartSite(r.epitopeStartSite[i]);
    
        runEpitopeEval(r, eStates, eStartSite, r.on+i);
        
    } }
    
	return 0;

}
