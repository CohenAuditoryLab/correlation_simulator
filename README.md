# correlation_simulator
Creating simulated spike trains with known amount of correlation and comparing their results to ACE fits.

TO RUN: Run IsingTestPipeline.m

  Within IsingTestPipeline.m, you can set parameters.
  
    Line 2: desired_c is an array of desired correlations between pairs of neurons. Its length must be <= numpairs.
    
    Line 3: numpairs is the number of pairs of neurons you want to generate.
    
    Line 8: fraction is the fraction of data from the generated neurons you want to run the ACE fitting algorithm on. This may need to be dependent on your computer's memory.
    
    Lines 32,33: fill in path to ACE repository on your computer.
    
    Lines 46, 49: fill in path to QEE repository on your computer.

The ACE Repository can be found here: https://github.com/johnbarton/ACE
  Within this repository, the scripts 'WriteCMSA' and 'WriteCMSAbin' can be found in the directory 'scripts'.
The QEE Repository was obtained from Cle Demu.
Examples running it on our data in Matlab can be found here: https://github.com/CohenAuditoryLab/auditory-ising/tree/master/ACE/running
