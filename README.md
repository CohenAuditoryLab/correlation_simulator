# correlation_simulator
Creating simulated spike trains with known amount of correlation and comparing their results to ACE fits.

Gen_STcorr_v3 first generates a set of pairs of neurons with known correlations and saves the spike times of these neurons as well as the 20ms bins in which they spike, their firing rates, and an upper triangular cross-correlation matrix across all generated neurons. This data is fed into convert_data_to_mat which takes the spike time data and formats it into a .dat file for fitting. Then, WriteCMSAbin is called, which writes the .p file (amongst others) for fitting. The ACE repository is called, which produces the .j file of parameters. Finally, the same .j file is fed into QEE, to see what data our fitted file outputs. Finally, the .dat file produced by QEE is compared to our original .dat file. 

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
