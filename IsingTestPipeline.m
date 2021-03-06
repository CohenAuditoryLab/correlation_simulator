%% Produce fake neurons 
desired_c = [.1];
numpairs = 2;
trainlength = Gen_STcorr_v3(desired_c, numpairs);

%% Convert to a dat file
file = [num2str(numpairs) '_pairs.mat'];
fraction = 0.1; % fraction of data you want to run ACE on - taken from beginning of array (t = 0)
convert_data_to_dat(file, fraction);

%% Convert them to a p file (input to ACE) 
%  also requires WeightCalculator

pieces = strsplit(file, '.');
% filetype = 'binary';
% filetype = 'numbers';
filetype = 'neuro';
filename = [pieces{1} '.dat']; 
theta = 0;
bin = 20;
redmethod = 'frequency';
redcut = 0;
gauge = 'least'; 
gapred = 0;
% WriteCMSA(filetype, filename, theta, redmethod, redcut, gauge, gapred);
WriteCMSAbin(filetype, filename, theta, bin) 

%% fit Ising to .dat file, produce .j file  

pieces = strsplit(file, '.');
inputfile = [pieces{1} '_po0_least'];
addpath(genpath('/Users/briannakarpowicz/Documents/CohenLab/ACE-master'));
cmd = ['/Users/briannakarpowicz/Documents/CohenLab/ACE-master/bin/ace -i ' inputfile ' -o out -g2 0.0002 -b 25000'];
status = system(cmd);
if status == 1
    disp('ACE failed to execute!')
elseif status == 0
    disp('ACE successfully executed.')
end 

%% Run QEE to reproduce original data set 
% ACE gives a .j file of params 
% use .j file of params to reproduce data set (feeding in proper power of
% 10 corresponding to length of original vectors)

addpath(genpath('/Users/briannakarpowicz/Documents/CohenLab/QEE'));
l = log10(trainlength);

status = system(['/Users/briannakarpowicz/Documents/CohenLab/QEE/qee.out -i out -o sim -mcr ' num2str(l)]);
if status == 1
    disp('QEE failed to execute!')
elseif status == 0
    disp('QEE successfully executed.');
end 

%% Compare sim--1 and original 1_pairs 

RMSE = compare_results(file)

