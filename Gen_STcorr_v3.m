function Gen_STcorr_v3(desired_c, Number_of_pairs)
%GEN_STcorr_v2 Generates inputed number of spike train pairs with desired
%pairwise cross-correlation.
%LAST EDITED: B. Karpowicz 3/18/19

% 3/18/19: if desired correlation cannot be reached within tolerance after 100
% iterations, return the maximum correlation reached. Fixed correlation
% bug.

% 3/14/19: rather than producing correlations between pairs,
% produce them between all neurons. Requires spike trains to be the same
% length, binned to the latest time of all units.

% INPUTS
%   desired_c : a vector of desired pairwise correlations for each pair of
%   neurons, such that desired_c(i) is the desired correlation for the ith
%   pair. If desired_c is not of length Number_of_pairs, the vector is
%   repeated using repmat until it reaches this length. Recommended
%   correlations ~0.1-0.2.
%
%   Number_of_pairs : an integer number of the desired pairs of neurons to
%   be generated

if length(desired_c) ~= Number_of_pairs
    num_reps = ceil(Number_of_pairs/length(desired_c));
    desired_c = repmat(desired_c, 1, num_reps);
    desired_c = desired_c(1:Number_of_pairs);
end

Train_in_group = 1;

contam_vec=[linspace(0,.1,Train_in_group)...
    linspace(.11,.20,Train_in_group)...
    linspace(.21,.30,Train_in_group)...
    linspace(.31,.40,Train_in_group)...
    linspace(.41,.50,Train_in_group)...
    linspace(.51,.60,Train_in_group)...
    linspace(.61,.70,Train_in_group)...
    linspace(.71,.80,Train_in_group)...
    linspace(.81,.90,Train_in_group)...
    linspace(.91,1,Train_in_group)]';

dat = load('simulation_parameters');
fs  = dat.fs; % sampling rate
wav = dat.waves; % mean waveforms for all neurons

%% Start Loop

pairs = cell(Number_of_pairs, 4);
pairs{1,1} = 'A';
pairs{1,2} = 'B';
pairs{1,3} = 'actual corr';
pairs{1,4} = 'desired corr';

spikes = cell(Number_of_pairs*2, 1);

for iter=1:Number_of_pairs
    
    ID=iter
    
    %% Parameters
    %generates contaminate spike train pairs, currently only works with single
    % pair of spike trains but could adjust further.
    %Some parameters pulled from eMouse
    
    nn        = 2; % number of simulated neurons (keep at 2 for now)
    fr_bounds = [.1 10]; % min and max of firing rates % reset to .1 to reflect lower average spiking
    t_record  = 1000; % duration in seconds of simulation. longer is better (and slower!) (1000)
    
    numclu=1:nn; % only will be 1 to 2 for now
    per2con = contam_vec(iter); % [0 1] percent of overlap (contamination) between two spikes; calculated as percent of smaller st
    %i.e.  A I B/ B where B is subset of A, 100% means perfect
    %subset; naturally since adding spikes without removing others
    %increase spike total, max percent is ratio of
    %numel(A) to numel(A I B)
    pull_reverse = 0;  %Normally take spikes from larger spike train to smaller spike train.  Set to 1 to reverse.
    %Makes max percent contamination percent of spikes B is
    %of A U B (e.g. numel(A)=100 numel(B)=50  means max is
    %50/150
    
    
    %% simulator pulled from eMouse - define st_A and st_B
    % UPDATE: 12/5/17 adding common drive parameter and thus changed for loop
    % for common drive
    
    wav = wav(:,:, randperm(size(wav,3), nn)); %select random wave forms for number of neurons chosen %step means must be across channels for same neuron
    
    
    NN = size(wav,3); % number of neurons %should be same as nn...curious why they did this again
    fr = fr_bounds(1) + (fr_bounds(2)-fr_bounds(1)) * rand(NN,1); % create variability in firing rates % uses uniform scale rather than shift 0 to diff between bounds
    
    spk_times = [];
    clu = [];
    
    for j = 1:length(fr) %% iterates through number of neurons; each is assigned a firing rate based on previous step
        
        dspks = int64(geornd(1/(fs/fr(j)), ceil(2*fr(j)*t_record),1)); %%int64 turns into signed 64 bit integer, geornd is pull from geometric distribution with p = input,
        % MATLAB geopdf is failures before 1 success.  above is p based on
        % sampling and firing rate and then making a column vector double the
        % length of the recording times the firing rate
        dspks(dspks<ceil(fs * 3/1000)) = [];  % remove ISIs below the refractory period %consequence: cant have below refract between neurons as well...
        res = cumsum(dspks); %%gives running calc of spike intervals in units of samples (e.g. 3000 samples to first 500 to second means 3500 from onset to second which is in cumsum)
        spk_times = cat(1, spk_times, res); %sets up spike times as these cumulative intervals (make sense essentially gives running time in samples)
        clu = cat(1, clu, j*ones(numel(res), 1)); %% %prevents need for matrix, just use cluster to index spk_times vectors vector as long as spike times of a single 1 and then all elements equal to number of neurons indicated
        
    end
    
    [spk_times, isort] = sort(spk_times);
    clu = clu(isort); %assigns sorted spikes randomly to cluster since cluster identity matched with UNsorted spikes.
    clu       = clu(spk_times<t_record*fs); %remove clusters ids that coincide with being out of the session
    spk_times = spk_times(spk_times<t_record*fs); %remove out of session spikes, hence above 2 * firing rate thing is arbitrary
    
    %we now should have two spike trains.  spk_times is a concatenation of all
    %spike times sorted (in units of sample number).  clu vector contains
    %maping for each value in spk_times
    
    %Define spike trains as sets where numel(A) > numel(B)
    ind_st_big=abs((numel(spk_times(clu==1))>=numel(spk_times(clu==2)))-1)+1;%not sure if faster than if or find, but maybe?  Is arithmatic logic gate
    
    st_A=spk_times(clu==ind_st_big); %A is always larger
    st_B=spk_times(clu~=ind_st_big); %B is train that will receive spikes from A hence is always subset
    
    %IS THIS NECESSARY? makes the initial correlations way higher
    %     st_A = st_A/fs * 1e3; %converting from sample # to ms
    %     st_B = st_B/fs * 1e3; %converting from sample # to ms
    
    fr_cm=cell(1,4); %maybe overkill but to make sure get correct firing rate for each spike train
    
    fr_cm{1}='A';
    fr_cm{3}='B';
    
    %just some logical noodling to avoid using find and getting the info
    %that I want
    fr_cm{2}=(ind_st_big(1)~=1)+1;
    fr_cm{4}=(numclu(numclu~=(ind_st_big(1)~=1)+1));
    
    %store spike times
    spikes{2*ID-1} = st_A;
    spikes{2*ID} = st_B;
end

%% need to rebin the spike times back into 1's and 0's

dt = 20; % ms (accepted ACE bin size)

curr_max = -1000;
for i = 1:length(spikes)
    m = max(spikes{i});
    if m > curr_max
        curr_max = m;
    end
end

bin_edges = 0:dt:double(curr_max);

for j = 1:Number_of_pairs
    
    st_A = spikes{2*j-1};
    st_B = spikes{2*j};
    
    binned_A = histc(double(st_A), bin_edges) >= 1;
    binned_B = histc(double(st_B), bin_edges) >= 1;
    
    %% correlation: 1) see overlap 2) copy over spikes 3) recalculate overlap to check
    
    %check if the two spike trains are overlaping in time at all...obviously if they don't kind of meh
    
    checking1=min(st_A)<max(st_B)&&min(st_B)<max(st_A);
    
    if checking1==0
        
        maxb_mina=max(st_B)-min(st_A);
        maxa_minb=max(st_A)-min(st_B);
        
    end
    
    initial_corr = xcorr(double(binned_A), double(binned_B), 0, 'coeff');
    
    %% copy over spikes
    
    if pull_reverse==1 % check if want to add spikes from smaller to larger set rather than larger to smaller
        
    else %standard protocol
        
        tolerance = 0.001;
        num_iters = 0;
        
        while initial_corr < desired_c(j) - tolerance & num_iters < 100
            
            disp(['Correlation: ' num2str(initial_corr)]);
            
            n2add=floor(abs((desired_c(ID)-initial_corr))*numel(st_B));% subtract out percent already overlap...usually so small this won' t matter but just in case
            % How many spikes are we adding over
            
            act_per=n2add/numel(st_B); %store actual percentage since rounding
            
            %uniform pull of spike times from larger set.  Can adjust to be only from certain fractions if you want contamination more concentrated
            pulls_ind=randperm(numel(st_A), n2add)'; %get indices to pull from st_A
            
            pulls=st_A(pulls_ind);%get values at those indices
            %the trick - keeping cluster assignment after sorting new spike train
            %Below maybe able to be implemented better...right now essentially
            %reclustering
            
            st_B=[st_B; pulls]; % add spikes from A to B
            new_spk_times=[st_A ; st_B]; % concat spikes unsorted in time but split by cluster
            
            new_clu=[ind_st_big*ones(numel(st_A),1);...
                numclu(numclu~=ind_st_big)*ones(numel(st_B),1)]; %redo clustering so that added spikes go to cluster B by first just assigning halfs
            
            %sort spike train and then redo cluster vector based on sort
            [new_spk_times, ind2match]=sort(new_spk_times);
            new_clu=new_clu(ind2match);
            
            [st_B, ind_as]=sort(st_B); %save original indicies of added spikes so that can be pulled out in visualize
            
            per_over_norm=sum(intersect(st_A,st_B)>0)/(numel(st_B)-numel(pulls)); %percent of st_B that are from A compared to original number of spikes in B
            
            per_over_total=sum(intersect(st_A,st_B)>0)/(numel(st_B)); % percent of st_B that are from A compared to new total spikes in st_B
            
            binned_A = histc(double(st_A), bin_edges) >= 1;
            binned_B = histc(double(st_B), bin_edges) >= 1;
            initial_corr = xcorr(double(binned_A), double(binned_B), 0, 'coeff');
            
            num_iters = num_iters + 1;
            
        end
        
        % update stored spike trains for later correlation computations 
        spikes{2*j-1} = st_A;
        spikes{2*j} = st_B;
    end
    
    disp(['Stopping Correlation: ' num2str(initial_corr)]);
    pairs{j+1, 1} = binned_A;
    pairs{j+1, 2} = binned_B;
    pairs{j+1, 3} = initial_corr;
    pairs{j+1, 4} = desired_c(j);
    
end

%% compute upper triangular correlations matrix for all pairs (since all same length now)

corrs = zeros(Number_of_pairs);

for i = 1:Number_of_pairs*2
    for j = i:Number_of_pairs*2
        st_A = spikes{i}; 
        st_B = spikes{j};
        binned_A = histc(double(st_A), bin_edges) >= 1;
        binned_B = histc(double(st_B), bin_edges) >= 1;
        corrs(i,j) = xcorr(double(binned_A), double(binned_B), 0, 'coeff');
    end 
end 
        
%%

directory = '/Users/briannakarpowicz/Documents/CohenLab/contamination_and_simulations_2017and18-master/Critical_Functions_and_Files/';

save([directory num2str(ID) '_pairs'], 'pairs', 'corrs');

end

