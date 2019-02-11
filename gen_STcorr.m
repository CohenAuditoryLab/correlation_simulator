function initial_corr = gen_STcorr(c_desired) %Number_of_pairs,
%   Divides numbers of trains by 10 and does sweep from 0-10, 11-20, 21-30
%   etc. all the way up to 100% contamination.  Saves each train into
%   folder for these to be analyzed for ISI violiations

%UPDATE: Should save firing rate too

% if mod(Number_of_pairs,10)~=0
%
%     Number_of_pairs= input('change to number divisible by 10')
%
% end

Train_in_group=1;%Number_of_pairs/10;

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

for iter=1:1 %Number_of_pairs
    
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
    fr(1) = fr(2); %%%% ADDED
    
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
    
    fr_cm=cell(1,4); %maybe overkill but to make sure get correct firing rate for each spike train
    
    fr_cm{1}='A';
    fr_cm{3}='B';
    
    %just some logical noodling to avoid using find and getting the info
    %that I want
    fr_cm{2}=(ind_st_big(1)~=1)+1;
    fr_cm{4}=(numclu(numclu~=(ind_st_big(1)~=1)+1));
    
    %%
    
    % need to rebin the spike times back into 1's and 0's 
    % use time bins of width 1/fs 
    
%     dt = 1/fs; %time step -- too small (makes array too large for matlab)
    dt = 5; % ms 
    bin_edges = 0:dt:double(max(max(st_A), max(st_B))); 
%     binned_A = -1*ones(1, length(bin_edges)-1); % initialize to -1 so can tell when updated with counts 
%     binned_B = -1*ones(1, length(bin_edges)-1);
%     for i = 1:length(bin_edges) - 1
%         lower = bin_edges(i);
%         upper = bin_edges(i+1);
%         binned_A(i) = sum(st_A > lower & st_A <= upper);
%         binned_B(i) = sum(st_B > lower & st_B <= upper);
%     end 
%     

    binned_A = histc(double(st_A), bin_edges);
    binned_B = histc(double(st_B), bin_edges);
    
    binned_A = binned_A >= 1; %turn back into 1, 0
    binned_B = binned_B >= 1;
    
    %% correlation: 1) compute initial correlation 2) add more 3) recalculate correlation to check
    
    % zeroth lag cross correlation - only one time bin
    % want to compute noise correlation or zero lag cross correlation
    
    initial_corr = xcorr(double(binned_A), double(binned_B), 0, 'coeff');
    
    % below is started but not finished nor debugged
%     while initial_corr ~= c_desired
%         if initial_corr < c_desired % add correlation
%             jitter = randn(1, length(st_B) - d);
%             jitter = floor(max(diff(st_B)) * jitter);
%             st_B = st_B + [jitter'; zeros(d)];
%         elseif initial_corr > c_desired % remove correlation
%             continue; % placeholder: TODO
%         else % do nothing
%             continue;
%         end
%         initial_corr = xcorr(double(st_A), double(st_B), 0, 'coeff');
%     end
%     % Save pair
    
end
end