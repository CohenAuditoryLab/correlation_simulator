function RMSE_all = compare_results(original_file)
%%
% ORIGINAL_FILE: the name of the .mat file containing the originally
% generated neurons.

% RMSE_all: contains the RMSE for each neuron compared to the newly
% generated data from QEE.

%%

original = load(original_file);
original = original.pairs;
reproduced = dlmread('sim--1.dat');

original_data = [];
for i = 2:size(original,1)
    neuron1 = original{i, 1};
    neuron2 = original{i, 2};
    original_data(end + 1, :) = neuron1;
    original_data(end + 1, :) = neuron2;
end

RMSE_all = zeros(size(pairs,2),1);
for j = 1:size(pairs,2)
    o = original_data(j, :);
    r = reproduced(j, :);
    RMSE_all(j) = sqrt(mean((o-r).^2));
end 

end