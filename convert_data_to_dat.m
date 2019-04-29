function convert_data_to_dat(file, fraction)
%% File is the name of the .mat file produced from Gen_STcorr_v3
    pieces = strsplit(file, '.');
    savename = pieces{1};

    load(file);

    data = [];
    smallest = min(cellfun('size',spikes,1));
    smallest = round(smallest * fraction);

%     for i = 2:size(pairs,1);
%         neuron1 = pairs{i, 1};
%         neuron2 = pairs{i, 2};
%         data(end + 1, :) = neuron1(1:100000);
%         data(end + 1, :) = neuron2(1:100000);
%     end 

    for i = 1:size(spikes,1)
        neuron = spikes{i};
        data(:,i) = neuron(1:smallest);
    end 

    disp(size(data));

    dlmwrite([savename '.dat'], data);
end 