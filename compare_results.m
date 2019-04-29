function RMSE = compare_results(original_file)

original = dlmread(original_file);
reproduced = dlmread('sim--1.dat');

RMSE = sqrt(mean((original-reproduced).^2));

end