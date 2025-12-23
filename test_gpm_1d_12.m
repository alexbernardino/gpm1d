function test_gpm_1d_12


% Author: Nicola Greggio
%         Instituto Superior Tecnico - Lisbon, Portugal
%         nicola.greggio@sssup.it - nicola.greggio@gmail.com
%

% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Init everything
%clear;

% find out which directory this m-file resides in
cd '../'
startupDirectory = pwd;

% add relative paths as absolute paths within the code directory
addpath( genpath( fullfile( startupDirectory, 'gpm1d-code' ) ) );

% back to the original directory
cd 'gpm1d-code'

init_global_variables;
folder_paths;

ouputFolder = '../Output Data/gpm_1d_output_data_12-5/';

num_points = [20, 100, 500, 2500];

weight_bins = [15,20,25];
mean_bins = [15,20,25];
var_bins = [15,20,25];

weight_bins = [30];
mean_bins = [30];
var_bins = [30];


num_reps = length(weight_bins) * length(mean_bins) * length(var_bins);

numRepetitions = 100; 
saveBoxPlot2File = 0;

for c = 1:length(num_points) % all cases of number of points
    counter = 0;
    tic;
    for i = 1:length(weight_bins)
        for j = 1:length(mean_bins)
            for k = 1:length(var_bins)
                outputFileName_base =  'gpm_1d_out_N_%06d_wb_%02d_mb_%02d_vb_%02d';
                outputFileName= sprintf(outputFileName_base, num_points(c), weight_bins(i), mean_bins(j), var_bins(k));
                %outputFileName= sprintf(outputFileName_base, mean_bins, std_bins);
                gpm_1d_12(numRepetitions, num_points(c), weight_bins(i), mean_bins(j), var_bins(k), saveBoxPlot2File, outputFileName, ouputFolder);
                toc
                counter = counter + 1;
                fprintf('numPoints = %d - counter = %d/%d\n', num_points(c), counter, num_reps);
            end
        end
    end
end
return;