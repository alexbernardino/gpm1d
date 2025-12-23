function test_gpm_1d_11


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





outputFolder = [OUTPUT_DATA_Fol 'gpm_1d_output_data_11/'];


mean_bins = [10,15,20,25];
std_bins = [3,6,9,12,15,18,21];
num_reps = length(mean_bins) * length(std_bins);

saveBoxPlot2File = 0;


% case 1
numRepetitions = 100; 
numPoints = 100;
exec_gpm1d(numRepetitions, numPoints, mean_bins, std_bins, num_reps, outputFolder, saveBoxPlot2File);

% case 2
numRepetitions = 100; 
numPoints = 1000; 
exec_gpm1d(numRepetitions, numPoints, mean_bins, std_bins, num_reps, outputFolder, saveBoxPlot2File);

% case 3
numRepetitions = 100; 
numPoints = 10000;
exec_gpm1d(numRepetitions, numPoints, mean_bins, std_bins, num_reps, outputFolder, saveBoxPlot2File);

% case 4
numRepetitions = 100; 
numPoints = 100000; 
exec_gpm1d(numRepetitions, numPoints, mean_bins, std_bins, num_reps, outputFolder, saveBoxPlot2File);

return;




function exec_gpm1d(numRepetitions, numPoints, mean_bins, std_bins, num_reps, outputFolder, saveBoxPlot2File)

    counter = 0;
    tic;
    for i = 1:length(mean_bins)
        for j = 1:length(std_bins)
            outputFileName_base =  'gpm_1d_out_N_%06d_mb_%02d_sb_%02d';
            outputFileName = sprintf(outputFileName_base, numPoints, mean_bins(i), std_bins(j));
            gpm_1d_11(numRepetitions, numPoints, mean_bins(i), std_bins(j), saveBoxPlot2File, outputFileName, outputFolder);
            toc
            counter = counter + 1;
            fprintf('numPoints = %d - counter = %d/%d\n', numPoints, counter, num_reps);
        end
    end

return;


