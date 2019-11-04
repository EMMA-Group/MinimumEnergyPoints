% loop through sets of files containing coordinates of ME points and, for
% each such file, compute the optimal Gaussian kernel parameter gamma for
% spherical interpolation of the constant one-function (same kernel width
% at all points). output the result in human-readable table format.

clearvars
close all
addpath(genpath('..'))

% files containing the point coordinates. assumed to contain points on
% hyperspheres, i.e. normalization is not checked. example: if the total
% file names are ~/point_directory/coordinate_file_X.txt with X a varying 
% alphanumerical value, then set the following.
% values
input_files_dir         = '~/point_directory/';
input_files_wildcard    = 'coordinate_file_*.txt';
input_files             = dir([input_files_dir,input_files_wildcard]);
num_files               = length(input_files);
if num_files == 0
    error(['No file found matching ',input_files_dir,input_files_wildcard]);
end

% create necessary containters and parameters
gamma_best  = zeros(num_files,1);
error_min   = zeros(num_files,1);
gamma_min   = 0.18;      % the more points, the larger this should be
gamma_max   = 0.8;      % the more points, the larger this should be
gamma_steps = 200;
gamma_all   = zeros(gamma_steps,1);
file_count  = 0;
N_points    = zeros(num_files,1);

% start the loop
for file = input_files'
    file_count = file_count + 1;
    current_file = [input_files_dir,file.name];
    fprintf(2,['reading ',current_file,'\n'])
    
    % import coordinates from file. here, it is assumed that the user wants
    % to read from .txt files. adapt to reading .mat files when appropriate
    X = dlmread(current_file);
    X = X'; % don't forget to check your case!
    N_points(file_count) = size(X,2);
    
    % do the gamma optimization
    [gamma_best(file_count), error_min(file_count), ~, gamma_all, ~] ...
        = SearchGamma(  X,      ... the support points of the Gaussian SBF
                        1e4,    ... number of Equal Area evaluation points
                        gamma_min, gamma_max, gamma_steps, pi, true );
end

fprintf(1,'    N       gamma \n')
sortrows([N_points,gamma_best],1)