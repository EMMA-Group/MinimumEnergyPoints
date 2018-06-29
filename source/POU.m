function [POU_RMSE, POU_diff] = POU(Xi_training, Xi_training_vali, gamma)
% POU  compute the error of partition of unity (POU) from provided geodesic
% distances and gamma.
% The measure of error can be adjusted within this function.
% 
% Inputs
% Xi_training       matrix of distances between training directions, size
%                   N-by-N
% Xi_training_vali  matrix of distances between training directions and
%                   validation directions, size N-by-Nvali
% gamma             kernel width parameter
%
% Outputs
% POU_RMSE          RMSE of POU for all Nvali validation directions
% POU_diff          pointwise difference of sum of weights and one
%
% See also SEARCHGAMMA, POI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% POU.m
% Copyright (C) 2018 by Felix Fritzen and Oliver Kunc
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% (the full license is distributed together with the software
% in a file name LICENSE)
%
%
% This program employs a modified version of the softwares
%
% 1) Recursive Zonal Equal Area (EQ) Sphere Partitioning Toolbox.
%    Release 1.10 2005-06-26
%
%    written by Paul Leopardi for the University of New South Wales.
% 
%    See COPYING in the subfolder eq_sphere_partitions for
%    licensing information regarding this software.
% 
%    See CHANGELOG in the subfolder eq_sphere_partitions for
%    a concise list of changes that were made to the original code.
%
% 2) UIGETVAR_PUTVAR
%    Release 1.0 2010-02-19
%
%    written by John D'Errico
%
%    See uigetvar.m in the subfolder uigetvar_putvar for more information.
%
%    See license.txt in the subfolder uigetvar_putvar for licensing
%    information regarding this software.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This software package is related to the research article
%
% Oliver Kunc and Felix Fritzen: ''
% JOURNAL NAME, Number/Volume, p. XX-YY, 2018
% DOI   ...
% URL   dx.doi.org/...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute the kernel matrix K
% size of K is N-by-N
K = exp( - gamma * Xi_training.^2 );

%% compute the kernel function zeta for all validation directions
% size of k is N-by-Nvali
zeta = exp( -gamma * Xi_training_vali .^ 2 );
     
%% determine weights w
% size of w is N-by-Nvali
% the i-th column of w contains the N weights
% of the i-th validation case
w = K \ zeta;

%% sum W of all weights
% for each validation case, i.e. sum the columns of w.
% this should be 1 in order to satisfy
% the partition of unity property (POU)
W = sum(w,1);

%% compute error of POU
% rooted mean square error w.r.t. 1-function
Nvali = size(Xi_training_vali,2);
POU_diff = W -  ones(1,Nvali);
POU_RMSE = norm( POU_diff ) / sqrt(Nvali);

%% some handy plotting routine. opens new figure window at every call!
% figure('Name', ['POU: diff. to one, gamma = ', num2str(gamma)], 'NumberTitle', 'off'); 
% plot(sort(POU_diff,'ascend'));
% xlabel('case of evaluation direction')
% ylabel('difference of interpolation to one')
% grid on

