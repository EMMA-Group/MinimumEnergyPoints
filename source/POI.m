function [POI_RMSE, POI_normdiff] = POI(Xi_training, Xi_training_vali, ...
    gamma, X_training, X_vali)
% POI  compute the error of partition of identity (POI) from provided
% geodesic distances and gamma.
% The measure of error can be adjusted within this function.
% 
% Inputs
% Xi_training       matrix of geodesic distances between training
%                   directions, size N-by-N
% Xi_training_vali  matrix of geodesic distances between training
%                   directions and validation directions, size N-by-Nvali
% gamma             kernel width parameter
% X_training        training directions, size D-by-N (D = dimension)
% X_vali            validation directions, size D-by-Nvali
%
% Outputs
% POI_RMSE          RMSE of POI for all Nvali validation directions
% POI_normdiff      pointwise frobenius norms of differences of validation
%                   directions and interpolated directions
%
% See also SEARCHGAMMA, POU

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% POI.m
% Copyright (C) 2019 by Felix Fritzen and Oliver Kunc
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
% JOURNAL NAME, Number/Volume, p. XX-YY, 2019
% DOI   ...
% URL   dx.doi.org/...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% consistency of inputs
if size(Xi_training_vali,2) ~= size(X_vali,2) || ...
   size(Xi_training,1) ~= size(Xi_training,2) || ...
   size(Xi_training,1) ~= size(Xi_training_vali,1)
    error('Dimension mismatch of input arguments.')
end
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

%% interpolate validation directions
% this should be X_vali in order to satisfy
% the partition of unity property (POI)
X_interpolated = X_training * w;

%% compute error of partition
% rooted mean square error w.r.t. identity on set of validation directions
temp = sum( (X_interpolated - X_vali).^2 );
POI_normdiff = sqrt( temp );
Nvali = size(Xi_training_vali,2);
POI_RMSE = sqrt( sum(temp)/Nvali );

%% some handy plotting routine. opens new figure window at every call!
% figure('Name', ['POI: diff. to identity, gamma = ', num2str(gamma)], 'NumberTitle', 'off'); 
% plot(sort(POU_normdiff,'ascend'));
% xlabel('case of evaluation direction')
% ylabel('POI_normdiff between vali directions and their interplants')
% grid on
