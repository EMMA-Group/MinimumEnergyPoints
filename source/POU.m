function [POU_L2error, POU_diff] = POU(Xi_training, Xi_training_vali, gamma)
% POU  compute the error of partition of unity (POU) from provided geodesic
% distances and gamma.
% The measure of error can be adjusted within this function.
%
% Inputs
% Xi_training       matrix of geodesic distances between training
%                   points, size N-by-N
% Xi_training_vali  matrix of geodesic distances between training
%                   points and validation points, size N-by-Nvali
% gamma             kernel width parameter
%
% Outputs
% POU_RMSE          RMSE of POU for all Nvali validation points
% POU_diff          pointwise difference of sum of weights and one
%
% See also SEARCHGAMMA, POI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% POU.m
% Copyright (C) 2018, Felix Fritzen and Oliver Kunc
% All rights reserved.
%
% This source code is licensed under the BSD 3-Clause License found in the
% LICENSE file in the root directory of this source tree.
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

%% compute the kernel matrix K
% size of K is N-by-N
K = exp( - gamma * Xi_training.^2 );
% fprintf(2,'WARNING: symmetry is hard-wired!\n')
% K = K + exp( - gamma * (pi-Xi_training).^2 );

%% compute the kernel function zeta for all validation points
% size of k is N-by-Nvali
zeta = exp( -gamma * Xi_training_vali .^ 2 );
% fprintf(2,'WARNING: symmetry is hard-wired!\n')
% zeta = zeta + exp( -gamma * (pi-Xi_training_vali) .^ 2 );

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
POU_L2error = norm( POU_diff ) / sqrt(Nvali); % norm() defaults to the 2-norm

%% some handy plotting routine. opens new figure window at every call!
% figure('Name', ['POU: diff. to one, gamma = ', num2str(gamma)], 'NumberTitle', 'off');
% plot(sort(POU_diff,'ascend'));
% xlabel('case of evaluation point')
% ylabel('difference of interpolation to one')
% grid on

