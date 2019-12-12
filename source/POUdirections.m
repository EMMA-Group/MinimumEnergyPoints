function [POU_err diff_one] = POUpoints(X_training, X_vali, gamma)
% POU  compute the error of partition of unity (POU) from provided
% points and gamma, calling POU.m
% The measure of error can be adjusted within this function.
%
% Inputs
% X_training    coordinate matrix of training points, size D-by-N
% X_vali        coordinate matrix of validation points, size D-by-Nvali
% gamma         kernel width parameter
%
% Outputs
% POU_err       RMSE of POU for all Nvali validation points
% diff_one      pointwise difference of sum of weights and one
%
% See also EQ_POINT_SET, ENERGY, DISTANCE,
% RENORMALIZECOLUMNS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% POUpoints.m
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
% Oliver Kunc and Felix Fritzen: 'Generation of energy-minimizing point sets on
% spheres and their application in mesh-free interpolation and differentiation',
% Advances in Computational Mathematics, Number/Volume, p. XX-YY, 2019
% DOI   10.1007/s10444-019-09726-5
% URL   doi.org/10.1007/s10444-019-09726-5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% compute the distance matrices
Xi_training = Distance(X_training);
Xi_training_vali = real(acos(min(1,max(-1,X_training'*X_vali))));

%% call the main POU function
[POU_err diff_one] = POU(Xi_training, Xi_training_vali, gamma);
