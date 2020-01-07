function [NN_MinMeanMax, meshnorm, meshratio, largestgap] = NN_and_mesh_statistics(X,sym_flag,Ny,Y)
%% NN_AND_MESH_STATISTICS compute nearest neighbor (NN) and mesh statistics
% of a spherical point set X. Needs very large point set Y for the
% mesh norm, e.g. size(Y,2)=1e5 is usually adequate. The
% largestgap_position is chosen among the columns of Y.
%
% The Euclidean distance function is employed atm. This function is modular
% w.r.t. the distance function.
%
% Inputs
% X          points ( D-by-N matrix )
% sym_flag   symmetry flag (i.e. consider point set X and -X if sym_flag == 1)
% Ny         [OPTIONAL] number of test points
% Y          [OPTIONAL] test points ( D-by-Ny matrix ). if not provided,
%            then will be initialized randomly.
%
% Outputs
% NN_MinMeanMax
% meshnorm          mesh norm (i.e. circumference at largestgap_position)
% meshratio         mesh ratio (2 * ratio of meshnorm and smallest nearest
%                   neighbor distance)
% largestgap_position        point coordinate within Y leading to the
%                   largest nearest neighbor distance w.r.t. X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% NN_and_mesh_statistics.m
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
% Advances in Computational Mathematics 45(5-6), pp. 3021-3056, 2019
% URL: doi.org/10.1007/s10444-019-09726-5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preparation
if ~exist('Ny','var')
    Ny  = max( 10 * size(X,2), 1e5 );  % the "full" sphere, should be very fine
elseif Ny <= 0
    Ny  = max( 10 * size(X,2), 1e5 );
end
if exist('Y','var') && size(Y,2)~=Ny
    error(['dimension mismatch: Ny = ',num2str(Ny),' size(Y) = ',num2str(size(Y))]);
end

if ~exist('Y','var')
    Y                    = RandomPoints(size(X,1),Ny);
end

%% mesh statistics
[~, NN]          = Distance(X,sym_flag);
NN_min           = min (NN);
NN_mean          = mean(NN);
NN_max           = max (NN);
NN_MinMeanMax    = [NN_min NN_mean NN_max];

%% mesh ratio, norm, and largest gap position
% mesh ratio = mesh norm / NN_min * 2
[meshnorm, largestgap]   = MeshNorm(X,Y,sym_flag);
meshratio                = 2.0 * meshnorm / NN_min;
