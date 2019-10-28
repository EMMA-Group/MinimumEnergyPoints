function [meshnorm, largestgap_position] = MeshNorm(X, Y, sym_flag)
%% MESHNORM returns the mesh norm of a spherical point set X, where the
% full sphere is approximated by another, very large point set Y. Also,
% optionally, the column of Y that has the maximum nearest neighbor
% distance to X is returned as largestgap_position.
%
% The Euclidean distance function is employed here.
%
% INPUTS
% X            matrix D-by-Nx, normalized columns
% Y            matrix D-by-Ny, normalized columns. Ny >> Nx
% sym_flag     symmetry flag, indicating whether or not to use symmetric
%              distance function
%
% OUTPUTS
% meshnorm     an approximation to the mesh norm
% largestgap_position   an approximation to the position of the largest gap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% MeshNorm.m
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

%% do some consistency checks
if size(X,1) ~= size(Y,1)
    error('size(X,1) must be the same as size(Y,1)');
end

if size(Y,2) <= size(X,2)*3 % TODO insert empirical factor depending on D and Nx
    warning('size(X,2) should be (significantly) smaller than size(Y,2)!');
end

if ~exist('sym_flag','var')
    sym_flag = 0;
elseif sym_flag ~= 0 && sym_flag ~= 1
    error('sym_flag must be either 0 or 1')
end

%% go
[XX]    = RenormalizeColumns(X);
[YY]    = RenormalizeColumns(Y);
M0      = XX'*YY;
% make sure that no values exceeding +/- 1 exist
M0      = max( min(  M0(:,:) , 1.0 ), -1.0 );
% compute d. first index: X, second index: Y
if sym_flag == 0
    d = sqrt(2*(1-M0));
else
    d = min(sqrt(2*(1-M0)),sqrt(2*(1+M0)));
end

% mesh norm definition: sup_Y min_X d
% now the task is: find the column for which the smallest entry is largest,
% and take that entry as meshnorm, and the corresponding column of Y as
% largest gap coordinates
min_d = min(d,[],1);   % row vector with nearest neighbor distances for all test points in Y
[meshnorm, largestgap_index] = max(min_d); % maximum nearest neighbor distance, and the corresponding point index of Y
if nargout > 1
    largestgap_position = YY(:,largestgap_index);
end

end
