function [ xi, xi_nn, Y, l ] = DistanceGeodesic( X, sym_flag )
%% DISTANCEGEODESIC  Compute the geodesic distances of the spherical points
% defined by the columns of X (with or without consideration of symmetry
% w.r.t. the origin). Operates on normalized point coordinate vectors.
%
% Inputs
% X           D x N matrix containing the point coordinates as columns
% sym_flag    [OPTIONAL] if 1, then xi_nn is of size 2*N x 2*N. The upper
%             left and lower right N x N blocks each contain the distances
%             of the explicitly provided points X, the off-diagonal blocks
%             contain the distances of X to their antipodes -X.
%
% Outputs
% xi          geodesic distance matrix. Size N x N if sym_flag==0 or not
%             provided, size 2*N x 2*N if sym_flag==1.
% xi_nn       nearest neighbor (=nn) distance of every point / within every
%             column of d. Size 1 x N if sym_flag==0 or not provided, size
%             1 x 2*N if sym_flag==1.
% Y           copy of X with normalized columns
% l           vector containing the length of each column of X
%
% See also RENORMALIZECOLUMNS, DISTANCE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% Distance.m
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


%% computes the geodesic distance in the form of a matrix
if exist('sym_flag','var')
    if (sym_flag ~= 1 && sym_flag ~= 0 )
        error('sym_flag must be either 0 or 1');
    end
else
    sym_flag = 0;
end

N       = size(X,2);
% normalize
[Y, l]  = RenormalizeColumns(X);
% pairwise scalar products
M0      = Y'*Y;
% make sure M0 is symmetric
M0      = 0.5*(M0+M0');
% make sure that no values exceeding +/- 1 exist
M0      = max( min(  M0(:,:) , 1.0 ), -1.0 );
% put ones on the diagonal
M0      = M0 - diag(diag(M0)) + eye(N);
% prevent spurious complex numbers
xi   = real( acos(M0) );

% if sym_flag==1, compute also the distances to the antipodes of X. Return
% the distance matrix d of the point coordinate matrix [X,-X]. Then d is of
% size 2*N x 2*N.
if( sym_flag == 1 )
	xi = [ xi, pi-xi; pi-xi, xi ];
end

if nargout>=2
    % NOTE: the diagonal must not be taken here as it contains a zero distance...
    if( sym_flag == 1 )
        xi_nn = min( xi + pi*eye(2*N) );
    else
        xi_nn = min( xi + pi*eye(N) );
    end
end
