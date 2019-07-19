function [ d, nn_d, Y, l ] = Distance( X, sym_flag )
%% DISTANCE  Compute the Euclidean distances of the directions
% defined by the columns of X (with or without consideration of point
% symmetry)
% 
% Inputs
% X           D x N matrix containing nonzero column vectors;
%             each column will be interpreted in terms of
%             one sampling direction
% sym_flag    [OPTIONAL]  1    --> use symmetrized distance function
%             [ i.e. Dist_ij = acos( abs (X_i . X_j)/(l_i * l_j) ) ]
%             otherwise: Dist_ij = acos( (X_i . X_j)/(l_i * l_j) ) ]
%
% Outputs
% d           Euclidean distance matrix, size N x N
% nn_d        nearest neighbor distance of every point / within every
%             column of d ( nn_d = min( d + 2*eye(N) ) ), size 1 x N
% Y           copy of X with normalized columns
% l           vector containing the length of each column of X
%
% See also RENORMALIZECOLUMNS, SEARCHGAMMAPOU
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% Distance.m
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


%% computes the Euclidean distance in the form of a matrix
if exist('sym_flag','var')
    if (sym_flag ~= 1 && sym_flag ~= 0 )
        error('sym_flag must be either 0 or 1');
    end
else
    sym_flag = 0;
end

N       = size(X,2);
[Y, l]  = RenormalizeColumns(X);
M0      = Y'*Y;
% make sure M0 is symmetric
M0      = 0.5*(M0+M0');
% make sure that no values exceeding +/- 1 exist
M0      = max( min(  M0(:,:) , 1.0 ), -1.0 );
% compute d
d = sqrt(2*(1-M0));

% if symmetrized points are provided, then the output nn_d has
% twice as many entries as number of points requested and
% d is of dimension (2*N)-by-(2*N)
if( sym_flag == 1 )
	d = [ d, sqrt(max(0,4-d.^2)); sqrt(max(0,4-d.^2)), d ];
end

if nargout>=2
    % NOTE: the diagonal must be modified as it contains a zero distance,
    % as far as the nearest neighbor statistics is concerned
    if( sym_flag == 1 )
        nn_d = min( d + 2*eye(2*N) );
    else
        nn_d = min( d + 2*eye(N) );
    end
end
