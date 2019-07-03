function [ xi, nn_xi, Y, l ] = GeodesicDistance( X, sym_flag )
%% GEODESICDISTANCE  Compute the geodesic distances of the directions
% defined by the columns of X (with or without consideration of point
% symmetry)
% 
% Inputs
% X           D x N matrix containing nonzero column vectors;
%             each column will be interpreted in terms of
%             one sampling direction
% sym_flag    [OPTIONAL]  1    --> use symmetrized distance function
%             [ i.e. Theta_ij = acos( abs (X_i . X_j)/(l_i * l_j) ) ]
%             otherwise: Theta_ij = acos( (X_i . X_j)/(l_i * l_j) ) ]
%
% Outputs
% xi          geodesic (i.e. angular) distance matrix xi_ij [rad]
% nn_xi       nearest neighbor distance of every point / within every
%             column of xi ( nn_xi = min( xi + pi*eye(N) ) ), size 1 x N
% Y           copy of X with normalized columns
% l           vector containing the length of each column of X
%
% See also RENORMALIZECOLUMNS, SEARCHGAMMAPOU
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% GeodesicDistance.m
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


%% computes the geodesic distance in the form of a matrix
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
% put ones on the diagonal
M0      = M0 - diag(diag(M0)) + eye(N);
%  if( sym_flag  == 1 )
%      % if symmetry conditions apply, then the absolute value must be taken
%      % NOTE: then the maximum output is pi/2, independent of the space
%      % dimension D!
%      M0  = abs(M0);
%  end
% prevent spurious complex numbers from popping up
xi   = real( acos(M0) );

% if symmetrized points are provided, then the output nn_xi has
% two times as many entries as number of points requested and
% xi is of dimension (2*N)-by-(2*N)
if( sym_flag == 1 )
	xi = [ xi, pi-xi; pi-xi, xi ];
end

if nargout>=2
    % NOTE: the diagonal must not be taken here as it contains a zero distance...
    if( sym_flag == 1 )
        nn_xi = min( xi + pi*eye(2*N) );
    else
        nn_xi = min( xi + pi*eye(N) );
    end
end
