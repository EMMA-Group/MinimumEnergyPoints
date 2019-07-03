function [ f, g ] = RepulsionEnergy( Xvec, D, energy_index, sym_flag )
%% REPULSIONENEREGY computes for a global vector of inputs X of size D * N
% (e.g. obtained via reshape of an D x N matrix) the overall energy f of
% the point set and [optional] the tangential part of the gradient of f.
% NOTE: This function is eligible for use with fminunc etc. via
% fminunc( @(X) RepulsionEnergy( X, D ), ... );
% where D is defined before.
%
% Inputs
% Xvec       matrix in vector shape of size D*N x 1
% D          dimension of the local data
% sym_flag   [OPTIONAL] 1 if symmetric, 0 if not
%
% Outputs
% f          global energy of the point set (scalar)
% g          [OPTIONAL] gradient of f (same dimension as X) projected onto the
%            plane, n.e. perpendicular to current direction
%
% See also FMINUNC, GEODESICDISTANCE, SEARCHGAMMAPOU, REPULSIONFORCE
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% RepulsionEnergy.m
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

if ~exist('sym_flag','var')
    sym_flag = 0;
end

%% step 1: compute geodesic distance map (matrix Theta, size N x N)
N = size(Xvec,1)/D;
% note: symmetry flag is considered in the energy function only
[ xi, ~, Y, l ] = GeodesicDistance( reshape( Xvec, D, N ) );

%% step 2a: evaluate energy
f = 0.;
for n=1:N

	% account for the symmetry of Theta, i.e. use only entries below diagonal
	f = f + PairPotential(xi((n+1):N,n), energy_index, sym_flag);
end
if( energy_index == 0 )
    if( sym_flag == 0 )
        f = ( 2 * f + 0* pi ) / (N*N); % xlogx specific
%         f = ( 2 * f + N * pi ) / N^2; % xlogx specific
    else
        f = ( 2 * f + N * pi ) / (2 * N^2 ); % xlogx specific
    end
else
    f = 2 * f / N^2; % Riesz s-energies only
end
if isnan(f)
    error('f is NaN')
end

if( nargout == 2 )
    %% step 2b: evaluate gradient
    % the chain rule is required here since Theta(n,j) is acos( Y_i . Y_j ),
    % with Y_k = X_k / || X_k ||
    g = zeros( D, N );
    for n=1:N
        % for each direction the gradient is a weighted
        % sum of all other directions
        G       = - dPairPotential(xi(:,n), energy_index, sym_flag) ./ sin(xi(:,n));
        % the point n does not induce any force on itself:
        G(n)    = 0; 
        y       = Y(:,n);
        g(:,n)  = Y*G;
        % the vector g is now projected into the tangent plane
        % defined by the unit vector y, n.e. the gradient of
        % W w.r.t. direction x_i is perpendicular to x_i by
        % definition
        g(:,n)  = (g(:,n) - y'*g(:,n) * y) /l(n);
    end
    % reshape g --> proper dimension for the use with fminunc etc.
    g = reshape(g, [],1) / N^2;
    if( sym_flag == 1)
        g = g/2;
    end
end
