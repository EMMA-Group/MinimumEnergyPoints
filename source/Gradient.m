function [ g ] = Gradient( X, D, energy_index, sym_flag )
%% GRADIENT computes the tangential part g of the gradient of the
% overall energy I, i.e. returns only the second argument of
% EnergyAndGradient
%
% Inputs
% X         matrix of size D*N x 1 (i.e. a vector)
% D         dimension of the local data
% sym_flag  0: use standard point set; 1: consider X and -X
%
% Outputs
% g         tangential part of gradient of I (same dimension as X)
%
% See also FMINUNC, DISTANCE, SEARCHGAMMAPOU, ENERGY
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% Gradient.m
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

%% step 1: compute Euclidean distance map (matrix Xi, size N x N)
N = size(X,1)/D;
% note: symmetry flag is considered in the energy function only
[ d, ~, Y, l ] = Distance( reshape( X, D, N ), 0 );


%% step 2b: evaluate gradient
% the chain rule is required here: derivative of I w.r.t. n-th
% direction is
%     1/N^2 * sum_(j=1,..,N) k'(d(n,j)) / d(n,j) * ( Y(:,n) - Y(:,j) )
% implemented in a vectorized manner
g = zeros( D, N );
for n=1:N
    % column vecotor of the kernel derivatives
    dPP     = dPointPotential(d(:,n), energy_index, sym_flag);
    % column vector containing the quotions of kernel derivatives and
    % corresponding distances d
    dPP_d   = dPP ./ d(:,n);
    % exclude undefined "self-gradient". physically: no force of a
    % point on itself
    dPP_d(n)= 0;
    % first part of the gradient ( Y(:,n)-part )
    g(:,n)  = sum(dPP_d) * Y(:,n);
    % second part of the gradient ( Y(:,j)-part )
    g(:,n)  = g(:,n) - Y * dPP_d;
    % the vector g is now projected into the tangent plane
    % defined by the unit vector y, n.e. the gradient of
    % I w.r.t. the direction x_i is perpendicular to x_i by
    % definition
    g(:,n)  = (g(:,n) - Y(:,n)'*g(:,n) * Y(:,n)) /l(n);
end
% reshape g --> proper dimension for the use with fminunc etc.
% and apply appropriate factor 1/N^2
g = reshape(g, [], 1) / N^2;
if( sym_flag == 1)
    g = g/2;
end
