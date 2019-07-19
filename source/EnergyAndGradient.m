function [ I, g ] = EnergyAndGradient( Xvec, D, energy_index, sym_flag )
%% ENERGYANDGRADIENT computes for a global vector of inputs X of size D * N
% (e.g. obtained via reshape of an D x N matrix) the overall energy I of
% the point set and [optional] the tangential part of the gradient g of I.
% NOTE: This function is eligible for use with fminunc etc. via
% fminunc( @(X) EnergyAndGradient( X, D ), ... );
% where D is defined before.
%
% Inputs
% Xvec       matrix in vector shape of size D*N x 1
% D          dimension of the local data
% sym_flag   [OPTIONAL] 1 if symmetric, 0 if not
%
% Outputs
% I          global energy of the point set (scalar)
% g          [OPTIONAL] gradient of I (same dimension as X) projected onto
%            the plane, i.e. perpendicular to current direction
%
% See also FMINUNC, DISTANCE, SEARCHGAMMAPOU, REPULSIONFORCE
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% EnergyAndGradient.m
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

%% step 1: compute Euclidean distance map (matrix Dist, size N x N)
N = size(Xvec,1)/D;
% note: symmetry flag is considered in PointPotential, dPoitPotential, and 
% in step 2, but NOT in the following call of Distance
[ d, ~, Y, l ] = Distance( reshape( Xvec, D, N ), 0 );

%% step 2a: evaluate energy
I = 0.;
for n=1:N

	% account for the symmetry of Dist, i.e. use only entries below diagonal
	I = I + PointPotential(d((n+1):N,n), energy_index, sym_flag);
end
if( energy_index == -2 )
    if( sym_flag == 0 )
        % add the diagonal part to the energy
%         I = ( 2 * I + 0 * PointPotential(0, energy_index, 0) ) / N^2;
        I = ( 2 * I + N * PointPotential(0, energy_index, 0) ) / N^2;
    else
        % if sym_flag==1, then I already contains the kernel values w.r.t.
        % antipodes, i.e. twice as many terms in the sum: N^2 -> 2*N^2
        % further, the diagonal of the antipodal part is zero -> only N
        % times the diagonal (i.e. N*2), as in the asymmetric case
        I = ( 2 * I + N * PointPotential(0, energy_index, 0) ) / (2 * N^2 );
    end
else % classical case (log, Riesz). ATTENTION: different normalization!
    if( sym_flag == 0 )
        I = 2 * I / (N * (N-1));
    else
        % if sym_flag==1, then consider also the diagonal of the antipodal
        % part, i.e. sum of k(x_i,-x_i) = N * k(2)
        I = (2 * I + N * PointPotential(2, energy_index, 0)) ...
            / (2 * N * (N-1) + N);
    end
end
if isnan(I)
    error('I is NaN')
end

if( nargout == 2 )
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
end
