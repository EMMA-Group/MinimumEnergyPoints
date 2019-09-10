function [ I, g ] = EnergyAndGradient( Xvec, D, energy_index, sym_flag )
%% ENERGYANDGRADIENT computes the s-energy I of a spherical point set.
% Optionally, the tangential part g of the gradient of I is computed, just
% as in the function Gradient. Also optionally, symmetrized kernels are
% employed (if sym_flag==1).
%
% "s" from the paper is energy_index here. The input vector Xvec contains
% the columns of the coordinate matrix X stacked on top of each other (e.g.
% via reshape). 
%
% For the Riesz cases and the log case, the normalization factor is
% N*(N-1). For the LOG case, the normalization factor N^2 is used (see
% paper Section 2.2).
%
% NOTE: This function is eligible for use with fminunc etc. via
% fminunc( @(x) EnergyAndGradient( x, D ), ... );
% This is the reason why it operates on Xvec and not on X.
%
% The gradient part of this function is identical to that of the function
% Gradient. Only reason for the inclusion of that code copy here is that,
% if known in advance, the Distance function needs to be called only once.
%
% Inputs
% Xvec          column vectors of point coordinates stacked on top of each
%               other, size D*N x 1
% D             dimension of the Euclidean space in which the points are
%               defined
% energy_index  index of the energy function to be used ("s" in the paper)
% sym_flag      [OPTIONAL] if this is 1, then symmetrized kernel functions
%               are employed. if this is 0 or omitted, then the kernel
%               functions are not modified.
%
% Outputs
% I             global energy of the point set (scalar)
% g             [OPTIONAL] gradient of I at each point. It is projected
%               onto the tangential plane of the sphere at each point and
%               is of size D x N.
%
% See also FMINUNC, DISTANCE, MINIMIZEENERGY, GRADIENT
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

%% preparations
if ~exist('sym_flag','var')
    sym_flag = 0;
end

N = size(Xvec,1)/D;

% Note: symmetry flag makes PointPotential and dPoitPotential employ
% symmetrized kernels. The following call of Distance does not require the
% symmetry flag for this purpose (in fact it MUST be called w/o symmetry).
[ d, ~, Y, l ] = Distance( reshape( Xvec, D, N ), 0 );

%% step 2a: evaluate energy
I = 0.;
for n=1:N
	% account for the symmetry of Dist, i.e. use only entries below diagonal
	I = I + PointPotential(d((n+1):N,n), energy_index, sym_flag);
end
if( energy_index == -2 )
    if( sym_flag == 0 )
        % add the diagonal part to the energy, i.e. "self-energy"
%         I = ( 2 * I + 0 * PointPotential(0, energy_index, 0) ) / N^2;
        I = ( 2 * I + N * PointPotential(0, energy_index, 0) ) / N^2;
    else
        % if sym_flag==1, then I already contains the kernel values w.r.t.
        % antipodes, i.e. twice as many terms in the sum: N^2 -> 2 * N^2
        % further, the diagonal of the antipodal part is zero -> only N
        % times the diagonal (i.e. N*2), as in the asymmetric case
        I = ( 2 * I + N * PointPotential(0, energy_index, 0) ) / (2 * N^2 );
    end
elseif( energy_index == -1 || energy_index > 0 ) % classical case (log, Riesz). ATTENTION: different normalization!
    if( sym_flag == 0 )
        I = 2 * I / (N * (N-1));
    else
        % if sym_flag==1, then consider also the diagonal of the antipodal
        % part, i.e. sum of k(x_i,-x_i) = N * k(2)
        I = (2 * I + N * PointPotential(2, energy_index, 0)) ...
            / (2 * N * (N-1) + N);
    end
else
    error(['EnergyAndGradient not yet implemented for energy_index = ',...
        num2str(energy_index)])
end
if isnan(I)
    error('I is NaN')
end

if( nargout == 2 )
    %% step 2b: evaluate gradient
    % the chain rule is required here: derivative of I w.r.t. coordinates
    % of the n-th point is
    %     1/N^2 * sum_(j=1,..,N) k'(d(j,n)) / d(n,j) * ( Y(:,n) - Y(:,j) )
    g = zeros( D, N );
    for n=1:N
        % column vector of the kernel derivatives evaluated at d(:,n)
        dPP     = dPointPotential(d(:,n), energy_index, sym_flag);
        % column vector containing the quotients of kernel derivatives and
        % corresponding distances d
        dPP_d   = dPP ./ d(:,n);
        % exclude undefined "self-gradient". physically: no force of a
        % point on itself
        dPP_d(n)= 0;
        % assemble the gradient (vectorized implementation)
        g(:,n)  = sum(dPP_d) * Y(:,n) - Y * dPP_d;
        % the vector g(:,n) is now projected onto the tangent plane of the
        % sphere at the n-th point Y(:,n)
        g(:,n)  = (g(:,n) - Y(:,n)'*g(:,n) * Y(:,n)) /l(n);
    end
    % reshape g --> proper dimension for the use with fminunc etc.
    g = reshape(g, [], 1);
    % apply appropriate normalization factor
    if( energy_index == -2 )
        g = g/N^2;
    else
        g = g/(N*(N-1));
    end
    if( sym_flag == 1)
        g = g/2;
    end
end