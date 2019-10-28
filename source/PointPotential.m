function U = PointPotential(d, energy_index, sym_flag)
%% POINTPOTENTIAL computes the potential of one point of a set of spherical
% points. Let this particular point be x_p. Then the return value is
% U = sum_(i=1,...,N) k_s(x_p, x_i)
% where s is the variable "s" in the paper, here s = energy_index.
%
% If sym_flag==1, then symmetrized kernels k_s are used. The antipodes are
% not considered explicitly, i.e. they are not regarded as individual
% points. Instead, symmetrized kernels are implemented.
%
% User-defined kernel functions have to be implemented here and in
% DPOINTPOTENTIAL.
%
% Input
% d             N x 1 vector of distances of all N points from a single
%               point of that set (i.e. including one zero self-distance).
%               In the symmetric case, d does NOT include the distances
%               to the antipodes, as these distances are computed here
%               explicitly and efficiently.
% energy_index  index of the energy function to be used ("s" in the paper)
% sym_flag      also include the contributions of the antipodes to the
%               potential
%
% Output
% U             potential of the point
%
% See also DPOINTPOTENTIAL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% PointPotential.m
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

%% NOTES ON USER-DEFINED KERNEL FUNCTIONS:
% - d contains the distances of all points to one point, hence the point
% potential is the sum of the kernel function evaluated at each entry of d.
% - also define the symmetrized kernel function
% - vectorized implementations are STRONGLY encouraged, as this is the
%   lowest level of the function hierarchy.

%% compute the point potential depending on the energy_index
if energy_index == -2
    %% LOG case
    U   = sum(  d .* (log(d/2)-1.0) ) + 2*length(d) ;
    if sym_flag == 1
        d_antipode = sqrt(4-d.^2);
        U   = U + sum(  d_antipode .* (log(d_antipode/2)-1.0) ) + 2*length(d);
    end
    % ATTENTION: NaN's are converted to 2's. This should (and does) work,
    % but is not entirely failsafe.
    U(isnan(U)) = 2;
elseif energy_index == -1
    %% log case
    U   = -sum( log(d) ) + log(2)*length(d) ;
    if sym_flag == 1
        U   = U - sum( log(sqrt(4-d.^2)) ) + log(2)*length(d);
    end
elseif energy_index > 0
    %% Riesz case
    U   = sum( 1./ (d.^energy_index) );
    if sym_flag == 1
        U   = U + sum( 1./ (sqrt(4-d.^2).^energy_index) );
    end
else
    error(['PointPotential not yet implemented for energy_index = ',...
        num2str(energy_index)])
end
