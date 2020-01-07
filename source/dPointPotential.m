function dU = dPointPotential(d, energy_index, sym_flag)
%% DPOINTPOTENTIAL computes the derivative of a point potential w.r.t. the
% distances of this point to all other points. Let this particular point be
% x_p. Then the i-th component of the return vector is
% dk_s(x_p, x_i)
% where s is the variable "s" in the paper (here s = energy_index), dk_s is
% the derivative of the kernel function k_s, and x_i is any of the N points
%
% If sym_flag==1, then symmetrized kernels k_s are used. The antipodes are
% not considered explicitly, i.e. they are not regarded as individual
% points. Instead, symmetrized kernels are implemented.
%
% User-defined kernel functions have to be implemented here and in
% POINTPOTENTIAL.
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
% dU            derivative of the potential of the current point w.r.t.
%               distances of this point to all points (including itself)
%
% See also POINTPOTENTIAL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% dPointPotential.m
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



%% NOTES ON USER-DEFINED KERNEL FUNCTIONS:
% - since the point potential is the sum of the kernel function evaluated
%   at each entry of d, the derivative with respect to d is the vector of
%   kernel derivatives evaluated at the entries of d.
% - also define the symmetrized kernel function's derivative
% - vectorized implementations are STRONGLY encouraged, as this is the
%   lowest level of the function hierarchy.

%% compute the derivative of the point potential
if energy_index == -2
    %% LOG case.
    dU   = log(d/2);
    if sym_flag == 1
        d_antipode = sqrt(4-d.^2);
        dU   = dU - log(d_antipode/2) .* d ./ d_antipode;
    end
elseif energy_index == -1
    %% log case.
    dU   = -1./d;
    if sym_flag == 1
        dU   = dU + d ./ (4 - d.*d);
    end
elseif energy_index > 0
    %% hyperbolic function as force (no cutoff)
    % dk = - s/d^(s+1)
    dU   = - energy_index ./ (d.^(energy_index+1));
    if sym_flag == 1
        % dk_sym = -s/d^(s+1) + s*d / sqrt(4-d^2)^(s+2)
        dU   = dU + energy_index * d ./ ( sqrt(4-d.^2).^(energy_index+2) );
    end
else
    error(['dPointPotential not yet implemented for energy_index = ',...
        num2str(energy_index)])
end
