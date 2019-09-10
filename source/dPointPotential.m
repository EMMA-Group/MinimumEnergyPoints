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