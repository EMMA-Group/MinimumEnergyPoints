function U = PointPotential(d, energy_index, sym_flag)

%% POINTPOTENTIAL computes the sum of repulsion energies corresponding to 
% given distances
%
% Input
% d             distances as vector
% energy_index  index of the energy function to be used
% sym_flag      use symmetrized point sets
%
% Output
% U             sum of kernel values of the distances in d
%
% See also DPOINTPOTENTIAL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% PointPotential.m
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

%% NOTE 1: the functions k and dk can be modified at will
% Hoewever, they must do the following:
% - given an input vector x, k(x) --> sum of k(x_i), with k(x_i) being
%   the value of point potential for a scalar distance x_i >= 0
% - dk(x) must return a vector, i.e. dk_i = diff( k(x), x_i )
% - vectorized implementations are STRONGLY encouraged
% See below for details
%% NOTE 2: the symmetry flag has already been checked for consistency in
% GENERATEDIRECTIONS, hence it is either 0 or 1.
%% begin choosing the energy function corresponding to energy_index

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
    error('PointPotential not yet implemented for energy_index !>= 0')
end

