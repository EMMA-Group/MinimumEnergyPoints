function dw = dPairPotential(xi, energy_index, sym_flag)

%% DPAIRPOTENTIAL computes the derivative of the repulsion energy that
% corresponds to input distance
%
% Input
% xi     geodesic distances as vector
% energy_index  index of the energy function to be used
% sym_flag      use symmetrized point sets
%
% Output
% dw     derivatives of repulsion energy corresponding to the entries of xi
%
% See also PAIRPOTENTIAL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% dPairPotential.m
% Copyright (C) 2018 by Felix Fritzen and Oliver Kunc
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
% JOURNAL NAME, Number/Volume, p. XX-YY, 2018
% DOI   ...
% URL   dx.doi.org/...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% NOTE 1: the functions W and dW can be modified at will
% Hoewever, they must do the following:
% - given an input vector x, W(x) --> sum of W(x_i), with W(x_i) being
%   the value of pair potential for a scalar distance x_i >= 0
% - dW(x) must return a vector, i.e. dW_i = diff( W(x), x_i )
% - vectorized implementations are STRONGLY encouraged
% See below for details
%% NOTE 2: the symmetry flag has already been checked for consistency in
% GENERATEDIRECTIONS, hence it is either 0 or 1.


if energy_index == 0
    %% Logarithmic energy function. Distinct symmetry cases for efficiency.
    dw   = log(xi/pi);
    if sym_flag == 1
        dw   = dw - log(1 - xi/pi);
    end
%     dw   = ( log( xi / pi ) + 1 ) /pi;
%     if sym_flag == 1
%         dw   = dw - ( log(1 - xi/pi) + 1 ) /pi;
%     end
%     if sym_flag == 0
%         dw   = - ( 1 - log(xi/pi) ) / pi ;
%     else
%         dw   = ( - log(xi/pi)  / pi +(1-xi/pi) .* log(1-xi/pi) ) / pi ;
% %     end
%     if sym_flag == 0
%         dw   = - 0.5 * cos(xi/2)./(sin(xi/2).^2);
%     else
%         dw   = - 0.5 * ( cos(xi/2)./(sin(xi/2).^2) - cos(pi/2 - xi/2)./(sin(pi/2 - xi/2).^2) ) ;
%     end

elseif energy_index > 0
    %% hyperbolic function as force (no cutoff)
    dw   = - energy_index ./ (xi.^(energy_index+1));
    if sym_flag == 1
        dw   = dw + energy_index ./ ( (pi-xi).^(energy_index+1) );
    end 
end


