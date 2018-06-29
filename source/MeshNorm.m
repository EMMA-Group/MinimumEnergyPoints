function [meshnorm, largestgap] = MeshNorm(X, Y, sym_flag)
%% MESHNORM returns the mesh norm of a spherical point set X, where the
% full sphere is approximated by another, very large point set Y. Also,
% optionally, the column of Y that has the maximum nearest neighbor
% distance to X is returned as largestgap.
%
% INPUTS
% X            matrix D-by-Nx, normalized columns
% Y            matrix D-by-Ny, normalized columns. Ny >> Nx
% sym_flag     symmetry flag, indicating whether or not to use symmetric
%              distance function (i.e. angular distance modulo pi)
%
% OUTPUTS
% meshnorm     an approximation to the mesh norm
% largestgap   an approximation to the position of the largest gap

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% MeshNorm.m
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

% firstly, do some consistency checks
if size(X,1) ~= size(Y,1)
    error('size(X,1) must be the same as size(Y,1)');
end

if size(Y,2) < size(X,2)
    warning('size(X,2) should be (significantly) small than size(Y,2) !!!');
end



if ~exist('sym_flag','var')
    sym_flag = 0;
elseif sym_flag ~= 0 && sym_flag ~= 1
    error('sym_flag must be either 0 or 1')
end

% then, get the distance matrix. first index: X, second index: Y
if sym_flag == 0
    D = real(acos(min(1,max(-1,X'*Y))));
else
    D = real(acos(abs(min(1,max(-1,X'*Y)))));
end


% recall mesh norm definition: sup_y min_x D
% now the task is: find the column for which the smallest entry is largest,
% and take that entry as meshnorm, and the corresponding column of Y as
% largest gap coordinates
minD = min(D,[],1);   % row vector with nearest neighbor distances for all test points in Y
[meshnorm, largestgap_index] = max(minD); % maximum nearest neighbor distance, and the corresponding point index of Y
if nargout > 1
    largestgap = Y(:,largestgap_index);
end

end
