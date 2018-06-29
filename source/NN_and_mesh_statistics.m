function [NN_xi_MinMeanMax, meshnorm, meshratio, largestgap] = NN_and_mesh_statistics(X,sym_flag,Ny)
%% NN_AND_MESH_STATISTICS compute nearest neighbor (NN) and mesh statistics
%
% Inputs
% X          directions ( D-by-N matrix )
% sym_flag   symmetry flag (i.e. consider point set X and -X if sym_flag == 1)
% Ny         number of test points
% 
% Outputs
% NN_xi_MinMeanMax  
% meshnorm          mesh norm
% meshratio         mesh ratio (2 * ratio of meshnorm and smallest nearest neighbor distance)
% largestgap        largest gap in mesh (point coordinate within Y leading to
%                   the largest distance w.r.t. X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% NN_and_mesh_statistics.m
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

% disp('Calculating NN statistics')
[~, NN_xi]          = GeodesicDistance(X,sym_flag);
NN_xi_min           = min(NN_xi);
NN_xi_mean          = mean(NN_xi);
NN_xi_max           = max(NN_xi);
NN_xi_MinMeanMax    = [NN_xi_min NN_xi_mean NN_xi_max];

% disp('Calculating "mesh" statistics')
if ~exist('Ny','var')
    Ny  = max( 10 * size(X,2), 1e5 );  % the "full" sphere, should be very fine
elseif Ny <= 0
    Ny  = max( 10 * size(X,2), 1e5 );
end

% mesh ratio = mesh norm / NN_xi_min * 2
Y                        = RandomDirections(size(X,1),Ny);
[meshnorm, largestgap]   = MeshNorm(X,Y,sym_flag);
meshratio                = 2.0 * meshnorm / NN_xi_min;
