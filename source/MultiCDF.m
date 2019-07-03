function [ fig ] = MultiCDF( X, sym_flag )
%MULTICDF plots empirical cumulative distribution functions for all of the
%directions in X. More precisely, X contains the directions as columns, and
%their distances are computed by calling GeodesicDistance.
%
% INPUT
% X             matrix containing directions as columns. will be normalized
% sym_flag      [optional] if 1, then plot CDF only for distances up to
%               pi/2
%
% OUTPUT
% fig           handle to the created figure
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% MultiCDF.m
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

% check inputs
if exist('sym_flag','var')
    if sym_flag ~= 0 && sym_flag ~= 1
        error('sym_flag must be either 0 or 1')
    end
else
    sym_flag = 0;
end

fig = figure;

% compute full geodesic distance matrix
GD = GeodesicDistance(X);

% plot all CDFs
for n=1:size(GD,1)
    % Plot
    cdfplot(GD(:,n))
    hold on
end

% rescale x-axis
if sym_flag == 1
    xlim([0 pi/2])
else
    xlim([0 pi])
end




end

