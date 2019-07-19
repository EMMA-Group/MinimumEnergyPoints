function [ xi, F ] = EmpiricalDistributionFunction( X, n )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% EmpiricalDistributionFunction.m
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

d = size(X,1);
N = size(X,2);
M = N;
if(nargin == 2)
    M = n;
end

R = RandomDirections(d, M);

Dist = R'*X;
Xi    = acos( min( max( Dist, -1 ), 1 ) );

figure;
subplot(1,2,1);
grid on;
hold on;
% xiref = linspace(0,pi,N);
% Fref  = 0.5*(1-cos(xiref));
Fref = linspace(0,1,N);
xiref  = acos(1-2*Fref);
err = 0.;
for j=1:M
    err = err + norm(sort(Xi(j,:),'ascend') - xiref)^2;
end
err = sqrt(err / norm(xiref)^2 / M);
fprintf('average error: %16.8f%%\n', err*100.);
for j=1:min(20,M)
    plot(sort(Xi(j,:),'ascend'),linspace(0,1,N), '-', 'linewidth',1, 'color', 'blue');
end
plot(xiref,Fref,'color','yellow');
xi = reshape(Xi, [],1);
xi = sort(xi, 'ascend');
F = linspace(0,1,M*N);
subplot(1,2,2);
plot(xi,F, '-', 'linewidth',3, 'color', 'blue');

x = linspace(0,pi,1000);
y = 0.5*(1-cos(x));
hold on;
grid on;
plot(x,y,'color', 'yellow');
y = 0.5*(1-cos(x));
