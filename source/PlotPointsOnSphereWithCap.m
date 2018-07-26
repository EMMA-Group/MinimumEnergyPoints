function [h] = PlotPointsOnSphereWithCap(X, capcenter, capradius)
% [] = PlotPointsOnSphere(X, Y, sym)
% Plot point set X on the surface of the unit sphere in 3 dimensional space
% and also a little sphere in the center, so symmetry can be judged on a
% visual basis. If point set Y is given, it is also plotted. If 'sym' is 1,
% then -X and -Y are also plotted.
%
%  Inputs
%  X           point set of size D-by-N (i.e. rows = coordinates)
%  Y           optional second point set
%  sym         if 0, then plot the coordinates X and Y as given. if 1, then
%              additionally plot -X and -Y in red color. default is 0.
%
%  Output
%  h           handle to created figure
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% PlotPointsOnSphereWithCap.m
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

if size(X,1)~=3
    error('size(X,1) must be 3');
end

figname = [num2str(size(X,2)) ' points'];

% if exist('sym','var')
%     if sym==1
%         figname = [figname ', symmetry on'];
%     elseif sym~=0
%         error('sym must be either 0 or 1');
%     end
% else
%     sym = 0;
% end

if exist('capcenter','var')
    if size(capcenter,1)~=3
        error('size(capcenter,1) must be 3');
    end
    capcenter = RenormalizeColumns(capcenter);
else
    capcenter = [0;0;1];
    disp('set cap center to default (north pole)')
end
if exist('capradius','var')
    if length(capradius) ~= 1
        error('length(capradius) must be 1, or capradius must not be given')
    elseif capradius > pi-0.01
        error('capradius must be pi-0.01 at most')
    elseif capradius < 0.01
        error('capradius must be 0.01 at moost')
    end
else
    capradius = pi/4;
    disp('set cap radius to default (pi/4)')
end




% define colors, sizes, and sphere's and transparency (alpha)
color_sphere=[0.78,0.81,0.83];
color_small_sphere=[0  1  1];
alpha_sphere=0.6;
alpha_small_sphere=alpha_sphere;
size_sphere =0.997;
size_small_sphere = size_sphere/60;
% color_points=[ 0    0.4470    0.7410];
color_X_on_cap=[ 1    0    0];
color_X_no_cap=[ 0    0    1];
% color_sym_points=[0.8500    0.3250    0.0980];
point_size = 6;



% create coordinates of sphere
[xx,yy,zz] = sphere(30);


% create figure
h = figure('Name',figname, ...
    'NumberTitle','off', ...
    'OuterPosition',[100 100 500 500]);

% plot main sphere
surf(xx*size_sphere, yy*size_sphere, zz*size_sphere, ...
    'LineStyle','none', ...
    'FaceColor',color_sphere, ...
    'FaceAlpha',alpha_sphere);
hold on;

% plot small center sphere for visual inspection of symmetry
surf(xx.*size_small_sphere, yy.*size_small_sphere, zz.*size_small_sphere, ...
    'LineStyle','none', ...
    'FaceColor',color_small_sphere, ...
    'FaceAlpha',alpha_small_sphere);


% plot spherical cap (before plotting points, s.t. points are on top)
if exist('capcenter','var')
    color_cap = [0 1 0];
    plot3(capcenter(1), capcenter(2), capcenter(3), '.', 'markersize', point_size+5, 'Color', color_cap);
    hold on
    
    % coordinates of cap border if capcenter was the north pole
    n_capborderpoints = 1000;
    phi = linspace(0, 2*pi, n_capborderpoints+1)';
    phi = phi(1:end-1);
    capborderpoints = zeros(n_capborderpoints,3); % num points x dimensions
    capborderpoints(:,1) = cos(phi)*sin(capradius);
    capborderpoints(:,2) = sin(phi)*sin(capradius);
    capborderpoints(:,3) = cos(capradius);
    
    % alpha = angular distance of capcenter to the north pole
    alpha = real(acos(min(1,max(-1,[0 0 1]*capcenter))));
    if alpha > 1e-3
        % assemble rotation matrix, mapping the north pole onto the capcenter
        %   direction of rotation axis
        r = RenormalizeColumns(cross([0;0;1],capcenter));

        %   cross product matrix of r
        cross_mat = [0 -r(3) r(2); r(3) 0 -r(1); -r(2) r(1) 0];
        %   rotation matrix
        R = cos(alpha)*eye(3) + sin(alpha)*cross_mat ...
            +  (1-cos(alpha))*(r*r');

        % rotate the cap
        capborderpoints = capborderpoints*R';
    end
    
    % plot the cap border
    plot3(capborderpoints(:,1), capborderpoints(:,2), capborderpoints(:,3), ...
        '.', 'markersize', point_size, 'color', color_cap);
end






% get points that are within cap, and those that aren't
idx = real( acos( min(1,max(-1,X'*capcenter)) ) ) < capradius;
X_on_cap = X(:,idx);
X_no_cap = X(:,~idx);

% plot X ...
% ... on cap
plot3(X_on_cap(1,:),X_on_cap(2,:),X_on_cap(3,:), ...
    '.','markersize',point_size,'Color',color_X_on_cap);
hold on
% ... no cap
plot3(X_no_cap(1,:),X_no_cap(2,:),X_no_cap(3,:), ...
    '.','markersize',point_size,'Color',color_X_no_cap);
hold off
axis equal
rotate3d


