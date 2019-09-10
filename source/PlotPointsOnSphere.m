function [h] = PlotPointsOnSphere(X, Y, sym_flag)
% PLOTPOINTSONSPHERE plots the point set X on the surface of the unit
% sphere in 3 dimensional space, i.e. size(X,1) must be 3.
%
% Also, a little sphere in the center is indicated, so symmetry can be
% judged on a visual basis (rotating the points such that symmetric pairs
% overlap when the front point covers the center.
%
% If a point set Y is given, it is also plotted (in red color).
%
% If sym_flag==1, then -X and -Y are also plotted.
%
% Examples
% 1) plot one point set X:
% >> PlotPointsOnSphere(X);
% 2) plot two point sets X and Y:
% >> PlotPointsOnSphere(X,Y);
% 3) plot two point sets X, Y and their mirror sets -X, -Y:
% >> PlotPointsOnSphere(X,Y,1);
% 4) plot one point set X and its mirror set -X:
% >> PlotPointsOnSphere(X,[],1)
%
% Inputs
% X           point set of size 3 x N
% Y           [OPTIONAL] second point set, might be empty (i.e. Y = [])
% sym_flag    [OPTIONAL] if 0, then plot the coordinates X (and Y) as
%             given. if 1, then additionally plot -X (and -Y) in red color.
%             defaults to 0.
%
% Output
% h           handle to created figure
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% PlotPointsOnSphere.m
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

%% preparation
if size(X,1)~=3
    error(['size(X,1) must be 3 but is ' num2str(size(X,1))]);
end
if exist('Y','var')
    if size(Y,1)~=3 && ~isempty(Y)
        error('if Y is provided, size(Y,1) must be 3, or Y=[]');
    end
    if size(Y,1)==3
        figname = [num2str(size(X,2)) ' points (X) and largest gap (Y)'];
    elseif isempty(Y)
        figname = [num2str(size(X,2)) ' points'];
    end
else
    Y = [];
    figname = [num2str(size(X,2)) ' points'];
end
if exist('sym_flag','var')
    if sym_flag==1
        figname = [figname ', symmetry on'];
    elseif sym_flag~=0
        error('sym must be either 0 or 1');
    end
else
    sym_flag = 0;
end

% create coordinates of sphere
[xx,yy,zz] = sphere(30);

% define colors, sizes, and sphere's and transparency (alpha)
color_sphere=[0.78,0.81,0.83];
color_small_sphere=[0  1  1];
alpha_sphere=0.7;
alpha_small_sphere=alpha_sphere;
size_sphere =0.997;
size_small_sphere = size_sphere/60;
% color_points=[ 0    0.4470    0.7410];
color_X=[ 0    0    1];
% color_sym_points=[0.8500    0.3250    0.0980];
color_Y=[1    0    0];
point_size = 10;

%% go
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

% plot X, and optionally Y
xhandle = plot3(X(1,:),X(2,:),X(3,:), ...
    '.','markersize',point_size+1,'Color',color_X);
if ~isempty(Y)
    yhandle = plot3(Y(1,:),Y(2,:),Y(3,:), ...
        'x','markersize',point_size, 'linewidth', 3,'Color',color_Y);
end

% if sym, then plot -X, and optionally -Y
if sym_flag==1
    xhandle_sym = plot3(-X(1,:),-X(2,:),-X(3,:), ...
        '.','markersize',point_size,'Color',...color_X);
        'red');
	if ~isempty(Y)
		yhandle = plot3(-Y(1,:),-Y(2,:),-Y(3,:), ...
			'x','markersize',point_size, 'linewidth', 3,'Color',color_Y);
	end
end

% some view options
axis equal vis3d
xlabel('x')
ylabel('y')
zlabel('z')

if ~isempty(inputname(1))
    Xname = inputname(1);
else
    Xname = 'X';
end
if nargin >= 2 && ~isempty(inputname(2))
    Yname = inputname(2);
else
    if sym_flag==0
        Yname = 'Y';
    else
        Yname = 'Y, -Y';
    end
end

if sym_flag==0
    if exist('yhandle','var')
        legend([xhandle yhandle], Xname, Yname);
    else
        legend(xhandle,'X');
    end
else
    if exist('yhandle','var')
        legend([xhandle yhandle xhandle_sym yhandle], Xname, Yname, ['-' Xname], ['-' Yname]);
        legend([xhandle xhandle_sym yhandle ], Xname, ['-' Xname], Yname);
    else
        legend([xhandle xhandle_sym], Xname, ['-' Xname]);
    end
end
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
ax=gca;
ax.XTick=[-1,0,1];
ax.YTick=[-1,0,1];
ax.ZTick=[-1,0,1];
hold off
rotate3d