%% Example program for the identification of point sets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% BatchGenerateDirections.m
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

%% free memory, close open windows, add parent directory and all of its subdirectories
close all;
clear all;
addpath(genpath('..'));

% options
sym_flag = 1;
energy_index=0;

for D=3
    nit     = 20;   % number of energy optimization iterations per cycle
    ncycle  = 200;  % number of energy optimization cycles (after which the points are projected onto the sphere)
    nitF    = 0;    % number of gradient optimization iterations per cycle
    ncycleF = 0;    % number of gradient optimization cycles (after which the points are projected onto the sphere)
    Nrange_two = [4 8 16 32 64 128 512 1024 2048 4096];
    Nrange_prime = [41 83 193 383 701 881 1319 1889 2383 3121 3797];
    Nrange = sort([Nrange_two, Nrange_prime]);
    for N=Nrange
        N
        if( D == 9 )  % decrease nit for large dimensions s.t. runtime becomes practical
            nit = 1.5*nit;
        end
        if N > 2000
            ncycle = 50;
        end
        %% first use fminunc ... then directly minimize gradient
        disp('Use fminunc to minimize energy...')
        tic
        [ X ]           = GenerateDirections( D, N, nit,  ncycle, [], energy_index, sym_flag );
        toc
%         PlotPointsOnSphere(X,X0,sym_flag);

        %% ... then directly minimize gradient
        disp('Minimize gradient by directly simulating repulsion...')
        tic
        [ X ] = GenerateDirectionsForce( D, N, nitF, ncycleF, X, energy_index, sym_flag );
        toc
        
        %% postprocessing
        % nearest neighbor (NN) and mesh statistics
        [NN_xi_MinMeanMax, meshnorm, meshratio] = NN_and_mesh_statistics(X,sym_flag);
        
        % search gamma: first on coarse grid, then on fine grid
        gamma_min = 0.1;
        gamma_max = 5;
        ngamma = 100;
        disp('coarse fit of gamma w.r.t. POU')
        [ gamma_best ] = SearchGamma( X, 4096, gamma_min, gamma_max, ngamma, 'U' );
        gamma_min = (gamma_min + 2*gamma_best) / 3;
        gamma_max = min(gamma_max, 2*gamma_best - gamma_min);
        disp('fine fit of gamma w.r.t. POU')
        [ gamma_best, POU_diff_to_1 ] = SearchGamma( X, 4096, gamma_min, gamma_max, ngamma, 'U' );
        
        % save .mat file
%         save(sprintf('exports/TESTdir_d%i_n%i_xlogx.mat', D, N), 'X', 'gamma_best', 'POU_diff_to_1', 'NN_xi_MinMeanMax', 'meshnorm', 'meshratio' );

        % plot against EQ points
%         Y = eq_point_set(D-1,N);
%         PlotPointsOnSphere(X,Y);
    end
end

