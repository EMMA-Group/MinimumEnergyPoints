%% Script_ME_Points This is a template-style example script call (e.g. for looping through parameters). For a GUI, see ../Start_ME_Points

%% Example script for generating multiple sets of minimum energy points
% loops over dimension D, number of points N, energy index s
% writes .txt and .mat files with the resulting point sets
% stores the parameters to the .mat file
% saves a log of all outputs
%
% For high-performance computations, make sure to comment non-essential
% parts of MinimizeEnergy and MinimizeGradient, immediately before the last
% last display output of the respective cycle.

addpath(genpath('..'));
c=clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% BatchGeneratePoints.m
% Copyright (C) 2018, Felix Fritzen and Oliver Kunc
% All rights reserved.
%
% This source code is licensed under the BSD 3-Clause License found in the
% LICENSE file in the root directory of this source tree.
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
% Oliver Kunc and Felix Fritzen: 'Generation of energy-minimizing point sets on
% spheres and their application in mesh-free interpolation and differentiation',
% Advances in Computational Mathematics 45(5-6), pp. 3021-3056, 2019
% URL: doi.org/10.1007/s10444-019-09726-5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% names of output files and of diary (i.e. logbook of all screen outputs)
filename_prefix = 'myprefix';
diary(['exports/log_BatchGenerateDirections_',date,'_',num2str(c(4)),'_',num2str(c(5))])

%% set the parameters

% define Euclidean dimensions, i.e. points will be on the (D-1)-sphere
all_D = [4,5,6];
% define number of points
all_N = [32];%,41,64,83,128,193,256,383,512,1024,1319,1889,2048];
% use symmetrized kernels?
sym_flag = 0;
% define cases of s to be considered: -2 = LOG, -1 = log, >0 = Riesz
all_s = [-2,-1,0.5,1,2,3];

% solver options
n_cycle_energy = 100; % number of cycles for energy minimization, unconstrained optimization
n_it_energy    = 100; % number of unconstrained energy iterations (per energy cycle) before points are projected onto the sphere again
n_cycle_gradient=  0; % number of cycles for unconstrained least square minimization of the energy's gradient. costly!
n_it_gradient  =   1; % number of unconstrained gradient iterations (per gradient cycle) before points are projected onto the sphere again

% additional information
n_fixpoints = 0; % fixpoints not yet implemented
creation_date = datestr(now);

%% loop
for D = all_D
    for N = all_N
        % initial points: Equal Area [Leopardi 2006] with perturbation
        % the same initial positions are used for all s
        Xstart_origin = 'Equal Area';
        Xstart_name = 'Equal Area';
        Xstart = RenormalizeColumns( 0.01*RandomPoints(D,N) + eq_point_set(D-1,N) );
        for s = all_s

            % something for the diary
            D
            N
            sym_flag
            s

            % store parameters in structure
            parameters = struct;
            parameters.D = D;
            parameters.N = N;
            parameters.sym_flag = sym_flag;
            parameters.energy_index = s;
            parameters.n_cycle_energy = n_cycle_energy;
            parameters.n_it_energy = n_it_energy;
            parameters.n_cycle_gradient = n_cycle_gradient;
            parameters.n_it_gradient = n_it_gradient;
            parameters.Xstart_origin = Xstart_origin;
            parameters.Xstart_name = Xstart_name;
            parameters.Xstart = Xstart;
            parameters.n_fixpoints = n_fixpoints;
            parameters.creation_date = creation_date;

            %% start the energy minimization by fminunc
            disp('------------------------------------')
            disp('Use fminunc to minimize energy...')
            tic
            [ X, nn_xi, X0, nn_xi0, cyc ] = MinimizeEnergy( D, N, n_it_energy,  n_cycle_energy, Xstart, s, sym_flag );
            toc
            %% start the explicit gradient minimization
            disp('------------------------------------')
            disp('Use lsqnonlin to minimize gradient...')
            tic
            X = MinimizeGradient( D, N, n_it_gradient, n_cycle_gradient, X, s, sym_flag );
            toc

            %% save to .mat and .txt file
            full_filename = ['exports/', filename_prefix, '_D', num2str(D), '_N', num2str(N), '_sym', num2str(sym_flag), '_s', num2str(s)];
            save([full_filename,'.mat'], 'parameters', 'X');
            dlmwrite([full_filename,'.txt'], X','delimiter','\t','precision','%30.23e'); % remember: X is transposed here!

            disp(['Finished gegnerating and exporting', num2str(N), ' directions in ', ...
                num2str(D), ' dimensions. filename:']);
            disp(full_filename);

        end % s
    end % N
end % D
