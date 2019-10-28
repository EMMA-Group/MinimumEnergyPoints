function [ X, d_nn, X0, d_nn0 ] = MinimizeGradient( D, N, nit, ncycle, Xstart, energy_index, sym_flag )
%% MINIMIZEGRADIENT produces a set of N minimum energy points on the
% (D-1)-dimensional hypersphere in the D-dimensional Euclidean space. The
% energy is the s-energy as in the paper, here symbolically s=energy_index.
% This function uses the Matlab/Octave built in function lsqnonlin to
% minimize the norm of the gradient, moving the points off the sphere.
% Therefore, the execution of lsqnonlin is limited to nit iterations, after
% which the point coordinates are re-normalized. Then, lsqnonlin is called
% again. Each execution of lsqnonlin is called a "cycle".
%
% The retry cycle loop (as in MinimizeEnergy) is not implemented atm
% because it doesn't appear to be necessary.
%
% Hint on usage: prefer MinimizeEnergy over MinimizeGradient! The former is
% much faster as the latter. Especially for large N, you should not use the
% latter at all. E.g. currently (2019), for N > 1000 it may take hours on a
% standard workstation to finish a single cycle with few iterations.
%
% Inputs
% D             dimension of the containing Euclidean space
% N             number of points
% nit           number of iterations per call of fminunc (i.e. per cycle)
% ncycle        number of fminunc calls (here: "cycles"). coordinates are
%               normalized after each cycle.
% Xstart        [OPTIONAL] initial point set, size D x N (or empty). If
%               empty or if not provided, then perturbed Equal Area Points
%               are used.
% energy_index  index of the energy function to be used, corresponding to
%               "s" in the paper.
%               energy_index == -2: LOG
%               energy_index == -1: log
%               energy_index  >  0: Riesz
% sym_flag      [OPTIONAL] use symmetrized kernel function if sym_flag==1
%
% Outputs
% X             matrix D x N containing the resulting points' Euclidean
%               coordinates as columns
% d_nn          nearest neighbor distance of resulting points
% X0            initial point set, D x N. equals Xstart if it was provided.
% d_nn0         nearest neighbor distance of initial points
%
% See also FMINUNC, RENORMALIZECOLUMNS, EQ_POINT_SET, MINIMIZEENERGY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% MinimizeGradient.m
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
% Oliver Kunc and Felix Fritzen: ''
% JOURNAL NAME, Number/Volume, p. XX-YY, 2019
% DOI   ...
% URL   dx.doi.org/...
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../eq_sphere_partitions'));

%% energy index
if ~isnumeric(energy_index)
    error('energy index must be numeric')
elseif ~(energy_index>0) && ~(energy_index==-1) && ~(energy_index==-2)
    error(['energy index ',num2str(energy_index),' is not implemented'])
end
if energy_index > 0
    disp(['Case: s = ',num2str(energy_index),' (Riesz)'])
elseif energy_index == -1
    disp('Case: s = log')
else
    disp('Case: s = LOG')
end

%% symmetry flag
if ~exist('sym_flag','var')
    sym_flag = 0;
elseif sym_flag ~= 0 && sym_flag ~= 1
    error('Symmetry flag must be either 0 or 1')
end

if sym_flag == 1
    disp('Symmetry flag: ON')
else
    disp('Symmetry flag: OFF')
end

%% initial guess: either provided, or equal area point set with random perturbation
if ~exist('Xstart','var') % this also means sym_flag == 0 because of ordering!!!
        X0          = 0.01 * RenormalizeColumns( randn( D, N ) );
        X0          = X0 + eq_point_set( D-1, N);
        X0          = RenormalizeColumns(X0);
else
    if prod(size(Xstart)==[D N])==1
        X0 = Xstart;
    else    % Xstart is not compatible with D and N
        if ~isempty(Xstart)
            error('Xstart must be either [] or of size [D N]')
        end % else: Xstart = []
        if sym_flag == 1
            N_X0 = N*2;         % double the amount of points...
        else
            N_X0 = N;
        end
        disp('creating initial points as perturbed EQ points')
        X0 = RenormalizeColumns( 0.01*randn( D, N_X0) + eq_point_set( D-1, N_X0) );
%         X0 = RandomPoints(D,N);
        if sym_flag == 1
            X0 = X0(:,1:end/2);  % ... and remove second half. EQ point sets are sorted nicely.
        end
    end
end
X = X0;

%% compute initial residual and information on the distribution
[ ~, d0 ]   = Distance(X0, sym_flag);
NN_min      = min(d0);  % minimum nearest neighbor distance
NN_mean     = mean(d0); % mean nearest neighbor distance
NN_max      = max(d0);  % maximum nearest neighbor distance
[I, g]      = EnergyAndGradient( reshape(X, [], 1), D, energy_index, sym_flag );

fprintf('initial:    I= %16.8e, |g| = %16.8e .. NN %6.4f / %6.4f / %6.4f ... NN range %7.5f; sqrt(Var(NN)) %10.6f\n', ...
    I, norm(g), NN_min, NN_mean, NN_max, NN_max - NN_min, sqrt(var(d0)) );

%% cycle loop
options = optimoptions(@lsqnonlin,      ...
    'algorithm',                'levenberg-marquardt',    ...
    'steptolerance',            1e-12, ...
    'functiontolerance',        1e-14, ...
    'optimalitytolerance',      1e-14, ...
    'MaxIterations',            nit, ...
    'display',                  'final');
dX          = X;    % difference of point coordinates between cycles
jcycle      = 0;    % cycle counter
sc          = 0;    % success indicator
while( jcycle<ncycle && sc==0)
    jcycle=jcycle+1;

    %%%% ACTUAL WORK OF THE CURRENT CYCLE IS HERE ...
    X   = reshape( lsqnonlin( @(x)Gradient(x,D, energy_index, sym_flag), reshape( X, [], 1 ), [],[],options ), D, N );
    X   = RenormalizeColumns(X);
    %%%%

    dX  = X - dX;
    if( norm(dX, 'fro')  < sqrt( D*N * 1e-12 ) )
		sc = 1;
		fprintf('# fix point found. stopping iteration.\n');
	else
		dX  = X; % backup current point
    end

    % TODO: implement retry cycle if necessary

    % this is just for outputting the evolution of the nearest neighbor
    % (NN) statistics. comment this for high-performance calls.
    [ I, g  ] = EnergyAndGradient( reshape(X, [], 1), D, energy_index, sym_flag );
    [ ~, d ] = Distance( X, sym_flag );
    NN_min      = min(d);
    NN_max      = max(d);
    NN_mean     = mean(d);
    fprintf('cycle %4i  I= %16.8e, |g| = %16.8e .. NN %6.4f / %6.4f / %6.4f ... NN range %7.5f; sqrt(Var(NN)) %10.6f\n', ...
        jcycle, I, norm(g), NN_min, NN_mean, NN_max, NN_max - NN_min, sqrt(var(d)) );
end

[~, d_nn ]  = Distance( X, sym_flag );
[~, d_nn0 ] = Distance( X0, sym_flag );
