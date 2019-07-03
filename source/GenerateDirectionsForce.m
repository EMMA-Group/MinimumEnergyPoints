function [ X, theta, X0, theta0 ] = GenerateDirectionsForce( D, N, nit, ncycle, Xstart, energy_index, sym_flag )
%% GENERATEDIRECTIONSFORCE produces a set of N directions of
%  dimension D such that the global energy computed from
%  the sum of pair potentials is minimized (an
%  approximate minimum is computed)
%
%  Inputs
%  D           dimension of the directions
%  N           number of requested directions
%  nit         number of iterations per cycle
%  ncycle      number of cycles
%  Xstart     initial point set
%  energy_index index of the energy function to be used
%  sym_flag   use symmetrized ernergy function
%
%  Outputs
%  X           matrix D-by-N containing the resulting directions
%  nn_xi       nearest neighbor distance [radian] of resulting directions
%  X0          initial directions, D-by-N
%  nn_xi0      nearest neighbor distance [radian] of initial directions
%
%
%  See also REPULSIONENERGY, GEODESICDISTANCE, FMINUNC,
%  RENORMALIZECOLUMNS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% GenerateDirectionsForce.m
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


addpath(genpath('../eq_sphere_partitions'));


%% settings for fminunc
%  NOTE: * after each nit iterations the point set is renormalized this
%          makes the abortion criterion more reliable (this is one cycle)
%        * this is repeated of a numbe rof cycles
%        * alternative settings are possible;
%        * for the trust region algorithm memory requirement is massive
%        * tolerances are usually not reached before the iteration count
%        * 'bfgs' instead of 'steepdesc' did not show better performance but
%          may be worth a try (then: less iterations <--> same computing time)
%  %  options = optimoptions('fminunc',      ...
%  %      'SpecifyObjectiveGradient', true,  ...
%  %      'algorithm',                'quasi-newton',    ...
%  %      'hessupdate',               'steepdesc',       ...
%  %      'steptolerance',            1e-6, ...
%  %      'functiontolerance',        1e-8, ...
%  %      'optimalitytolerance',      1e-8, ...
%  %      'MaxIterations',            nit, ...
%  %      'display',                  'none');
options = optimoptions(@lsqnonlin,      ...
    'algorithm',                'levenberg-marquardt',    ...
    'steptolerance',            1e-12, ...
    'functiontolerance',        1e-14, ...
    'optimalitytolerance',      1e-14, ...
    'MaxIterations',            nit, ...
    'display',                  'final');


%% equal area point set with random perturbation as initial guess
if size(Xstart,1)==D && size(Xstart,2)==N
    disp('Using provided directions as initial values')
    X0 = Xstart;
else
    disp(['Creating ' num2str(N) ' perturbed EQ point set as initial positions. Ignoring Xstart.'])
    X0 = RenormalizeColumns(0.01 * randn( D, N ) + eq_point_set( D-1, N));
end
X      = X0;

%% compute initial residual and information on the distribution
%% compute distance maps
[ D0, xi0 ] = GeodesicDistance(X0);
NN_min      = min(xi0);
NN_max      = max(xi0);
NN_mean     = mean(xi0);
[ff, gg]    = RepulsionEnergy( reshape(X, [], 1), D, energy_index, sym_flag );


fprintf('initial:    f, |g| %16.8e, %16.8e .. NN %6.4f / %6.4f / %6.4f ... range %7.5f; sqrt(V(NN)) %10.6f\n', ...
    ff, norm(gg), NN_min, NN_mean, NN_max, NN_max - NN_min, sqrt(var(xi0)) );

%% loop of the cycles
nretry = 0;
dX = X;
jcycle=0;
sc=0;
while(jcycle<ncycle && sc == 0)
    jcycle=jcycle+1;
    X   = reshape( lsqnonlin( @(x)RepulsionForce(x,D, energy_index, sym_flag), reshape( X, [], 1 ), [],[],options ), D, N );
    X   = RenormalizeColumns(X);
    dX  = X - dX;
	if( norm(dX, 'fro')  < sqrt( D*N * 1e-12 ) )
		sc = 1;
		fprintf('# fixed point found. stopping iteration.\n');
	else
		dX  = X; % backup current direction
	end	
%  	else
%  		nretry = 0;
%  	end
    [ D1, xi1 ] = GeodesicDistance( X );
    % additional function calls for evaluation of the direction set after the N-th iteration
    % not needed except for displaying the quality of the iterate (can be commented)
    [ ff, gg  ] = RepulsionEnergy( reshape(X, [], 1), D, energy_index, sym_flag );
    NN_min      = min(xi1);
    NN_max      = max(xi1);
    NN_mean     = mean(xi1);
    fprintf('cycle %4i  f, |g| %16.8e, %16.8e .. NN %6.4f / %6.4f / %6.4f ... range %7.5f; sqrt(V(NN)) %10.6f\n', ...
        jcycle, ff, norm(gg), NN_min, NN_mean, NN_max, NN_max - NN_min, sqrt(var(xi1)) );
end

[~, theta ]  = GeodesicDistance( X,  sym_flag );
[~, theta0 ] = GeodesicDistance( X0, sym_flag );
