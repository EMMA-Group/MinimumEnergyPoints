function [ gamma_best, RMSE, V ] = SearchGamma( X, vali, ...
    gamma_min, gamma_max, ngamma, partition_case, xitrunc, do_plot )
% SearchGamma  Find best gamma in a provided range w.r.t. to the rooted 
% mean square error of the partition of unity (POU). Uses Gaussian kernel.
% 
% Inputs
% X           D x N matrix containing normalized columns; each column is
%             one sampling/training direction
% vali        encodes the directions on which the POU will be evaluated
%              if integer -> equal area partition with random perturbation
%              if matrix of size(vali,1)=D -> use the columns as directions
%              if string, load directions from textfile via dlmread
% gamma_min   minimum gamma to be considered
% gamma_max   maximum gamma to be considered
% ngamma      number of gammas to investigate
% partition_case    if 'U': partition of unity, i.e. the objective is to
%                   optimize gamma s.t. the (scalar) one function is
%                   approximated best
%                   if 'I': partition of identity, i.e. the objective is to
%                   optimize gamma s.t. the (vectorial) identity function
%                   (i.e. mapping directions onto themselves) is
%                   approximated best
% xitrunc     [OPTIONAL] max. angle at which the kernel is truncated
% do_plot     [OPTIONAL] if true, then plot
%
% Outputs
% gamma_best        gamma minimizing the average POU error
% RMSE              rooted mean square error (mean w.r.t. all validation
%                   directions)
% V                 the validation directions ( D x Nvali )
%
% See also POU, POI, GEODESICDISTANCE
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% SearchGamma.m
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

disp('Entering SearchGamma')
addpath(genpath('../eq_sphere_partitions'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% determine the size of the input data
N = size(X,2);
D = size(X,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% choose or generate validation directions
if ischar(vali)
    try
        Y=dlmread(vali);
    catch
        msgbox(['ERROR in SearchGamma: could not read file ',vali]);
        return
    end
    V=RenormalizeColumns(Y);
    Nvali = size(Y,2);
    disp([' using ' num2str(Nvali) ' validaton directions from file ' vali])
elseif length(vali)==1
    disp([' using ' num2str(vali) ' random validation directions'])
    Nvali = vali;
%     V = RenormalizeColumns( eq_point_set(D-1,Nvali) + 0.01* RenormalizeColumns(randn( D, Nvali )) );
    V = RandomDirections(D,Nvali);
elseif size(vali,1)==D
    V = vali;
    Nvali = size(V,2);
    disp([' using ' num2str(Nvali) ' directly provided validation directions'])
else
    error(['argument vali (inputname: ' inputname(2) ') does not meet requirements'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% geodesic distance between the validation data and input data
% size(xi_input_vali) = [N   Nvali]

xi_input_vali = real( acos( min( 1, max( -1, X'*V) ) ) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute distance matrix for the provided point set
[xi_input, NN_input] = GeodesicDistance( X );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% output some information on the input data, and do input checks
fprintf('###################################################################\n');
fprintf('# information on the provided sampling directions:\n');
fprintf('# min nearest neighbor distance ............ %8.5f rad\n', min(NN_input) );
fprintf('# mean nearest neighbor distance ........... %8.5f rad\n', mean(NN_input) );
fprintf('# max nearest neighbor distance ............ %8.5f rad\n', max(NN_input) );
fprintf('# std. dev. of nearest neighbor distance ... %8.5f rad\n', sqrt(var(NN_input)) );
fprintf('# normalized third moment .................. %8.5f\n',  mean( (NN_input-mean(NN_input)).^3)/(var(NN_input).^1.5) );
fprintf('###################################################################\n');
fprintf('# try finding gamma minimizing the\n# average ')
if strcmp(partition_case,'I')
    fprintf('POI');
elseif strcmp(partition_case,'U')
    fprintf('POU');
else
    msgbox('ERROR in SearchGamma: partition_case must be either ''I'' or ''U''')
end
fprintf(' error over the validation set\n');
fprintf('#\n# search parameters:\n');
fprintf('# gamma_min ................................ %8.5f\n', gamma_min );
fprintf('# gamma_min ................................ %8.5f\n', gamma_max );
fprintf('# n_gamma .................................. %8i\n', ngamma );
fprintf('###################################################################\n');
fprintf('#  gamma        RMSE\n');
RMSE_best = -1;
gamma_best = gamma_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% truncation
if exist('xitrunc','var') && ~isempty(xitrunc)
    ximax = xitrunc;
    ndel = 0;
    for n=1:N
        for j=1:Nvali
            if( xi_input_vali(n,j) > ximax )
                xi_input_vali(n,j) = 100*pi;
                ndel = ndel+1;
            end
        end
    end
    fprintf('# remaining fraction of validation points: %10.5f\n', 1 - ndel/N/Nvali);
    pattern = ones(N,N);
    for n1=1:N
        for n2=1:N
            if( xi_input(n1,n2) > ximax )
                xi_input(n1,n2) = 100*pi;
                xi_input(n2,n1) = 100*pi;
                pattern(n1,n2) = 0;
                pattern(n2,n1) = 0;
            end
        end
    end
    fprintf('# sparse kernel with xi_max = %10.5f\n', ximax);
    fprintf('# remaining: %10.5f\n', sum(sum(pattern))/N^2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop for finding gamma
RMSE_all=zeros(ngamma,1);
i_gamma=1;
gamma_with_warning = [];
lastwarn('');            % empty lastwarn s.t. new warnings can be detected
gammaspace = linspace( gamma_min, gamma_max, ngamma );
for gamma = gammaspace
% compute RMSE. case consistency of partition_case has been checked above
    if strcmp(partition_case,'U')
        RMSE = POU(xi_input, xi_input_vali, gamma);
    else
        RMSE = POI(xi_input, xi_input_vali, gamma, X, V);
    end
    if ~strcmp(lastwarn,'')
        gamma_with_warning(end+1) = gamma;  % track gammas which issued warnings
        lastwarn('');
    end
    RMSE_all(i_gamma)  = RMSE;
    if ( RMSE < RMSE_best || RMSE_best < 0 )   % RMSE_best was initialized  as -1
        gamma_best  = gamma;
        RMSE_best    = RMSE;
    end
    fprintf('%10.5f    %16.8e\n', gamma, RMSE );
    i_gamma=i_gamma+1;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% outputs
fprintf('best fit: %10.6f, err: %10.5e\n', gamma_best, RMSE_best );
% figure;
% title('rel. LSQ POU/POI error\n');

if exist('do_plot','var') && do_plot==true
    semilogy( linspace(gamma_min,gamma_max,ngamma), RMSE_all, 'linewidth', 3 )
    hold on
    % plot gammas with warnings
    plot(gamma_with_warning,ones(size(gamma_with_warning))*min(RMSE_all)/5,'r.')
    grid minor
    xlabel('gamma')
    ylabel('RMSE')
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% development / debugging code
% show all gammas at which a warning occured:
% format longeng
% gamma_with_warning
% format short

% write gamma vs error plot to gnuplottable txt file:
% dlmwrite('some_file_name.txt', [gammaspace' RMSE_all], 'delimiter', '\t', ...
%     'precision', '%28.20e');

