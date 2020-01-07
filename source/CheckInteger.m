function [ x ] = CheckInteger( x, lower_bound )
%CHECKINTEGER Auxiliary function
%   Checks whether the input argument x is
%     1) of integer type (and, first, of numerical type)
%     2) greater than or equal to lower_bound
% returns x in any case, but may open a dialogue box, informing the user
% that the check was unsuccessful.
%
% Inputs
% x                 scalar
% lower_bound       scalar
%
% Outputs
% x                 unmodified input
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% CheckInteger.m
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

if ~isnumeric(lower_bound)
    msgbox(['WARNING: lower_bound in CheckInteger is not of numerical type, but is ', ...
        class(lower_bound), ' instead']);
end
if ~isnumeric(x)
    msgbox(['WARNING: input is not of numerical type, but is ', ...
        class(x), ' instead']);
elseif isnan(x) || isinf(x) || x < lower_bound || floor(x)~=x
    msgbox(['WARNING: input is ', num2str(x), ...
        ' but should be an integer >= ', num2str(lower_bound)]);
end

end

