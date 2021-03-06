% scratch function!
function Y = RotateToPoles( X, X_eq, sym_flag )
%% Rotate any point set s.t. its some point coincides with
% that of the corresponding EQ point set.
% ATTENTION: NO FAILSAFE CHECKS YET!
%
% IN:
% X         arbitrary point set
% X_eq      EQ point set
% sym_flag  0: take X as it is. 1: symmetrize X by [X,-X]
%
% OUT:
% Y         rotated point set

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% RotateToPoles.m
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

s = 0;
if( nargin > 2)
    s = sym_flag;
end
if( s ~= 0 )
    X = [ X, -X ];
end
[~, imax]=max(X(3,:));
n = [0;0;1] + X(:,imax);
n = n/norm(n);
P = eye(3)-2*n*n';
Y = P*X;

% now find phi such that Y rotate around [0;0;1] by phi is closest to X
% F=@(phi) norm( [ cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0,0,1] * Y-X_eq, 'fro');
xmin = inf;
F=[];
for phi=0:pi/200:pi/2
    [~,nn] = Distance( [[ cos(phi), sin(phi), 0; -sin(phi), cos(phi), 0; 0,0,1] * Y, X_eq]);
    f=norm(nn);
    F=[F;f];
    if( f < xmin )
        xmin = f;
        phimin = phi;
    end
end
Y= [ cos(phimin), sin(phimin), 0; -sin(phimin), cos(phimin), 0; 0,0,1] * Y;
[~,nn] = Distance( [Y, X_eq]);
f=norm(nn);
% F
% max(F)
% min(F)
% xmin
