function PlotEDF(X, sym_flag, ax)
%PLOTEDF plots the empirical cumulative distribution function of a given
% set of points X (size D x N, where D is the Euclidean dimension and N is
% is the number of points) to the provided axes ax, considering the
% antipodes _explicitly_ if sym_flag==1.
%
% See also PLOTPOINTSONSPHERE, DISTANCE, ECDF
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPYRIGHT NOTES
%
% PlotEDF.m
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
% Advances in Computational Mathematics, Number/Volume, p. XX-YY, 2019
% DOI   10.1007/s10444-019-09726-5
% URL   doi.org/10.1007/s10444-019-09726-5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preparation
% compute Euclidean distances
d = Distance(X,sym_flag);
% plot vertical line with LaTeX label "\sqrt{2}"
hold(ax,'on')
plot(ax,[sqrt(2),sqrt(2)],[0,1],'k-.')
text(1.4142,0.95,'x=$\sqrt{2}\rightarrow$','HorizontalAlignment','right','Interpreter','latex')

%% plot EDF's
xlim([0,2])
for i_plot = 1:size(X,2)
    ecdf(ax, d(:,i_plot));
end
grid(ax,'on')
xlabel('x [-] (Euclidean distance)')
ylabel('P( |X-Y| < x ) [-]')
hold(ax,'off')

end
