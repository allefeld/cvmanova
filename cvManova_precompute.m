function [XXs, betas, xis] = cvManova_precompute(Xrun, Yrun)

% cross-validated MANOVA, GLM precomputation
%
% [XXs, betas, xis] = cvManova_precompute(Xrun, Yrun)
%
% Xrun:     cell array of per-run design matrices X
% Yrun:     cell array of per-run data matrices Y
% XXs:      cell array of per-run X' X
% betas:    cell array of per-run GLM parameter estimates
% xis:      cell array of per-run GLM residuals
%
%
% Copyright (C) 2013 Carsten Allefeld


Xrun = Xrun(:);
Yrun = Yrun(:);
nRuns = size(Xrun, 1);

betas = cell(nRuns, 1);
xis = cell(nRuns, 1);
XXs = cell(nRuns, 1);

% for each run
for ri = 1 : nRuns
    % estimate GLM
    beta = pinv(Xrun{ri}) * Yrun{ri};
    xi = Yrun{ri} - Xrun{ri} * beta;
    % store precomputed matrices, possibly truncating
    betas{ri} = beta;
    xis{ri} = xi;
    XXs{ri} = Xrun{ri}' * Xrun{ri};
end


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

