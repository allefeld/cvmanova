function lambda = estReg(xis, fE)

% lambda = estReg(xis, fE)
%
% Implements the method of
% J. Schäfer and K. Strimmer, A shrinkage approach to large-scale covariance
% estimation and implications for functional genomics, Statistical Applications
% in Genetics and Molecular Biology, vol. 4, no. 1, 2005.
%
%
% Copyright (C) 2015 Carsten Allefeld


fprintf('\nestimating optimal regularization parameter lambda\n')

nRuns = numel(xis);
df = (nRuns - 1) * fE;

lambda = nan(nRuns, 1);
for ri = 1 : nRuns
    
    x = cat(1, xis{(1 : nRuns) ~= ri});

    [n, p] = size(x);
    
    % the empirical covariance matrix and its variance
    
    % inelegant and time-consuming, but necessary to allow larger n*N^2
    S = zeros(p, p);
    VarS = zeros(p, p);
    for i = 1 : p
        for j = i : p
            w = x(:, i) .* x(:, j);
            mw = sum(w) / n;
            S(i, j) = mw * n / df;
            VarS(i, j) = sum((w - mw) .^ 2) * n / df .^ 3;
            S(j, i) = S(i, j);
            VarS(j, i) = VarS(i, j);
        end
    end
    
    % helper quantities
    sumVarS = sum(VarS(:));
    sumDiagVarS = sum(diag(VarS));
    sumOffdiagVarS = sumVarS - sumDiagVarS;
    sumSS = sum(S(:) .^ 2);
    sumDiagSS = sum(diag(S .^ 2));
    sumOffdiagSS = sumSS - sumDiagSS;
    
    % shrinkage intensity for target "diagonal, unequal variance"
    l = sumOffdiagVarS / sumOffdiagSS;
    lambda(ri) = max(0, min(1, l));
    
end

fprintf('  %d voxels, %d runs, %g df per run\n', p, nRuns, fE)

fprintf('  simple approximation: %g\n', p / (p + df));

fprintf('  estimated from data per cv fold:\n  ')
fprintf('  %g', lambda)
fprintf('\n')

lambda = mean(lambda);
fprintf('  used value (mean across cv folds): %g\n', lambda);
fprintf('\n')


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

