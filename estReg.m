function lambda = estReg(xis, fE)

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

