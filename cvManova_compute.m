function mDl = cvManova_compute(vi, XXs, betas, xis, Cs, permute, lambda)

% cross-validated MANOVA, main computation
%
% mDl = cvManova_compute(vi, XXs, betas, xis, Cs, permute, lambda)
%
% vi:       indices of voxels to use
% XXs:      cell array of per-run precomputed X' X
% betas:    cell array of per-run precomputed GLM parameter estimates
% xis:      cell array of per-run precomputed GLM residuals
% Cs:       cell array of contrast matrices
% permute:  whether to compute permutations
% lambda:   regularization parameter (0–1)
% mDl:      cross-validated MANOVA statistic, raw version
%
% Copyright (C) 2013-2014 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

nRuns = size(betas, 1);


persistent sp                       % sign permutations
persistent CCs                      % canonical contrast matrices

if isempty(vi)
    % initialization call
    
    % generate sign permutations
    sp = signPermutations(nRuns);
    nPerms = size(sp, 2) / 2;       % the two halves are equivalent
    if ~permute
        nPerms = 1;                 % identity permutation only
    end
    sp = sp(:, 1 : nPerms);

    % prepare canonical contrast matrices
    nContrasts = numel(Cs);
    CCs = cell(nContrasts, 1);
    for ci = 1 : nContrasts
        % compute canonical contrast matrix
        CCs{ci} = pinv(Cs{ci}') * Cs{ci}';
    end
    
    % answer initialization call
    mDl = nan(1, nContrasts * nPerms);
    return
else
    nPerms = size(sp, 2);
    nContrasts = size(CCs, 1);
end

if nargin < 7, lambda = 0; end

% pre-compute per-run E
Es = cell(nRuns, 1);
for k = 1 : nRuns
    x = xis{k}(:, vi);          % weirdly, this is faster than
    Es{k} = x' * x;             % Es{k} = xis{k}(:, vi)' * xis{k}(:, vi);
end
clear x

% pre-compute inverse of per-fold summed E
iEls = cell(nRuns, 1);
for l = 1 : nRuns
    ks = 1 : nRuns;
    ks(l) = [];
    
    El = sum(cat(3, Es{ks}), 3);

    % shrinkage regularization towards diagonal
    El = (1 - lambda) * El + lambda * diag(diag(El));
    
    iEls{l} = eye(size(El)) / El;   % faster than inv
end
clear Es


mDl = zeros(nContrasts, nPerms);

% for each contrast
for ci = 1 : nContrasts
    % number of regressors involved in contrast
    nConReg = size(CCs{ci}, 1);
    
    % pre-compute per-run betaDelta
    betaDelta = cell(nRuns, 1);
    for k = 1 : nRuns
        betaDelta{k} = CCs{ci} * betas{k}(1 : nConReg, vi);
    end
    
    % pre-compute per-run H
    Hs = cell(nRuns, nRuns);
    for k = 1 : nRuns
        for l = 1 : nRuns
            if l == k, continue, end
            
            Hs{k, l} = betaDelta{k}' * XXs{l}(1 : nConReg, 1 : nConReg) * betaDelta{l};
        end
    end
    clear betaDelta
    
    % for each permutation
    for pi = 1 : nPerms
        
        % for each cross-validation fold
        for l = 1 : nRuns
            ks = 1 : nRuns;
            ks(l) = [];
            
            % sign-permuted, summed H
            Hl = sum(bsxfun(@times, ...
                cat(3, Hs{ks, l}), ...
                reshape(sp(ks, pi) * sp(l, pi)', 1, 1, nRuns - 1) ...
                ), 3);
            
            % fold-wise D
            Dl = sum(sum(Hl' .* iEls{l}));  % faster than trace(Hl * iEls{l})
            
            % sum across cross-validation folds
            mDl(ci, pi) = mDl(ci, pi) + Dl;
        end
    end
end

% sum -> mean
mDl = mDl / nRuns;

% return row vector
mDl = mDl(:)';

