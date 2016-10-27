function D = cvManovaCore(vi, varargin)

% cross-validated MANOVA
% core implementation of the method proposed by Allefeld and Haynes (2014)
%
% cvManovaCore([], Ys, Xs, Cs, fE, permute = false, lambda = 0);
%
% Ys:       cell array of per-session data matrices
% Xs:       cell array of per-session design matrices
% Cs:       cell array of contrast vectors or matrices
% fE:       error degrees of freedom
% permute:  whether to compute permutation values
% lambda:   regularization parameter (0–1)
%
% D = cvManovaCore(vi, ...);
%
% vi:       voxels, indices into columns of Ys{:}
% D:        pattern distinctness, contrasts × permutations as a row vector
%
% The function has two call syntaxes. The first form initializes internal
% (persistent) variables based on the parameters, and has to be used once
% for each data set to be analyzed. The second form performs the actual
% analysis on the specified voxels, and can be repeated for each set of
% voxels of interest. This two-step approach provides improved performance
% for searchlight analysis. To release the memory occupied by the internal
% variables, use 'clear cvManovaCore'.
%
% It is assumed that the data and design matrices have been whitened and
% possibly filtered. fE is the residual number of degrees of freedom,
% i.e. the number of scans per session minus the rank of the design matrix
% and minus further loss of dfs due to filtering. If this number varies
% across sessions, use its mean.
%
% See also cvManovaSearchlight, cvManovaRegion
%
%
% Copyright (C) 2016 Carsten Allefeld

% check consistent dimensions between Ys, Xs
% check that Cs does not exceed minimum number of regressors
% rank deficiency should not be error


persistent m n betas xis XXs nContrasts CCs nPerms sp fE lambda


%% initialization

if isempty(vi)
    % extract arguments
    Ys = varargin{1};
    Xs = varargin{2};
    Cs = varargin{3};
    fE = varargin{4};
    if numel(varargin) >= 5
        permute = varargin{5};
    else
        permute = false;
    end
    if numel(varargin) >= 6
        lambda = varargin{6};
    else
        lambda = 0;
    end
    m = numel(Ys);
    n = cellfun(@(x) size(x, 1), Ys);
    nContrasts = numel(Cs);

    
    % check input
    assert(isequal(cellfun(@(x) size(x, 1), Ys), cellfun(@(x) size(x, 1), Xs)), ...
        'inconsistent number of scans between data and design!')
    assert(all(diff(cellfun(@(x) size(x, 2), Ys)) == 0), ...
        'inconsistent number of voxels within data!')

    % check contrasts
    qMin = min(cellfun(@(x) size(x, 2), Xs));
    for ci = 1 : nContrasts
        qC = find(all(Cs{ci} == 0, 2) == false, 1, 'last');
        Cs{ci} = Cs{ci}(1 : qC, :); % trim trailing all-zero rows
        assert(qC <= qMin, ...
            'contrast %d exceeds the %d common regressors!', ci, qMin)
        for si = 1 : m
            assert(inestimability(Cs{ci}, Xs{si}) <= 1e-6, ...
                'contrast %d is not estimable in session %d!', ci, si)
        end
    end
    
    % estimate GLM parameters and errors, and prepare design inner products
    betas = cell(m, 1);
    xis = cell(m, 1);
    XXs = cell(m, 1);
    for si = 1 : m
        betas{si} = pinv(Xs{si}) * Ys{si};
        xis{si} = Ys{si} - Xs{si} * betas{si};
        XXs{si} = Xs{si}' * Xs{si};
    end
    
    % prepare contrast projectors
    CCs = cell(nContrasts, 1);
    for ci = 1 : nContrasts
        CCs{ci} = pinv(Cs{ci}') * Cs{ci}';
    end
    
    % generate sign permutations
    sp = signPermutations(m);
    nPerms = size(sp, 2) / 2;       % the two halves are equivalent
    if ~permute
        nPerms = 1;                 % neutral permutation only
    end
    sp = sp(:, 1 : nPerms);
    
    % answer initialization call by runSearchlight
    D = nan(1, nContrasts * nPerms);
    return
end


%% computation

% precompute per-session E
Es = cell(m, 1);
for k = 1 : m
    x = xis{k}(:, vi);          % weirdly, this is faster than
    Es{k} = x' * x;             % Es{k} = xis{k}(:, vi)' * xis{k}(:, vi);
end
clear x

% precompute inverse of per-fold summed E
iEls = cell(m, 1);
for l = 1 : m
    ks = 1 : m;
    ks(l) = [];
    
    El = sum(cat(3, Es{ks}), 3);

    % shrinkage regularization towards diagonal
    El = (1 - lambda) * El + lambda * diag(diag(El));
    
    iEls{l} = eye(size(El)) / El;   % faster than inv
end
clear Es


p = numel(vi);
D = zeros(nContrasts, nPerms);
% for each contrast
for ci = 1 : nContrasts
    % number of regressors involved in contrast
    qCC = size(CCs{ci}, 1);
    
    % precompute per-session betaDelta
    betaDelta = cell(m, 1);
    for k = 1 : m
        betaDelta{k} = CCs{ci} * betas{k}(1 : qCC, vi);
    end
    
    % precompute per-session H
    Hs = cell(m, m);
    for k = 1 : m
        for l = 1 : m
            if l == k, continue, end
            
            Hs{k, l} = betaDelta{k}' * XXs{l}(1 : qCC, 1 : qCC) * betaDelta{l};
        end
    end
    clear betaDelta
    
    % for each permutation
    for pi = 1 : nPerms
        
        % for each cross-validation fold
        for l = 1 : m
            ks = 1 : m;
            ks(l) = [];
            
            % sign-permuted, summed H
            Hl = sum(bsxfun(@times, ...
                cat(3, Hs{ks, l}), ...
                    reshape(sp(ks, pi) * sp(l, pi)', 1, 1, m - 1) ...
                ), 3);
            
            % fold-wise D
            Dl = sum(sum(Hl' .* iEls{l}));  % faster than trace(Hl * iEls{l})
            
            % bias correction (fold-specific!)
            Dl = (sum(fE(ks)) - p - 1) / sum(n(ks)) * Dl;
            
            % sum across cross-validation folds
            D(ci, pi) = D(ci, pi) + Dl;
        end
    end
end

% mean across folds
D = D / m;

% return row vector
D = reshape(D, 1, []);


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
