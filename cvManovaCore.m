function D = cvManovaCore(vi, varargin)

% move separation into runs into loadDataSPM
% include here precompute, existing initialization, and bias correction

% assumes whitened and filtered data and design matrices
% misc.n from mean size(Ys, 1)

% cvManovaCore([], Ys, Xs, Cs, fE, permute = false, lambda = 0);
% D = cvManovaCore(vi, ...);

% check consistent dimensions between Ys, Xs
% check that Cs does not exceed minimum number of regressors


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
    n = mean(cellfun(@(x) size(x, 1), Ys));

    % check contrasts
    for ci = 1 : nContrasts
        if size(Cs{ci}, 2) > rank(Cs{ci})
            error('contrast %d is misspecified!', ci)
        end
        for ri = 1 : m
            if inestimability(Cs{ci}, Xs{ri}) > 1e-6
                error('contrast %d is not estimable in run %d!', ci, ri)
            end
        end
    end
    
    % estimate GLM parameters and errors, and prepare design inner products
    betas = cell(m, 1);
    xis = cell(m, 1);
    XXs = cell(m, 1);
    for ri = 1 : m
        betas{ri} = pinv(Xs{ri}) * Ys{ri};
        xis{ri} = Ys{ri} - Xs{ri} * betas{ri};
        XXs{ri} = Xs{ri}' * Xs{ri};
    end
    
    % prepare contrast projectors
    nContrasts = numel(Cs);
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
    
    % answer initialization call
    D = nan(1, nContrasts * nPerms);
    return
end


%% computation

% precompute per-run E
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


D = zeros(nContrasts, nPerms);
% for each contrast
for ci = 1 : nContrasts
    % number of regressors involved in contrast
    nConReg = size(CCs{ci}, 1);
    
    % precompute per-run betaDelta
    betaDelta = cell(m, 1);
    for k = 1 : m
        betaDelta{k} = CCs{ci} * betas{k}(1 : nConReg, vi);
    end
    
    % precompute per-run H
    Hs = cell(m, m);
    for k = 1 : m
        for l = 1 : m
            if l == k, continue, end
            
            Hs{k, l} = betaDelta{k}' * XXs{l}(1 : nConReg, 1 : nConReg) * betaDelta{l};
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
            
            % sum across cross-validation folds
            D(ci, pi) = D(ci, pi) + Dl;
        end
    end
end

% mean across folds, bias correction
p = numel(vi);
D = ((m - 1) * fE - p - 1) / ((m - 1) * n) * D / m;

% return row vector
D = reshape(D, 1, []);


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

