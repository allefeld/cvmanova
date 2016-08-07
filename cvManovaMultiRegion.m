function [Ds, ps] = cvManovaMultiRegion(dirName, region, Cs, lambda, permute, highpass)

% cross-validated MANOVA on region
%
% [D, p] = cvManovaRegion(dirName, region, Cs, lambda, permute = false)
%
% dirName:  directory where the SPM.mat file referring to an estimated
%           model is located
% region:   region mask, logical 3d-volume
% Cs:       cell array of contrast matrices
% lambda:   regularization parameter (0–1)
% permute:  whether to compute permutation values of the test statistic
% D:        pattern distinctness, contrasts x permutations
% p:        number of voxels in the region
%
% If region is empty, all in-mask voxels are used.
% If lambda is not specified or empty, its optimal value is estimated using
% the method of
%   J. Schäfer and K. Strimmer, A shrinkage approach to large-scale covariance
%   estimation and implications for functional genomics, Statistical Applications
%   in Genetics and Molecular Biology, vol. 4, no. 1, 2005.
% implemented in estReg.
%
% Copyright (C) 2015 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.


fprintf('\n\ncvManovaRegion\n\n')

if nargin < 5
    permute = false;
end
if nargin < 4
    lambda = [];
end    

nContrasts = numel(Cs);

if dirName(end) ~= filesep, dirName = [dirName filesep]; end


% load data, design matrix etc.
[Ys, X, region, misc] = loadDataSPM_Region(dirName, region,1,highpass);

% determine per-run design and data matrices
nRuns = misc.m;
nRegion = numel(region);
Xrun = cell(nRuns,1);
Yrun = cell(nRuns,nRegion);
% for each run
for nr = 1:nRegion
    for ri = 1 : nRuns    
        Yrun{ri,nr} = Ys{nr}(misc.sRow{ri}, :);
    end;
    Ys(nr) = cell(1,1);

end
for ri = 1 : nRuns
    Xrun{ri} = X(misc.sRow{ri}, misc.sCol{ri});
end;
clear X Ys

% check contrasts
for ci = 1 : nContrasts
   
    if size(Cs{ci}, 2) > rank(Cs{ci})
        keyboard
        error('contrast %d is misspecified!', ci)
    end
    for ri = 1 : nRuns
        if inestimability(Cs{ci}, Xrun{ri}) > 1e-6
            error('contrast %d is not estimable in run %d!', ci, ri)
        end
    end
end
    

lambda_input = lambda;
nContrasts = numel(Cs);
Ds = zeros(nContrasts,nRegion);
for nr = 1:nRegion    
    % precomputation
    fprintf('\nprecomputing GLM runwise\n')
    [XXs, betas, xis] = cvManova_precompute(Xrun, Yrun(:,nr));
    Yrun(:,nr) = cell(nRuns,1);
    
    % estimate regularization parameter
    if isempty(lambda_input)
        lambda = estReg(xis, misc.fE);        
    end
    
    % bias correction factor
    p = sum(region{nr}(:));
    factor = ((misc.m - 1) * misc.fE - p - 1) / ((misc.m - 1) * misc.n);
    
    % run searchlight
    fprintf('\ncomputing cross-validated MANOVA on region\n')
    cvManova_compute(nan(0, 1), XXs, betas, xis, Cs, permute, lambda);
    mDl = cvManova_compute(1 : p, XXs, betas, xis, Cs, permute, lambda);
    
    % separate contrast and permutation dimensions
    
    nPerms = numel(mDl) / nContrasts;
    mDl = reshape(mDl, nContrasts, nPerms);
    
    % compute the unbiased estimate of the pattern discriminability D
    D = factor * mDl;
    Ds(1:nContrasts,nr) = D;
    ps(nr) = p;
end;
