function cvManovaSearchlight(dirName, slRadius, Cs, permute)

% cross-validated MANOVA on searchlight
%
% cvManovaSearchlight(dirName, slRadius, Cs, permute = false)
%
% dirName:  directory where the SPM.mat file referring to an estimated
%           model is located
% slRadius: radius of the searchlight sphere, in voxels
% Cs:       cell array of contrast matrices
% permute:  whether to compute permutation values of the test statistic
%
% Output files are written to the same directory:
% spmD_C####_P####.nii:
%           images of the pattern discriminability D
%           contrast and permutation are identified by numbers
% VPSL.img  an image of the number of voxels within each searchlight
% cms.mat   a record of the analysis parameters
%
% Copyright (C) 2013 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.


fprintf('\n\ncvManovaSearchlight\n\n')

if nargin < 4
    permute = false;
end

nContrasts = numel(Cs);
outpattern = 'spmD_C%04d_P%04d.nii';

if dirName(end) ~= '/', dirName = [dirName '/']; end

% simplify saving images
wd = cd;
cd(dirName)

% load data, design matrix etc.
[Y, X, mask, misc] = loadDataSPM(dirName);

% determine per-run design and data matrices
nRuns = misc.m;
Xrun = cell(nRuns, 1);
Yrun = cell(nRuns, 1);
% for each run
for ri = 1 : nRuns
    Yrun{ri} = Y(misc.sRow{ri}, :);
    Xrun{ri} = X(misc.sRow{ri}, misc.sCol{ri});
end
clear Y X

% check contrasts
for ci = 1 : nContrasts
    if size(Cs{ci}, 2) > rank(Cs{ci})
        error('contrast %d is misspecified!', ci)
    end
    for ri = 1 : nRuns
        if inestimability(Cs{ci}, Xrun{ri}) > 1e-6
            error('contrast %d is not estimable in run %d!', ci, ri)
        end
    end
end
    
% precomputation
fprintf('\nprecomputing GLM runwise\n')
[XXs, betas, xis] = cvManova_precompute(Xrun, Yrun);
clear Xrun Yrun

% determine voxels per searchlight, and save as image
fprintf('\ncomputing voxels per searchlight image\n')
p = runSearchlight(slRadius, mask, @(vi)(size(vi, 1)));
saveMRImage(reshape(p, size(mask)), 'VPSL.nii', misc.v2mm, ...
    'voxels per searchlight');

% error check
pMax = max(p(:));
if pMax > 0.9 * (misc.m - 1) * misc.fE
    error('insufficient amount of data for searchlight size %d!', pMax)
end

% bias correction factor
factor = ((misc.m - 1) * misc.fE - p - 1) / ((misc.m - 1) * misc.n);
clear p

% run searchlight
fprintf('\ncomputing cross-validated MANOVA on searchlight\n')
mDl = runSearchlight(slRadius, mask, @cvManova_compute, ...
    XXs, betas, xis, Cs, permute);

% separate contrast and permutation dimensions
nContrasts = numel(Cs);
nPerms = size(mDl, 2) / nContrasts;
mDl = reshape(mDl, [], nContrasts, nPerms);

% compute the unbiased estimate of the pattern discriminability D
D = bsxfun(@times, factor, mDl);
clear mDl

% save results
for ci = 1 : nContrasts
    for pi = 1 : nPerms
        fn = sprintf(outpattern, ci, pi);
        saveMRImage(reshape(D(:, ci, pi), size(mask)), fn, misc.v2mm, ...
            'pattern discriminability');
    end
end

% analysis parameters
save cms.mat slRadius Cs permute misc nPerms

cd(wd)