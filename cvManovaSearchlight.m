function cvManovaSearchlight(dirName, slRadius, Cs, permute, lambda)

% cross-validated MANOVA on searchlight
%
% cvManovaSearchlight(dirName, slRadius, Cs, permute = false, lambda = 0)
%
% dirName:  directory where the SPM.mat file referring to an estimated
%           model is located
% slRadius: radius of the searchlight sphere, in voxels
% Cs:       cell array of contrast matrices
% permute:  whether to compute permutation values of the test statistic
% lambda:   regularization parameter (0–1)
%
% Output files are written to the same directory:
% spmD_C####_P####.nii:
%           images of the pattern discriminability D
%           contrast and permutation are identified by numbers
% VPSL.img  an image of the number of voxels within each searchlight
% cms.mat   a record of the analysis parameters
%
%
% Copyright (C) 2013-2016 Carsten Allefeld


fprintf('\n\ncvManovaSearchlight\n\n')

if nargin < 4
    permute = false;
end
if nargin < 5
    lambda = 0;
end

outpattern = 'spmD_C%04d_P%04d.nii';

if dirName(end) ~= filesep
    dirName = [dirName filesep];
end

% load data, design matrix etc.
[Ys, Xs, mask, misc] = loadDataSPM(dirName);

% for checkpointing, compute unique ID that encodes parameters:
%   SPM.mat & referenced data -> <timestamp of SPM.mat>,
%   slRadius, Cs, permute, lambda
% encode as string
if ~exist('gencode.m', 'file')
    addpath([spm('dir') filesep 'matlabbatch'])
end
uid = gencode({['SPM.mat of ' getfield(dir([dirName 'SPM.mat']), 'date')], ...
    slRadius, Cs, permute, lambda});
uid = sprintf('%s\n', uid{:});
% compute Fletcher-16 checksum
uid = dec2hex(fletcher16(uid), 4);

% run searchlight
fprintf('\ncomputing cross-validated MANOVA on searchlight\n')
[D, p] = runSearchlight(['cmsCheckpoint' uid '.mat'], slRadius, mask, ...
    @cvManovaCore, Ys, Xs, Cs, misc.fE, permute, lambda);
 
% separate contrast and permutation dimensions
nContrasts = numel(Cs);
nPerms = size(D, 2) / nContrasts;
D = reshape(D, [], nContrasts, nPerms);

% simplify saving results by changing to directory
wd = cd(dirName);
% ensure change back on exit
cleanupObj = onCleanup(@() cd(wd));     

% save results
for ci = 1 : nContrasts
    for pi = 1 : nPerms
        fn = sprintf(outpattern, ci, pi);
        spmWriteImage(reshape(D(:, ci, pi), size(mask)), fn, misc.mat, ...
            'descrip', 'pattern discriminability')
    end
end

% save voxels per searchlight as image
spmWriteImage(reshape(p, size(mask)), 'VPSL.nii', misc.mat, ...
    'descrip', 'voxels per searchlight')

% save analysis parameters
save cmsParameters.mat slRadius Cs permute misc nPerms


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

