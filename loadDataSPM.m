function [Y, X, mask, misc] = loadDataSPM(dirName)

% load fMRI data via SPM.mat
%
% [Y, X, mask, misc] = loadDataSPM(dirName)
%
% dirName:  name of directory that contains SPM.mat
% Y:        MR data (within mask), samples x voxels
% X:        design matrix, samples x regressors
% mask:     brain mask, logical 3d-volume
% misc:     struct with additional data:
%     v2mm  voxels to mm transformation matrix
%     sRow  rows for each session
%     sCol  columns for each session
%     m     number of sessions
%     n     number of images per session
%     fE    residual degrees of freedom per session
%
% Y & X and are high-pass filtered and whitened.
% Y includes only those voxels selected through mask.
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

SPMname = [dirName 'SPM.mat'];

fprintf('loading data\n')
fprintf(' via %s\n', SPMname)

% load SPM.mat
load(SPMname, 'SPM');

% recreate volumes, for backwards compatibility (with SPM2?)
fnames = {SPM.xY.VY.fname};
if ~exist(fnames{1}, 'file')
    % and to enable reading moved data files
    fnames = patchPath(fnames, dirName);
end
fprintf(' reading volume information\n')
VY = cell2mat(spm_vol(fnames));
% but copy scaling information
[VY(:).pinfo] = SPM.xY.VY(:).pinfo;
nImages = numel(VY);
nVoxels = prod(VY(1).dim);

% check dimensions and orientations of all images
spm_check_orientations(VY);

% read mask image
if isfield(SPM, 'VM')
    VM = spm_vol([dirName SPM.VM.fname]);
else
    VM = spm_vol([dirName 'mask.img']);
end
mask = logical(spm_read_vols(VM));
nMaskVoxels = sum(mask(:));

% check memory
memNeed = nMaskVoxels * nImages * 8 / 1024 / 1024;  % assuming double precision
memFree = systemFree / 1024;
fprintf(' memory needed  %7.1f MiB\n', memNeed)
fprintf('        free    %7.1f MiB\n', memFree)
if memFree < memNeed, warning('not enough memory!'), end

% read and mask data
Y = nan(nImages, nMaskVoxels);
fprintf(' reading images\n')
for i = 1 : nImages             % slow, but saves memory
    y = spm_read_vols(VY(i));
    Y(i, :) = y(mask)';                             % apply mask
    if mod(i, 50) == 0
        fprintf(' %d of %d images read\n', i, nImages)
    end
end
fprintf('\n')

% get whitening matrix
if isfield(SPM.xX, 'W')
    W = SPM.xX.W;
else
    fprintf(' * SPM.mat does not define whitening matrix!\n')
    W = 1;
end

% whiten and high-pass filter data
fprintf(' filtering and whitening data\n')
Y = spm_filter(SPM.xX.K, W * Y);

% get whitened and high-pass filtered design matrix
if isfield(SPM.xX, 'xKXs')
    X = SPM.xX.xKXs.X;
else
    fprintf(' * SPM.mat does not contain whitened and high-pass filtered design matrix!\n')
    X = spm_filter(SPM.xX.K, W * SPM.xX.X);
end

% degrees of freedom
Tdf = sum(SPM.nscan);                               % total
Kdf = sum(arrayfun(@(x)(size(x.X0, 2)), SPM.xX.K)); % loss from high-pass filter
Xdf = rank(X);                                      % loss from regressors
Rdf = Tdf - Kdf - Xdf;                              % residual
fprintf(' df: %d - %d - %d = %d\n', Tdf, Kdf, Xdf, Rdf);

% output
misc.v2mm = VY(1).mat;
misc.m = size(SPM.nscan, 2);
misc.sRow = {SPM.Sess.row};
misc.sCol = {SPM.Sess.col};
for si = 1 : misc.m                                 % add constant regressors
    misc.sCol{si} = [misc.sCol{si}, SPM.xX.iB(si)];
    % NOT SURE WHETHER THIS IS GENERALLY CORRECT
end
misc.fE = Rdf / misc.m;                             % if not consistent across sessions,
misc.n = Tdf / misc.m;                              % then this is an approximation
