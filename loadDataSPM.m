function [Y, X, mask, misc] = loadDataSPM(dirName, region)

% load fMRI data via SPM.mat
%
% [Y, X, mask, misc] = loadDataSPM(dirName, region = [])
%
% dirName:  name of directory that contains SPM.mat
% region:   optional additional region mask, logical 3D volume
% Y:        MR data (within mask), scans x voxels
% X:        design matrix, scans x regressors
% mask:     analysis brain mask, logical 3D volume; possibly combined with region mask
% misc:     struct with additional data:
%     mat   voxels to mm transformation matrix
%     sRow  rows for each session
%     sCol  columns for each session
%     m     number of sessions
%     n     number of images per session
%     fE    residual degrees of freedom per session
%
% Y & X and are high-pass filtered and whitened.
% Y includes only those voxels selected through mask.
%
%
% Copyright (C) 2013-2016 Carsten Allefeld


% default argument values
if nargin < 2
    region = [];
end

% load SPM.mat
SPMname = fullfile(dirName, 'SPM.mat');
fprintf('loading data\n')
fprintf(' via %s\n', SPMname)
load(SPMname, 'SPM');

% get data volumes and check matching voxel grid
VY = SPM.xY.VY;
nImages = numel(VY);
spm_check_orientations(VY);

% check whether data might have been moved
if ~exist(VY(1).fname, 'file')
    SPMold = fullfile(SPM.swd, 'SPM.mat');
    comLen = min(numel(SPMname), numel(SPMold));
    comPart = comLen - find(diff(...
        [SPMname(end - comLen + 1 : end) == SPMold(end - comLen + 1 : end) 1] ...
        ), 1, 'last');
    fprintf(2, ' if analysis and data folders were moved together, try\n');
    fprintf(2, '   spm_changepath(''%s'', ''%s'', ''%s'')\n', ...
        dirName, SPMold(1 : end - comPart), SPMname(1 : end - comPart));
end

% read analysis brain mask image
if isfield(SPM, 'VM')
    try
        mask = logical(spm_read_vols(SPM.VM));
    catch
        % SPM8 stores the filename without the path
        VM = spm_vol([dirName SPM.VM.fname]);
        mask = logical(spm_read_vols(VM));
    end
    fprintf(' using analysis brain mask from SPM.VM\n');
else
    mask = true(VY(1).dim);
    fprintf(' no analysis brain mask!\n')        % should not happen
end

% possibly apply region mask
if isempty(region)
    region = true(size(mask));
    fprintf(' no region mask!\n')
end
if all(size(region) == size(mask))
    mask = mask & logical(region);
else
    error('region mask doesn''t match\n')
end
fprintf(' %d voxels\n', sum(mask(:)));

% read and mask data
fprintf(' reading images\n')
pattern = SPM.xY.P(1, :);
pattern(~all(diff(SPM.xY.P) == 0)) = '?';
fprintf('  from %s\n', pattern)
[Y, mask] = spmReadVolsMasked(VY, mask);

% get design matrix
X = SPM.xX.X;

% whiten data and design matrix
if isfield(SPM.xX, 'W')
    fprintf(' whitening\n')
    W = SPM.xX.W;
    Y = W * Y;
    X = W * X;
else
    fprintf(' * SPM.mat does not define whitening matrix!\n')
end

% high-pass filter data and design matrix
fprintf(' high-pass-filtering\n')
Y = spm_filter(SPM.xX.K, Y);
X = spm_filter(SPM.xX.K, X);

% degrees of freedom
Tdf = sum(SPM.nscan);                               % total
Kdf = sum(arrayfun(@(x)(size(x.X0, 2)), SPM.xX.K)); % loss from filter
Xdf = rank(X);                                      % loss from regressors
Rdf = Tdf - Kdf - Xdf;                              % residual
fprintf(' df: %d - %d - %d = %d', Tdf, Kdf, Xdf, Rdf);
% other than SPM, we assume that whitening is perfect; for comparison
fprintf('   [SPM: trRV = %g  erdf = %g]\n', SPM.xX.trRV, SPM.xX.erdf)

% miscellaneous output
misc.mat = VY(1).mat;                               % voxels to mm transformation
misc.m = size(SPM.nscan, 2);                        % number of sessions
misc.sRow = {SPM.Sess.row};                         % scans of each session
% regressors for each session
% X has two parts, corresponding to SPM.xX.iC and SPM.xX.iB, each of
% which is block diagonal.
% The SPM.Sess.col for all sessions together are identical to
% SPM.xX.iC, but do not include the constant regressors in SPM.xX.iB.
misc.sCol = {SPM.Sess.col};
for si = 1 : misc.m
    % add constant regressors
    misc.sCol{si} = [misc.sCol{si}, SPM.xX.iB(si)];
end
% number of scans and degrees of freedom per session
% if not consistent across sessions, then this is an approximation
misc.fE = Rdf / misc.m;
misc.n = Tdf / misc.m;


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

