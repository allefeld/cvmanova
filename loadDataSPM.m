function [Y, X, mask, misc] = loadDataSPM(dirName, region, whiten, highpass)

% load fMRI data via SPM.mat
%
% [Y, X, mask, misc] = loadDataSPM(dirName, region = [], whiten = true, highpass = true)
%
% dirName:  name of directory that contains SPM.mat
% region:   optional additional region mask, logical 3d-volume
% whiten:   whether to whiten data and design matrix
% highpass: whether to high-pass filter data and design matrix
% Y:        MR data (within mask), volumes x voxels
% X:        design matrix, volumes x regressors
% mask:     brain mask, logical 3d-volume, possibly combined with region mask
% misc:     struct with additional data:
%     v2mm  voxels to mm transformation matrix
%     sRow  rows for each session
%     sCol  columns for each session
%     m     number of sessions
%     n     number of images per session
%     fE    residual degrees of freedom per session
%
% If not otherwise specified, Y & X and are high-pass filtered and whitened.
% Y includes only those voxels selected through mask.
%
%
% Copyright (C) 2013-2016 Carsten Allefeld


% change "v2mm" to mat


% default argument values
if nargin < 2
    region = [];
end
if nargin < 3
    whiten = true;
end
if nargin < 4
    highpass = true;
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

% read brain mask image
if isfield(SPM, 'VM')
    try
        mask = logical(spm_read_vols(SPM.VM));
    catch
        % SPM8 stores the filename without the path
        VM = spm_vol([dirName SPM.VM.fname]);
        mask = logical(spm_read_vols(VM));
    end
    fprintf(' using brain mask from SPM.VM\n');
else
    mask = true(VY(1).dim);
    fprintf(' no brain mask!\n')        % should not happen
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

% check memory
% memory needed, assuming double precision
memNeed = sum(mask(:)) * nImages * 8 / 1024 / 1024;
% during whitening and filtering temporarily twice as much is needed
memNeed = memNeed * 2;
% available system memory
memFree = systemFree / 1024;
fprintf(' memory needed  %7.1f MiB\n', memNeed)
fprintf('        free    %7.1f MiB\n', memFree)
if memFree < memNeed, fprintf(2, 'not enough memory!\n'); end

% read and mask data
fprintf(' reading images\n')
pattern = SPM.xY.P(1, :);
pattern(~all(diff(SPM.xY.P) == 0)) = '?';
fprintf('  from %s\n', pattern)
[Y, mask] = spmReadVolsMasked(VY, mask);

% get design matrix
X = SPM.xX.X;

% whiten data and design matrix
if whiten
    if isfield(SPM.xX, 'W')
        fprintf(' whitening\n')
        W = SPM.xX.W;
        Y = W * Y;
        X = W * X;
    else
        fprintf(' * SPM.mat does not define whitening matrix!\n')
    end
end

% high-pass filter data and design matrix
if highpass
    fprintf(' high-pass-filtering\n')
    Y = spm_filter(SPM.xX.K, Y);
    X = spm_filter(SPM.xX.K, X);
end

% degrees of freedom
Tdf = sum(SPM.nscan);                               % total
if highpass                                         % loss from high-pass filter
    Kdf = sum(arrayfun(@(x)(size(x.X0, 2)), SPM.xX.K));
else
    Kdf = 0;
end
Xdf = rank(X);                                      % loss from regressors
Rdf = Tdf - Kdf - Xdf;                              % residual
fprintf(' df: %d - %d - %d = %d\n', Tdf, Kdf, Xdf, Rdf);

% miscellaneous output
misc.v2mm = VY(1).mat;                              % voxels to mm transformation
misc.m = size(SPM.nscan, 2);                        % number of sessions
misc.sRow = {SPM.Sess.row};                         % volumes for each session
misc.sCol = {SPM.Sess.col};                         % regressors for each session
for si = 1 : misc.m                                 % add constant regressors
    misc.sCol{si} = [misc.sCol{si}, SPM.xX.iB(si)];
    % NOT SURE WHETHER THIS IS GENERALLY CORRECT
end
misc.fE = Rdf / misc.m;                             % if not consistent across sessions,
misc.n = Tdf / misc.m;                              % then this is an approximation


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

