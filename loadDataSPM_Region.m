function [Ys, X, region, misc] = loadDataSPM_Region(dirName, region, whiten, highpass)

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
% Copyright (C) 2013-2016 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

% modularize pure data reading code -> spmReadVols
% change "v2mm" to mat

if nargin < 2
    region = [];
end
if nargin < 3
    whiten = true;
end
if nargin < 4
    highpass = false;    
end

SPMname = fullfile(dirName, 'SPM.mat');

fprintf('loading data\n')
fprintf(' via %s\n', SPMname)

% load SPM.mat
load(SPMname, 'SPM');

% recreate volumes, for backwards compatibility (with SPM2?)
fnames = SPM.xY.P;
fnames = patchPath(fnames, dirName);     % and to enable reading moved data files

fprintf(' reading volume information\n')
VY = spm_vol(fnames);
% but copy scaling information
[VY(:).pinfo] = SPM.xY.VY(:).pinfo;
nImages = numel(VY);

% check dimensions and orientations of all images
spm_check_orientations(VY);

% read brain mask image
if isfield(SPM, 'VM')
    fprintf(' using brain mask from SPM.VM\n');
    VM = spm_vol(fullfile(dirName, SPM.VM.fname));
else
    fprintf(' using brain mask from mask.img\n');
    VM = spm_vol(fullfile(dirName, 'mask.img'));
end
mask = logical(spm_read_vols(VM));
% possibly apply region mask
if isempty(region)
    error(' NO REGION')
end


for nr = 1:numel(region)
    if all(size(region{nr}) == size(mask))
        region{nr} = mask & logical(region{nr});
    else
        error('region mask doesn''t match\n')
    end
    nMaskVoxels(nr) = sum(region{nr}(:));
end;

% check memory
% memory needed, assuming double precision
memNeed = sum(nMaskVoxels(nr)) * nImages * 8 / 1024 / 1024;
% during whitening and filtering temporarily twice as much is needed
memNeed = memNeed * 2;
% available system memory
memFree = systemFree / 1024;
fprintf(' memory needed  %7.1f MiB\n', memNeed)
fprintf('        free    %7.1f MiB\n', memFree)
if memFree < memNeed, warning('not enough memory!'), end

Ys = cell(1,numel(region));
for nr = 1:numel(region)
    Ys{nr} = nan(nImages, nMaskVoxels(nr));
end;




% read and mask data
fprintf(' reading images\n')

for i = 1 : nImages
    V = VY(i);
    
    % read data directly (faster)
    y = V.private.dat(:, :, :, V.n(1));
    y = reshape(y, [], V.dim(3));
    y = bsxfun(@plus, bsxfun(@times, y, V.pinfo(1, :)), V.pinfo(2, :));
    for nr = 1:numel(region)
        Ys{nr}(i, :) = y(region{nr});
    end;
    
    % NOTE: if directly accessing memory-mapped data doesn't work,
    % comment the block above and uncomment the block below.
    
    %     % read data using SPM routines
    %     y = spm_read_vols(V);
    %     Y(i, :) = y(mask);
    
    if mod(i, 50) == 0
        fprintf(' %d of %d images read\n', i, nImages)
    end
end
fprintf('\n')


% get design matrix
X = SPM.xX.X;

% whiten data and design matrix
if whiten
    if isfield(SPM.xX, 'W')
        fprintf(' whitening\n')
        W = SPM.xX.W;
        for nr = 1:numel(region)
            Ys{nr} = W * Ys{nr};
        end;
        X = W * X;
    else
        fprintf(' * SPM.mat does not define whitening matrix!\n')
    end
end

% high-pass filter data and design matrix
if highpass
    fprintf(' high-pass-filtering\n')
    if highpass ~= 1
        for ii = 1:4
            K(ii) = struct('HParam', highpass,'row',    SPM.Sess(ii).row,'RT',     SPM.xY.RT);
        end;
        K = spm_filter(K);
        K(1).HParam
    else
        K = SPM.xX.K;
    end;
    
    
    for nr = 1:numel(region)
        Ys{nr} = spm_filter(K, Ys{nr});
    end;
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
end;
misc.fE = Rdf / misc.m;                             % if not consistent across sessions,
misc.n = Tdf / misc.m;                              % then this is an approximation

