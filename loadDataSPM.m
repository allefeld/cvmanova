function [Ys, Xs, mask, misc] = loadDataSPM(dirName, regions)

% load fMRI data via SPM.mat
%
% [Ys, Xs, mask, misc] = loadDataSPM(dirName, regions = {})
%
% dirName:  name of directory that contains SPM.mat
% regions:  optional additional region mask(s),
%           (cell array of) logical 3D volume(s) or filename(s)
% Ys:       MR data (within mask) for each session, scans x voxels
% Xs:       design matrix for each session, scans x regressors
% mask:     analysis brain mask, logical 3D volume;
%           possibly combined with region mask
% misc:     struct with additional data:
%     mat   voxels to mm transformation matrix
%     fE    residual degrees of freedom per session
%     rmvi  cell array of mask voxel indices for each region
%
% Y & X and are high-pass filtered and whitened.
% Y includes only those voxels selected through mask.
%
%
% Copyright (C) 2013-2016 Carsten Allefeld


% default argument values
if nargin < 2
    regions = {};
end

% load SPM.mat
SPMname = fullfile(dirName, 'SPM.mat');
fprintf('loading data\n')
fprintf(' via %s\n', SPMname)
load(SPMname, 'SPM');

% get data volumes and check matching voxel grid
VY = SPM.xY.VY;
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
assert(isfield(SPM, 'VM'), ' no analysis brain mask in SPM.VM!')
try
    mask = logical(spm_read_vols(SPM.VM));
catch
    % SPM8 stores the filename without the path
    VM = spm_vol([dirName SPM.VM.fname]);
    mask = logical(spm_read_vols(VM));
end
fprintf(' %d in-mask voxels\n', sum(mask(:)));

% possibly apply region mask(s)
if isempty(regions)
    fprintf(' no region mask\n')
    rmvi = {};
else
    if ~iscell(regions)
        regions = {regions};
    end
    nRegions = numel(regions);
    for i = 1 : nRegions
        if ~isnumeric(regions{i})
            regions{i} = logical(spmReadVolMatched(regions{i}, VY(1)));
        end
    end
    try
        regions = logical(cat(4, regions{:}));
    catch
        error('region masks don''t match!')
    end
    assert(isequal(size(regions(:, :, :, 1)), size(mask)), ...
        'region masks don''t match!')
    % restrict brain mask to conjunction of regions
    mask = mask & any(regions, 4);
    % determine mask voxel indices for each region
    regions = reshape(regions, [], nRegions);
    rmvi = cell(nRegions, 1);
    for i = 1 : nRegions
        rmvi{i} = find(regions(mask(:), i));
        fprintf('  %d in-mask voxels in region %d\n', numel(rmvi{i}), i)
    end
    fprintf(' %d in-mask voxels in regions\n', sum(mask(:)));
end

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

% separate Y and X into session blocks
m = numel(SPM.nscan);
Xs = cell(m, 1);
Ys = cell(m, 1);
for si = 1 : m
    Ys{si} = Y(SPM.Sess(si).row, :);
    % SPM.Sess(:).col does not include constant regressors,
    % get those from SPM.xX.iB
    Xs{si} = X(SPM.Sess(si).row, [SPM.Sess(si).col, SPM.xX.iB(si)]);
end
clear Y X

% degrees of freedom for each session
Tdf = nan(m, 1);
Kdf = nan(m, 1);
Xdf = nan(m, 1);
Rdf = nan(m, 1);
for si = 1 : m
    Tdf(si) = SPM.nscan(si);                    % total
    Kdf(si) = rank(SPM.xX.K(si).X0);            % loss from filter
    Xdf(si) = rank(Xs{si});                     % loss from regressors
    Rdf(si) = Tdf(si) - Kdf(si) - Xdf(si);      % residual
end
fprintf(' df: %d - %d - %d = %d', sum(Tdf), sum(Kdf), sum(Xdf), sum(Rdf));
% other than SPM, we assume that whitening is perfect; for comparison
fprintf('   [SPM: trRV = %g  erdf = %g]\n', SPM.xX.trRV, SPM.xX.erdf)

% miscellaneous output
% voxels to mm transformation
misc.mat = VY(1).mat;
% residual degrees of freedom for each session
misc.fE = Rdf;
% mask voxel indices for each region
misc.rmvi = rmvi;


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
