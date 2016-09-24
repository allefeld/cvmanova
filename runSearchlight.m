function res = runSearchlight(checkpoint, slRadius, mask, fun, varargin)

% general purpose searchlight:
% apply function to data contained in a sliding spherical window
%
% res = runSearchlight(checkpoint, slRadius, mask, fun, ...)
%
% checkpoint:   name of checkpoint file ([] to disable checkpointing)
% slRadius:     radius of searchlight, in voxels
% mask:         3-dimensional logical array indicating which voxels to use
% fun:          function to call with voxel indices within window
% ...:          additional arguments are passed through to the function
% res:      results, 2-dimensional matrix, voxels x output values
%
% The function has to be of the form r = fun(mvi, ...)
% mvi:      column vector of linear indices into the mask voxels
% r:        row vector of results
% The function is called with mvi = [] first, and the output is used to
% initialize res.
%
% A voxel is included in the searchlight if its distance from the center is
% *smaller than or equal to* the radius. Note that fractional values are
% possible. Run slSizes for a table of meaningful values and resulting
% searchlight sizes
%
% Intermediate results are saved at regular intervals to the checkpoint
% file, if given. On a subsequent run, if the checkpoint file exists, its
% contents are loaded and the computation continues from that point. 
%
% See also slSizes
%
%
% Copyright (C) 2013-2016 Carsten Allefeld


% normalize checkpoint file name, preserving emptiness
if ~isempty(checkpoint)
    [~, ~, ext] = fileparts(checkpoint);
    if ~strcmp(ext, '.mat')
        checkpoint = [checkpoint '.mat'];
    end
end

% volume dimensions
dim = size(mask);
nVolumeVoxels = prod(dim);
nMaskVoxels = sum(mask(:));

% determine searchlight voxel offsets relative to center voxel
% prototype searchlight
[dxi, dyi, dzi] = ndgrid(-ceil(slRadius) : ceil(slRadius));
PSL = (dxi .^ 2 + dyi .^ 2 + dzi .^ 2 <= slRadius .^ 2);
% spatial offsets
dxi = dxi(PSL);
dyi = dyi(PSL);
dzi = dzi(PSL);
% index offsets
PSL(dim(1), dim(2), dim(3)) = 0;
di = find(PSL);
cInd = find((dxi == 0) & (dyi == 0) & (dzi == 0));
di = di - di(cInd);                                                         %#ok<FNDSB>
clear PSL cInd

% mapping from volume to mask voxel indices
vvi2mvi = nan(nVolumeVoxels, 1);
vvi2mvi(mask) = 1 : nMaskVoxels;

% initialize result volume(s)
res = fun(nan(0, 1), varargin{:});
res = repmat(reshape(res, 1, []), nVolumeVoxels, 1);

fprintf(' running searchlight\n')
fprintf('  searchlight size: %d\n', size(di, 1))
tic
t = 0;
cmvi = 0;
for cvvi = 1 : nVolumeVoxels    % searchlight center *volume voxel index*
    if (cvvi == 1) && ~isempty(checkpoint)
        % load checkpoint file, if existing
        if exist(checkpoint, 'file')
            load(checkpoint)
            fprintf('  *restart*  %6d voxels  %5.1f %%\n', ...
                cmvi, cmvi / nMaskVoxels * 100)
        end
    end
    
    % process only if center is within mask
    if mask(cvvi)
        cmvi = cmvi + 1;        % searchlight center *mask voxel index*
        
        % searchlight center coordinates
        [xi, yi, zi] = ind2sub(dim, cvvi);
        % searchlight voxel coordinates; limit to volume boundaries
        ind = (xi + dxi >= 1) & (xi + dxi <= dim(1)) & ...
            (yi + dyi >= 1) & (yi + dyi <= dim(2)) & ...
            (zi + dzi >= 1) & (zi + dzi <= dim(3));
        % searchlight voxel volume indices
        vvi = cvvi + di(ind);
        % discard out-of-mask voxels
        vvi = vvi(mask(vvi) == 1);
        % translate to mask voxel indices
        mvi = vvi2mvi(vvi);
        
        % call function and store output
        res(cvvi, :) = fun(mvi, varargin{:});
    end
    
    
    nt = toc;
    if (nt - t > 30) || (cvvi == nVolumeVoxels)
        % progress report
        t = nt;
        fprintf(' %6.1f min  %6d voxels  %5.1f %%\n', ...
            t / 60, cmvi, cmvi / nMaskVoxels * 100)
        
        if (cvvi < nVolumeVoxels) && ~isempty(checkpoint)
            % save checkpoint file
            save(checkpoint, 'res', 'cvvi', 'cmvi')
        end
    end
end
% delete checkpoint file after completion
if ~isempty(checkpoint)
    spm_unlink(checkpoint)
end


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
