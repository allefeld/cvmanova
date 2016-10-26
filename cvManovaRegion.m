function [D, p] = cvManovaRegion(dirName, regions, Cs, permute, lambda)

% cross-validated MANOVA on region
%
% [D, p] = cvManovaRegion(dirName, regions, Cs, permute = false, lambda = 0)
%
% dirName:  directory where the SPM.mat file referring to an estimated
%           model is located
% regions:  region mask(s), (cell array of) logical 3D volume(s)
% Cs:       cell array of contrast vectors or matrices
% lambda:   regularization parameter (0–1)
% permute:  whether to compute permutation values
% D:        pattern distinctness, contrasts × permutations × regions
% p:        number of voxels in the region(s)
%
%
% Copyright (C) 2015-2016 Carsten Allefeld


fprintf('\n\ncvManovaRegion\n\n')

if nargin < 4
    permute = false;
end
if nargin < 5
    lambda = 0;
end    
if islogical(lambda)
    error('incorrect order of parameters?')
end
if isempty(regions)
    error('no region mask specified!')
end

if dirName(end) ~= filesep
    dirName = [dirName filesep];
end

% load data, design matrix etc.
[Ys, Xs, ~, misc] = loadDataSPM(dirName, regions);
nRegions = numel(misc.rmvi);

% compute on regions
fprintf('\ncomputing cross-validated MANOVA on region\n')
cvManovaCore(nan(0, 1), Ys, Xs, Cs, misc.fE, permute, lambda);
clear Ys Xs
D = [];
for i = 1 : nRegions
    D = cat(3, D, cvManovaCore(misc.rmvi{i}));
end
p = cellfun(@numel, misc.rmvi);
clear cvManovaCore          % free memory

% separate contrast and permutation dimensions of result
D = reshape(D, numel(Cs), [], nRegions);


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
