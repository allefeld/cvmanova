function [D, p] = cvManovaRegion(dirName, region, Cs, lambda, permute)

% cross-validated MANOVA on region
%
% [D, p] = cvManovaRegion(dirName, region, Cs, lambda = 0, permute = false)
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
% If region is [], all in-mask voxels are used.
%
%
% Copyright (C) 2015 Carsten Allefeld


fprintf('\n\ncvManovaRegion\n\n')

if nargin < 4
    lambda = 0;
end    
if nargin < 5
    permute = false;
end
% *** sequence of parameters different between Region and Searchlight!

if dirName(end) ~= filesep, dirName = [dirName filesep]; end

% load data, design matrix etc.
[Ys, Xs, ~, misc] = loadDataSPM(dirName, region);

% compute on region
fprintf('\ncomputing cross-validated MANOVA on region\n')
cvManovaCore(nan(0, 1), Ys, Xs, Cs, misc.fE, permute, lambda);
p = size(Ys{1}, 2);
clear Ys Xs
D = cvManovaCore(1 : p);

% separate contrast and permutation dimensions
D = reshape(D, numel(Cs), []);


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

