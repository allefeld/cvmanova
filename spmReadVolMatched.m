function Y = spmReadVolMatched(V, Vtemplate, hold)

% read MR image data such that it matches the voxel grid of a template
%
% Y = spmReadVolMatched(V, Vtemplate, hold = 0)
%
% V:            volume struct to read from (or filename)
% Vtemplate:    volume struct to match (or filename)
% hold:         interpolation method (see spm_slice_vol)
%                    0         : nearest neighbour
%                    1         : trilinear interpolation
%                    2 – 127   : higher-order Lagrange interpolation
%                    -127 – -1 : sinc interpolation
% Y:            data (3D array)
%
% The image data are interpolated while reading in order to match the voxel
% grid of the template image. The returned data Y has the same dimensions
% and the same physical voxel locations as data read from the template, so
% that a one-to-one correspondence between voxels exists.
%
% The function can be used e.g. to read an ROI mask that has been
% inverse-normalized from MNI to subject space; such an image is in subject
% space, but not necessarily with the same voxel grid as the original
% subject-space data. In this case it is best to use nearest-neighbor
% interpolation (hold = 0) on a mask image that has been inverse normalized
% using spm_write_sn with flags.interp = 0, flags.preserve = 0, and a small
% voxel size (1 mm or smaller).
%
% Vtemplate may also be an array of volume structs, if the voxel grid is
% consistent. Only the fields mat (4 × 4 matrix) and dim (1 × 3 vector)
% have to be present.
%
% See also spmCoords, spm_write_sn, spm_slice_vol
%
%
% This file is part of v3 of cvmanova, see
% https://github.com/allefeld/cvmanova/releases
%
% Copyright (C) 2016 Carsten Allefeld
% adapted from spm_read_vols and John Ashburner's reslice.m, v1.42
% see https://www2.warwick.ac.uk/fac/sci/statistics/staff/academic-research/nichols/scripts/spm/johnsgems/#Gem8


% support filename
if ischar(V)
    V = spm_vol(V);
end
if ischar(Vtemplate)
    Vtemplate = spm_vol(Vtemplate);
end

if nargin < 3
    hold = 0;
end

if numel(V) > 1
    error('can only read from single volume!')
end

% in case the template is an array of volume structs,
% make sure that mat and dim are consistent
spm_check_orientations(Vtemplate);

% read data
Y = nan(Vtemplate(1).dim);
for k = 1 : Vtemplate(1).dim(3)    % for each plane
    % affine transform between voxel grids (inv(V.mat) * Vtemplate.mat)
    % plus z-direction offset to read the correct plane
    A = V.mat \ Vtemplate(1).mat * spm_matrix([0 0 k]);
    % read via spm_slice_vol
    Y(:, :, k) = spm_slice_vol(V, A, Vtemplate(1).dim(1 : 2), hold);
end


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
