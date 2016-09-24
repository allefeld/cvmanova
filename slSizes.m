function slSizes(rMax)

% tabulate searchlight radii that lead to different searchlight sizes
%
% slSizes(rMax = 5)
%
% rMax: maximum radius
%
%
% Copyright (C) 2016 Carsten Allefeld


if nargin == 0
    rMax = 5;
end

% distances from center voxel on grid
[dxi, dyi, dzi] = ndgrid(-ceil(rMax) : ceil(rMax));
d = sqrt(dxi .^ 2 + dyi .^ 2 + dzi .^ 2);

% occurring distances
r = unique(d(d <= rMax));

p = nan(size(r));
prec = 1 - floor(log10(min(diff(r))));      % necessary precision
fprintf('slRadius  pMax\n--------  ----\n')
for i = 1 : numel(r)
    % number of voxels within radius
    p(i) = numel(find(d <= r(i)));
    rd = ceil(r(i) * 10^prec) / 10^prec;    % round up
    fprintf('  %-5g   % 4d\n', rd, p(i))
    assert(numel(find(d <= rd)) == p(i))    % sufficient precision?
end


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
