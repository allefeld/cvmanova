function pMax = slSize(slRadius)

% searchlight size as a function of searchlight radius
% 
% single query:
% pMax = slSize(slRadius)
% slRadius:     searchlight radius
% pMax:         searchlight size
%
% tabulate radii and sizes:
% slSizes(slRadius = 5)
% slRadius:     maximum searchlight radius
%
%
% Copyright (C) 2016 Carsten Allefeld


if nargin == 0
    slRadius = 5;
end

% distances from center voxel on grid
[dxi, dyi, dzi] = ndgrid(-ceil(slRadius) : ceil(slRadius));
d = sqrt(dxi .^ 2 + dyi .^ 2 + dzi .^ 2);

% single query
if nargout > 0
    pMax = nnz(d <= slRadius);
    return
end

% tabulate radii and sizes
r = unique(d(d <= slRadius));
pMax = nan(size(r));
prec = 1 - floor(log10(min(diff(r))));      % necessary precision
fprintf('slRadius  pMax\n--------  ----\n')
for i = 1 : numel(r)
    % number of voxels within radius
    pMax(i) = nnz(d <= r(i));
    rd = ceil(r(i) * 10^prec) / 10^prec;    % round up
    fprintf('  %-5g   % 4d\n', rd, pMax(i))
    assert(numel(find(d <= rd)) == pMax(i)) % sufficient precision?
end
clear pMax


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
