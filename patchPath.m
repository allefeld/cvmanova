function fnames = patchPath(fnames, dirName)

% patch the path portion of filenames to access them after being moved
%
% fnames = patchPath(fnames, dirname)
%
% Assumption: SPM.mat and image files resided in subfolders of some common
% folder and were moved together. In moving, at least one element of the
% folder hierarchy common to SPM.mat and image files was preserved.
%
% Copyright (C) 2013 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

% find common beginning of all filenames
f = cell2mat(fnames(:));
ind = find(var(f) > 0, 1, 'first');
b = f(1, 1 : ind - 1);

% decompose into folder hierarchy
bp = regexp(b, '[^/]*/', 'match');
dp = regexp(dirName, '[^/]*/', 'match');

% use last common element to determine from and to root folders
[~, bpi, dpi] = intersect(bp, dp);
if ~isempty(bpi)
    from = [bp{1 : max(bpi)}];
    to = [dp{1 : max(dpi)}];
else
    from = [bp{:}];
    to = [dp{:}];
end

% patch folders
fprintf(' patching paths\n  from %s\n  to   %s\n', from, to)
f = [repmat(to, size(f, 1), 1), f(:, size(from, 2) + 1 : end)];
fnames = cellstr(f);

