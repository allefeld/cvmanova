function [cMatrix, cName] = contrasts(fLevel, fName)

% generate contrasts for main effects and interactions of a factorial design
% 
% [cMatrix, cName] = contrasts(fLevel)
% [cMatrix, cName] = contrasts(fLevel, fName)
%
% fLevel:   number of levels for each factor (row vector)
% fName:    names of factors (cell array)
% cMatrix:  contrast matrices (cell array)
% cName:    contrast names (cell array)
%
% the first factor is the one being enumerated slowest.
% contrasts are not orthonormalized!
%
%
% This file is part of v3 of cvmanova, see
% https://github.com/allefeld/cvmanova/releases
%
% Copyright (C) 2013–2016 Carsten Allefeld


% number of factors
nf = size(fLevel, 2);

if nargin < 2
    % generate generic names for factors
    fName = num2cell(char((1 : nf) + '@'));
end

% generate sorted list of contrast signatures
nc = 2 ^ nf - 1;
cs = bitget(ones(nf, 1) * (1 : nc), (1 : nf)' * ones(1, nc))';
[~, ind] = sort(sum(cs, 2));
cs = cs(ind, :);
cs = logical(cs);

% compute contrast elements
e = cell(2, nf);
for fi = 1 : nf
    e{1, fi} = ones(fLevel(fi), 1);
    e{2, fi} = -(diff(eye(fLevel(fi)))');
end

% compute contrasts
cMatrix = cell(nf, 1);
cName = cell(nf, 1);
for ci = 1 : nc
    % contrast matrix
    cMatrix{ci} = 1;
    for fi = 1 : nf
        cMatrix{ci} = kron(cMatrix{ci}, e{cs(ci, fi) + 1, fi});
    end
    % contrast name
    l = sprintf('×%s', fName{cs(ci, :)});
    cName{ci} = l(2:end);
end

% if no output, print contrasts
if nargout == 0
    for ci = 1 : nc
        fprintf('%s:\n', cName{ci})
        disp(cMatrix{ci})
        fprintf('\n')
    end
    clear cMatrix cName
end


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
