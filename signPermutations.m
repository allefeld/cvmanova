function [perms, nPerms] = signPermutations(n, maxPerms)

% generate sign permutations
%
% [perms, nPerms] = signPermutations(n, maxPerms = 5000)
%
% n:         number of data points
% maxPerms:  maximum number of permutations
% perms:     permutations, n × nPerms array of ±1
% nPerms:    number of permutations
%
% Permutations are randomly selected if the full enumeration
% is larger than maxPerms. The first permutation is always the neutral
% permutation, ones(n, 1).
%
%
% This file is part of v3 of cvmanova, see
% https://github.com/allefeld/cvmanova/releases
%
% Copyright (C) 2013–2016 Carsten Allefeld


if nargin < 2
    maxPerms = 5000;    % supports z-values of up to +-3
end

% determine permutations; first one is always the neutral permutation
if 2 ^ n <= maxPerms
    % full enumeration of permutations
    nPerms = 2 ^ n;
    perms = bitget(ones(n, 1) * (0 : nPerms - 1), (1 : n)' * ones(1, nPerms));
else
    % random (Monte Carlo) selection of permutations
    nPerms = maxPerms;  
    perms = [zeros(n, 1), (rand(n, nPerms - 1) > 0.5)];
end

perms = (-1) .^ perms;    % permuted signs


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
