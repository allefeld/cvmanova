function [perms, nPerms] = signPermutations(n, maxPerms)

% generate sign permutations
%
% perms = signPermutations(n, maxPerms = 5000)
%
% Permutations are randomly selected if the full enumeration
% is larger than maxPerms.
%
%
% Copyright (C) 2013 Carsten Allefeld


if nargin < 2
    maxPerms = 5000;    % adapted for z-values up to +-3
end

% determine permutations; first one is always the identity permutation
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

