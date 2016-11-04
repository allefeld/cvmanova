function checksum = fletcher16(data)

% compute Fletcher-16 checksum of byte array
%
% checksum = fletcher16(data)
%
% data:     data to be checksummed,
%           array of integers in the range 0 to 255 (0xFF)
% checksum: Fletcher-16 checksum,
%           16-bit integer where each of the bytes is in the range 0 to 254 (0xFE)
%
% https://en.wikipedia.org/wiki/Fletcher%27s_checksum
%
%
% This file is part of v3 of cvmanova, see
% https://github.com/allefeld/cvmanova/releases
%
% Copyright (C) 2016 Carsten Allefeld


if ~isreal(data) || any(data < 0) || any(data > 255) || any(data ~= fix(data))
    error('invalid data!')
end

% vectorized version of
% https://en.wikipedia.org/wiki/Fletcher%27s_checksum#Straightforward
%
% Computations are done in double precision, because sum(..., 'native') is
% slower and the 53 bits of integer precision of a double are not much less
% than those of a uint64.
% Modulo operations are deferred to the end for efficiency, with the
% drawback that the algorithm will break down for too large data arrays due
% to loss of least significant bits. This can however be detected.

data = double(data(:));
sum1 = mod(sum(data), 255);
sum2 = sum(cumsum(data));
assert(sum2 < 2^53, 'insufficient precision of implementation!')
sum2 = mod(sum2, 255);
checksum = sum2 * 256 + sum1;


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
