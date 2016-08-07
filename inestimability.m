function ie = inestimability(C, X)

% degree of inestimability of a contrast wrt a design matrix
%
% ie = inestimability(C, X)
%
% C:    contrast matrix
% X:    design matrix
%
% For an estimable contrast, ie should be 0 save for numerical error.
% A number smaller than 1 indicates that the contrast has an estimable
% part. For a completely inestimable contrast, ie = 1.
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

% maximum "0" observed so far: 4.79*eps

if size(C, 1) < size(X, 2)
    C(size(X, 2), end) = 0;
end

% determine orthonormal basis of the range of the contrast matrix
RC = orth(C);

% determine orthonormal basis of the null space of the design matrix
NX = null(X);

% NX' * RC is a matrix of the inner products ("correlations")
% of all the C-range vectors with all the X-null vectors.
% Squared values indicate the proportion of variance of an effect described
% by a C-range vector that falls within the null-space of X.
% The worst-case proportion across the range of C is given by the matrix norm:
ie = norm(NX' * RC);
