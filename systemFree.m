function kbytes = systemFree

% returns the amount of free memory in units of 1024 bytes
%
% kbytes = systemFree()
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

if ispc
    % Windows: largest possible size for a single new array.
    user = memory;
    kbytes = user.MaxPossibleArrayBytes / 1024;
else
    % Works on Linux and probably other Unices. Macs, too?
    % Memory that is available, possibly by flushing cache & buffers.
    [~, out] = system('free -k');
    kbytes = str2double(out(178:188));
end
