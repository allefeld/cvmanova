function kbytes = systemFree

% returns the amount of available physical memory in units of 1024 bytes
%
% kbytes = systemFree()
%
% Copyright (C) 2014 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.


switch computer
    
    case {'PCWIN', 'PCWIN64'}
        % Windows: largest possible size for a single new array.
        user = memory;
        kbytes = user.MaxPossibleArrayBytes / 1024;
        
    case {'GLNX86', 'GLNXA64'}
        % read and parse /proc/meminfo
        fid = fopen('/proc/meminfo','r');
        d = textscan(fid, '%s%d%s');
        fclose(fid);
        names = d(1);
        names = names{1};
        names = strrep(names, ':', '');
        names = strrep(names, '(', '_');
        names = strrep(names, ')', '');
        data = d(2);
        data = data{1};
        data = num2cell(data);
        meminfo = cell2struct(data, names);
        
        kbytes = double(meminfo.MemFree);   % but see below...
        
    otherwise
        fprintf('systemFree is not implemented for this platform\n')
        kbytes = nan;
end


% Old version for Linux, gives unreliable results
%     % Works on Linux and probably other Unices. Macs, too?
%     % Memory that is available, possibly by flushing cache & buffers.
%     [~, out] = system('free -k');
%     kbytes = str2double(out(178:188));


% Not-that-old version that's not reliable either. WTF.
%         if isfield(meminfo, 'MemAvailable')
%             kbytes = meminfo.MemAvailable;
%         else
%             % low memory watermark:
%             % Not sure whether this can be determined simply as 4 x the sum of all
%             % "low" lines in /proc/zoneinfo. Instead use a value that's hopefully on
%             % the safe side.
%             wmark_low = 200000;
%             % recreate computations in Linux/fs/proc/meminfo.c
%             available = meminfo.MemFree - wmark_low;
%             pagecache = meminfo.Active_file + meminfo.Inactive_file;
%             pagecache = pagecache - min(pagecache / 2, wmark_low);
%             available = available + pagecache;
%             available = available + meminfo.SReclaimable ...
%                 - min(meminfo.SReclaimable / 2, wmark_low);
%             if available < 0
%                 available = 0;
%             end
%             kbytes = double(available);
%         end
