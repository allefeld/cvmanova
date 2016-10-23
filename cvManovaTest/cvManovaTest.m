% test the implementation of 'cross-validated MANOVA'
%
% This script uses the data of subject 1 from Haxby et al. (2001), which
% are downloaded automatically.
% Source: http://dev.pymvpa.org/datadb/haxby2001.html
%
% It can also be taken as an example for how to use the implementation.
%
% With Matlab 8.5.0 (R2015a) and SPM12 r6685, this script
% produces the following pattern distinctness values on a region:
%   5.443, 4.421
% and generates image files with the following MD5 checksums:
%   bf96acecfe2499e19e141e3baaaaacc2  spmD_C0001_P0001.nii
%   57be2e6a77f6d477ce25d62b168af77b  spmD_C0002_P0001.nii
%
%
% Copyright (C) 2016 Carsten Allefeld

clear

% select subject; 6 is missing anatomy, 5 has corrupted data
sub = 'subj1';
fprintf('analyzing data of %s from Haxby et al. (2001)\n', sub)

% init
spm defaults fmri
spm_jobman initcfg

% prerequisites of cvMANOVA
cvManovaTest_getdata
cvManovaTest_preprocess
cvManovaTest_model

% set up contrasts
% notice the transpose operators!
Cs = {};
% 1) main effect of stimulus
Cs{1} = [ 1 -1  0  0  0  0  0  0
          0  1 -1  0  0  0  0  0
          0  0  1 -1  0  0  0  0
          0  0  0  1 -1  0  0  0
          0  0  0  0  1 -1  0  0
          0  0  0  0  0  1 -1  0
          0  0  0  0  0  0  1 -1]';
% 2) main effect, meaningful stimuli only
Cs{2} = [ 1 -1  0  0  0  0  0  0
          0  1 -1  0  0  0  0  0
          0  0  1 -1  0  0  0  0
          0  0  0  1 -1  0  0  0
          0  0  0  0  1 -1  0  0
          0  0  0  0  0  1 -1  0]';

% run cvMANOVA on region
region = logical(spm_read_vols(spm_vol(fnRegion)));
[D, p] = cvManovaRegion(modelDir, region, Cs);

% run cvMANOVA on searchlight of radius 3
cvManovaSearchlight(modelDir, 3, Cs)

% compute checksums of searchlight results
fprintf('\ncvManovaRegion: D = %.3f, %.3f\n', D)
fprintf('\ncvManovaSearchlight, MD5 checksums of results:\n')
d = dir([modelDir filesep 'spmD_*.nii']);
for l = 1 : numel(d)
    fid = fopen([modelDir filesep d(l).name], 'rb');
    bytes = fread(fid, 'uint8=>uint8');
    fclose(fid);
    md5Ins = java.security.MessageDigest.getInstance('MD5');
    md5Str = lower(reshape(dec2hex(typecast(...
        md5Ins.digest(bytes), 'uint8'))', 1, []));
    fprintf('%s  %s\n', md5Str, d(l).name)
end
fprintf('\nconsider deleting the directory %s and its contents\n', modelDir)


% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
