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
%   5.443427, 4.421474
% and generates image files with the following MD5 checksums:
%   03adb4e589c9e1da8f08829c839b26d9  spmD_C0001_P0001.nii
%   70b2d9cb8839b502578ab0f9c1ffbe55  spmD_C0002_P0001.nii

clear

% select subject
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
fprintf('\ncvManovaRegion: D = %.6f, %.6f\n', D)
fprintf('\ncvManovaSearchlight, MD5 checksums of results:\n')
d = dir([modelDir filesep 'spmD_*.nii']);
for l = 1 : numel(d)
    % use only data
    Y = spm_read_vols(spm_vol([modelDir filesep d(l).name]));
    % use only most significant 15 bits and sign
    Y = int16(round(Y(:) / max(abs(Y(:))) * 2^15));
    % compute and print MD5
    md5Ins = java.security.MessageDigest.getInstance('MD5');
    md5 = md5Ins.digest(typecast(Y, 'uint8'));
    md5 = lower(reshape(dec2hex(typecast(md5, 'uint8'))', 1, []));
    fprintf('%s  %s\n', md5, d(l).name)
end
fprintf('\nconsider deleting the directory %s and its contents\n', modelDir)
