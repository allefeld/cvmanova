% test the implementation of 'cross-validated MANOVA'
%
% This script uses the data of subject 1 from Haxby et al. (2001), which
% are downloaded automatically.
% Source: http://dev.pymvpa.org/datadb/haxby2001.html
%
% It can also be taken as an example for how to use the implementation.
%
% With Matlab 8.5.0 (R2015a) and SPM12 r6685, this script
% produces the following pattern distinctness values on regions:
%   region 1, contrast 1:  D = 5.443427
%   region 1, contrast 2:  D = 1.021870
%   region 2, contrast 1:  D = 0.314915
%   region 2, contrast 2:  D = 0.021717
%   region 3, contrast 1:  D = 1.711423
%   region 3, contrast 2:  D = 0.241187
% and generates images with the following checksums:
%   03adb4e589c9e1da8f08829c839b26d9  spmD_C0001_P0001.nii
%   7a8f0d5918363c213e0d749a1bfdd665  spmD_C0002_P0001.nii
%   8bfe2b4261920127b2fcf5fe5358a340  spmDs_C0001_P0001.nii
%   e7d2c583c5159feb671dea7ff2b72570  spmDs_C0002_P0001.nii


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
% conditions: face, house, cat, bottle, scissors, shoe, chair, scrambledpix
Cs = {};
% 1) main effect of stimulus
Cs{1} = [ 1 -1  0  0  0  0  0  0
          0  1 -1  0  0  0  0  0
          0  0  1 -1  0  0  0  0
          0  0  0  1 -1  0  0  0
          0  0  0  0  1 -1  0  0
          0  0  0  0  0  1 -1  0
          0  0  0  0  0  0  1 -1]';
% 2) main effect of category within object
Cs{2} = [ 0  0  0  1 -1  0  0  0
          0  0  0  0  1 -1  0  0
          0  0  0  0  0  1 -1  0]';
% notice the transpose operators!

% run cvMANOVA on regions
[D, p] = cvManovaRegion(modelDir, fnRegions, Cs);

% run cvMANOVA on searchlight of radius 3
cvManovaSearchlight(modelDir, 3, Cs)

% check results
cvManovaTest_check
fprintf('\nconsider deleting the directory %s and its contents\n', modelDir)
