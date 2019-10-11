% slice time information is not available

% make sure realigned BOLD data are there
fnrBOLD = [sub filesep 'rbold.nii'];
if ~exist(fnrBOLD, 'file')
    fprintf('realigning volumes\n')
    matlabbatch = {};
    for ri = 1 : nRuns
        vi = (ri - 1) * nVolsPerRun + (1 : nVolsPerRun);
        vn = arrayfun(@(i) sprintf('%s,%d', fnBOLD, i), vi, 'UniformOutput', false);
        matlabbatch{1}.spm.spatial.realign.estwrite.data{ri} = vn';
    end
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    spm_jobman('run', matlabbatch)
end

