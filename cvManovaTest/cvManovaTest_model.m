% make sure specified model is there
modelDir = [sub filesep 'model'];
fnSPM = [modelDir filesep 'SPM.mat'];
if ~exist(fnSPM, 'file')
    fprintf('specifying model\n')
    load(fnDesign);
    V = spm_vol(fnBOLD);    % for motion parameters
    mat = cat(4, V.mat);
    mkdir(modelDir)
    
    fmri_spec = struct;
    fmri_spec.dir = {modelDir};
    fmri_spec.timing.units = 'secs';
    fmri_spec.timing.RT = TR;
    fmri_spec.timing.fmri_t = 25;
    fmri_spec.timing.fmri_t0 = 12;
    for ri = 1 : nRuns
        vi = (ri - 1) * nVolsPerRun + (1 : nVolsPerRun);
        vn = arrayfun(@(i) sprintf('%s,%d', [pwd filesep fnrBOLD], i), vi, 'UniformOutput', false);
        Q = nan(nVolsPerRun, 6);
        for i = 1 : nVolsPerRun
            qq = spm_imatrix(mat(:, :, vi(i)) / mat(:, :, vi(1)));
            Q(i, :) = qq(1 : 6);
        end
        
        fmri_spec.sess(ri).scans = vn;
        for ci = 1 : numel(conditions)
            fmri_spec.sess(ri).cond(ci).name = conditions{ci};
            fmri_spec.sess(ri).cond(ci).onset = onsets{ri, ci};
            fmri_spec.sess(ri).cond(ci).duration = durations{ri, ci};
            fmri_spec.sess(ri).cond(ci).tmod = 0;
            fmri_spec.sess(ri).cond(ci).orth = 1;
        end
        fmri_spec.sess(ri).multi = {''};
        for mi = 1 : 6
            fmri_spec.sess(ri).regress(mi).name = sprintf('motion%d', mi);
            fmri_spec.sess(ri).regress(mi).val = Q(:, mi);
        end
        fmri_spec.sess(ri).multi_reg = {''};
        fmri_spec.sess(ri).hpf = 128;
    end
%     fmri_spec.fact(1).name = 'stim';
%     fmri_spec.fact(1).levels = numel(conditions);
    fmri_spec.bases.hrf.derivs = [0 0];
    fmri_spec.volt = 1;
    fmri_spec.global = 'None';
    fmri_spec.mthresh = 0.8;
    fmri_spec.mask = {''};
    fmri_spec.cvi = 'AR(1)';
    
    matlabbatch = {};
    matlabbatch{1}.spm.stats.fmri_spec = fmri_spec;
    spm_jobman('run', matlabbatch(1))
end

% make sure estimated model is there
if numel(dir([modelDir filesep 'beta_*.*'])) == 0
    fprintf('estimating model\n')
    matlabbatch = {};
    matlabbatch{1}.spm.stats.fmri_est.spmmat = {fnSPM};
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch(1))
end
