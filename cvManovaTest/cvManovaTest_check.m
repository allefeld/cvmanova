% print ROI results
fprintf('\ncvManovaRegion:\n')
for i = 1 : size(D, 3)
    for j = 1 : size(D, 1)
        fprintf('  region %d, contrast %d:  D = %.6f\n', i, j, D(j, 1, i))
    end
end

% compute checksums of searchlight results
fprintf('\ncvManovaSearchlight, MD5 checksums of results:\n')
d = dir([modelDir filesep 'spmD*.nii']);
for l = 1 : numel(d)
    % use only image data
    Y = spm_read_vols(spm_vol([modelDir filesep d(l).name]));
    % use only most significant 15 bits and sign
    Y = int16(round(Y(:) / max(abs(Y(:))) * 2^15));
    % compute and print MD5
    md5Ins = java.security.MessageDigest.getInstance('MD5');
    md5 = md5Ins.digest(typecast(Y, 'uint8'));
    md5 = lower(reshape(dec2hex(typecast(md5, 'uint8'))', 1, []));
    fprintf('  %s  %s\n', md5, d(l).name)
end
