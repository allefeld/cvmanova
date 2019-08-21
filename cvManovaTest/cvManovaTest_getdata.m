% make sure subdirectory with unpacked data is there
if ~exist(sub, 'dir')
    fn = [sub '-2010.01.14.tar.gz'];
    url = ['http://data.pymvpa.org/datasets/haxby2001/' fn];
    
    if ~exist(fn, 'file')
        fprintf('downloading %s\n', fn)
        urlwrite(url, fn);
    end
    
    fprintf('unpacking %s\n', fn)
    untar(fn)
    delete(fn)
end

% make sure gunzipped BOLD data are there
fnBOLD = [sub filesep 'bold.nii'];
if ~exist(fnBOLD, 'file')
    fprintf('unpacking %s.gz\n', fnBOLD)
    gunzip([fnBOLD '.gz'])
    delete([fnBOLD '.gz'])
end

% make sure design information is there
fnDesign = [sub filesep 'design.mat'];
if ~exist(fnDesign, 'file')
    % other design information
    TR = 2.5;
    nRuns = 12;
    nVolsPerRun = 121;
    % extract stimulus information
    fnLabels = [sub filesep 'labels.txt'];
    fprintf('extracting stimulus timing information from %s\n', fnLabels)
    % volume labels
    labels = readtable(fnLabels, 'Delimiter', ' ');
    %   sequence within run: (12 s rest, 24 s block) × 8, 12 s rest = 300 s
    % extract condition for each block (from its middle)
    conditions = {'face', 'house', 'cat', 'bottle', ...
        'scissors', 'shoe', 'chair', 'scrambledpix'};   % + 'rest'
    nConds = numel(conditions);
    [~, vsi] = ismember(labels{:, 1}, conditions);  % stimulus of volumes
    vsi = reshape(vsi, nVolsPerRun, nRuns);
    bsi = vsi(round((24 :36: 300) / TR) + 1, :);    % blocks × runs
    % prepare stimulus onset and duration information (in seconds)
    %   each block has 12 stimulus presentations of the same condition
    %   each stimulus was shown for 500 ms and was followed by a 1500 ms ISI
    figure
    col = hsv(8);
    onsets = cell(nRuns, nConds);
    durations = cell(nRuns, nConds);
    for ri = 1 : nRuns
        subplot(4, 3, ri)
        % plot volume indicators
        plot(1 : nVolsPerRun, ones(1, nVolsPerRun), 'k.')
        axis([0, nVolsPerRun + 1, 0.01, 1.2])
        title(sprintf('run %d', ri))
        hold all
        for ci = 1 : nConds
            % block in which condition is realized
            bi = find(bsi(:, ri) == ci);
            % onset of block
            bo = (bi - 1) * 36 + 12;
            % onsets of stimuli in block
            onsets{ri, ci} = bo + (0 : 11) * 2;
            durations{ri, ci} = 0.5 * ones(1, 12);
            % plot stimulus type for each volume,
            % stimulus onsets for each condition
            plot(vsi(:, ri) == ci, 'Color', col(ci, :))
            plot(onsets{ri, ci} / TR + 1, 1.1 * ones(1, 12), ...
                '.', 'Color', col(ci, :))
        end
    end
    % save design information
    save(fnDesign, 'TR', 'nRuns', 'nVolsPerRun', ...
        'conditions', 'onsets', 'durations')
end

% make sure gunzipped region mask data are there
fnRegions = {[sub filesep 'mask4_vt.nii'], ...
             [sub filesep 'mask8_face_vt.nii'], ...
             [sub filesep 'mask8_house_vt.nii']};
for i = 1 : numel(fnRegions)
    if ~exist(fnRegions{i}, 'file')
        fprintf('unpacking %s.gz\n', fnRegions{i})
        gunzip([fnRegions{i} '.gz'])
        delete([fnRegions{i} '.gz'])
    end
end
