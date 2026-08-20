function ttestStruct = extract_ttest_features(sub, searchWord, dataDir, varargin)

%% ============================= VARARGIN MANAGER =========================
if mod(length(varargin), 2) ~= 0
    error('Check input variables, mismatch field - value')
end

for i = 1:2:length(varargin)
    name  = varargin{i};
    value = varargin{i+1};

    switch name
        case 'freqlimits'
            freqLimits = value;
        case 'type' % Sync or Desync
            type = value;
        case 'cond1'
            cond1 = value;
        case 'cond2'
            cond2 = value;
        case 'isContra1'
            isContra1 = value;
        case 'isContra2'
            isContra2 = value;
        case 'searchWord1'
            searchWord1 = value;
        case 'searchWord2'
            searchWord2 = value;
        case 'timelimitsttest'
            timeLimitsTtest = value;
        case 'rename_ttest'
            renameTtest = value;
        case 'tthreshold'
            thr = value;
        case 'bandname'
            bandNames = value;
    end
end

%% ============================ MAIN ======================================

numConds = numel(type);
% Check if file was already computed
ttestFile = find_presults_bst(sub, renameTtest, 1,'timefreq');
if isempty(ttestFile)
    
    fprintf("Found %d conditions, starting %d loops\n", numConds, numConds);
    
    % Compute ttest file
    sFiles1 = bst_process('CallProcess', 'process_select_files_timefreq', [], [], ...
        'subjectname',   sub, ...
        'condition',     cond1, ...
        'tag',           searchWord1, ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0, ...
        'outprocesstab', 'process2a');  
    sFiles1(~startsWith({sFiles1.Comment}, searchWord1)) = [];

    sFiles2 = bst_process('CallProcess', 'process_select_files_timefreq', [], [], ...
        'subjectname',   sub, ...
        'condition',     cond2, ...
        'tag',           searchWord2, ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0, ...
        'outprocesstab', 'process2b');  
    sFiles2(~startsWith({sFiles2.Comment}, searchWord2)) = [];
    
    % Process: t-test equal [70ms,148ms 1-30Hz]          H0:(A=B), H1:(A<>B)
    ttestFile = bst_process('CallProcess', 'process_test_parametric2', sFiles1, sFiles2, ...
        'timewindow',    timeLimitsTtest', ...
        'freqrange',     [1,30], ...
        'rows',          '', ...
        'isabs',         0, ...
        'avgtime',       0, ...
        'avgrow',        0, ...
        'avgfreq',       0, ...
        'matchrows',     1, ...
        'iszerobad',     1, ...
        'Comment',       renameTtest, ...
        'test_type',     'ttest_unequal', ...  % Student's t-test   (equal variance)        A,B~N(m,v)t = (mean(A)-mean(B)) / (Sx * sqrt(1/nA + 1/nB))Sx = sqrt(((nA-1)*var(A) + (nB-1)*var(B)) / (nA+nB-2)) df = nA + nB - 2
        'tail',          'two');  % Two-tailed
end

if isempty(ttestFile)
    ttestStruct = struct();
    return
end
try
    ttest = in_bst_data(ttestFile.FileName);
catch
    ttest = in_bst_data(ttestFile{1});
end
roiNames = ttest.RowNames;
numRois = numel(roiNames);

ttestStruct = struct();
ttestStruct.ttest = ttest;

for c = 1:numConds
    currBand = bandNames{c};
    currType = type{c};
    isType = false(numRois,1);
    meanTypeVal = zeros(numRois,1);
    currF = ttest.Freqs;
    freqsIdx = currF >= freqLimits(c,1) & currF <= freqLimits(c,2);
    for r = 1:numRois
        currP = squeeze(ttest.pmap(r,:,:));
        currT = squeeze(ttest.tmap(r,:,:));
        pmask = currP < thr;
        tmasked = currT.*pmask;
        tmasked(tmasked == 0) = nan;
        tmasked = tmasked(:,freqsIdx);
        if length(find(tmasked)) > 1
            if strcmpi(currType, 'Sync') && ~all(isnan(tmasked(:))) && any(tmasked(:) > 0)
                isType(r) = true;
                meanTypeVal(r) = mean(tmasked(tmasked > 0), 'omitmissing');
            elseif strcmpi(currType, 'Desync') && ~all(isnan(tmasked(:))) && any(tmasked(:) < 0)
                isType(r) = true;
                meanTypeVal(r) = mean(tmasked(tmasked < 0), 'omitmissing');
            end
        end
    end
    ttestStruct.(currBand) = struct();
    ttestStruct.(currBand).Type = currType;
    ttestStruct.(currBand).Times = ttest.Time;
    ttestStruct.(currBand).isActive = isType;
    ttestStruct.(currBand).Freqs = currF(freqsIdx);
    ttestStruct.(currBand).MeanValue = meanTypeVal;
    ttestStruct.(currBand).Thr = thr;
    
end
end