function ITPC = PPA_ITPC(sub, mainDir, varargin)
% Calculate Inter Trial Phase Coherence following FieldTrip method
% (https://www.fieldtriptoolbox.org/faq/spectral/itc/)


%% ====================== VARARGIN MANAGER ================================
minT = 0.07; maxT = 0.15; useAllTime = 0;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'abilitate', abilitate = varargin{i+1};
        case 'cond', cond = varargin{i+1};
        case 'searchQuery', searchQuery = varargin{i+1};
        case 'minT', minT = varargin{i+1};
        case 'maxT', maxT = varargin{i+1};
        case 'useAllTime', useAllTime = varargin{i+1};
    end
end

if ~abilitate
    fprintf("Function PPA_ITPC skipped\n")
    return
end
assert(~isempty(cond), 'cond parameter is mandatory for PPA_ITPC');
assert(~isempty(searchQuery), 'searchQuery parameter is mandatory for PPA_ITPC');

%% ============================= INPUT FILES ==============================

% Find Time Frequency files
files = bst_process('CallProcess', 'process_select_files_timefreq', [], [], ...
    'subjectname',   sub, ...
    'condition',     cond, ...
    'tag',           searchQuery, ...
    'includebad',    0, ...
    'includeintra',  0, ...
    'includecommon', 0, ...
    'outprocesstab', 'no');

AllTF = [];

for f = 1:numel(files)
    tmpF = in_bst_data(files(f).FileName);

    if f == 1
        times = tmpF.Time;
        timeIdx = zeros(length(times),1);
        if useAllTime == 1 && isempty(minT) && isempty(maxT)
            minT = times(1);
            maxT = times(end);
        end

        timeIdx(times >= minT & times <= maxT) = 1;
        timeIdx = logical(timeIdx);
        times = times(timeIdx);
        freqs = tmpF.Freqs;
        nTime = sum(timeIdx);
        nFreq = length(tmpF.Freqs);
        nTrials = numel(files);
        roiNames = {};
        for r = 1:numel(tmpF.RowNames)
            currRow = tmpF.RowNames{r};
            roiName = extractBefore(currRow, ' @');
            roiNames{end+1} = roiName;
        end
        numRois = numel(roiNames);
        AllTF = zeros(nTime, nFreq, numRois, nTrials);

    end

    for r = 1:numRois
        phaseMap = angle(squeeze(tmpF.TF(r,:,:)))';   % freq × time
        phaseMap = phaseMap(:, timeIdx);             % freq × selected time

        % time × freq
        phaseMap = phaseMap.';                       % time × freq

        AllTF(:,:,r,f) = exp(1i * phaseMap);
    end
end

ITPCvals = abs(sum(AllTF, 4)) / nTrials; 

ITPC = struct();
ITPC.ITPC  = ITPCvals;
ITPC.times = times;
ITPC.freqs = freqs;
ITPC.roiNames = roiNames;

end



