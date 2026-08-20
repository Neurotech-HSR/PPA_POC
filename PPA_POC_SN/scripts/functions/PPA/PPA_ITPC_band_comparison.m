function itpcComp = PPA_ITPC_band_comparison(subName, itpcStruct, varargin)
% PPA_ITPC_BAND_COMPARISON  Contrasti TF per banda tra Task vs Rest e
%   Task vs Contralateral, condizionati al superamento di una soglia ITPC.
%
% LOGICA:
%   1. Gate: la ROI deve avere ITPC >= minValue in task nella banda
%      (usa itpcStruct già calcolato da PPA_ITPC + PPA_extract_significant_ITPCs)
%   2. Se gate passa: carica TF per-trial (task, rest, contra), filtra per
%      banda e time window, media per trial → [nTrials x 1]
%   3. Wilcoxon paired task vs rest  e  task vs contra
%
% INPUT
%   subName     : nome soggetto in Brainstorm
%   itpcStruct  : output di PPA_ITPC + PPA_extract_significant_ITPCs
%                 (deve contenere i campi .Theta, .Alpha, .Beta con .RelRois)
%
% VARARGIN
%   'taskCond'        : condizione task BST          (default: 'Immagine')
%   'restCond'        : condizione rest BST          (default: 'Rest')
%   'taskSearch'      : tag file TF task ipsi        (default: 'for_ITPC')
%   'restSearch'      : tag file TF rest             (default: 'for_ITPC')
%   'contraSearch'    : tag file TF contra           (default: 'for_ITPC')
%   'tWin'            : [1x2] time window secondi   (default: [0.07, 0.20])
%   'alpha'           : soglia significatività       (default: 0.05)
%   'minITPCvalue'    : soglia gate ITPC             (default: 0.3)
%   'bands'           : struct con bandName→[fLow fHigh]
%                       default: Theta[4,8] Alpha[8,13] Beta[13,30]
%
% OUTPUT
%   itpcComp.(safeName).(bandName).gatePass       logical
%   itpcComp.(safeName).(bandName).VsRest.p       double
%   itpcComp.(safeName).(bandName).VsRest.sig     logical
%   itpcComp.(safeName).(bandName).VsRest.medDiff double
%   itpcComp.(safeName).(bandName).VsRest.obsDir  +1/-1/0
%   itpcComp.(safeName).(bandName).VsContra.*     (stesso schema)
%
% NOTA: funzione standalone, richiamabile da script di plotting.

%% ====================== DEFAULTS ========================================
taskCond     = 'Immagine';
restCond     = 'Rest';
taskSearch   = 'for_ITPC';
restSearch   = 'for_ITPC';
contraSearch = 'for_ITPC';
tWin         = [0.07, 0.20];
alpha        = 0.05;
minITPCvalue = 0.3;

% Bande default
bands = struct();
bands.Theta = [4,  8];
bands.Alpha = [8,  13];
bands.Beta  = [13, 30];

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'taskCond',     taskCond     = varargin{i+1};
        case 'restCond',     restCond     = varargin{i+1};
        case 'taskSearch',   taskSearch   = varargin{i+1};
        case 'restSearch',   restSearch   = varargin{i+1};
        case 'contraSearch', contraSearch = varargin{i+1};
        case 'tWin',         tWin         = varargin{i+1};
        case 'alpha',        alpha        = varargin{i+1};
        case 'minITPCvalue', minITPCvalue = varargin{i+1};
        case 'bands',        bands        = varargin{i+1};
    end
end

bandNames = fieldnames(bands);
roiNames  = itpcStruct.roiNames;
nROI      = numel(roiNames);

%% ====================== CARICA FILE TF ==================================
fprintf('[ITPC_band_comparison] Carico file TF per %s...\n', subName)

tfTask  = load_tf_files(subName, taskCond,  taskSearch);
tfRest  = load_tf_files(subName, restCond,  restSearch);
tfContra= load_tf_files(subName, taskCond,  contraSearch, true);  % contra flag

%% ====================== MAIN LOOP =======================================
itpcComp = struct();

for r = 1:nROI
    rName  = roiNames{r};
    sName  = matlab.lang.makeValidName(rName);

    for b = 1:numel(bandNames)
        bName  = bandNames{b};
        fLims  = bands.(bName);

        %% STEP 1: Gate ITPC
        gatePass = false;
        if isfield(itpcStruct, bName) && isfield(itpcStruct.(bName), 'RelRois')
            gatePass = itpcStruct.(bName).RelRois(r);
        end

        itpcComp.(sName).(bName).gatePass     = gatePass;
        itpcComp.(sName).(bName).VsRest       = empty_result();
        itpcComp.(sName).(bName).VsContra     = empty_result();

        if ~gatePass
            continue
        end

        %% STEP 2: Estrai magnitudine TF per-trial nella banda e time window
        magTask  = extract_band_magnitude(tfTask,   rName, fLims, tWin);
        magRest  = extract_band_magnitude(tfRest,   rName, fLims, tWin);
        magContra= extract_band_magnitude(tfContra, rName, fLims, tWin);

        %% STEP 3: Wilcoxon paired
        itpcComp.(sName).(bName).VsRest   = wilcoxon_comparison(magTask, magRest,   alpha);
        itpcComp.(sName).(bName).VsContra = wilcoxon_comparison(magTask, magContra, alpha);
    end
end

end


%% ========================================================================
%                           FUNZIONI LOCALI
%% ========================================================================

function tfData = load_tf_files(subName, cond, searchTag, isContra)
% Carica tutti i file TF per-trial e li organizza per ROI.
% Output struct:
%   tfData.roiNames   {nROI x 1}
%   tfData.freqs      [1 x nFreq]
%   tfData.times      [1 x nTime]
%   tfData.(safeName) [nTrials x nFreq x nTime]  magnitudine

    if nargin < 4, isContra = false; end

    files = bst_process('CallProcess', 'process_select_files_timefreq', [], [], ...
        'subjectname',   subName, ...
        'condition',     cond, ...
        'tag',           searchTag, ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0, ...
        'outprocesstab', 'no');

    % Filtra per verso: escludi contralaterali se isContra=false e viceversa
    hasContraComment = startsWith({files.Comment}, 'Contra', 'IgnoreCase', true);
    if isContra
        files(~hasContraComment) = [];
    else
        files(hasContraComment)  = [];
    end

    if isempty(files)
        tfData = struct('roiNames', {{}}, 'freqs', [], 'times', []);
        return
    end

    % Leggi primo file per struttura
    tmp      = in_bst_data(files(1).FileName);
    freqs    = tmp.Freqs;
    times    = tmp.Time;
    nROI     = size(tmp.TF, 1);
    nFreq    = numel(freqs);
    nTime    = numel(times);
    nTrials  = numel(files);

    % Estrai nomi ROI (formato 'ROIname @ scout')
    roiNames = cell(nROI, 1);
    for r = 1:nROI
        raw = tmp.RowNames{r};
        if contains(raw, ' @')
            roiNames{r} = strtrim(extractBefore(raw, ' @'));
        else
            roiNames{r} = strtrim(raw);
        end
    end

    tfData.roiNames = roiNames;
    tfData.freqs    = freqs;
    tfData.times    = times;

    % Inizializza matrici per-ROI: [nTrials x nFreq x nTime]
    for r = 1:nROI
        sName           = matlab.lang.makeValidName(roiNames{r});
        tfData.(sName)  = nan(nTrials, nFreq, nTime);
    end

    % Riempi per ogni trial
    for t = 1:nTrials
        d = in_bst_data(files(t).FileName);
        for r = 1:nROI
            sName  = matlab.lang.makeValidName(roiNames{r});
            % d.TF: [nROI x nTime x nFreq] in BST → permuta a [1 x nFreq x nTime]
            slice  = abs(permute(d.TF(r, :, :), [1, 3, 2]));  % [1 x nFreq x nTime]
            tfData.(sName)(t, :, :) = slice;
        end
    end
end


function magVec = extract_band_magnitude(tfData, roiName, fLims, tWin)
% Estrae la magnitudine media nella banda [fLims] e time window [tWin]
% per ogni trial → [nTrials x 1].
% Restituisce [] se la ROI non è presente o tfData è vuoto.

    if isempty(tfData.roiNames)
        magVec = [];
        return
    end

    sName = matlab.lang.makeValidName(roiName);
    if ~isfield(tfData, sName)
        magVec = [];
        return
    end

    freqs = tfData.freqs;
    times = tfData.times;

    fIdx  = freqs >= fLims(1) & freqs <= fLims(2);
    tIdx  = times >= tWin(1)  & times <= tWin(2);

    if ~any(fIdx) || ~any(tIdx)
        magVec = [];
        return
    end

    % tfData.(sName): [nTrials x nFreq x nTime]
    slice  = tfData.(sName)(:, fIdx, :);      % [nTrials x nFband x nTime]
    slice  = slice(:, :, tIdx);               % [nTrials x nFband x nWin]
    magVec = mean(slice, [2, 3], 'omitnan');  % [nTrials x 1]
end


function res = wilcoxon_comparison(a, b, alpha)
% Wilcoxon signed-rank paired tra due vettori di magnitudini TF.
    res = empty_result();
    if isempty(a) || isempty(b)
        return
    end

    n = min(numel(a), numel(b));
    a = a(1:n);
    b = b(1:n);
    validIdx = ~isnan(a) & ~isnan(b);

    if sum(validIdx) < 5
        return
    end

    [p, ~]       = signrank(a(validIdx), b(validIdx));
    medDiff      = median(a(validIdx) - b(validIdx));
    obsDir       = sign(medDiff);

    res.p       = p;
    res.sig     = p < alpha;
    res.medDiff = medDiff;
    res.obsDir  = obsDir;
end


function res = empty_result()
    res = struct('p', NaN, 'sig', false, 'medDiff', NaN, 'obsDir', 0);
end