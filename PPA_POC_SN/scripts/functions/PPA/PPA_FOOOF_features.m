function fooofStruct = PPA_FOOOF_features(subName, condName, varargin)
% PPA_FOOOF_FEATURES  Estrae exponent e offset FOOOF per-trial dagli scout.
%
% Legge i file Brainstorm 'Scouts_FOOOF_Exponent_{cond}' e
% 'Scouts_FOOOF_Offset_{cond}' (o i corrispettivi Contralateral_Scouts_*)
% e restituisce una struct per-ROI con i vettori trial-by-trial.
%
% INPUT
%   subName   : nome soggetto in Brainstorm
%   condName  : cartella condizione (es. 'Immagine', 'Rest')
%
% VARARGIN
%   'prefix'      : prefisso nome file BST (default: 'Scouts')
%                   usare 'Contralateral_Scouts' per il controlaterale
%   'expTag'      : tag completo per exponent (override automatico)
%   'offTag'      : tag completo per offset   (override automatico)
%
% OUTPUT
%   fooofStruct.(safeName).exponent  [nTrials x 1]
%   fooofStruct.(safeName).offset    [nTrials x 1]
%   fooofStruct.roiNames             {nROI x 1}
%
% NOTA: funzione senza side-effects, richiamabile da script di plotting.

%% ====================== DEFAULTS ========================================
prefix = 'Scouts';
expTag = '';
offTag = '';

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'prefix', prefix = varargin{i+1};
        case 'exptag', expTag = varargin{i+1};
        case 'offtag', offTag = varargin{i+1};
    end
end

if isempty(expTag), expTag = sprintf('%s_FOOOF_Exponent_%s', prefix, condName); end
if isempty(offTag), offTag = sprintf('%s_FOOOF_Offset_%s',   prefix, condName); end
rename = split(condName,'_');
rename = rename{end};

%% ====================== CARICA FILE BST =================================
expFiles = bst_find_tf(subName, condName, expTag);
offFiles = bst_find_tf(subName, condName, offTag);

if isempty(expFiles)
    error('PPA_FOOOF_features: file exponent non trovato — tag: %s', expTag);
end
if isempty(offFiles)
    error('PPA_FOOOF_features: file offset non trovato — tag: %s', offTag);
end

%% ====================== ESTRAI PER-TRIAL ================================
% I file FOOOF estratti da process_extract_fooof sono timefreq con una
% singola "frequenza" (il parametro); dimensioni: [nROI x 1 x nTrials]
% oppure, per come BST salva i risultati scalari: F è [nROI x nTrials]

expVals = load_fooof_param(expFiles);  % [nROI x nTrials]
offVals = load_fooof_param(offFiles);  % [nROI x nTrials]

roiNames = expFiles.roiNames;
if size(roiNames,1) == 1; roiNames = roiNames'; end;
roiNames = rearrange_rowNames(roiNames);
nROI     = numel(roiNames);

%% ====================== ASSEMBLA STRUCT OUTPUT ==========================
fooofStruct          = struct();
fooofStruct.roiNames = roiNames;

for r = 1:nROI
    safeName = matlab.lang.makeValidName(roiNames{r});
    fooofStruct.(safeName).exponent = expVals(r, :)';  % [nTrials x 1]
    fooofStruct.(safeName).offset   = offVals(r, :)';
end

end


%% ====================== HELPER PRIVATI ==================================

function tfFiles = bst_find_tf(subName, cond, tag)
% Cerca file timefreq in BST; aggiunge campo .roiNames dal primo file.
    raw = bst_process('CallProcess', 'process_select_files_timefreq', [], [], ...
        'subjectname',   subName, ...
        'condition',     cond, ...
        'tag',           tag, ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0, ...
        'outprocesstab', 'no');

    if isempty(raw)
        tfFiles = [];
        return
    end

    % Filtra per match esatto del comment
    raw(~cellfun(@(c) strncmp(c, [tag ' '], length(tag)+1), {raw.Comment})) = [];

    if isempty(raw)
        tfFiles = [];
        return
    end

    % Leggi il primo file per estrarre i nomi delle ROI
    tmp = in_bst_data(raw(1).FileName);
    if isfield(tmp, 'RowNames')
        raw(1).roiNames = tmp.RowNames;
    else
        raw(1).roiNames = arrayfun(@(i) sprintf('ROI_%02d', i), ...
            1:size(tmp.TF,1), 'UniformOutput', false);
    end

    tfFiles      = raw;
end


function vals = load_fooof_param(tfFiles)
% Carica tutti i trial e li impila: output [nROI x nTrials]
    nTrials = numel(tfFiles);
    tmp0    = in_bst_data(tfFiles(1).FileName);
    nROI    = size(tmp0.TF, 1);
    vals    = nan(nROI, nTrials);

    for t = 1:nTrials
        d = in_bst_data(tfFiles(t).FileName);
        % TF è [nROI x 1 x 1] dopo process_extract_fooof → squeeze
        vals(:, t) = squeeze(d.TF(:, 1, 1));
    end
end
function roiNames = rearrange_rowNames(roiNames)
% Normalizza i nomi ROI: rimuove ' @ scout', suffissi ' L' e ' R'.
    roiNames = split(roiNames, ' @');
    roiNames = roiNames(:, 1);
    for r = 1:numel(roiNames)
        currR = roiNames{r};
        if endsWith(currR, ' R')
            tmp = split(currR, ' R'); roiNames{r} = tmp{1};
        elseif endsWith(currR, ' L')
            tmp = split(currR, ' L'); roiNames{r} = tmp{1};
        end
    end
end
