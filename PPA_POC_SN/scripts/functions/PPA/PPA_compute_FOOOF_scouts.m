function PPA_compute_FOOOF_scouts(subName, mainDir, varargin)
%PPA_COMPUTE_FOOOF_SCOUTS  Computes FOOOF (specparam) on scout sources and
%   extracts the aperiodic parameters (exponent and offset).
%
%   Pipeline, per condition:
%     1. Look for a per-trial Welch PSD (Output:'all', not normalized)
%        -> name: Scouts_PSD_Welch_{cond}  or  Contralateral_Scouts_PSD_Welch_{cond}
%     2. If it doesn't exist, compute it from the scout files identified
%        by scoutsName
%     3. Apply process_fooof -> Scouts_FOOOF_{cond}
%     4. Extract exponent -> Scouts_FOOOF_Exponent_{cond}
%     5. Extract offset   -> Scouts_FOOOF_Offset_{cond}
%
%   Parameters (varargin, name-value pairs):
%     abilitate   - bool, if false the function is skipped
%     cond        - string, Brainstorm condition folder (e.g. 'Image', 'Rest')
%     scoutsName  - search tag for the source scout files
%     prefix      - (optional) prefix for output file names
%                   default = 'Scouts'; use 'Contralateral_Scouts' for the
%                   contralateral hemisphere

%% ====================== VARARGIN MANAGER ================================
abilitate  = false;
cond       = [];
scoutsName = [];
prefix     = 'Scouts';

for i = 1:2:length(varargin)
    switch varargin{i}
        case 'abilitate',  abilitate  = varargin{i+1};
        case 'cond',       cond       = varargin{i+1};
        case 'scoutsName', scoutsName = varargin{i+1};
        case 'prefix',     prefix     = varargin{i+1};
    end
end

if ~abilitate
    fprintf('Function PPA_compute_FOOOF_scouts skipped\n')
    return
end

assert(~isempty(cond),       'Parameter cond is required by PPA_compute_FOOOF_scouts');
assert(~isempty(scoutsName), 'Parameter scoutsName is required by PPA_compute_FOOOF_scouts');

%% ====================== FILE NAMES =======================================
% Build the search/output tags once, so every step below can check for
% (and reuse) an already-computed result under the same name.
welchName    = sprintf('%s_PSD_Welch_%s',      prefix, cond);
fooofName    = sprintf('%s_FOOOF_%s',          prefix, cond);
expName      = sprintf('%s_FOOOF_Exponent_%s', prefix, cond);
offName      = sprintf('%s_FOOOF_Offset_%s',   prefix, cond);

%% ====================== STEP 1: WELCH PSD ===============================
welchFile = find_bst_tf('timefreq', subName, cond, welchName);

if isempty(welchFile)
    fprintf('  [FOOOF] Welch PSD not found for %s - %s, computing...\n', subName, cond)

    scoutFiles = bst_process('CallProcess', 'process_select_files_matrix', [], [], ...
        'subjectname',   subName, ...
        'condition',     cond, ...
        'tag',           scoutsName, ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0, ...
        'outprocesstab', 'no');
    % Tag-based search can return partial matches; keep only files whose
    % comment starts exactly with "<scoutsName> " to avoid unrelated hits.
    scoutFiles(~cellfun(@(c) strncmp(c, [scoutsName ' '], length(scoutsName)+1), ...
        {scoutFiles.Comment})) = [];

    if isempty(scoutFiles)
        error('PPA_compute_FOOOF_scouts: no scout file found for %s in %s (tag: %s)', ...
            subName, cond, scoutsName)
    end

    % Welch PSD, per-trial (Output 'all'), used as input for FOOOF/specparam
    % fitting (which expects raw, non band-limited power spectra).
    welchFile = bst_process('CallProcess', 'process_psd', scoutFiles, [], ...
        'timewindow',  [], ...
        'win_length',  1, ...
        'win_overlap', 50, ...
        'units',       'physical', ...
        'win_std',     0, ...
        'edit',        struct( ...
            'Comment',         'Power', ...
            'TimeBands',       [], ...
            'Freqs',           [], ...
            'ClusterFuncTime', 'none', ...
            'Measure',         'power', ...
            'Output',          'all', ...
            'SaveKernel',      0));

    welchFile = bst_process('CallProcess', 'process_set_comment', welchFile, [], ...
        'tag',     welchName, ...
        'isindex', 1);

    fprintf('  [FOOOF] Welch PSD computed: %s\n', welchName)
else
    fprintf('  [FOOOF] Welch PSD already present: %s\n', welchName)
end

%% ====================== STEP 2: FOOOF ===================================
fooofFile = find_bst_tf('timefreq', subName, cond, fooofName);

if isempty(fooofFile)
    fprintf('  [FOOOF] Applying specparam on %s...\n', welchName)

    % Fit range 1-45 Hz, up to 3 peaks; 'fixed' aperiodic mode assumes no
    % spectral knee (appropriate for this frequency range).
    fooofFile = bst_process('CallProcess', 'process_fooof', welchFile, [], ...
        'implementation', 'matlab', ...
        'freqrange',      [1, 45], ...
        'powerline',      '-5', ...
        'method',         'leastsquare', ...
        'peakwidth',      [4, 30], ...
        'maxpeaks',       3, ...
        'minpeakheight',  3, ...
        'proxthresh',     2, ...
        'apermode',       'fixed', ...
        'guessweight',    'none', ...
        'sorttype',       'param', ...
        'sortparam',      'frequency', ...
        'sortbands',      {});

    fooofFile = bst_process('CallProcess', 'process_set_comment', fooofFile, [], ...
        'tag',     fooofName, ...
        'isindex', 1);

    fprintf('  [FOOOF] FOOOF computed: %s\n', fooofName)
else
    fprintf('  [FOOOF] FOOOF already present: %s\n', fooofName)
end

%% ====================== STEP 3: EXTRACT EXPONENT ========================
expFile = find_bst_tf('timefreq', subName, cond, expName);

if isempty(expFile)
    expFile = bst_process('CallProcess', 'process_extract_fooof', fooofFile, [], ...
        'fooof', 6);  % 6 = Exponent

    expFile = bst_process('CallProcess', 'process_set_comment', expFile, [], ...
        'tag',     expName, ...
        'isindex', 1);

    fprintf('  [FOOOF] Exponent extracted: %s\n', expName)
else
    fprintf('  [FOOOF] Exponent already present: %s\n', expName)
end

%% ====================== STEP 4: EXTRACT OFFSET ==========================
offFile = find_bst_tf('timefreq', subName, cond, offName);

if isempty(offFile)
    offFile = bst_process('CallProcess', 'process_extract_fooof', fooofFile, [], ...
        'fooof', 7);  % 7 = Offset

    offFile = bst_process('CallProcess', 'process_set_comment', offFile, [], ...
        'tag',     offName, ...
        'isindex', 1);

    fprintf('  [FOOOF] Offset extracted: %s\n', offName)
else
    fprintf('  [FOOOF] Offset already present: %s\n', offName)
end

end


%% ====================== HELPER ==========================================
function tfFile = find_bst_tf(fileType, subName, cond, tag)
%FIND_BST_TF  Look up a timefreq file in Brainstorm; returns [] if absent.
%   Centralizes the "search by tag, then strictly filter by exact tag
%   prefix" pattern reused across all steps above, so results aren't
%   polluted by files that merely contain the tag as a substring.
    tfFile = bst_process('CallProcess', ...
        sprintf('process_select_files_%s', fileType), [], [], ...
        'subjectname',   subName, ...
        'condition',     cond, ...
        'tag',           tag, ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0, ...
        'outprocesstab', 'no');

    if ~isempty(tfFile)
        tfFile(~cellfun(@(c) strncmp(c, [tag ' '], length(tag)+1), ...
            {tfFile.Comment})) = [];
    end
end
