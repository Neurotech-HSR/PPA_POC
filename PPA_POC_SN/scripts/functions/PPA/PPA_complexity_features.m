function featStruct = PPA_complexity_features(trialMatrix, fs, varargin)
% PPA_COMPLEXITY_FEATURES  Calcola SampEn, HFD, LZ e PE per-trial su una ROI.
%
% Strategia (allineata alla versione electrode-level):
%   SampEn, HFD  — sull'intero segnale (stime più stabili con più campioni)
%   LZ, PE       — nella time window [tWin(1), tWin(2)] (fenomeno time-locked)
%
% INPUT
%   trialMatrix  : [nTrials x nSamples] — timeseries per-trial di una ROI
%   fs           : frequenza di campionamento (Hz)
%
% VARARGIN
%   'tWin'       : [1x2] finestra temporale in secondi per LZ e PE
%                  riferita all'onset (t=0 = primo campione se non specificato)
%                  default: intera trial (fallback se non fornita)
%   'tVec'       : [1 x nSamples] vettore temporale reale dei campioni
%                  se non fornito, costruito da 0 con passo 1/fs
%   'm'          : embedding dim SampEn      (default: 2)
%   'r'          : tolerance SampEn          (default: 0.2)
%   'kmax'       : kmax HFD                  (default: 8)
%   'pe_order'   : ordine PE                 (default: 3)
%   'pe_delay'   : delay PE                  (default: 1)
%
% OUTPUT
%   featStruct.(roiName) con campi:
%     .sampen  [nTrials x 1]
%     .hfd     [nTrials x 1]
%     .lz      [nTrials x 1]
%     .pe      [nTrials x 1]
%     .params  struct
%
% DIPENDENZE: sampen.m, HigFracDim.m, calc_lz_complexity.m, PE.m

%% ====================== DEFAULTS ========================================
tWin     = [];
tVec     = [];
m        = 2;
r        = 0.2;
kmax     = 8;
pe_order = 3;
pe_delay = 1;

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'twin',     tWin     = varargin{i+1};
        case 'tvec',     tVec     = varargin{i+1};
        case 'm',        m        = varargin{i+1};
        case 'r',        r        = varargin{i+1};
        case 'kmax',     kmax     = varargin{i+1};
        case 'pe_order', pe_order = varargin{i+1};
        case 'pe_delay', pe_delay = varargin{i+1};
    end
end

%% ====================== VALIDAZIONE =====================================
assert(ismatrix(trialMatrix) && size(trialMatrix,2) > 1, ...
    'trialMatrix deve essere [nTrials x nSamples]');
nTrials  = size(trialMatrix, 1);
nSamples = size(trialMatrix, 2);

% Costruisci vettore temporale se non fornito
if isempty(tVec)
    tVec = (0:nSamples-1) / fs;
end
assert(numel(tVec) == nSamples, 'tVec deve avere lunghezza nSamples');

% Maschera temporale per LZ e PE
if ~isempty(tWin)
    tMask = tVec >= tWin(1) & tVec <= tWin(2);
else
    tMask = true(1, nSamples);  % fallback: intera trial
end
nWin = sum(tMask);

%% ====================== CALCOLO PER-TRIAL ================================
sampenVals = nan(nTrials, 1);
hfdVals    = nan(nTrials, 1);
peVals     = nan(nTrials, 1);

for t = 1:nTrials
    sig_full   = trialMatrix(t, :);
    sig_window = sig_full(tMask);

    % --- SampEn — intero segnale ---
    try
        sampenVals(t) = sampen(sig_full, m, r, 'chebychev');
    catch ME
        warning('PPA_complexity_features: SampEn fallita trial %d — %s', t, ME.message);
    end

    % --- HFD — intero segnale ---
    try
        hfdVals(t) = HigFracDim(sig_full, kmax);
    catch ME
        warning('PPA_complexity_features: HFD fallita trial %d — %s', t, ME.message);
    end

    % --- PE — time window ---
    if nWin > 20
        try
            windowSize = nWin - pe_order * pe_delay;
            if windowSize >= pe_order + 1
                pe_out     = PE(sig_window, pe_delay, pe_order, windowSize);
                peVals(t)  = pe_out(end);
            end
        catch ME
            warning('PPA_complexity_features: PE fallita trial %d — %s', t, ME.message);
        end
    end
end

%% ====================== OUTPUT ==========================================
featStruct.sampen = sampenVals;
featStruct.hfd    = hfdVals;
featStruct.pe     = peVals;
featStruct.params = struct('m', m, 'r', r, 'kmax', kmax, ...
                           'pe_order', pe_order, 'pe_delay', pe_delay, ...
                           'fs', fs, 'tWin', tWin, 'nWin', nWin, ...
                           'nTrials', nTrials, 'nSamples', nSamples);
end