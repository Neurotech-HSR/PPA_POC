function itpcStruct = PPA_extract_significant_ITPCs(itpcStruct, bandLim, bandName, minValue)

% Find ROIs with a minimum value of ITPC in a specified band.
% INPUT:    
% itpcStruct --> Output struct from function calculate_ITPC_2
% bandLim    --> Band limits (e.g. [8, 13])
% bandName   --> Band name field in final struct
% minValue   --> Minimum ITPC value to find in time-frequence
% OUTPUT
% itpcStruct --> itpcStruct with ITPC related fields

%% ========================= VARARGIN MANAGER =============================
if nargin < 1
    error("Lanciata la funzione senza parametri!")
end
if nargin < 2
    error("Selezionare la banda di frequenza!")
end
if nargin < 3
    error("Nome della banda non inserito!")
end
if nargin < 4
    fprintf("<strong>Nessun valore di ITPC minimo inserito, uso il default (0.5)</strong>\n")
    minValue = 0.5;
end

%% ============================== MAIN ====================================
% For each ROI, filter in the band, find if at least one value is above
% threshold
roiNames = itpcStruct.roiNames;
numRois = numel(roiNames);
roiBool = false(numRois,1);
freqs = itpcStruct.freqs;
freqsIdx = freqs >= bandLim(1) & freqs <= bandLim(2);
values = zeros(numRois,1);
for r = 1:numRois
    currItpc = squeeze(itpcStruct.ITPC(:,freqsIdx,r));
    isRelevant = length(find(currItpc >= minValue)) > 0;
    if isRelevant
        mask = currItpc >= minValue;
        currItpc = currItpc.*mask;
        currItpc(currItpc == 0) = nan;
        roiBool(r) = true;
        values(r) = mean(currItpc, 'all','omitmissing');
    end

end

itpcStruct.(bandName) = struct();
itpcStruct.(bandName).Thr = minValue;
itpcStruct.(bandName).Value = values;
itpcStruct.(bandName).RelRois = roiBool;

end