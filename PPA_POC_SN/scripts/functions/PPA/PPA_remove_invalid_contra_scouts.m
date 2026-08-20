function PPA_remove_invalid_contra_scouts(subName, mainDir, varargin)
% Default values
atlasName = '';
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'atlasName', atlasName = varargin{i+1};
        case 'abilitate', abilitate = varargin{i+1};
    end
end

if ~abilitate
    sprintf("Function PPA_compute_ERSD_scouts skipped")
    return
end
% Obbligatori
assert(~isempty(atlasName), 'atlasName parameter is mandatory for extract_band_relative_powers');

%% Carica la superficie
[subStruct, iSubject] = bst_get('Subject', subName);
% Find central_low cortex
surfIdx = find(strcmp({subStruct.Surface.Comment}, 'central_15002V'));
surfStruct = in_bst_data(subStruct.Surface(surfIdx).FileName);
[invalidScouts, invalidIdx, report] = check_contra_projectable(subStruct.Surface(surfIdx).FileName, atlasName);
% Obtain atlas
atlas = surfStruct.Atlas(strcmp({surfStruct.Atlas.Name}, atlasName));
atlas.Scouts(invalidIdx) = [];
% Rename scouts
numScouts = numel(atlas.Scouts);
for s = 1:numScouts
    prevName = atlas.Scouts(s).Label;
    atlas.Scouts(s).Label = regexprep(prevName, '\d+', num2str(s));
end
surfStruct.Atlas(strcmp({surfStruct.Atlas.Name}, atlasName)) = atlas;
bst_save(file_fullpath(subStruct.Surface(surfIdx).FileName), surfStruct, 'v7')
fillBsSeeds(subStruct.Surface(surfIdx).FileName)
end