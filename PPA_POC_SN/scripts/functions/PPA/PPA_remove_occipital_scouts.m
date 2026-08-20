function PPA_remove_occipital_scouts(subName, mainDir, varargin)
% Default values
atlasName = '';
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'atlasName', atlasName = varargin{i+1};
        case 'abilitate', abilitate = varargin{i+1};
    end
end

if ~abilitate
    sprintf("Funzione PPA_compute_ERSD_scouts skipped")
    return
end
% Mandatory
assert(~isempty(atlasName), 'atlasName parameter is mandatory for extract_band_relative_powers');

%% Carica la superficie
[subStruct, iSubject] = bst_get('Subject', subName);
% Find central_low cortex
surfIdx = find(strcmp({subStruct.Surface.Comment}, 'central_15002V'));
surfStruct = in_bst_data(subStruct.Surface(surfIdx).FileName);
atlas = surfStruct.Atlas(strcmp({surfStruct.Atlas.Name}, atlasName));
removeIdx = [];
for s = 1:numel(atlas.Scouts)
    if endsWith(atlas.Scouts(s).Region, 'O')
        removeIdx = [removeIdx; s];
    end
end
if isempty(removeIdx)
    fprintf("No occipital scout found\n")
    return
else
    atlas.Scouts(removeIdx) = [];
    surfStruct.Atlas(strcmp({surfStruct.Atlas.Name}, atlasName)) = atlas;
    bst_save(file_fullpath(subStruct.Surface(surfIdx).FileName), surfStruct, 'v7')
    fillBsSeeds(subStruct.Surface(surfIdx).FileName)
end

end