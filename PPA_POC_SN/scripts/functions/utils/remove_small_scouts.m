function remove_small_scouts(surfPath, atlasName, minVerts, minArea)
% =========================================================================
% FUNZIONE PER RIMUOVERE DEGLI SCOUT DA UN ATLANTE SOTTO UNA CERTA SOGLIA
% DI AREA E VERTICI
% Input:
%       - surfPath --> fullpath alla superficie in cui è presente l'atlante
%       - atlasName --> nome dell'atlante. Può essere il nome completo, o
%       una parte dell'atlante (deve essere univoca la ricerca)
%       - minVerts --> numero minimo di vertici. Default: 3
%       - minArea --> superficie minima dello scout. Default: 5 cm^2
% =========================================================================
% =========================================================================
% Remove all scouts under a certan amount of vertices and size (cm^2)
% INPUT
%       - surfPath --> surface fullpath of the atlas
%       - atlasName --> atlas full or partial name. Must be unique
%       - minVerts --> numero minimo di vertici. Default: 3
%       - minArea --> superficie minima dello scout. Default: 5 cm^2
% =========================================================================
% ======================= INPUT MANAGER ===================================
    if nargin < 1
        error('Set surface path!')
    end
    
    if nargin < 2
        error('Set the name of the atlas!')
    end
    
    if nargin < 3 || ~isnumeric(minVerts)
        minVerts = 3; % at least 3 vertices
        if ~isnumeric(minVerts)
            sprintf("minVerts not set, " + ...
                "continue with default value (%d)", minVerts)
        else
            sprintf("minVerts not set, " + ...
                "continue with default value (%d)", minVerts)
        end
    end
    
    if nargin < 4 || ~isnumeric(minArea)
        minArea = 5; % cm^2
        if ~isnumeric(minArea)
            sprintf("minArea not set, " + ...
                "continue with default value (%d)", minArea)
        else
            sprintf("minArea not set, " + ...
                "continue with default value (%d)", minArea)
        end
    end
    
% ======================== LOADING ========================================
    surfStruct = in_bst_data(surfPath);
    
    atlases = surfStruct.Atlas;
    availableAtlases = {atlases.Name};
    atlasIdx = find(strcmpi(availableAtlases, atlasName));
    
    if isempty(atlasIdx)
        atlasIdx = find(contains(availableAtlases, atlasName));
    
        if isempty(atlasIdx)
            error('Atlas non trovato')
        end
    end
    
    % Check scouts field area-related 
    if ~isfield(surfStruct, 'VertArea')
        [~, VertArea] = tess_area(surfStruct.Vertices, surfStruct.Faces);
    else
        VertArea = surfStruct.VertArea;
    end
    
    scouts = surfStruct.Atlas(atlasIdx).Scouts;
    numScouts = numel(scouts);
    newScouts = repmat(db_template('Scout'), 0);
    
    for s = 1:numScouts
        currScout = surfStruct.Atlas(atlasIdx).Scouts(s);
        scoutArea = scout_area(currScout, surfStruct,VertArea);
        numVerts = length(currScout.Vertices);
    
        if ~(numVerts <= minVerts || scoutArea <= minArea)
            newScouts(end+1) = db_template('Scout');
            newScouts(end).Vertices = currScout.Vertices;
            newScouts(end).Label    = currScout.Label;
            newScouts(end).Function = currScout.Function;
            newScouts(end).Region   = currScout.Region;
            newScouts(end).Color    = currScout.Color;
        else
            fprintf("<strong>One scout does not satisfies requisites:\n</strong>")
            sprintf("Area --> %.2f, threshold: %.2f", scoutArea, minArea)
            sprintf("Number of vertices --> %d, threshold %d", numVerts, minVerts)
        end
    end
    
    surfStruct.Atlas(atlasIdx).Scouts = newScouts;
    bst_save(surfPath, surfStruct, 'v7');
end