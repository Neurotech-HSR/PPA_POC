function PPA_fmri_import(subName,mainDir, varargin)

    %% Varagin manager
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'abilitate'
                abiltiate = varargin{i+1};
            case 'num_augmentation'
                num_augmentation = varargin{i+1};
            case 'area'
                area = varargin{i+1};
            case 'scoutsName'
                scoutsName = varargin{i+1};
            case 'atlasName'
                atlasName = varargin{i+1};
            case 'mniSearch'
                mniSearch = varargin{i+1};
            case 'roiAreaRange'
                roiAreaRange = varargin{i+1};
            case 'isContralateral'
                isContralateral = varargin{i+1};
        end
    end
    
    if ~abiltiate
        fprintf("<strong>Function PPA_fmri_import skipped</strong>\n")
        return
    end

    %% Get original folders and data folder
    origSubname = split(subName, '_');
    try
        Timepoint = join({origSubname{end-1:end}}, '_');
        Timepoint = Timepoint{1};
    catch
        Timepoint = 'Task_T0';
    end
    origSubname = join({origSubname{1:end-2}}, '_');
    origSubname = origSubname{1};
    %% fMRI check
    [subStruct, iSubject] = bst_get('Subject', subName);
    fmriFiles = subStruct.Anatomy(contains({subStruct.Anatomy.Comment},mniSearch));
    if isempty(fmriFiles)
        subjFMRIDir = fullfile(mainDir, origSubname, Timepoint, 'func');
        fprintf("<strong>No functional data detected, starting loading...</strong>\n")
        fmriFile = dir(fullfile(subjFMRIDir, sprintf('*%s*.nii',mniSearch)));
        if isempty(fmriFile)
            fprintf("<strong>No functional data found!</strong>\n")
            keyboard
        end
        fmriFile = fullfile(fmriFile.folder, fmriFile.name);
        [MriFileMni, sMriMni] = import_mri(iSubject, fmriFile, 'ALL-MNI-ATLAS', 0);
        [subStruct, iSubject] = bst_get('Subject', subName);
        fmriFiles = subStruct.Anatomy(contains({subStruct.Anatomy.Comment},mniSearch));
    end

    %% Surface Atlas check
    surfaceComment = 'central_15002V';
    surfPath = subStruct.Surface(strcmpi({subStruct.Surface.Comment}, surfaceComment)).FileName;
    surfStruct = in_bst_data(surfPath);
    % Atlas check
    isAtlas = any(strcmp({surfStruct.Atlas.Name}, atlasName));
    if ~isAtlas
        fprintf("<strong>No scout related to functional data detected, computing...</strong>\n")
        view_surface(surfPath)
        panel_scout('LoadScouts', file_fullpath(fmriFiles.FileName), 1)
        close(gcf);
        % Rename atlas
        surfStruct = in_bst_data(surfPath);
        rename_atlas(surfStruct, surfPath, atlasName);
        surfStruct = in_bst_data(surfPath);
        fmriAtlas = surfStruct.Atlas(end);
    else
        fmriAtlas = surfStruct.Atlas(strcmpi({surfStruct.Atlas.Name},atlasName));
    end
    scoutNames = {fmriAtlas.Scouts.Label};

    %% ROI expansion
    allAreas = [];
    [~, VertArea] = tess_area(surfStruct.Vertices, surfStruct.Faces);
    for s = 1:numel(scoutNames)
        allAreas = [allAreas; scout_area(fmriAtlas.Scouts(s),[],VertArea)];
    end
    perc10cm2 = sum(allAreas>=roiAreaRange(1) & allAreas<=roiAreaRange(2)) / length(allAreas) * 100;
    if perc10cm2 < 50
        fmriAtlas = grow_all_scouts(fmriAtlas, surfStruct, num_augmentation);
        fmriAtlas = subdivide_scouts_by_area(fmriAtlas, surfPath, area);
        % Remove occipital scouts
        fmriAtlas = remove_occipital_scouts(fmriAtlas);
        % Rename scouts
        fmriAtlas = rename_scouts(fmriAtlas, scoutsName);
        surfStruct.Atlas(strcmpi({surfStruct.Atlas.Name},atlasName)) = fmriAtlas;
        bst_save(file_fullpath(surfPath), surfStruct, 'v7')
        fillBsSeeds(file_fullpath(surfPath))
    end

    %% Additional occipatal scouts removal
    PPA_remove_occipital_scouts(subName, mainDir, 'abilitate', 1, 'atlasName', atlasName);
    bst_progress('Stop');

    %% Medial scouts removal
    PPA_remove_medial_int_scouts(surfPath, atlasName);
    %% Remove scouts that could not be projected contralaterally
    PPA_remove_invalid_contra_scouts(subName, mainDir, 'abilitate', 1, 'atlasName', atlasName);
    bst_progress('Stop');

    %% Ricarica la superficie per precauzione
    surfStruct = in_bst_data(surfPath);
    %% Check contralateral atlas
    if isContralateral
        surfStruct = in_bst_data(surfPath);
        contraName = sprintf("Contralateral_%s", atlasName);
        contraIdx = find(strcmp({surfStruct.Atlas.Name},contraName));
        if isempty(contraIdx)
            fmriAtlas = surfStruct.Atlas(strcmp({surfStruct.Atlas.Name}, atlasName));
            project_contralateral(surfStruct, surfPath, fmriAtlas, file_fullpath(surfPath));
            bst_progress('Stop');
        end
    end

    %% Sub-functions

    function rename_atlas(surfStruct, surfFilePath, atlasName)
        surfStruct.Atlas(end).Name = atlasName;
        bst_save(file_fullpath(surfFilePath), surfStruct, 'v7');
    end

    function SOatlas = grow_all_scouts(SOatlas, surfStruct, nRepeat)
        % GROW_ALL_SCOUTS  Expand all scouts on surface
        %
        %   SOatlas = grow_all_scouts(SOatlas, surfStruct)
        %
        %   Input:
        %     SOatlas  - .Scouts atlas struct (fields like Vertices, Seed, ...)
        %     surfStruct - struct surface with VertConn field (sparse, NxN) 
        %                  and Vertices (Nx3 coordinates)
        %     nRepeat    - Number of increments to apply for scouts growth
        %
        %   Output:
        %     SOatlas  - updated atlas
        
        if nargin<3
            nRepeat = 5;
        end
        patchVertices = surfStruct.Vertices;
        VertConn      = surfStruct.VertConn;
        
        for r = 1:nRepeat
            for i = 1:length(SOatlas.Scouts)
                vi      = SOatlas.Scouts(i).Vertices;
                seedXYZ = patchVertices(SOatlas.Scouts(i).Seed, :);
            
                % Expand through connectivity
                viNew = tess_scout_swell(vi, VertConn);
            
                if ~isempty(viNew)
                    % Seeds distance fro, candidates
                    newPts       = patchVertices(viNew, :);
                    distFromSeed = sqrt(sum((newPts - seedXYZ).^2, 2))';
            
                    % Spherical threshold: mean + 1.5*std of candidates
                    % distance
                    sphereRadius = mean(distFromSeed) + 1.5 * std(distFromSeed);
            
                    % Keep candiadtes within sphere
                    vi = union(vi, viNew(distFromSeed <= sphereRadius));
            
                    % Add isolated vertices surrounded by scout
                    iOut   = setdiff(1:size(VertConn,1), vi);
                    iOut   = setdiff(iOut, find(~any(VertConn)));   % exlude disconnected
                    iAlone = iOut(~any(VertConn(iOut, iOut)));
                    if ~isempty(iAlone)
                        iConn = find(VertConn(iAlone, :));
                        if all(ismember(iConn, vi))
                            vi = union(vi, iAlone);
                        end
                    end
                end
            
                SOatlas.Scouts(i).Vertices = vi;
            
                % If Scout is empty, use only the Seed
                if isempty(SOatlas.Scouts(i).Vertices)
                    SOatlas.Scouts(i).Vertices = SOatlas.Scouts(i).Seed;
                end
            end
        end
    
        % Merge overlapping scouts
        while true
            merged = false; 
        
            for s = 1:numel(SOatlas.Scouts)
                if s == numel(SOatlas.Scouts)
                    break
                end
                currScout    = SOatlas.Scouts(s);
                currVerts    = currScout.Vertices;
                joinIdxs     = false(numel(SOatlas.Scouts), 1);
                joinIdxs(s)  = true;
                joinVertices = currVerts;
                joinLabels   = currScout.Label;
        
                for j = s+1:numel(SOatlas.Scouts)
                    nextScout = SOatlas.Scouts(j);
                    nextVerts = nextScout.Vertices;
                    if ~isempty(intersect(joinVertices, nextVerts))
                        joinVertices = unique([joinVertices, nextVerts]);
                        joinIdxs(j)  = true;
                        joinLabels   = sprintf("%s & %s", joinLabels, nextScout.Label);
                    end
                end
        
                if sum(joinIdxs) > 1
                    newScout          = db_template('Scout');
                    newScout.Vertices = sort(joinVertices);
                    newScout.Label    = char(joinLabels);
                    newScout.Color    = rand(1,3);
                    center            = mean(patchVertices(newScout.Vertices,:), 1);
                    [~, imin]         = min(sum(bst_bsxfun(@minus, patchVertices(newScout.Vertices,:), center).^2, 2));
                    newScout.Seed     = newScout.Vertices(imin(1));
        
                    SOatlas.Scouts(joinIdxs) = [];
                    SOatlas.Scouts(end+1)    = newScout;
        
                    merged = true;
                    break
                end
            end
        
            if ~merged
                break
            end
        end
    end
    
    function SOatlas = subdivide_scouts_by_area(SOatlas, surfPath, targetArea_cm2)
    % SUBDIVIDE_SCOUTS_BY_AREA  Subdivide scouts by area threshold.
    %
    %   SOatlas = subdivide_scouts_by_area(SOatlas, surfStruct, targetArea_cm2)
    %
    %   Input:
    %     SOatlas       - atlas struct (.Scouts) (array con .Vertices, .Seed, .Label, ...)
    %     surfStruct      - Surface struct with .Vertices (Nx3), .VertConn (sparse NxN), .VertArea (Nx1, in m2)
    %     targetArea_cm2  - squared centimeters area target
    %
    %   Output:
    %     SOatlas  - updated atlas
    
    if nargin < 3
        targetArea_cm2 = 10;
    end
    surfStruct = in_bst_data(surfPath);
    sNewScouts = [];
    [~, VertArea] = tess_area(surfStruct.Vertices, surfStruct.Faces);
    
    for i = 1:length(SOatlas.Scouts)
        vi       = SOatlas.Scouts(i).Vertices;
        nVertices = length(vi);
    
        % total Area in squared metres
        totalArea = sum(VertArea(vi)) * 100 * 100;
        nClust    = round(totalArea / targetArea_cm2);
    
        % At least 5 vertices per scout
        if nClust > nVertices / 5
            nClust = floor(nVertices / 5);
        end
    
        % If scout cannot be subdivided, return it as it is
        if nClust <= 1
            sNewScouts = [sNewScouts, SOatlas.Scouts(i)];
            continue;
        end
    
        % Spectral clustering based on scout local connectivity
        ScoutConn = surfStruct.VertConn(vi, vi);
        Labels    = tess_cluster(ScoutConn, nClust);
        uniqueLabels = unique(Labels);
    
        % Generate new scout for each cluster
        for iLabel = 1:length(uniqueLabels)
            s            = db_template('Scout');
            s.Vertices   = vi(Labels == uniqueLabels(iLabel));
            s.Label      = sprintf('%s.%d', SOatlas.Scouts(i).Label, iLabel);
            s.Function   = SOatlas.Scouts(i).Function;
            s.Region     = SOatlas.Scouts(i).Region;
            s.Color      = SOatlas.Scouts(i).Color .* (1 - iLabel/length(uniqueLabels)/2);
            sNewScouts   = [sNewScouts, s];
        end
    end
    
    % Recompute Seeds
    SOatlas.Scouts = sNewScouts;
    fillBsSeeds(surfPath);
    SOatlas.Scouts = bst_setScoutRegion(SOatlas.Scouts, surfStruct.Atlas(strcmp({surfStruct.Atlas.Name},'Desikan-Killiany')).Scouts);
    
    end
    function SOatlas = remove_occipital_scouts(SOatlas)
        sScouts = SOatlas.Scouts;
        removeIdx = false(numel(sScouts),1);
        for s = 1:numel(sScouts)
            region = sScouts(s).Region;
            if endsWith(region, 'O')
                removeIdx(s) = true;
            end
        end
        SOatlas.Scouts(removeIdx) = [];
    end
    
    function SOatlas = rename_scouts(SOatlas, scoutsName)
        sScouts = SOatlas.Scouts;
        for s = 1:numel(sScouts)
            sScouts(s).Label = char(sprintf("%s_%d",scoutsName, s));
        end
        SOatlas.Scouts = sScouts;
    end
end