function [invalidScouts, invalidIdx, report] = check_contra_projectable(srcSurfFile, atlasName)

%% Load surface file
sSurf = in_tess_bst(srcSurfFile);

%% Find atlas
atlasIdx = strcmp({sSurf.Atlas.Name}, atlasName);
assert(any(atlasIdx), 'Atlas "%s" not found.', atlasName);
sAtlasStruct.Name   = sSurf.Atlas(atlasIdx).Name;
sAtlasStruct.Scouts = sSurf.Atlas(atlasIdx).Scouts;
sAtlas = sAtlasStruct;

%% bst_project_scouts_contra replica
hasCached = isfield(sSurf, 'tess2tess_interp') && ...
            all(isfield(sSurf.tess2tess_interp, {'Wmat_LR','Wmat_RL'})) && ...
            ~isempty(sSurf.tess2tess_interp.Wmat_LR);

hasSphereLR = isfield(sSurf, 'Reg') && isfield(sSurf.Reg, 'SphereLR') && ...
              isfield(sSurf.Reg.SphereLR, 'Vertices') && ...
              size(sSurf.Reg.SphereLR.Vertices,1) == size(sSurf.Vertices,1);

needsTemplate = ~hasCached && ~hasSphereLR;

if needsTemplate
    % Find template surface
    sSubjectDef = bst_get('Subject', 0);
    if ~isempty(strfind(srcSurfFile, '_low'))
        defCortexFile = 'tess_cortex_pial_low.mat';
    else
        defCortexFile = 'tess_cortex_pial_high.mat';
    end
    iCortexDef = find(~cellfun(@(c)isempty(strfind(c, defCortexFile)), {sSubjectDef.Surface.FileName}));
    assert(~isempty(iCortexDef), 'Superficie default %s non trovata.', defCortexFile);
    defSurfFile = sSubjectDef.Surface(iCortexDef(1)).FileName;

    % Project on template
    [nScoutProj, sSurfDef, sAtlas] = bst_project_scouts(srcSurfFile, defSurfFile, sAtlas, 0, 0);
    assert(nScoutProj > 0, 'Proiezione scout sul template fallita.');

    % lH/rH on template
    [rH, lH] = tess_hemisplit(sSurfDef);
else
    % No template: lH/rH on subject surface
    [rH, lH] = tess_hemisplit(sSurf);
    defSurfFile = [];
end

%% Check scouts
invalidScouts = {};
invalidIdx    = [];
report = struct('Label',{}, 'nVertices',{}, 'nLeft',{}, 'nRight',{}, ...
                'fracLeft',{}, 'fracRight',{}, 'isInvalid',{});

scouts = sAtlas.Scouts;
for i = 1:numel(scouts)
    verts  = scouts(i).Vertices(:);
    nLeft  = numel(intersect(verts, lH));
    nRight = numel(intersect(verts, rH));
    nTot   = numel(verts);

    isInvalid = ((nLeft > 0) && (nRight > 0)) | ((nLeft == 0) && (nRight == 0));

    report(i).Label     = scouts(i).Label;
    report(i).nVertices = nTot;
    report(i).nLeft     = nLeft;
    report(i).nRight    = nRight;
    report(i).fracLeft  = nLeft  / nTot;
    report(i).fracRight = nRight / nTot;
    report(i).isInvalid = isInvalid;

    if isInvalid
        invalidScouts{end+1} = scouts(i).Label; %#ok
        invalidIdx(end+1)    = i;               %#ok
    end
end

%% Report
fprintf('\n--- check_contra_projectable: atlas "%s" ---\n', atlasName);
fprintf('Total scouts: %d | Not projectable: %d\n\n', numel(scouts), numel(invalidScouts));
if ~isempty(invalidScouts)
    fprintf('%-40s  %6s  %5s  %5s\n', 'Label', 'nVert', '%L', '%R');
    fprintf('%s\n', repmat('-',1,62));
    for i = 1:numel(report)
        if report(i).isInvalid
            fprintf('%-40s  %6d  %4.1f%%  %4.1f%%\n', ...
                report(i).Label, report(i).nVertices, ...
                report(i).fracLeft*100, report(i).fracRight*100);
        end
    end
else
    fprintf('All scouts can be projected contralaterally.\n');
end
fprintf('\n');