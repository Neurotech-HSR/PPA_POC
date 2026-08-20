function sourceFile = get_source_model(cond, sub,searchQuery)

%% ========================== VARAGIN MANAGER =============================
if nargin < 1
    error('get_source_model with no inputs')
end

if isempty(cond)
    error('No condition selected')
end

if isempty(sub)
    subStruct = bst_get('Subject');
    sub = subStruct.Name;
end

%% ========================== MAIN ========================================
if nargin < 3
    searchQuery = 'dSPM';
    protocolInfo = bst_get('ProtocolInfo');
    dataDir = protocolInfo.STUDIES;
    subFolder = fullfile(dataDir, sub, cond);
    sFiles = dir(subFolder);
    sFiles = sFiles(contains({sFiles.name}, searchQuery));
    if numel(sFiles) > 1
        sprintf(['Found many models for subject %s. ' ...
            'Starting computation of a new model, with personalized name... '], sub)
        sourceFile = [];
    elseif numel(sFiles) == 1
        sourceFile = fullfile(sFiles.folder, sFiles.name);
    else
        sourceFile = [];
    end

else
    sFiles = bst_process('CallProcess', 'process_select_files_results', [], [], ...
        'subjectname',   sub, ...
        'condition',     cond, ...
        'tag',           searchQuery, ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0, ...
        'outprocesstab', 'no');  % No
    if numel(sFiles) > 1
        sprintf(['Found many models with same searchQuery for subject %s. ' ...
            'Starting computation of a new model, with personalized name...'], sub)
        sourceFile = [];
    elseif numer(sFiles) == 1
        sourceFile = sFiles.FileName;
    else
        sprintf(['No source models with selected parameters. ' ...
            'Starting computation of a new model, with personalized name...'])
        sourceFile = [];
    end
end

end