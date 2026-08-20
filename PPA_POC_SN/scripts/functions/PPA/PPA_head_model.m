% function PPA_head_model(subName, cond, headModelStruct)
function PPA_head_model(subName,mainDir, varargin)

    %% Varargin manager

    % Initialize
    getFrom = [];

    % Varargin
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'abilitate', abilitate = varargin{i+1};
            case 'cond', cond = varargin{i+1};
            case 'headModelStruct', headModelStruct = varargin{i+1};
            case 'getFrom', getFrom = varargin{i+1};
            case 'possibleFolders', possibleFolders = varargin{i+1};
        end
    end
    if ~abilitate
        fprintf("<strong>Function PPA_head_model skipped</strong>\n")
    return
    end

    assert(~isempty(cond), 'PPA_head_model: cond is mandatory')
    assert(~isempty(headModelStruct), 'PPA_head_model: headModelStruct is mandatory')
    [subStruct, iSubject] = bst_get('Subject', subName);
    dataDir = bst_get('ProtocolInfo').STUDIES;

    %% Check head model existence
    subDataDir = fullfile(dataDir, subName);
    if ~isempty(possibleFolders)
        headFiles = [];
        for f = 1:numel(possibleFolders)
            tmp = dir(fullfile(subDataDir, possibleFolders{f}, 'headmodel*.mat'));
            headFiles = [headFiles; tmp];
        end
    else
        headFiles = dir(fullfile(subDataDir, '*/headmodel*.mat'));
    end
    if isempty(headFiles)
        fprintf("No head model found, proceed to compute\n")
        PPA_compute_head_model(subName, cond, headModelStruct);
    else
        % Search file in current condition
        if any(endsWith({headFiles.folder}, cond))
            fprintf("Head model found. Skipping computation\n")
            return
        else
            fprintf("Head model found in other folder, copying in current directory\n")
            if isempty(getFrom)
                error("Folder from which obtain head model file not specified")
            end
            sourceFolder = fullfile(dataDir,subName, getFrom);
            targetFolder = fullfile(dataDir,subName, cond);
            PPA_copy_head_model(sourceFolder, targetFolder);
        end
    end
end
