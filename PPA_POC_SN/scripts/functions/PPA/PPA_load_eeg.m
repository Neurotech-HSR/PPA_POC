function chanRelFile = PPA_load_eeg(subName,mainDir, varargin)
% PPA_load_eeg - Carica e verifica i file EEG di task e rest in Brainstorm.
%
% Sintassi:
%   chanRelFile = PPA_load_eeg(subName, Name, Value, ...)
%
% Input:
%   subName         - Nome del soggetto in Brainstorm
%
% Name-Value pairs:
%   'mainDir'       - Path della cartella principale del soggetto
%   'taskName'      - Tag identificativo del task di attività
%   'restName'      - Tag identificativo del task di rest
%   'searchWordTask'- Pattern glob per trovare il file EEG del task
%   'searchWordRest'- Pattern glob per trovare il file EEG del rest
%   'toRefine'      - Boolean, esegui refine registration (default: false)
%   'searchChan'    - Nome cartella da cui copiare il channel file (default: [])
%
% Output:
%   chanRelFile     - Path del channel file EEG relativo al task

% PPA_load_eeg - Load and check Task and Rest EEG files in Brainstorm
%
% Sintax:
%   chanRelFile = PPA_load_eeg(subName, Name, Value, ...)
%
% Input:
%   subName         - Subject name in Brainstorm
%
% Name-Value pairs:
%   'mainDir'       - Subject folder path
%   'taskName'      - Task tag name
%   'restName'      - Rest tag name
%   'searchWordTask'- Pattern to find Task EEG file
%   'searchWordRest'- Pattern to find Rest EEG file
%   'toRefine'      - Boolean, perform refine registration (default: false)
%   'searchChan'    - folder from which get channel.mat file (default: [])
%
% Output:
%   chanRelFile     - Channel path file

    %% Parsing varargin
    toRefine       = false;
    taskName       = [];
    restName       = [];
    searchWordTask = [];
    searchWordRest = [];
    searchChan     = [];

    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'toRefine',       toRefine       = varargin{i+1};
            case 'taskName',       taskName       = varargin{i+1};
            case 'restName',       restName       = varargin{i+1};
            case 'searchWordTask', searchWordTask = varargin{i+1};
            case 'searchWordRest', searchWordRest = varargin{i+1};
            case 'searchChan',     searchChan     = varargin{i+1};
            case 'abilitate',     abilitate     = varargin{i+1};
            otherwise
                warning('PPA_load_eeg: unexpected parameter: %s', varargin{i});
        end
    end
    
    if ~abilitate
        fprintf("<strong>Function PPA_load_eeg skipped</strong>\n")
        return
    end
    %% Mandatory input validation
    assert(~isempty(mainDir),        'PPA_load_eeg: ''mainDir'' is mandatory.');
    assert(~isempty(taskName),       'PPA_load_eeg: ''taskName'' is mandatory.');
    assert(~isempty(restName),       'PPA_load_eeg: ''restName'' is mandatory.');
    assert(~isempty(searchWordTask), 'PPA_load_eeg: ''searchWordTask'' is mandatory.');
    assert(~isempty(searchWordRest), 'PPA_load_eeg: ''searchWordRest'' is mandatory.');
    
    %% Get original name and data folder
    origSubname = split(subName, '_');
    try
        Timepoint = join({origSubname{end-1:end}}, '_');
        Timepoint = Timepoint{1};
    catch
        Timepoint = 'Task_T0';
    end
    origSubname = join({origSubname{1:end-2}}, '_');
    origSubname = origSubname{1};

    %% Find EEG files
    subjEegDir = fullfile(mainDir, origSubname,Timepoint, 'eeg');

    imgEeglab = dir(fullfile(subjEegDir, sprintf("*%s*.set",searchWordTask)));
    if numel(imgEeglab) ~= 1
        error('PPA_load_eeg: found %d EEG files for task (expected1). Pattern: %s', ...
              numel(imgEeglab), searchWordTask);
    end
    imgFile = fullfile(imgEeglab.folder, imgEeglab.name);

    restEeglab = dir(fullfile(subjEegDir, sprintf("*%s*.set",searchWordRest)));
    if numel(restEeglab) ~= 1
        error('PPA_load_eeg: found %d EEG fiels for rest (expected 1). Pattern: %s', ...
              numel(restEeglab), searchWordRest);
    end
    restFile = fullfile(restEeglab.folder, restEeglab.name);

    %% Loading
    chanRelFile = load_eeg_files(subName, imgFile, restFile, toRefine);

    %% -----------------------------------------------------------------------
    function chanRelFile = load_eeg_files(subName, imgFile, restFile, toRefine)
    % Load EEG files in Brainstorm if not present. Rename condition
    % folders. Manage refine registration

        dataDir = bst_get('ProtocolInfo').STUDIES;

        %% Check files presence in Brainstorm
        loadedActivity = bst_process('CallProcess', 'process_select_files_data', [], [], ...
            'subjectname',   subName, ...
            'condition',     '', ...
            'tag',           taskName, ...
            'includebad',    0, ...
            'includeintra',  0, ...
            'includecommon', 0, ...
            'outprocesstab', 'no');
        if ~isempty(loadedActivity)
            loadedActivity(~cellfun(@(c) strncmp(c, [taskName ' '], length(taskName)+1), {loadedActivity.Comment})) = [];
        end

        loadedRest = bst_process('CallProcess', 'process_select_files_data', [], [], ...
            'subjectname',   subName, ...
            'condition',     '', ...
            'tag',           restName, ...
            'includebad',    0, ...
            'includeintra',  0, ...
            'includecommon', 0, ...
            'outprocesstab', 'no');
        if ~isempty(loadedRest)
            loadedRest(~cellfun(@(c) strncmp(c, [restName ' '], length(restName)+1), {loadedRest.Comment})) = [];
        end

        isActivity = ~isempty(loadedActivity);
        isRest     = ~isempty(loadedRest);

        %% Import missing
        if ~isActivity || ~isRest
            toRefine = true;
            fprintf('<strong>EEG files not found. Loading...</strong>\n');

            if ~isActivity
                fprintf('Loading Task data: %s\n', taskName);
                activityFiles = bst_process('CallProcess', 'process_import_data_epoch', [], [], ...
                    'subjectname',   subName, ...
                    'condition',     '', ...
                    'datafile',      {{imgFile}, 'EEG-EEGLAB'}, ...
                    'iepochs',       [], ...
                    'eventtypes',    taskName, ...
                    'createcond',    0, ...
                    'channelalign',  1, ...
                    'usectfcomp',    1, ...
                    'usessp',        1, ...
                    'freq',          [], ...
                    'baseline',      [], ...
                    'blsensortypes', 'EEG');                        %#ok<NASGU>
                bst_process('CallProcess', 'process_set_comment', activityFiles, [], ...
                    'tag',     taskName, ...
                    'isindex', 1);
            end

            if ~isRest
                fprintf('Loading Rest data: %s\n', restName);
                restFiles = bst_process('CallProcess', 'process_import_data_epoch', [], [], ...
                    'subjectname',   subName, ...
                    'condition',     '', ...
                    'datafile',      {{restFile}, 'EEG-EEGLAB'}, ...
                    'iepochs',       [], ...
                    'eventtypes',    restName, ...
                    'createcond',    0, ...
                    'channelalign',  1, ...
                    'usectfcomp',    1, ...
                    'usessp',        1, ...
                    'freq',          [], ...
                    'baseline',      [], ...
                    'blsensortypes', 'EEG');                        %#ok<NASGU>
                bst_process('CallProcess', 'process_set_comment', restFiles, [], ...
                    'tag',     restName, ...
                    'isindex', 1);
            end

        else
            fprintf('<strong>EEG files already present. Did you execute refine registration previously?</strong>\n');
            toRefine = askYesNo('Perform refine registration? (y/n): ');
        end

        %% Rename folders if necessary
        renameConditionFolder(subName, taskName, dataDir);
        renameConditionFolder(subName, restName, dataDir);

        %% Refine registration
        if toRefine
            fprintf('Starting refine registration for task: %s\n', taskName);
            tmp = bst_process('CallProcess', 'process_select_files_data', [], [], ...
                'subjectname',   subName, ...
                'condition',     '', ...
                'tag',           taskName, ...
                'includebad',    0, ...
                'includeintra',  0, ...
                'includecommon', 0, ...
                'outprocesstab', 'no');
            chanRelFile = tmp(1).ChannelFile;
            channel_align_manual(chanRelFile, 'EEG', 1);
            waitfor(gcf);
        else
            chanRelFile = [];
        end

        %% Copy channel file (optional)
        if ~isempty(searchChan)
            copyChanFile(subName, taskName, restName, searchChan, dataDir);
        end
    end

    %% -----------------------------------------------------------------------
    function renameConditionFolder(subName, expectedName, dataDir)

        files = bst_process('CallProcess', 'process_select_files_data', [], [], ...
            'subjectname',   subName, ...
            'condition',     '', ...
            'tag',           expectedName, ...
            'includebad',    0, ...
            'includeintra',  0, ...
            'includecommon', 0, ...
            'outprocesstab', 'no');

        if isempty(files)
            warning('renameConditionFolder: no files found for tag %s.', expectedName);
            return
        end
        % Keep files with expected name
        mask = cellfun(@(x) strncmp(x, [expectedName ' ('], length(expectedName)+2), {files.Comment});
        files = files(mask);

        currentName = files(1).Condition;
        if ~strcmp(currentName, expectedName)
            src = fullfile(dataDir, subName, currentName);
            dst = fullfile(dataDir, subName, expectedName);
            movefile(src, dst);
            fprintf('Changed folder name: %s -> %s\n', currentName, expectedName);
        end
        db_reload_subjects(files(1).iItem);
        % Reload all conditions for the subject
        [~, iSubject] = bst_get('Subject',subName);
        db_reload_conditions(iSubject);
    end 

    %% -----------------------------------------------------------------------
    function copyChanFile(subName, taskName, restName, searchChan, dataDir)

        srcEntry = dir(fullfile(dataDir, subName, searchChan, 'channel.mat'));
        if isempty(srcEntry)
            error('copyChanFile: channel.mat not found in in %s', ...
                  fullfile(dataDir, subName, searchChan));
        end
        srcFile = fullfile(srcEntry.folder, srcEntry.name);

        for condName = {taskName, restName}
            dstFile = fullfile(dataDir, subName, condName{1}, 'channel.mat');
            if isfile(dstFile), delete(dstFile); end
            copyfile(srcFile, dstFile);
            fprintf('Channel file copied in: %s\n', condName{1});
        end
    end 

    %% -----------------------------------------------------------------------
    function result = askYesNo(prompt)
        while true
            answer = input(prompt, 's');
            if strcmpi(answer, 'y')
                result = true;
                return
            elseif strcmpi(answer, 'n')
                result = false;
                return
            else
                disp('Risposta non valida. Digitare y o n.');
            end
        end
    end 

end 