function sFiles = create_link_files_brainstorm(files, model)
    sFiles = {};
    numFiles = numel(files);

    % Prepare Brainstorm base folder (e.g. subj/stimName)
    if ~isfield(files, 'FileName')
        lastFolder = files(1).folder;
        bsBaseFolder = split(lastFolder, '/');
        bsBaseFolder = [bsBaseFolder{end-1}, '/', bsBaseFolder{end}];
    
        for f = 1:numFiles
            unionString = ['link|', file_short(model), '|', bsBaseFolder, '/', files(f).name];
            sFiles{end+1} = unionString;
        end
    else    
        for f = 1:numFiles
            unionString = ['link|', file_short(model) '|', files(f).FileName];
            sFiles{end+1} = unionString;
        end
    end

end