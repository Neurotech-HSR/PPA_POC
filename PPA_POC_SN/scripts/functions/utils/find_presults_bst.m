function allFiles = find_presults_bst(subName, comment, isIntra, subtype)
    % Search for comparison files in current Brainstorm protocol.
    % Comparison files have in the filename p + file tipology (i.e presults,
    % ptimefreq, pdata, pmatrix)
    % INPUT
    %   -subName      --------> Name of the subject in Brainstorm. If
    %           empty, comparison files will be searche il all directories
    %   - comment     --------> Display name of the comparison in
%               Brainstorm
    %   - isIntra     --------> boolean, check (True) in the @intra folder
    %   - subType     --------> file type, results/data/matrix/timefreq
    dataDir = bst_get('ProtocolInfo').STUDIES;
    
    if isempty(subName)
        presultFiles = dir(fullfile(dataDir, '*/*',sprintf('p%s*.mat', subtype)));
    else
        presultFiles = dir(fullfile(dataDir,subName, '*',sprintf('p%s*.mat', subtype)));
    end
    
    allFiles = {};
    for f = 1:numel(presultFiles)
        if ~isIntra
            if endsWith(presultFiles(f).folder, 'intra')
                continue
            end
            tmp = in_bst_data(fullfile(presultFiles(f).folder,presultFiles(f).name));
            if contains(tmp.Comment, comment)
                allFiles{end+1} = fullfile(presultFiles(f).folder,presultFiles(f).name);
            end
        else
        tmp = in_bst_data(fullfile(presultFiles(f).folder,presultFiles(f).name));
            if contains(tmp.Comment, comment)
                allFiles{end+1} = fullfile(presultFiles(f).folder,presultFiles(f).name);
            end
        end
    end