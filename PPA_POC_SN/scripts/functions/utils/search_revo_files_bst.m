function revoFiles = search_revo_files_bst(subName)
    % Get protocolInfo
    iProtocol = bst_get('iProtocol');
    anatDir = bst_get('ProtocolInfo', iProtocol).SUBJECTS;
    subDir = fullfile(anatDir, subName);
    if ~exist(subDir, 'dir')
        error(sprintf("Folder: %s does not exist!", subDir))
    end
    revoFiles = dir(fullfile(subDir, '*tex.mat'));

end