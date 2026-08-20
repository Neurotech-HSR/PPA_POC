function chanFile = get_eeg_chanFile_bst(subName, eegName)
    iProtocol = bst_get('iProtocol');
    dataDir = bst_get('ProtocolInfo', iProtocol).STUDIES;
    subDir = fullfile(dataDir, subName, eegName);
    if ~exist(subDir, 'dir')
        error(sprintf("Folder: %s does not exist!", subDir))
    end
    chanFile = dir(fullfile(subDir, 'channel.mat'));
end