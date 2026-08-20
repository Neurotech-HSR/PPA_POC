function PPA_load_channelmat(subName, mainDir, varargin)
% Load Brainstorm channelmat file.
% INPUT:
%   - subName: Name of the subject in Brainstorm
%   - mainDir: ''
%   - eegName:    --->  name of the eeg files to search in the database
%   - cond:       --->  name of the eeg condition (folder) in the database
%   - channelFile --->  full filepath of the channel.mat file presented in 
%                       the github folder


%% Varargin manager
cond = [];
eegName = [];
chFile = [];
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'eegName', eegName = varargin{i+1};
        case 'cond', cond = varargin{i+1};
        case 'channelFile', chFile = varargin{i+1};
        case 'abilitate', abilitate = varargin{i+1};
    end
end

if ~abilitate
    sprintf("Function PPA_compute_ERSD_scouts skipped")
    return
end
% Required parameters
assert(~isempty(cond), 'Parameter cond is required by extract_band_relative_powers');
assert(~isempty(eegName), 'Parameter eegName is required by extract_band_relative_powers');
assert(~isempty(chFile), 'Parameter chFile is required by extract_band_relative_powers');

%% Get EEG files
eegFiles = bst_process('CallProcess', 'process_select_files_data', [], [], ...
                'subjectname',   subName, ...
                'condition',     cond, ...
                'tag',           eegName, ...
                'includebad',    0, ...
                'includeintra',  0, ...
                'includecommon', 0, ...
                'outprocesstab', 'no');  % No

% Load channelFile

% Process: Set channel file
bst_process('CallProcess', 'process_import_channel', eegFiles, [], ...
    'channelfile',  {chFile, 'BST'}, ...
    'usedefault',   '', ...  % 
    'channelalign', 0, ...
    'fixunits',     0, ...
    'vox2ras',      0);
end