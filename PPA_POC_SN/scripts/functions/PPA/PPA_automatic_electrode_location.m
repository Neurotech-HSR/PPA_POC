function PPA_automatic_electrode_location(subName, mainDir, varargin)
% PPA_AUTOMATIC_ELECTRODE_LOCATION  Assigns/imports electrode positions
%   for a Brainstorm subject/condition, or copies an already-registered
%   channel file from another condition instead of recomputing it.
%
%   Parameters (varargin, name-value pairs):
%     condFolder  - (required) Brainstorm condition folder containing the
%                   channel file to create/update
%     abilitate   - bool, if false the function is skipped entirely
%                   (used to toggle pipeline steps on/off from a config)
%     getFrom     - (optional) condition folder to copy an existing,
%                   already-located channel.mat from, instead of running
%                   electrode assignment again (saves time when multiple
%                   conditions share the same cap/montage)

%% Varargin manager
% Defaults
condFolder = [];
getFrom = [];
abilitate = true;
for i = 1:2:length(varargin)
    switch varargin{i}
        case 'condFolder'
            condFolder = varargin{i+1};
        case 'abilitate'
            abilitate = varargin{i+1};
        case 'getFrom'
            getFrom = varargin{i+1};
    end
end

if ~abilitate
    fprintf("<strong>Function PPA_automatic_electrode_location skipped</strong>\n")
    return
end
assert(~isempty(condFolder), 'PPA_automatic_electrode_location: ''condFolder'' is required.');

%% Get channel file
% Brainstorm stores per-subject/condition data under STUDIES; build the
% expected path to this condition's channel file.
dataDir = bst_get('ProtocolInfo').STUDIES;
chanFile = fullfile(dataDir, subName, condFolder, 'channel.mat');
if ~isfile(chanFile)
    error("Path %s is not a file!", chanFile)
end

%% Main
if isempty(getFrom)
    % No source condition specified: run electrode assignment from
    % scratch on this condition's channel file.
    % assign_electrodes_2(subName, chanFile)
    assign_electrodes_3(subName, chanFile)
    % assign_electrodes(subName, chanFile)
else
    % getFrom specified: reuse electrode positions already computed for
    % another condition (e.g. same recording session/montage) instead of
    % rerunning the (slower) assignment routine.
    copyChan = fullfile(dataDir, subName, getFrom, 'channel.mat');
    targetChan = chanFile;
    copyfile(copyChan, targetChan)
end
end
