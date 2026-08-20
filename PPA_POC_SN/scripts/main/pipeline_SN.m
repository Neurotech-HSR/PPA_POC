clear
close all
clc

%% Brainstorm setup

basepath = pwd;
scriptFolder = fileparts(basepath);
% Add Brainstorm and utilities
addpath(fullfile(fileDir,'..','..','Brainstorm_related/Brainstorm_current_version/brainstorm3'))
addpath(genpath(fullfile(fileDir,'..','functions')))
if ~brainstorm('status')
    brainstorm
end
% Check Braisntorm active protocol
protocolStruct = bst_get('ProtocolInfo');
% Check if Protocol Name is the expected one
protocolName = 'SN_POC';
if ~strcmp(protocolStruct.Comment, protocolName)
    error("Wrong protocol selected in Brainstorm")
end
protocolSubjects = bst_get('ProtocolSubjects');
numSubjects = numel(protocolSubjects.Subject);
dataDir = protocolStruct.STUDIES;
anatDir = protocolStruct.SUBJECTS;
bsSubjs = {protocolSubjects.Subject.Name};
bsSubjs(ismember(bsSubjs, {'Group_analysis'})) = [];

%% Configuration
[pipelineDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
pipelineType = 'SN'; 
configFile = fullfile(pipelineDir, sprintf('Config_%s.json',pipelineType));
config = jsondecode(fileread(configFile));
selProcess = config.Pipeline;
SUBJECT = config.Subject;
Timepoint = config.TimePoint;
pipelineStruct = config.(selProcess);
MainDir = config.MainDir;
Processes = pipelineStruct.Processes;
BsSubName = [SUBJECT, '_', Timepoint];
CONSTANTS = {sprintf('%s_%s', config.Subject, config.TimePoint), MainDir};

%% Single subject pipeline

src = pipelineStruct.Processes;
if ~iscell(src), src = num2cell(src); end
subName = config.Subject;
fprintf('\n=== Subject: %s ===\n', subName);

for p = 1:numel(src)
    proc = src{p};

    if isfield(proc, 'function')
        runProcess(proc, CONSTANTS);
    end
end

ff_fanfare;

%% Functions
function runProcess(proc, CONSTANTS)
    funcName   = proc.function;
    params     = proc.params;
    fh = str2func(funcName);
    fprintf('  -> %s\n', funcName);

    if isempty(params)
        fh(CONSTANTS{:});
    else
        extraArgs = struct2namevalue(params);
        fh(CONSTANTS{:}, extraArgs{:});
    end
end

function nv = struct2namevalue(s)
    fields = fieldnames(s);
    nv = cell(1, numel(fields) * 2);
    for i = 1:numel(fields)
        nv{2*i-1} = fields{i};
        nv{2*i}   = s.(fields{i});
    end
end
