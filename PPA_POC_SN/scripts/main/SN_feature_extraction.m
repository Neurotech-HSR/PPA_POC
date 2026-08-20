clear; close all; clc;

%% Addpaths

[fileDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
% Brainstorm addpath
addpath(fullfile(fileDir,'..','..','Brainstorm_related/Brainstorm_current_version/brainstorm3'))
addpath(genpath(fullfile(fileDir,'..','functions')))
if ~brainstorm('status')
    brainstorm
end

protocolName = 'SN_POC';
protocolStruct = bst_get('ProtocolInfo');
if ~strcmp(protocolStruct.Comment, protocolName)
    error('Wrong Protocol set in Brainstorm')
end

%% Subjects
subjs = bst_get('ProtocolSubjects').Subject;
subjs(ismember({subjs.Name}, {'Group_analysis'})) = [];
numSubjs = numel(subjs);

%% Results dir
resultsDir = fullfile(fileDir, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir)
end

%% Compute features
for s = 1:numSubjs
    subName = subjs(s).Name;
    PPA_SN_ROI_features_v2(subName, resultsDir);
end