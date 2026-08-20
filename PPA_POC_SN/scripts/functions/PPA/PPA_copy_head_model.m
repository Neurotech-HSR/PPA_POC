function PPA_copy_head_model(sourceFolder, targetFolder)


[~,cond,~] = fileparts(targetFolder);
% Get kernel from sourceFolder
kernel = dir(fullfile(sourceFolder, '*headmodel*.mat'));
if numel(kernel) ~=1
    error("Anomalous number of head model found")
end
kernel = fullfile(kernel.folder, kernel.name);

[sSrcStudy, ~] = bst_get('AnyFile', file_short(kernel));
[sDestStudies, iDestStudies] = bst_get('StudyWithSubject', sSrcStudy.BrainStormSubject);
% Find study condition
iTarget = iDestStudies(strcmp([sDestStudies.Condition], cond));
db_set_headmodel(file_short(kernel), iTarget)

end