function PPA_load_mri(subName,mainDir, varargin)

    %% Varargin manager
    abilitate = varargin{end};
    if ~abilitate
        fprintf("<strong>Funzione PPA_fmri_import saltata</strong>\n")
        return
    end
    
    [subStruct, iSubject] = bst_get('Subject', subName);

    %% Add subject if missing
    if isempty(subStruct)
        db_add_subject(subName, [], 0, 0);
        [subStruct, iSubject] = bst_get('Subject', subName);
    end

    %% Get original name and file lookup folder
    origSubname = split(subName, '_');
    try
        Timepoint = join({origSubname{end-1:end}}, '_');
        Timepoint = Timepoint{1};
    catch
        Timepoint = 'Task_T0';
    end
    origSubname = join({origSubname{1:end-2}}, '_');
    origSubname = origSubname{1};
    %% T1 check
    t1Idx = subStruct.iAnatomy;
    if isempty(t1Idx)
        % Get MRI file from data folder
        mriFile = dir(fullfile(mainDir, origSubname,Timepoint, 'anat', '*.nii'));
        if isempty(mriFile)
            error("No nifti files for subject %s", subName)
        elseif numel(mriFile) > 1
            error("Multiple files detected for subject %s", subName)
        end
        mriFile = fullfile(mriFile.folder, mriFile.name);
        [~,~] = import_mri(iSubject, mriFile, 'Nifti1', 1, [], 'mri');
        fprintf("<strong>MRI successfully imported. Adjust Nasion, LPA and RPA points\n</strong>\n")
        disp("Close figure to continue")
        waitfor(gcf);
        % Reload the subject
        [subStruct, iSubject] = bst_get('Subject', subName);
        iAnatomy = subStruct.iAnatomy;
        if isempty(iAnatomy)
            disp("Couldn't find correct index for MRI file. Set it up, then press F5");
            keyboard
        end
        % If T1 was loaded now, perform CAT12
        bst_call(@process_segment_cat12, 'ComputeInteractive', iSubject, iAnatomy)
    else
        % Check if CAT12 was computed: Schafer atlas should be present
        numSchafers = sum(contains({subStruct.Anatomy.Comment}, 'Schaefer'));
        if numSchafers == 0
            % Compute CAT12
            bst_call(@process_segment_cat12, 'ComputeInteractive', iSubject, subStruct.iAnatomy)
        end
    end

    %% Check MNI normalization
    MriFile = subStruct.Anatomy(subStruct.iAnatomy).FileName;
    mri = in_bst_data(MriFile);
    if isfield(mri.NCS, 'R')
        if isempty(mri.NCS.R)
            process_mni_normalize('ComputeInteractive', MriFile)
        end
    end
    clear mri;

    %% Check BEM surfaces
    bemSurfs = subStruct.Surface(startsWith({subStruct.Surface.Comment},'bem'));
    if isempty(bemSurfs)
        iAnatomy = subStruct.iAnatomy;
        bst_call(@process_generate_bem, 'ComputeInteractive', iSubject, iAnatomy)
    elseif numel(bemSurfs) == 3
        % Check all 3 BEMs are present
        if ~(sum(contains({bemSurfs.Comment}, 'head') == 1) && sum(contains({bemSurfs.Comment}, 'innerskull') == 1) && sum(contains({bemSurfs.Comment}, 'outerskull') == 1))
           fprintf("<strong>Anomalous BEMs detected, check\n</strong>")
           keyboard
        end
    else
        fprintf("<strong>Numero anomalo di superfici BEM rilevate. Controllare</strong>\n")
        keyboard
    end
    
    %% Check Revoscan
    texSurf = subStruct.Surface(endsWith({subStruct.Surface.FileName},'tex.mat') & strcmp({subStruct.Surface.SurfaceType},'Scalp'));
    if isempty(texSurf)
        fprintf("No Revoscan tex found, starting importation...\n")
        revoFile = dir(fullfile(mainDir, origSubname, Timepoint, 'revopoint', '*.obj'));
        if isempty(revoFile)
            error('No Revoscan surface found in folder %s', fullfile(mainDir, subName, Timepoint, 'revopoint'))
        end
        if numel(revoFile) > 1
            error("Found many Revoscan surfaces in folder %s", fullfile(mainDir, subName, Timepoint, 'revopoint'))
        end
        revoFile = fullfile(revoFile.folder, revoFile.name);
        bst_call(@import_surfaces, iSubject, revoFile, 'WFTOBJ', 0) 
        [subStruct, iSubject] = bst_get('Subject', subName);
        revoSurface = subStruct.Surface(strcmp({subStruct.Surface.SurfaceType}, 'Other')).FileName;
        newFileName = db_surface_type(revoSurface, 'Scalp');
        panel_protocols('SelectNode', [], newFileName);
        panel_protocols('UpdateNode',          'Subject', iSubject)
        tess_align_fiducials(newFileName)
        f = gcf;
        waitfor(f);
        % Revoscan files just loaded does not have compiled fields like
        % Curvature. To do so, open and close figure
        view_surface(newFileName)
        f = gcf;
        waitfor(f);
    end

end