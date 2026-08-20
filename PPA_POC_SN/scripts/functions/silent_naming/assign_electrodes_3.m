function assign_electrodes_3(subName, chanRelFile)
% ASSIGN_ELECTRODES_3
% Identifica gli elettrodi EEG sulla superficie Revopoint del soggetto.
%
% INPUT:
%   subName     : [str/char] nome del soggetto nel database Brainstorm
%   chanRelFile : path (full o relativo) del file canali EEG
%
% PIPELINE:
%   1. Rileva elettrodo azzurro (FCz) e centroidi verdi
%   2. Identifica 8 landmark verdi tramite PCA/assi geometrici
%   3. Procrustes sui landmark -> warpa template verde
%   4. Hungarian globale -> assegna tutti i verdi
%   5. Graph search (stima da vicini verdi) -> searchPos per i bianchi
%   6. Per ogni bianco: cubo locale -> DBSCAN su punti scuri -> cluster
%      più sferico -> posizione finale
%   7. Plot interattivo con correzione manuale (Ctrl+A + Shift+click)
%   8. Salvataggio su Brainstorm
% INPUT:
%   subName     : [str/char] subject name in Braisntorm database
%   chanRelFile : relative of full path of channel file
%
% PIPELINE:
%   1. Find FCz (cyan electrode, reference) and green electrodes
%   2. Identifiy 8 landmarks through PCA and geometrical axis
%   3. Procrustes on landmark -> warp green electrodes template
%   4. Global Hungarian -> assign green electrodes
%   5. Graph search (estimate from neighbours) -> estimate searchPos for white
%                   electrodes
%   6. For each white candidate: local cube --> DBSCAN on dark points --->
%      sèherical cluster -> final estimated position
%   7. Interactive plot for manual correction 
%   8. Saving in Brainstorm new electrodes position

% Find electrode on the Revoscan surface loaded in Brainstorm.

    if ~contains(path, '/gea/home3/dati/Brainstorm_Revopoint/brainstorm_250902_src/brainstorm3')
        addpath /gea/home3/dati/Brainstorm_Revopoint/brainstorm_250902_src/brainstorm3
    end
    if ~contains(path, '/gea/home3/dati/PPA_POC/matscripts')
        addpath(genpath('/gea/home3/dati/PPA_POC/matscripts'))
    end
    if ~brainstorm('status')
        brainstorm
    end

    revoFiles = search_revo_files_bst(subName);
    if numel(revoFiles) > 1
        disp('Many Revoscan surfaces found')
        keyboard
    end
    revoFilepath = fullfile(revoFiles(1).folder, revoFiles(1).name);
    revoStruct   = in_bst_data(revoFilepath);
    ChannelMat   = in_bst_data(chanRelFile);
    G = build_cap_graph_robust(ChannelMat.Channel);
    Vertices = revoStruct.Vertices;
    Faces    = revoStruct.Faces;

    % =========================================================
    %% 1. Cyan electrode (FCz) and green centroids
    % =========================================================
    cyanCenter  = get_cyan_eletctrode(revoStruct);
    greenCenters = get_green_eletctrodes(revoStruct);

    % Assign axis indexes
    xAx = 1; yAx = 2; zAx = 3;

    gX = greenCenters(:, xAx);
    gY = greenCenters(:, yAx);
    gZ = greenCenters(:, zAx);

    % =========================================================
    %% 2. Green landmarks
    % =========================================================

    % Cz: closest green centroid to the intersection of positive Z axis and mesh
    distFromZax = sqrt(gX.^2 + gY.^2);
    meshDistFromZax = sqrt(Vertices(:,xAx).^2 + Vertices(:,yAx).^2);
    meshZvals       = Vertices(:,zAx);
    % Considera solo la metà superiore della mesh (Z > mediana)
    zMeshMed = median(meshZvals);
    topMeshMask = meshZvals > zMeshMed;
    meshDistFromZax_top = meshDistFromZax;
    meshDistFromZax_top(~topMeshMask) = inf;
    [~, meshCzVtx] = min(meshDistFromZax_top);
    meshCzPos = Vertices(meshCzVtx, :);  % posizione [x,y,z] dell'intersezione

    % Closest green centroid --> Cz
    distToCzMesh = vecnorm(greenCenters - meshCzPos, 2, 2);
    [~, czIdx]   = min(distToCzMesh);

    % Fp1/Fp2: green centroids with maxima X, distinguished by Y values
    [~, fpCandidates] = sort(gX, 'descend');
    fpCandidates = fpCandidates(1:2);
    if gY(fpCandidates(1)) > gY(fpCandidates(2))
        fp1Idx = fpCandidates(1);
        fp2Idx = fpCandidates(2);
    else
        fp1Idx = fpCandidates(2);
        fp2Idx = fpCandidates(1);
    end

    % O1/Oz/O2: minimal X, reordered by Y
    gX_norm = (gX - min(gX)) / (max(gX) - min(gX));
    gZ_norm = (gZ - min(gZ)) / (max(gZ) - min(gZ));
    gY_norm = (gY - min(gY)) / (max(gY) - min(gY));  % 0=right, 1=left, 0.5=midline
    gY_mid  = abs(gY_norm - 0.5);                     % 0=midline, 0.5= lateral extreme
    
    % low X + low Z + near midline --> occipital
    occScore = gX_norm + gZ_norm + gY_mid;
    [~, iSortOcc] = sort(occScore, 'ascend');
    oIdxs = iSortOcc(1:3);
    [~, iSortY] = sort(gY(oIdxs), 'ascend');
    o2Idx = oIdxs(iSortY(1));   % Y minima  -> right
    ozIdx = oIdxs(iSortY(2));   % Y median -> midline
    o1Idx = oIdxs(iSortY(3));   % Y maxima -> left

    % T7/T8: X close to 0, Z > 15° percentile (esclude TP9/TP10)
    lateralMask    = abs(gX) < prctile(abs(gX), 40);
    zMedianLateral = prctile(abs(gZ), 15);
    t7Mask         = lateralMask & gZ > zMedianLateral;
    [~, t7Idx] = max(gY .* double(t7Mask) - 1e6*double(~t7Mask));
    [~, t8Idx] = min(gY .* double(t7Mask) + 1e6*double(~t7Mask));

    landmarkFoundIdxs = [fp1Idx, fp2Idx, o1Idx, ozIdx, o2Idx, czIdx, t7Idx, t8Idx];
    seedNamesRaw      = {'Fp1','Fp2','O1','Oz','O2','Cz','T7','T8'};
    foundLandmarks    = greenCenters(landmarkFoundIdxs, :);

    % Verify landmark plot
    figure
    patch('Vertices', Vertices, 'Faces', Faces, 'FaceColor', 'interp', ...
          'EdgeColor', 'none', 'FaceVertexCData', revoStruct.Color, 'FaceAlpha', 0.6)
    axis equal; hold on
    lmColors = {'y','g','c','m','w','b','r',[1 0.5 0]};
    for k = 1:numel(seedNamesRaw)
        scatter3(foundLandmarks(k,1), foundLandmarks(k,2), foundLandmarks(k,3), ...
                 75, lmColors{k}, 'filled')
        text(foundLandmarks(k,1), foundLandmarks(k,2), foundLandmarks(k,3), ...
             seedNamesRaw{k}, 'Color', lmColors{k}, 'FontSize', 14)
    end
    scatter3(cyanCenter(1), cyanCenter(2), cyanCenter(3), 150, 'k', 'filled')
    text(cyanCenter(1), cyanCenter(2), cyanCenter(3), 'FCz', 'Color','k','FontSize',10)
    title('Green Landmark Identified')

    % =========================================================
    %% 3. Setup channels and template
    % =========================================================
    channels = ChannelMat.Channel;

    greenChans = {'Fp1','Fp2','F7','FC5','F3','T7','TP9','P7','CP5','O1','Oz','O2', ...
                  'P3','C3','Cz','FC1','FC2','Fz','F4','F8','FC6','T8','C4','CP6',  ...
                  'P8','TP10','P4','CP2','POz','Pz','CP1'};

    % Channel colors
    for k = 1:numel(channels)
        if any(strcmpi(channels(k).Name, greenChans))
            channels(k).Color = 'green';
        else
            channels(k).Color = 'white';
        end
    end

    % Seed map
    seed_map = containers.Map('KeyType','char','ValueType','any');
    for k = 1:numel(seedNamesRaw)
        seed_map(seedNamesRaw{k}) = foundLandmarks(k,:);
    end

    % Template positions (green)
    template_pos = containers.Map('KeyType','char','ValueType','any');
    for k = 1:numel(channels)
        lbl = channels(k).Name;
        if any(strcmpi(lbl, greenChans))
            p = channels(k).Loc(:,1);
            template_pos(lbl) = p';
        end
    end

    % Template posizions (all channels)
    template_pos_all = containers.Map('KeyType','char','ValueType','any');
    for k = 1:numel(channels)
        lbl = channels(k).Name;
        p   = channels(k).Loc(:,1);
        if norm(p) > 1e-10
            template_pos_all(lbl) = p';
        end
    end

    % =========================================================
    %% 4. Procrustes on seeds -> warpa green template
    % =========================================================
    src = zeros(numel(seedNamesRaw), 3);
    dst = zeros(numel(seedNamesRaw), 3);
    valid_seed = true(numel(seedNamesRaw), 1);
    for k = 1:numel(seedNamesRaw)
        lbl = seedNamesRaw{k};
        if isKey(template_pos, lbl)
            src(k,:) = template_pos(lbl);
            dst(k,:) = foundLandmarks(k,:);
        else
            valid_seed(k) = false;
            warning('Seed %s not found in template.', lbl);
        end
    end
    src = src(valid_seed,:);
    dst = dst(valid_seed,:);

    [R, t] = rigid_transform_3D(src, dst);

    src_w = (R * src' + t * ones(1,size(src,1)))';
    res   = vecnorm(src_w - dst, 2, 2);
    fprintf('\nSeed Residuals Procrustes [m]:\n');
    validIdx = find(valid_seed);
    for k = 1:numel(validIdx)
        fprintf('  %s: %.4f\n', seedNamesRaw{validIdx(k)}, res(k));
    end

    % Warp template
    template_w = containers.Map('KeyType','char','ValueType','any');
    lbls_tmpl  = keys(template_pos);
    for k = 1:numel(lbls_tmpl)
        lbl = lbls_tmpl{k};
        pw  = R * template_pos(lbl)' + t;
        template_w(lbl) = pw';
    end

    % =========================================================
    %% 5. Global Hungarian -> assign green electrodes
    % =========================================================
    seed_idxs     = landmarkFoundIdxs;
    all_idxs      = 1:size(greenCenters,1);
    free_idxs     = setdiff(all_idxs, seed_idxs);
    detected_free = greenCenters(free_idxs, :);

    G_green  = build_green_graph(greenChans);
    assigned = graph_propagation(G_green, seed_map, template_w, detected_free);

    % Update green electrodes Loc in ChannelMat
    for k = 1:numel(channels)
        lbl = channels(k).Name;
        if ~any(strcmpi(lbl, greenChans)), continue; end
        if isKey(seed_map, lbl)
            channels(k).Loc(:,1) = seed_map(lbl)';
        elseif isKey(assigned, lbl)
            channels(k).Loc(:,1) = assigned(lbl)';
        else
            warning('Green electrode %s not assigned.', lbl);
        end
    end

    % Gather all greens
    all_green_labels = {};
    src_all = []; dst_all = [];
    for k = 1:numel(channels)
        lbl = channels(k).Name;
        if ~any(strcmpi(lbl, greenChans)), continue; end
        if ~isKey(template_pos, lbl), continue; end
        if isKey(seed_map, lbl)
            meas = seed_map(lbl);
        elseif isKey(assigned, lbl)
            meas = assigned(lbl);
        else
            continue
        end
        src_all = [src_all; template_pos(lbl)];
        dst_all = [dst_all; meas];
        all_green_labels{end+1} = lbl;
    end

    % Verify green by plot
    figure; hold on
    patch('Vertices', Vertices, 'Faces', Faces, 'FaceColor', 'interp', ...
          'EdgeColor', 'none', 'FaceVertexCData', revoStruct.Color, 'FaceAlpha', 0.4)
    axis equal
    assignedLabels = keys(assigned);
    for k = 1:numel(assignedLabels)
        lbl = assignedLabels{k};
        pos = assigned(lbl);
        scatter3(pos(1),pos(2),pos(3),60,'g','filled')
        text(pos(1),pos(2),pos(3),lbl,'Color','k','FontSize',8)
    end
    for k = 1:numel(seedNamesRaw)
        pos = foundLandmarks(k,:);
        scatter3(pos(1),pos(2),pos(3),100,'y','filled')
        text(pos(1),pos(2),pos(3),seedNamesRaw{k},'Color','y','FontSize',10)
    end
    scatter3(cyanCenter(1),cyanCenter(2),cyanCenter(3),150,'c','filled')
    text(cyanCenter(1),cyanCenter(2),cyanCenter(3),'FCz','Color','k','FontSize',10)
    title('Assegnamento elettrodi verdi'); xlabel('X'); ylabel('Y'); zlabel('Z')

    % =========================================================
    %% 5b. CHECK: if cyan electrode between Cz, Fz, FC1 and FC2?
    % =========================================================
    cyanCheckLabels = {'Cz','Fz','FC1','FC2'};
    cyanCheckPos    = zeros(0,3);
    for ci = 1:numel(cyanCheckLabels)
        lbl = cyanCheckLabels{ci};
        if isKey(seed_map, lbl)
            cyanCheckPos(end+1,:) = seed_map(lbl);
        elseif isKey(assigned, lbl)
            cyanCheckPos(end+1,:) = assigned(lbl);
        end
    end
    if ~isempty(cyanCheckPos)
        checkCenter = mean(cyanCheckPos, 1);
        checkRadius = 0.030;  % 3 cm intorno al centro
        inSphere = vecnorm(Vertices - checkCenter, 2, 2) < checkRadius;
        colorsInSphere = revoStruct.Color(inSphere, :);
        % Punti azzurri: B alto, G medio-alto, R basso
        cyanMask = colorsInSphere(:,1) < 0.4 & colorsInSphere(:,2) > 0.4 & colorsInSphere(:,3) > 0.5;
        nCyanPts = sum(cyanMask);
        if nCyanPts > 10
            fprintf('[CHECK 5] OK: found %d cyan point around Cz/Fz/FC1/FC2.\n', nCyanPts);
        else
            fprintf('[CHECK 5] WARN: few cyan points (%d) around Cz/Fz/FC1/FC2. Check assignation.\n', nCyanPts);
        end
    else
        fprintf('[CHECK 5] WARN: not every electrode among (Cz/Fz/FC1/FC2) have been assigned.\n');
    end

    % =========================================================
    %% 5c. GEOMETRICAL CHECK: FP1/FP2, T7/T8, O1/O2/Oz
    % =========================================================
    % Gather reference positions
    checkMap = containers.Map('KeyType','char','ValueType','any');
    allAssigned = [seed_map; assigned];   
    refLabels = {'Fp1','Fp2','T7','T8','O1','O2','Oz'};
    for ci = 1:numel(refLabels)
        lbl = refLabels{ci};
        if isKey(seed_map, lbl)
            checkMap(lbl) = seed_map(lbl);
        elseif isKey(assigned, lbl)
            checkMap(lbl) = assigned(lbl);
        end
    end

    fprintf('\n--- CHECK GEOMETRICI CANALI VERDI ---\n');
    % Raccoglie X di tutti i centroidi verdi per determinare "range massimo"
    % Gather X of all green centroids to determine maximum range
    allGX = greenCenters(:, xAx);
    xMax  = max(allGX);
    xMin  = min(allGX);

    % FP1 ans FP2: positive X and bewteen maximums
    for lbl = {'Fp1','Fp2'}
        lbl = lbl{1};
        if ~isKey(checkMap, lbl), fprintf('[CHECK 6] WARN: %s not assigned.\n', lbl); continue; end
        p = checkMap(lbl);
        px = p(xAx); py = p(yAx);
        if px > 0 && px > (xMax - 0.3*(xMax-xMin))
            fprintf('[CHECK 6] OK  : %s positive X=%.4f close to maxmimum (xMax=%.4f).\n', lbl, px, xMax);
        else
            fprintf('[CHECK 6] WARN: %s X=%.4f unlikely positive or close to maximum (xMax=%.4f).\n', lbl, px, xMax);
        end
    end

    % T7: positive Y and near maximum
    if isKey(checkMap,'T7')
        p  = checkMap('T7');
        py = p(yAx);
        allGY = greenCenters(:, yAx);
        yMax  = max(allGY);
        if py > 0 && py > (yMax - 0.25*(yMax - min(allGY)))
            fprintf('[CHECK 6] OK  : T7  positive Y=%.4f and close to maximum(yMax=%.4f).\n', py, yMax);
        else
            fprintf('[CHECK 6] WARN: T7  Y=%.4f unlikely positive or close to maximum (yMax=%.4f).\n', py, yMax);
        end
    else
        fprintf('[CHECK 6] WARN: T7 not assigned.\n');
    end

    % T8: negative Y and near minima
    if isKey(checkMap,'T8')
        p  = checkMap('T8');
        py = p(yAx);
        allGY = greenCenters(:, yAx);
        yMin  = min(allGY);
        if py < 0 && py < (yMin + 0.25*(max(allGY) - yMin))
            fprintf('[CHECK 6] OK  : T8  negative Y=%.4f and close to minima (yMin=%.4f).\n', py, yMin);
        else
            fprintf('[CHECK 6] WARN: T8  Y=%.4f unlikely to be negative or close to minima (yMin=%.4f).\n', py, yMin);
        end
    else
        fprintf('[CHECK 6] WARN: T8 not assigned.\n');
    end

    % O1, O2, Oz: X negative and similar, Z similar; O1.Y > Oz.Y > O2.Y
    occLabels = {'O1','Oz','O2'};
    occOk = true;
    occPos = zeros(3,3);
    for ci = 1:3
        lbl = occLabels{ci};
        if ~isKey(checkMap, lbl)
            fprintf('[CHECK 6] WARN: %s not assigned.\n', lbl);
            occOk = false;
        else
            occPos(ci,:) = checkMap(lbl);
        end
    end
    if occOk
        ox1 = occPos(1,xAx); ox2 = occPos(2,xAx); ox3 = occPos(3,xAx);
        oz1 = occPos(1,zAx); oz2 = occPos(2,zAx); oz3 = occPos(3,zAx);
        oy1 = occPos(1,yAx); oy2 = occPos(2,yAx); oy3 = occPos(3,yAx);
        xNeg  = all([ox1,ox2,ox3] < 0);
        xSim  = (max([ox1,ox2,ox3]) - min([ox1,ox2,ox3])) < 0.015;
        zSim  = (max([oz1,oz2,oz3]) - min([oz1,oz2,oz3])) < 0.020;
        yOrd  = (oy1 > oy2) && (oy2 > oy3);  % O1.Y > Oz.Y > O2.Y
        if xNeg,  fprintf('[CHECK 6] OK  : O1/Oz/O2 hanno X negativa.\n');
        else,     fprintf('[CHECK 6] WARN: O1/Oz/O2 non hanno tutte X negativa (X=[%.3f %.3f %.3f]).\n', ox1,ox2,ox3); end
        if xSim,  fprintf('[CHECK 6] OK  : O1/Oz/O2 hanno X simile (range=%.4f).\n', max([ox1,ox2,ox3])-min([ox1,ox2,ox3]));
        else,     fprintf('[CHECK 6] WARN: O1/Oz/O2 X non simile (range=%.4f).\n',   max([ox1,ox2,ox3])-min([ox1,ox2,ox3])); end
        if zSim,  fprintf('[CHECK 6] OK  : O1/Oz/O2 hanno Z simile (range=%.4f).\n', max([oz1,oz2,oz3])-min([oz1,oz2,oz3]));
        else,     fprintf('[CHECK 6] WARN: O1/Oz/O2 Z non simile (range=%.4f).\n',   max([oz1,oz2,oz3])-min([oz1,oz2,oz3])); end
        if yOrd,  fprintf('[CHECK 6] OK  : Ordine Y corretto: O1(%.3f) > Oz(%.3f) > O2(%.3f).\n', oy1,oy2,oy3);
        else,     fprintf('[CHECK 6] WARN: Ordine Y non corretto: O1(%.3f), Oz(%.3f), O2(%.3f). Atteso O1>Oz>O2.\n', oy1,oy2,oy3); end
    end
    fprintf('--- FINE CHECK GEOMETRICI ---\n\n');

    % =========================================================
    %% 6. Bianchi: stima posizione da vicini verdi + DBSCAN locale
    % =========================================================
    % whiteChans = setdiff({channels.Name}, greenChans, 'stable');
    whiteChans = setdiff({channels.Name}, greenChans);
    whiteChans = whiteChans(~cellfun(@isempty, whiteChans));

    % Grafo completo verde+bianco
    % G_full = build_full_graph(greenChans, whiteChans);
    G_full = G;

    % Mappa verdi assegnati
    green_assigned = containers.Map('KeyType','char','ValueType','any');
    for k = 1:numel(all_green_labels)
        green_assigned(all_green_labels{k}) = dst_all(k,:);
    end

    % Stima searchPos per ogni bianco tramite vicini verdi
    whiteChans_ordered = {};
    expectedPos_white  = [];
    for ki = 1:numel(whiteChans)
        lbl = whiteChans{ki};
        exp = estimate_position_from_green_neighbors(lbl, G_full, green_assigned, template_pos_all);
        if isempty(exp)
            fprintf('  SKIP %s: nessun vicino verde\n', lbl);
            continue
        end
        whiteChans_ordered{end+1} = lbl;
        expectedPos_white = [expectedPos_white; exp];
    end
    fprintf('Elettrodi bianchi da cercare: %d\n', numel(whiteChans_ordered));

    % Costruisci array channels_white con searchPos aggiornata
    channels_white = [];
    for k = 1:numel(whiteChans_ordered)
        lbl = whiteChans_ordered{k};
        idx = find(strcmpi({channels.Name}, lbl));
        if isempty(idx), continue; end
        ch = channels(idx);
        ch.Color = 'white';
        channels_white = [channels_white, ch];
    end

    % Chiama assign_electrodes_local con le nuove searchPos
    % Raccoglie posizioni verdi assegnate per il vincolo di prossimità
    green_lbl_list   = keys(green_assigned);
    green_pos_matrix = zeros(numel(green_lbl_list), 3);
    for gi = 1:numel(green_lbl_list)
        green_pos_matrix(gi,:) = green_assigned(green_lbl_list{gi});
    end

    assignedElec_white = assign_electrodes_local(revoStruct, channels_white, ...
        expectedPos_white, 'cubeSpan', 0.020, 'colorThr', 0.25, 'doPlot', false, ...
        'forbiddenPos', green_pos_matrix, 'minDistElec', 0.015);

    % Aggiorna ChannelMat bianchi
    for k = 1:numel(assignedElec_white)
        lbl = assignedElec_white(k).Name;
        disp(lbl)
        idx = find(strcmpi({channels.Name}, lbl));
        if isempty(idx)
            continue; 
        end
        channels(idx).Loc(:,1) = assignedElec_white(k).Pos(:);
    end
    ChannelMat.Channel = channels;

    % =========================================================
    %% 7. Raccogli tutti gli elettrodi per plot interattivo
    % =========================================================
    % Costruisci assignedElec completo (verdi + bianchi + FCz)
    assignedElec = struct('Name',{},'Pos',{},'Sphericity',{},'UsedThr',{},'Color',{});

    % Verdi
    for k = 1:numel(all_green_labels)
        lbl = all_green_labels{k};
        assignedElec(end+1).Name       = lbl;
        assignedElec(end).Pos          = dst_all(k,:);
        assignedElec(end).Sphericity   = NaN;
        assignedElec(end).UsedThr      = NaN;
        assignedElec(end).Color        = 'green';
    end

    % Bianchi
    for k = 1:numel(assignedElec_white)
        assignedElec(end+1) = assignedElec_white(k);
    end

    % FCz (azzurro)
    idx_fcz = find(strcmpi({channels.Name}, 'FCz'));
    if ~isempty(idx_fcz)
        assignedElec(end+1).Name       = 'FCz';
        assignedElec(end).Pos          = cyanCenter;
        assignedElec(end).Sphericity   = NaN;
        assignedElec(end).UsedThr      = NaN;
        assignedElec(end).Color        = 'cyan';
    end

    % =========================================================
    %% Plot interattivo con correzione manuale
    % =========================================================
    close all
    hFig = figure('Name', 'Verifica posizioni elettrodi', 'WindowState', 'maximized');
    hAx  = axes('Parent', hFig);
    patch('Parent', hAx, 'Vertices', Vertices, 'Faces', Faces, ...
          'FaceColor', 'interp', 'EdgeColor', 'none', ...
          'FaceVertexCData', revoStruct.Color, 'FaceAlpha', 0.9)
    axis(hAx, 'equal'); hold(hAx, 'on')

    for a = 1:numel(assignedElec)
        switch assignedElec(a).Color
            case 'green'; c = 'green';
            case 'cyan';  c = 'cyan';
            otherwise;    c = 'blue';
        end
        scatter3(assignedElec(a).Pos(1), assignedElec(a).Pos(2), assignedElec(a).Pos(3), ...
                 50, c, 'filled', 'Parent', hAx)
        text(assignedElec(a).Pos(1), assignedElec(a).Pos(2), assignedElec(a).Pos(3)+0.005, ...
             assignedElec(a).Name, 'FontSize', 10, 'Color', c, 'Parent', hAx)
    end
    title(hAx, 'Ruota la figura. Ctrl+A per correggere un elettrodo, R per rotate3D.')

    hFig.UserData.elecToFix    = [];
    hFig.UserData.assignedElec = assignedElec;
    hFig.UserData.Vertices     = Vertices;
    hFig.UserData.Faces        = Faces;
    hFig.UserData.hAx          = hAx;
    hFig.UserData.Color        = revoStruct.Color;
    hFig.UserData.rotateEnabled = false;
    rotate3d(hFig, 'off')

    uicontrol('Parent', hFig, 'Style', 'pushbutton', 'Units', 'normalized', ...
              'Position', [0.8 0.05 0.1 0.10], 'String', 'Save electrodes', ...
              'Callback', @(src,evt) returnElecCallback(src))
    set(hFig, 'WindowButtonDownFcn', @(src,event) clickCallback(src,event))
    set(hFig, 'KeyPressFcn',         @(src,event) keyPressCallback(src,event))
    waitfor(hFig)

    % =========================================================
    %% 8. Aggiorna ChannelMat finale e salva
    % =========================================================
    for c = 1:numel(channels)
        chanName    = channels(c).Name;
        assignedIdx = find(strcmpi({assignedElec.Name}, chanName));
        if isempty(assignedIdx), continue; end
        channels(c).Loc(:,1) = assignedElec(assignedIdx(1)).Pos(:);
    end
    ChannelMat.Channel = channels;

    fprintf('\nChannelMat pronto. Decommentare bst_save per salvare.\n')
    % keyboard
    bst_save(file_fullpath(chanRelFile), ChannelMat, 'v6')

    %% Chiedi per conferma che sia tutto corretto
    fprintf("<strong>Se è tutto corretto premi F5, altrimenti sposta manualmente gli elettrodi mancanti</strong>\n")
    keyboard
    %% Proietta sulla BEM
    % Assicurati che sia selezionata la BEM
    % Controlla ci siano le 3 BEM
    [subStruct, iSubject] = bst_get('Subject', subName);
    if (sum(contains({subStruct.Surface.FileName}, 'outerskull') == 1) && ...
        sum(contains({subStruct.Surface.FileName}, 'innerskull') == 1) && ...
        sum(contains({subStruct.Surface.FileName}, 'head_bem') == 1))
        % Accertati che siano selezionate
        % Innerskull
        panel_protocols('SelectNode', [], subStruct.Surface(contains({subStruct.Surface.FileName}, 'innerskull')).FileName);
        iSurface = find(contains({subStruct.Surface.FileName}, 'innerskull'));
        db_surface_default(iSubject, 'InnerSkull', iSurface);
        panel_protocols('RepaintTree');
        % Outerskull
        panel_protocols('SelectNode', [], subStruct.Surface(contains({subStruct.Surface.FileName}, 'outerskull')).FileName);
        iSurface = find(contains({subStruct.Surface.FileName}, 'outerskull'));
        db_surface_default(iSubject, 'OuterSkull', iSurface);
        panel_protocols('RepaintTree');
        % Head
        panel_protocols('SelectNode', [], subStruct.Surface(contains({subStruct.Surface.FileName}, 'head_bem')).FileName);
        iSurface = find(contains({subStruct.Surface.FileName}, 'head_bem'));
        db_surface_default(iSubject, 'Scalp', iSurface);
        panel_protocols('RepaintTree');
        panel_protocols('UpdateNode',          'Subject', iSubject)
    else
        error('Una o più BEM surfaces non trovate')
    end

    % Proietta i canali sulla superficie
    process_channel_project('Compute', chanRelFile, 'EEG');
    % =============================================================
    %% FUNZIONI INTERNE
    % =============================================================

    % ---------------------------------------------------------
    function G = build_green_graph(greenChans)
        edges = {
            'Fp1','F7';   'Fp1','F3';   'Fp2','F8';   'Fp2','F4';
            'F7','FC5';   'FC5','C3';   'FC5','T7';   'F3','FC1';
            'FC1','C3';   'FC1','Fz';   'F8','FC6';   'FC6','C4';
            'FC6','T8';   'F4','FC2';   'FC2','C4';   'FC2','Fz';
            'Fz','Cz';
            'T7','TP9';   'T7','CP5';   'TP9','P7';   'CP5','P7';
            'CP5','CP1';  'CP5','C3';   'C3','Cz';    'C3','CP1';
            'T8','TP10';  'T8','CP6';   'TP10','P8';  'CP6','P8';
            'CP6','CP2';  'CP6','C4';   'C4','Cz';    'C4','CP2';
            'Cz','CP1';   'Cz','CP2';
            'P7','O1';    'P7','P3';    'P8','O2';    'P8','P4';
            'P3','CP1';   'P3','Pz';    'P3','POz';   'P4','CP2';
            'P4','Pz';    'P4','POz';
            'O1','Oz';    'O2','Oz';    'Oz','POz';
            'Pz','POz';   'CP1','Pz';   'CP2','Pz';
            'Fp1','Fp2';  'F3','F4';    'F7','F8';
        };
        keep = cellfun(@(a,b) any(strcmpi(a,greenChans)) && any(strcmpi(b,greenChans)), ...
                       edges(:,1), edges(:,2));
        edges = edges(keep,:);
        G.edges = edges;
        G.nodes = greenChans;
        G.adj   = containers.Map('KeyType','char','ValueType','any');
        for i = 1:size(edges,1)
            a = edges{i,1}; b = edges{i,2};
            if ~isKey(G.adj,a), G.adj(a) = {}; end
            if ~isKey(G.adj,b), G.adj(b) = {}; end
            G.adj(a) = [G.adj(a), {b}];
            G.adj(b) = [G.adj(b), {a}];
        end
    end

    % ---------------------------------------------------------
    function G = build_full_graph(greenChans, whiteChans)
        allChans = [greenChans, whiteChans];
        edges = {
            'Fp1','AF7';  'Fp2','AF8';  'Fp1','AF3';  'Fp2','AF4';
            'AF7','F7';   'AF8','F8';   'AF3','F3';   'AF4','F4';
            'AF7','AF3';  'AF8','AF4';  'AF3','AF4';
            'F7','FT7';   'F8','FT8';
            'FT7','FC5';  'FT8','FC6';  'FT7','T7';   'FT8','T8';
            'T7','TP7';   'T8','TP8';
            'TP7','TP9';  'TP8','TP10'; 'TP7','CP5';  'TP8','CP6';
            'TP7','P7';   'TP8','P8';
            'P7','PO7';   'P8','PO8';
            'PO7','PO3';  'PO8','PO4'; 'PO7','O1';   'PO8','O2';
            'PO3','O1';   'PO4','O2';  'PO3','POz';  'PO4','POz';
            'PO3','P3';   'PO4','P4';
            'Fp1','Fpz';  'Fp2','Fpz'; 'Fpz','Fz';
            'Fz','FCz';   'FCz','Cz';  'Cz','CPz';
            'CPz','Pz';   'Pz','POz';
            'FC1','FCz';  'FC2','FCz';
            'CP1','CPz';  'CP2','CPz';
            % archi verdi (ridondanti ma necessari per connettività)
            'Fp1','F7';   'Fp1','F3';  'Fp2','F8';   'Fp2','F4';
            'F7','FC5';   'F8','FC6';  'F3','FC1';   'F4','FC2';
            'FC1','C3';   'FC2','C4';  'FC1','Fz';   'FC2','Fz';
            'T7','TP9';   'T8','TP10'; 'T7','CP5';   'T8','CP6';
            'C3','Cz';    'C4','Cz';   'C3','CP1';   'C4','CP2';
            'CP5','P7';   'CP6','P8';  'CP1','P3';   'CP2','P4';
            'P3','Pz';    'P4','Pz';   'P7','O1';    'P8','O2';
            'P3','POz';   'P4','POz';  'CP1','Pz';   'CP2','Pz';
            'Fz','Cz';    'Cz','Pz';
        };
        keep = cellfun(@(a,b) any(strcmpi(a,allChans)) && any(strcmpi(b,allChans)), ...
                       edges(:,1), edges(:,2));
        edges = edges(keep,:);
        G.edges    = edges;
        G.nodes    = allChans;
        % G.adj      = containers.Map('KeyType','char','ValueType','any');
        G.adj = containers.Map('KeyType','char','ValueType','any');
        G.is_green = containers.Map('KeyType','char','ValueType','logical');
        for i = 1:numel(allChans)
            G.adj(allChans{i}) = {};
        end
        for i = 1:size(edges,1)
            a = edges{i,1}; b = edges{i,2};
        
            G.adj(a) = [G.adj(a), {b}];
            G.adj(b) = [G.adj(b), {a}];
        end
        % for i = 1:numel(allChans)
        %     G.is_green(allChans{i}) = any(strcmpi(allChans{i}, greenChans));
        % end
        % for i = 1:size(edges,1)
        %     a = edges{i,1}; b = edges{i,2};
        %     if ~isKey(G.adj,a), G.adj(a) = {}; end
        %     if ~isKey(G.adj,b), G.adj(b) = {}; end
        %     G.adj(a) = [G.adj(a), {b}];
        %     G.adj(b) = [G.adj(b), {a}];
        % end
    end

    % ---------------------------------------------------------
    function [R, t] = rigid_transform_3D(A, B)
        centroid_A = mean(A,1);
        centroid_B = mean(B,1);
        AA = A - centroid_A;
        BB = B - centroid_B;
        H  = AA' * BB;
        [U,~,V] = svd(H);
        R = V * U';
        if det(R) < 0
            V(:,3) = -V(:,3);
            R = V * U';
        end
        t = (centroid_B - centroid_A * R')';  % 3x1
    end

    % ---------------------------------------------------------
    function assigned = graph_propagation(G, seed_map, template_w, detected)
        assigned = containers.Map(seed_map.keys, seed_map.values);
        in_queue = containers.Map(seed_map.keys, true(1, seed_map.Count));
        queue    = keys(seed_map);
        max_iter = numel(G.nodes) * 3;
        iter     = 0;

        % Fase 1: BFS per espandere la coda (ordine di visita)
        while ~isempty(queue) && iter < max_iter
            iter = iter + 1;
            % Nodo più vincolato (max vicini già assegnati)
            best_node = ''; best_count = -1;
            for qi = 1:numel(queue)
                nd = queue{qi};
                if ~isKey(G.adj, nd), continue; end
                n_a = sum(cellfun(@(nb) isKey(assigned,nb), G.adj(nd)));
                if n_a > best_count
                    best_count = n_a;
                    best_node  = nd;
                end
            end
            if isempty(best_node), break; end
            queue = queue(~strcmp(queue, best_node));
            for nb = G.adj(best_node)
                nb = nb{1};
                if isKey(assigned,nb) || isKey(in_queue,nb), continue; end
                if ~isKey(template_w,nb), continue; end
                queue{end+1} = nb;
                in_queue(nb) = true;
            end
        end

        % Fase 2: Hungarian globale
        unassigned_labels = {};
        for k = 1:numel(G.nodes)
            lbl = G.nodes{k};
            if isKey(assigned,lbl), continue; end
            if ~isKey(template_w,lbl), continue; end
            unassigned_labels{end+1} = lbl;
        end

        n_lbl = numel(unassigned_labels);
        n_det = size(detected,1);
        if n_lbl == 0 || n_det == 0
            warning('graph_propagation: nulla da assegnare in fase Hungarian.');
            return
        end

        cost = zeros(n_lbl, n_det);
        for i = 1:n_lbl
            ep = template_w(unassigned_labels{i});
            for j = 1:n_det
                cost(i,j) = norm(ep - detected(j,:));
            end
        end

        [row_ind, col_ind] = hungarian(cost);
        for i = 1:numel(row_ind)
            assigned(unassigned_labels{row_ind(i)}) = detected(col_ind(i),:);
        end
    end

    % ---------------------------------------------------------
    function [row_ind, col_ind] = hungarian(costMatrix)
        [nR, nC] = size(costMatrix);
        n = max(nR, nC);
        C = inf(n,n);
        C(1:nR,1:nC) = costMatrix;
        assignment = munkres_simple(C);
        row_ind = []; col_ind = [];
        for i = 1:nR
            j = assignment(i);
            if j > 0 && j <= nC
                row_ind(end+1) = i;
                col_ind(end+1) = j;
            end
        end
    end

    % ---------------------------------------------------------
    function assignment = munkres_simple(C)
        n = size(C,1);
        C = C - min(C,[],2);
        C = C - min(C,[],1);
        assignment = zeros(1,n);
        for iter = 1:100
            [rowCov, colCov, assignment] = cover_zeros(C, n);
            if sum(colCov) == n, break; end
            % mask_uncov = ~rowCov' & ~colCov;
            % uncov = C(~rowCov, mask_uncov);
            rowCov = rowCov(:)';   % forza 1×n
            colCov = colCov(:)';   % forza 1×n
            mask_uncov = ~colCov & ~rowCov;
            uncov = C(~rowCov, mask_uncov);
            minVal = min(uncov(:));
            if isinf(minVal), break; end
            C(~rowCov,:) = C(~rowCov,:) - minVal;
            C(:,colCov)  = C(:,colCov)  + minVal;
        end
    end

    % ---------------------------------------------------------
    function [rowCov, colCov, assignment] = cover_zeros(C, n)
        assignment = zeros(1,n);
        rowCov = false(1,n); colCov = false(1,n);
        for i = 1:n
            for j = 1:n
                if C(i,j)==0 && ~rowCov(i) && ~colCov(j)
                    assignment(i) = j;
                    rowCov(i) = true;
                    colCov(j) = true;
                end
            end
        end
        rowCov(:) = false;
        colCov(assignment(assignment>0)) = true;
    end

    % ---------------------------------------------------------
    function expected = estimate_position_from_green_neighbors(lbl, G, green_assigned, template_pos)
        expected = [];
        if ~isKey(G.adj, lbl), return; end
        if ~isKey(template_pos, lbl), return; end

        nbs       = G.adj(lbl);
        green_nbs = nbs(cellfun(@(nb) isKey(green_assigned,nb), nbs));
        if isempty(green_nbs), return; end

        p_tmpl    = template_pos(lbl);
        weights   = zeros(1, numel(green_nbs));
        positions = zeros(numel(green_nbs), 3);

        for i = 1:numel(green_nbs)
            nb        = green_nbs{i};
            vec_tmpl  = p_tmpl - template_pos(nb);
            positions(i,:) = green_assigned(nb) + vec_tmpl;
            weights(i)     = 1 / (norm(vec_tmpl) + 1e-9);
        end
        weights  = weights / sum(weights);
        expected = weights * positions;
    end

    % ---------------------------------------------------------
    function assignedElec = assign_electrodes_local(revoStruct, channels, newElecPos, varargin)
    % Assegna ogni elettrodo bianco: cerca cluster scuri (DBSCAN) nel cubo
    % centrato sulla searchPos stimata dal graph search.

        cubeSpan     = 0.015;
        colorThr     = 0.2;
        dbscanEps    = 0.002;
        dbscanMinPts = 5;
        maxColorThr  = 0.8;
        colorThrStep = 0.1;
        doPlot       = false;
        forbiddenPos = zeros(0,3);  % posizioni già assegnate (verdi + bianchi precedenti)
        minDistElec  = 0.015;       % distanza minima da elettrodi già assegnati [m]

        for i = 1:2:numel(varargin)
            switch varargin{i}
                case 'cubeSpan';     cubeSpan     = varargin{i+1};
                case 'colorThr';     colorThr     = varargin{i+1};
                case 'dbscanEps';    dbscanEps    = varargin{i+1};
                case 'dbscanMinPts'; dbscanMinPts = varargin{i+1};
                case 'maxColorThr';  maxColorThr  = varargin{i+1};
                case 'colorThrStep'; colorThrStep = varargin{i+1};
                case 'doPlot';       doPlot       = varargin{i+1};
                case 'forbiddenPos'; forbiddenPos = varargin{i+1};
                case 'minDistElec';  minDistElec  = varargin{i+1};
            end
        end

        Verts     = revoStruct.Vertices;
        Faces_loc = revoStruct.Faces;
        elecNames = {channels.Name}';
        assignedElec = struct('Name',{},'Pos',{},'Sphericity',{},'UsedThr',{},'Color',{});

        if doPlot
            hFigLoc = figure('Name','Assegnazione bianchi','WindowState','maximized');
        end

        for e = 1:numel(elecNames)
            elecName  = elecNames{e};
            searchPos = newElecPos(e,:);

            % Cubo di ricerca
            cubeIdx = find( ...
                Verts(:,1) > searchPos(1)-cubeSpan & Verts(:,1) < searchPos(1)+cubeSpan & ...
                Verts(:,2) > searchPos(2)-cubeSpan & Verts(:,2) < searchPos(2)+cubeSpan & ...
                Verts(:,3) > searchPos(3)-cubeSpan & Verts(:,3) < searchPos(3)+cubeSpan);

            % Fallback: proietta su superficie se cubo vuoto
            if isempty(cubeIdx)
                dirVec = searchPos / norm(searchPos);
                V_c    = Verts;
                proj   = V_c * dirVec';
                dist   = vecnorm(V_c - proj .* dirVec, 2, 2);
                dist(proj <= 0) = inf;
                [~, idx_nn] = min(dist);
                searchPos = Verts(idx_nn,:);
                cubeIdx = find( ...
                    Verts(:,1) > searchPos(1)-cubeSpan & Verts(:,1) < searchPos(1)+cubeSpan & ...
                    Verts(:,2) > searchPos(2)-cubeSpan & Verts(:,2) < searchPos(2)+cubeSpan & ...
                    Verts(:,3) > searchPos(3)-cubeSpan & Verts(:,3) < searchPos(3)+cubeSpan);
            end

            bestCenter = searchPos;
            bestCirc   = -Inf;
            usedThr    = colorThr;

            % Itera su TUTTI i threshold e raccoglie il cluster più sferico
            % tra tutti le combinazioni threshold/cluster trovate.
            % Questo previene che fili neri nell'intorno (colorThr basso)
            % vengano preferiti rispetto al cerchietto scuro corretto.
            allCandidates = struct('center',{},'sphericity',{},'thr',{},'nPts',{});

            for thr = colorThr : colorThrStep : maxColorThr
                colors_loc = revoStruct.Color(cubeIdx,:);
                colorMask  = all(colors_loc < thr, 2);
                verts_thr  = Vertices(cubeIdx(colorMask), :);
                if size(verts_thr,1) < dbscanMinPts, continue; end

                lbls_thr       = dbscan(verts_thr, dbscanEps, dbscanMinPts);
                uniqueLbls_thr = unique(lbls_thr);
                validLbls_thr  = uniqueLbls_thr(uniqueLbls_thr ~= -1);
                if isempty(validLbls_thr), continue; end

                for g = 1:numel(validLbls_thr)
                    clust = verts_thr(lbls_thr == validLbls_thr(g), :);
                    if size(clust,1) < 5, continue; end
                    [~, ~, latent] = pca(clust);
                    latent = latent / sum(latent);
                    % Circolarità 2D nel piano tangente alla superficie.
                    % latent è in ordine decrescente (PCA MATLAB).
                    % latent(3) ≈ 0 su superficie curva → min/max è sempre
                    % ~0 sia per ring circolari sia per fili. Usiamo invece
                    % il rapporto dei due autovalori nel piano (latent(2)/latent(1)):
                    %   cerchio: [.5 .5 0] → sph = 1
                    %   filo:    [1  0  0] → sph = 0
                    if latent(1) > 1e-9
                        sph = latent(2) / latent(1);
                    else
                        sph = 0;
                    end
                    cand.center     = mean(clust, 1);
                    cand.sphericity = sph;
                    cand.thr        = thr;
                    cand.nPts       = size(clust,1);
                    allCandidates(end+1) = cand;
                end
            end

            % Scegli il candidato con sfericità massima tra tutti i threshold
            if ~isempty(allCandidates)
                % --- Vincolo di prossimità -----------------------------------
                % Scarta candidati entro minDistElec da elettrodi già assegnati
                % (verdi passati via forbiddenPos + bianchi già processati).
                if ~isempty(forbiddenPos)
                    validCand = true(1, numel(allCandidates));
                    for ci2 = 1:numel(allCandidates)
                        dists = vecnorm(forbiddenPos - allCandidates(ci2).center, 2, 2);
                        if min(dists) < minDistElec
                            validCand(ci2) = false;
                        end
                    end
                    nFilt = sum(~validCand);
                    if nFilt > 0
                        fprintf('  [prox] %s: scartati %d candidati troppo vicini a elettrodi esistenti\n', ...
                            elecName, nFilt);
                    end
                    allCandidates = allCandidates(validCand);
                end
                % -------------------------------------------------------------
                if ~isempty(allCandidates)
                    [~, bestIdx] = max([allCandidates.sphericity]);
                    bestCenter   = allCandidates(bestIdx).center;
                    bestCirc     = allCandidates(bestIdx).sphericity;
                    usedThr      = allCandidates(bestIdx).thr;
                end
            end

            assignedElec(e).Name       = elecName;
            assignedElec(e).Pos        = bestCenter;
            assignedElec(e).Sphericity = bestCirc;
            assignedElec(e).UsedThr    = usedThr;
            assignedElec(e).Color      = 'white';

            % Aggiorna forbiddenPos con la posizione appena assegnata,
            % così i bianchi successivi non vi si sovrappongono.
            forbiddenPos = [forbiddenPos; bestCenter]; %#ok<AGROW>

            fprintf('%s → sph=%.3f thr=%.2f | pos=[%.4f %.4f %.4f]\n', ...
                elecName, bestCirc, usedThr, bestCenter(1), bestCenter(2), bestCenter(3))

            if doPlot
                clf(hFigLoc)
                patch('Vertices', Verts, 'Faces', Faces_loc, 'FaceColor','interp', ...
                      'EdgeColor','none','FaceVertexCData',revoStruct.Color,'FaceAlpha',.6)
                axis equal; hold on
                xlim(searchPos(1)+[-cubeSpan cubeSpan])
                ylim(searchPos(2)+[-cubeSpan cubeSpan])
                zlim(searchPos(3)+[-cubeSpan cubeSpan])
                % Mostra tutti i centroidi candidati trovati, colorati per threshold
                if ~isempty(allCandidates)
                    thrVals  = [allCandidates.thr];
                    thrUniq  = unique(thrVals);
                    cmap_plt = lines(numel(thrUniq));
                    for ci2 = 1:numel(allCandidates)
                        cIdx = find(thrUniq == allCandidates(ci2).thr, 1);
                        scatter3(allCandidates(ci2).center(1), ...
                                 allCandidates(ci2).center(2), ...
                                 allCandidates(ci2).center(3), ...
                                 40, cmap_plt(cIdx,:), 'filled')
                    end
                end
                scatter3(bestCenter(1),bestCenter(2),bestCenter(3),200,'r','filled')
                scatter3(searchPos(1),searchPos(2),searchPos(3),150,'y','filled')
                text(searchPos(1),searchPos(2),searchPos(3),elecName,'Color','y','FontSize',11)
                title(sprintf('%s | sph=%.3f | thr=%.2f', elecName, bestCirc, usedThr))
                drawnow; pause(0.3)
            end
        end
    end

    % ---------------------------------------------------------
    function clickCallback(src, event)
        if ~(ismember('shift', src.CurrentModifier) && strcmp(src.SelectionType,'extend'))
            return
        end
        if isempty(src.UserData.elecToFix), return; end
        clickPos = event.IntersectionPoint;
        elecIdx  = find(strcmp({src.UserData.assignedElec.Name}, src.UserData.elecToFix));
        if isempty(elecIdx), return; end
        src.UserData.assignedElec(elecIdx).Pos = clickPos;
        redraw_electrodes(src)
        src.UserData.elecToFix = [];
        title(src.UserData.hAx, 'Ruota la figura. Ctrl+A per correggere, R per rotate3D.')
    end

    % ---------------------------------------------------------
    function keyPressCallback(src, event)
        if strcmp(event.Key, 'r')
            src.UserData.rotateEnabled = ~src.UserData.rotateEnabled;
            if src.UserData.rotateEnabled
                rotate3d(src,'on');  disp('Rotate3D ON')
            else
                rotate3d(src,'off'); disp('Rotate3D OFF')
            end
            return
        end
        if ~isempty(event.Modifier) && any(strcmp(event.Modifier,'control')) && strcmp(event.Key,'a')
            elecToFix = inputdlg('Nome elettrodo da correggere:','Correzione',1);
            if isempty(elecToFix), return; end
            src.UserData.elecToFix = elecToFix{1};
            title(src.UserData.hAx, sprintf('Shift+click sulla nuova posizione per %s', elecToFix{1}))
        end
    end

    % ---------------------------------------------------------
    function redraw_electrodes(src)
        cla(src.UserData.hAx)
        patch('Parent', src.UserData.hAx, ...
              'Vertices', src.UserData.Vertices, 'Faces', src.UserData.Faces, ...
              'FaceColor','interp','EdgeColor','none', ...
              'FaceVertexCData',src.UserData.Color,'FaceAlpha',0.9)
        axis(src.UserData.hAx,'equal'); hold(src.UserData.hAx,'on')
        for a = 1:numel(src.UserData.assignedElec)
            switch src.UserData.assignedElec(a).Color
                case 'green'; c = 'green';
                case 'cyan';  c = 'cyan';
                otherwise;    c = 'blue';
            end
            scatter3(src.UserData.assignedElec(a).Pos(1), ...
                     src.UserData.assignedElec(a).Pos(2), ...
                     src.UserData.assignedElec(a).Pos(3), 50, c, 'filled', ...
                     'Parent', src.UserData.hAx)
            text(src.UserData.assignedElec(a).Pos(1), ...
                 src.UserData.assignedElec(a).Pos(2), ...
                 src.UserData.assignedElec(a).Pos(3)+0.005, ...
                 src.UserData.assignedElec(a).Name, 'FontSize',10,'Color',c, ...
                 'Parent', src.UserData.hAx)
        end
    end

    % ---------------------------------------------------------
    function returnElecCallback(src)
        fig = ancestor(src,'figure');
        assignedElec = fig.UserData.assignedElec;
        assignin('base','assignedElec',assignedElec);
        disp('assignedElec salvata nel workspace')
        close(fig)
    end


function G = build_cap_graph_robust(Channel)
% BUILD_CAP_GRAPH_ROBUST
% Grafo robusto per cuffia EEG 63 canali basato su coordinate reali.
%
% INPUT:
%   Channel : struct array Brainstorm con campi
%       .Name
%       .Loc   (3x1 o 3xN)
%
% OUTPUT:
%   G.nodes = labels
%   G.edges = cell Nx2
%   G.adj   = containers.Map
%
% Uso:
%   G = build_cap_graph_robust(Channel);

% ---------------------------------------------------------
% PARAMETRI
% ---------------------------------------------------------
kNN          = 6;      % vicini locali per nodo
maxDistScale = 1.75;   % tolleranza distanza locale
useSymmetry  = true;   % aggiunge specchi L/R

% ---------------------------------------------------------
% LABELS + POSIZIONI
% ---------------------------------------------------------
N = numel(Channel);

labels = cell(1,N);
XYZ    = zeros(N,3);

for i = 1:N
    labels{i} = strtrim(Channel(i).Name);

    p = Channel(i).Loc(:,1);
    p = p(:)';
    XYZ(i,:) = p ./ norm(p);   % normalizza su sfera
end

% ---------------------------------------------------------
% DISTANZE GEODESICHE
% ---------------------------------------------------------
D = zeros(N,N);

for i = 1:N
    for j = i+1:N
        c = dot(XYZ(i,:), XYZ(j,:));
        c = max(-1,min(1,c));
        ang = acos(c);
        D(i,j) = ang;
        D(j,i) = ang;
    end
end

% ---------------------------------------------------------
% BUILD EDGES KNN ROBUSTO
% ---------------------------------------------------------
E = {};

for i = 1:N

    di = D(i,:);
    di(i) = inf;

    [vals,idx] = sort(di,'ascend');

    localThr = vals(kNN) * maxDistScale;

    keep = idx(vals <= localThr);

    keep = keep(1:min(numel(keep),kNN));

    for j = keep
        E(end+1,:) = {labels{i}, labels{j}}; %#ok<AGROW>
    end
end

% ---------------------------------------------------------
% SIMMETRIA L/R
% ---------------------------------------------------------
if useSymmetry

    for i = 1:N
        a = labels{i};

        b = mirror_label(a);

        if isempty(b), continue; end

        j = find(strcmpi(labels,b),1);
        if isempty(j), continue; end

        % collega i due emisferi solo se simili in Z/X
        if abs(XYZ(i,1)-XYZ(j,1)) < 0.25
            E(end+1,:) = {a,b}; %#ok<AGROW>
        end
    end
end

% ---------------------------------------------------------
% RENDI UNDIRECTED + UNIQUE
% ---------------------------------------------------------
for i = 1:size(E,1)
    a = E{i,1};
    b = E{i,2};

    if strcmpi(a,b)
        E{i,1} = '';
        E{i,2} = '';
        continue
    end

    if strcmpi(a,b)
        continue
    end


    if strcmp(a,b) > 0
        E(i,:) = {b,a};
    end
end

E = E(~cellfun(@isempty,E(:,1)),:);
[~,ia] = unique(strcat(E(:,1),"__",E(:,2)),'stable');
E = E(ia,:);

% ---------------------------------------------------------
% ADJ MAP
% ---------------------------------------------------------
adj = containers.Map('KeyType','char','ValueType','any');

for i = 1:N
    adj(labels{i}) = {};
end

for i = 1:size(E,1)
    a = E{i,1};
    b = E{i,2};

    adj(a) = unique([adj(a), {b}],'stable');
    adj(b) = unique([adj(b), {a}],'stable');
end

% ---------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------
G.nodes = labels;
G.edges = E;
G.adj   = adj;

fprintf('Grafo robusto creato: %d nodi, %d archi\n', N, size(E,1));

end

% =========================================================
function out = mirror_label(lbl)

out = '';

tok = regexp(lbl,'^([A-Za-z]+)(\d+)$','tokens','once');

if isempty(tok)
    return
end

prefix = tok{1};
num = str2double(tok{2});

if mod(num,2)==1
    out = sprintf('%s%d',prefix,num+1);
else
    out = sprintf('%s%d',prefix,num-1);
end
end
end