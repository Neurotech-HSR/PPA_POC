%% === Interactive features GUI ===
% Visualizza le top-10 ROI sulla superficie corticale con features v2.
% Click su una ROI → pannello interattivo con plot per tipo di feature.
%
% Features supportate:
%   SampEn / HFD / PE      → boxplot tri-condizione (Task | Rest | Contra)
%   FOOOF_Exp / FOOOF_Off  → boxplot tri-condizione + spettro 1/f medio
%   ERSD_*_VsRest/Contra   → t-map banda×tempo (p-masked)
%   ITPC_*_VsRest/Contra   → TF time×freq con bande evidenziate
%
% Pannello destro (SEMPRE visibile, indipendente dalla feature selezionata):
%   - ERP medio della ROI  → media sui file matrix 'SN_Scout (#xx)' / Immagine
%   - ERSD media della ROI → media sui file timefreq 'Scouts_ERSD_Immagine (#xx)'
%   I file controlaterali ('Contralateral_*') sono esclusi.

clear; close all; clc

%% ====================== SETUP ===========================================
[fileDir, ~, ~] = fileparts(matlab.desktop.editor.getActiveFilename);
% Brainstorm addpath
addpath(fullfile(fileDir,'..','..','Brainstorm_related/Brainstorm_current_version/brainstorm3'))
addpath(fullfile(fileDir,'..','functions'))
if ~brainstorm('status')
    brainstorm
end
protocolInfo = bst_get('ProtocolInfo');
if ~strcmp(protocolInfo.Comment, 'SN_POC')
    error('Wrong protocol loaded in Brainstorm')
end

resultsDir = fullfile(fileDir,'results');
subName    = ''; % Input Brainstorm subject name here
atlasName  = 'SN_atlas';

%% ====================== CARICA DATI =====================================
subResultsDir = fullfile(resultsDir, subName);
assert(exist(subResultsDir, 'dir') == 7, 'Subject folder not found: %s', subResultsDir);

rankTable   = readtable(fullfile(subResultsDir, 'ROI_ranked.xlsx'));
boolTable   = readtable(fullfile(subResultsDir, 'ROI_boolean_features.xlsx'));
detailTable = readtable(fullfile(subResultsDir, 'ROI_detailed_features.xlsx'));
% roiRankTable = readtable(fullfile(subResultsDir, 'ROI_ranked.csv'));
% Load allResults from .mat 
matFile = fullfile(subResultsDir, 'ROI_features.mat');
assert(exist(matFile, 'file') == 2, 'ROI_features.mat not found: %s', matFile);
loaded = load(matFile, 'allResults');
allRes = loaded.allResults;

% ROI names from rankTable (sorted by TotalScore)
roiNames = rankTable.ROI;
numRois  = numel(roiNames);

if ~ismember('Row', boolTable.Properties.VariableNames)
    % Prima colonna è RowNames
    boolTable.Properties.RowNames = boolTable{:,1};
    boolTable(:,1) = [];
    detailTable.Properties.RowNames = detailTable{:,1};
    detailTable(:,1) = [];
else
    boolTable.Properties.RowNames = boolTable.Row;
    boolTable.Row = [];
    detailTable.Properties.RowNames = detailTable.Row;
    detailTable.Row = [];
end

%% ======================== CORTICAL SURFACE ==============================
subStruct    = bst_get('Subject', subName);
mriStruct    = in_bst_data(subStruct.Anatomy(subStruct.iAnatomy).FileName);
surfFileName = subStruct.Surface(subStruct.iCortex).FileName;
surfStruct   = in_bst_data(surfFileName);
atlasIdx     = find(strcmpi({surfStruct.Atlas.Name}, atlasName));
assert(~isempty(atlasIdx), 'Atlas %s not found for %s', atlasName, subName);

Vertices = surfStruct.Vertices;
Faces    = surfStruct.Faces;

% DLPFC MNI → SCS
dlpfc_mni = [-38, 30, 30] / 1000;
dlpfc_scs = cs_convert(mriStruct, 'mni', 'scs', dlpfc_mni);

%% ========================== MAIN FIGURE =================================
roiList = roiNames(1:min(10, numRois));
nPlot   = numel(roiList);
colors  = lines(nPlot);

nCols = ceil(sqrt(nPlot));
nRows = ceil(nPlot / nCols);
 
hFig  = figure('WindowState', 'maximized', ...
               'Name', sprintf('Top ROI — %s', subName), 'NumberTitle', 'off');
hAxes = gobjects(nPlot, 1);
 
for i = 1:nPlot
    hAxes(i) = subplot(nRows, nCols, i);
 
    patch('Vertices', Vertices, 'Faces', Faces, ...
          'FaceColor', [0.28 0.28 0.28], 'EdgeColor', 'none');
    camlight(0, 40, 'infinite');   camlight(180, 40, 'infinite');
    camlight(0, -90, 'infinite');  camlight(90, 0, 'infinite');
    camlight(-90, 0, 'infinite');
    light('Style', 'infinite', 'Color', [0.4 0.4 0.4], 'Parent', hAxes(i));
    camlight('headlight');
    axis equal off;  lighting gouraud;  hold on; material dull;
 
    % ROIs
    roiIdx_atlas = find(strcmpi({surfStruct.Atlas(atlasIdx).Scouts.Label}, roiList{i}));
    if isempty(roiIdx_atlas)
        title(hAxes(i), sprintf('[?] %s', roiList{i}), 'Interpreter', 'none', 'FontSize', 7);
        continue
    end
    ROIverts     = surfStruct.Atlas(atlasIdx).Scouts(roiIdx_atlas).Vertices';
    roiVertColor = NaN(size(Vertices, 1), 3);
    roiVertColor(ROIverts, :) = repmat([0.9 0.15 0.1], numel(ROIverts), 1);
    patch('Vertices', Vertices, 'Faces', Faces, ...
          'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceVertexCData', roiVertColor);
 
    if ~isempty(surfStruct.Atlas(atlasIdx).Scouts(roiIdx_atlas).Seed)
        roiCentre = Vertices(surfStruct.Atlas(atlasIdx).Scouts(roiIdx_atlas).Seed, :);
    else
        roiCentre = mean(Vertices(ROIverts, :));
    end
    dlpfc_dist = norm(roiCentre - dlpfc_scs) * 1000;
 
    % DLPFC center (blu) on surface
    d      = dlpfc_scs / norm(dlpfc_scs);
    pts    = dlpfc_scs + linspace(0, 0.1, 20)' * d;
    Dists  = sqrt(sum((reshape(Vertices, [], 1, 3) - reshape(pts, 1, [], 3)).^2, 3));
    [minD, vIdx]   = min(Dists, [], 1);
    [~, bestSamp]  = max(minD);
    pSurf = Vertices(vIdx(bestSamp), :);
    scatter3(pSurf(1), pSurf(2), pSurf(3), 80, [0.1 0.3 0.9], 'filled');
 
    % Get camera on ROIs
    campos(roiCentre * 1.5);
    camlight('headlight');
 
    % Number of active features
    roiSafe  = matlab.lang.makeValidName(roiList{i});
    nFeatROI = 0;
    if ismember(roiSafe, boolTable.Properties.RowNames)
        nFeatROI = sum(table2array(boolTable(roiSafe, :)) == 1, 'omitnan');
    end
 
    title(hAxes(i), ...
          {strrep(roiList{i}, '_', ' '), ...
           sprintf('%.1f mm  |  %d feat', dlpfc_dist, nFeatROI)}, ...
          'Interpreter', 'none', 'FontSize', 7);
    % rotate3d on
end
 
sgtitle(sprintf('Top %d ROI — %s', nPlot, subName), ...
        'FontWeight', 'bold', 'Interpreter', 'none');
 
% Callback click
hFig.WindowButtonDownFcn = @(src, ~) exploreFeatures(src, hAxes, roiList, ...
    boolTable, detailTable, subName, allRes);



%% ========================================================================
%                        INTERACTIVE FUNCTIONS
%% ========================================================================

function exploreFeatures(hFig, hAxes, roiList, boolTable, detailTable, subName, allRes)
    % Right click --> activate rotation on figure
    if strcmp(hFig.SelectionType, 'alt')
        tmp = rotate3d;
        if strcmpi(tmp.Enable, 'on')
            rotate3d off
        else
            rotate3d on
        end
        
        return
    end
    clickedAx = gca;
    axIdx     = find(hAxes == clickedAx);
    if isempty(axIdx), return; end
    roiName = roiList{axIdx};
    roiSafe = matlab.lang.makeValidName(roiName);

    if ~ismember(roiSafe, boolTable.Properties.RowNames)
        msgbox(sprintf('ROI %s not found in table', roiName), 'Info');
        return
    end

    featNames  = boolTable.Properties.VariableNames;
    boolRow    = table2array(boolTable(roiSafe, :));
    activeIdx  = find(boolRow == 1);
    activeFeats= featNames(activeIdx);

    % Interactive figure
    prefixName = 'ROI Explorer';

    % Close previously opened figures with same name
    prevFigs = findobj( 'Type', 'Figure');
    prevFigs(~startsWith({prevFigs.Name}, prefixName)) = [];
    if ~isempty(prevFigs)
        for f = 1:numel(prevFigs)
            close(prevFigs(f))
        end
    end
    
    figName = sprintf('%s — %s', prefixName, roiName);
    hInteract = figure('Name', figName, ...
                       'Units', 'normalized', 'Position', [0.03 0.08 0.94 0.80], ...
                       'NumberTitle', 'off');
    set(gca, 'Visible', 'off');

    hPanLeft = uipanel('Parent', hInteract, ...
                       'Title', sprintf('Active features (%d)', numel(activeFeats)), ...
                       'Units', 'normalized', 'Position', [0.0 0.0 0.17 1.0], ...
                       'FontSize', 9, 'FontWeight', 'bold');

    hPanRight = uipanel('Parent', hInteract, ...
                        'Title', 'Feature details', ...
                        'Units', 'normalized', 'Position', [0.175 0.0 0.50 1.0], ...
                        'FontSize', 9, 'FontWeight', 'bold');

    % ERP + ERSD always visible panel
    hPanSig = uipanel('Parent', hInteract, ...
                      'Title', sprintf('ROI — %s', roiName), ...
                      'Units', 'normalized', 'Position', [0.68 0.0 0.32 1.0], ...
                      'FontSize', 9, 'FontWeight', 'bold');

    nFeat = numel(activeFeats);
    btnH  = min(0.065, 0.93 / nFeat);
    hBtns = gobjects(nFeat, 1);

    for f = 1:nFeat
        yPos     = 0.97 - f * (btnH + 0.005);
        featDisp = strrep(activeFeats{f}, '_', ' ');
        % Coloured labels
        if startsWith(activeFeats{f}, 'SampEn') || startsWith(activeFeats{f}, 'HFD') || startsWith(activeFeats{f}, 'PE')
            bgcol = [0.85 0.92 0.98];
        elseif startsWith(activeFeats{f}, 'FOOOF')
            bgcol = [0.92 0.88 0.98];
        elseif startsWith(activeFeats{f}, 'ERSD')
            bgcol = [0.88 0.96 0.88];
        elseif startsWith(activeFeats{f}, 'ITPC')
            bgcol = [0.98 0.93 0.82];
        else
            bgcol = [0.93 0.93 0.93];
        end
        hBtns(f) = uicontrol('Parent', hPanLeft, 'Style', 'radiobutton', ...
                              'String', featDisp, 'Units', 'normalized', ...
                              'Position', [0.03 yPos 0.94 btnH], ...
                              'Value', f == 1, 'FontSize', 7.5, ...
                              'BackgroundColor', bgcol, ...
                              'Tag', activeFeats{f});
    end

    for f = 1:nFeat
        set(hBtns(f), 'Callback', @(s, ~) onFeatureSelect(s, hBtns, ...
            roiName, hPanRight, subName, allRes, detailTable, featNames, ...
            table2array(boolTable)));
    end

    % Initial features plot
    if nFeat > 0
        plotFeature(activeFeats{1}, roiName, hPanRight, subName, allRes, ...
                    detailTable, featNames, table2array(boolTable));
    else
        ax = axes('Parent', hPanRight);
        text(ax, 0.5, 0.5, sprintf('No active feature for:\n%s', roiName), ...
             'HorizontalAlignment','center','Units','normalized', ...
             'FontSize',11,'Interpreter','none');
        axis(ax,'off');
    end

    drawnow limitrate
    plot_roi_signals(roiName, subName, hPanSig);
end


function onFeatureSelect(src, hBtns, roiName, hPanRight, subName, allRes, detailTable, featNames, boolMat)
    set(hBtns, 'Value', 0);
    set(src,   'Value', 1);
    plotFeature(src.Tag, roiName, hPanRight, subName, allRes, detailTable, featNames, boolMat);
end

function plotFeature(featName, roiName, hPan, subName, allRes, detailTable, featNames, boolMat)
    delete(hPan.Children);

    if startsWith(featName, 'SampEn') || startsWith(featName, 'HFD') || startsWith(featName, 'PE')
        plot_complexity(featName, roiName, hPan, allRes, subName);

    elseif startsWith(featName, 'FOOOF_Exp') || startsWith(featName, 'FOOOF_Off')
        plot_fooof(featName, roiName, hPan, allRes, subName);

    elseif startsWith(featName, 'ERSD')
        if contains(featName, 'VsRest')
            compFile = 'Comparison_Scouts_ERSD_freqs_Immagine_vs_Rest';
        else
            compFile = 'Comparison_Scouts_ERSD_freqs_Immagine_vs_contralateral';
        end
        plot_ersd(featName, roiName, compFile, subName, hPan);

    elseif startsWith(featName, 'ITPC')
        plot_itpc(featName, roiName, hPan, allRes);

    else
        ax = axes('Parent', hPan);
        text(ax, 0.5, 0.5, sprintf('Unavailable plot for:\n%s', strrep(featName,'_',' ')), ...
             'HorizontalAlignment','center','FontSize',11,'Units','normalized','Interpreter','none');
        axis(ax,'off');
    end
end


%% ========================================================================
%               PLOT COMPLEXITY (SampEn, HFD, PE)
%% ========================================================================

function plot_complexity(featName, roiName, hPan, allRes, subName)
 
    sName = matlab.lang.makeValidName(roiName);
    if startsWith(featName, 'SampEn'), field = 'sampen';
    elseif startsWith(featName, 'HFD'), field = 'hfd';
    else,                               field = 'pe';
    end
 
    ax1 = axes('Parent', hPan, 'Position', [0.08 0.35 0.88 0.58]);
    ax2 = axes('Parent', hPan, 'Position', [0.08 0.04 0.88 0.24]);
 
    try
        % Read parameters from allRes if available
        fs   = 256;
        tWin = [0.07, 0.2];
        if isfield(allRes, 'params')
            if isfield(allRes.params, 'fs'),   fs   = allRes.params.fs;   end
            if isfield(allRes.params, 'tWin'), tWin = allRes.params.tWin; end
        end
        dTask  = compute_complexity_from_bst(subName, 'Immagine', 'SN_scout', roiName, field, fs, tWin);
        dRest  = compute_complexity_from_bst(subName, 'Rest',     'SN_rest_scout',            roiName, field, fs, tWin);
        dContra= compute_complexity_from_bst(subName, 'Immagine', 'Contralateral_SN_scout',  roiName, field, fs, tWin);
    catch ME
        plot_wilcoxon_summary(ax1, allRes, sName, field, featName);
        axis(ax2, 'off');
        title(ax1, sprintf('Data unavailable: %s', ME.message), ...
              'Interpreter','none','FontSize',8,'Color',[0.7 0 0]);
        return
    end
 
    % --- Boxplot + jittered dots ---
    cols = [0.15 0.45 0.8; 0.5 0.5 0.5; 0.85 0.35 0.1];
    lbls = {'Task ipsi','Rest','Task contra'};
    nT = numel(dTask); nR = numel(dRest); nC = numel(dContra);
    allD = [dTask; dRest; dContra];
    grp  = [ones(nT,1); 2*ones(nR,1); 3*ones(nC,1)];
 
    boxplot(ax1, allD, grp, 'Labels', lbls, ...
            'Colors', cols, 'Widths', 0.45, 'Symbol', '');
    set(findobj(ax1,'Tag','Box'),     'LineWidth', 1.8);
    set(findobj(ax1,'Tag','Median'),  'LineWidth', 2.0, 'Color', 'k');
    set(findobj(ax1,'Tag','Whisker'), 'LineWidth', 1.2);
    hold(ax1,'on');
 
    jitterW  = 0.15;
    dataSets = {dTask, dRest, dContra};
    for g = 1:3
        d    = dataSets{g};
        xPos = g + (rand(numel(d),1)-0.5) * 2 * jitterW;
        scatter(ax1, xPos, d, 22, cols(g,:), 'filled', 'MarkerFaceAlpha', 0.55);
    end
 
    grid(ax1,'on'); box(ax1,'off');
    annotate_wilcoxon(ax1, allRes, sName, field, max([nT nR nC]));
    hold(ax1,'off');
 
    ylabel(ax1, strrep(field,'_',' '), 'FontSize', 9);
    title(ax1, sprintf('%s — %s', strrep(featName,'_',' '), strrep(roiName,'_',' ')), ...
          'Interpreter','none','FontSize',10,'FontWeight','bold');
 
    % --- Stats summary below ---
    plot_wilcoxon_summary(ax2, allRes, sName, field, '');
    axis(ax2,'off');
end

function vals = compute_complexity_from_bst(subName, cond, scoutTag, roiName, field, fs, tWin)
% Recompute complexity features from scouts timeseries
    files = bst_process('CallProcess', 'process_select_files_matrix', [], [], ...
        'subjectname',   subName, ...
        'condition',     cond, ...
        'tag',           scoutTag, ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0, ...
        'outprocesstab', 'no');
    files(~cellfun(@(c) strncmp(c, [scoutTag ' '], length(scoutTag)+1), {files.Comment})) = [];
    if isempty(files)
        error('Nessun file scout trovato: %s / %s / %s', subName, cond, scoutTag)
    end
 
    % ROI index
    tmp    = in_bst_data(files(1).FileName);
    rNames = tmp.Description;
    rNames = cellfun(@(s) strtrim(regexp(s,' @','split','once')), rNames, 'UniformOutput', false);
    rNames = cellfun(@(c) c{1}, rNames, 'UniformOutput', false);
    rIdx   = find(strcmpi(rNames, roiName), 1);
    if isempty(rIdx)
        rIdx = find(startsWith(rNames, roiName), 1);
    end
    if isempty(rIdx)
        error('ROI "%s" non trovata negli scout', roiName)
    end
 
    nTrials  = numel(files);
    nSamples = size(tmp.Value, 2);
 
    % Build tVec
    if isfield(tmp, 'Time') && ~isempty(tmp.Time)
        tVec = tmp.Time;
    else
        tVec = (0:nSamples-1) / fs;
    end
 
    % Assemble matrix [nTrials x nSamples]
    trialMat = nan(nTrials, nSamples);
    for t = 1:nTrials
        d = in_bst_data(files(t).FileName);
        trialMat(t,:) = d.Value(rIdx, :);
    end
 
    % Compute feature with PPA_complexity_features
    if isempty(tWin)
        cx = PPA_complexity_features(trialMat, fs);
    else
        cx = PPA_complexity_features(trialMat, fs, 'tWin', tWin, 'tVec', tVec);
    end
    vals = cx.(field);
end


function annotate_wilcoxon(ax, allRes, sName, field, n)
    % Significance bar in the plot
    yLims  = ylim(ax);
    yRange = yLims(2) - yLims(1);
    yTop   = yLims(2) + 0.05 * yRange;
    hold(ax, 'on');

    pairs = {1, 2, 'compResults_vsRest';   % Task vs Rest
             1, 3, 'compResults_vsContra'}; % Task vs Contra

    offsets = [0.10, 0.20];
    for p = 1:size(pairs,1)
        x1     = pairs{p,1};  x2 = pairs{p,2};
        resKey = pairs{p,3};
        if ~isfield(allRes, resKey) || ~isfield(allRes.(resKey), sName), continue; end
        entry = allRes.(resKey).(sName).(field);
        if ~entry.sig, continue; end

        yLine = yTop + offsets(p) * yRange;
        plot(ax, [x1 x1 x2 x2], [yLine-0.01*yRange, yLine, yLine, yLine-0.01*yRange], ...
             'k-', 'LineWidth', 1);
        pStr = format_pval(entry.p);
        text(ax, (x1+x2)/2, yLine + 0.01*yRange, pStr, ...
             'HorizontalAlignment','center','FontSize',8,'FontWeight','bold');
    end
    hold(ax, 'off');
    ylim(ax, [yLims(1), yTop + 0.35*yRange]);
end


function plot_wilcoxon_summary(ax, allRes, sName, field, featName)
    lines_txt = {strrep(featName,'_',' '), '', 'Statistiche Wilcoxon:'};
    pairs = {'compResults_vsRest','Task vs Rest'; 'compResults_vsContra','Task vs Contra'};
    for p = 1:size(pairs,1)
        key   = pairs{p,1};  lbl = pairs{p,2};
        if ~isfield(allRes, key) || ~isfield(allRes.(key), sName), continue; end
        e = allRes.(key).(sName).(field);
        if isnan(e.p)
            lines_txt{end+1} = sprintf('%s: n.d.', lbl);
        else
            lines_txt{end+1} = sprintf('%s: p=%.3f  medΔ=%.3f  sig=%d', ...
                lbl, e.p, e.medianDiff, double(e.sig));
        end
    end
    text(ax, 0.05, 0.5, lines_txt, 'Units','normalized','FontSize',10, ...
         'VerticalAlignment','middle','Interpreter','none');
    axis(ax, 'off');
end


%% ========================================================================
%               PLOT FOOOF (Exponent, Offset)
%% ========================================================================

function plot_fooof(featName, roiName, hPan, allRes, subName)
 
    sName = matlab.lang.makeValidName(roiName);
    if contains(featName, 'Exp'), field = 'fooof_exp'; srcField = 'exponent';
    else,                          field = 'fooof_off'; srcField = 'offset';
    end
 
    switch field
        case 'fooof_exp'
            scoutName = 'Scouts_FOOOF_Exponent_Immagine';
            contraName= 'Contralateral_Scouts_FOOOF_Exponent_Immagine';
            restName  = 'Scouts_FOOOF_Exponent_Rest';
        case 'fooof_off'
            scoutName = 'Scouts_FOOOF_Offset_Immagine';
            contraName= 'Contralateral_Scouts_FOOOF_Offset_Immagine';
            restName  = 'Scouts_FOOOF_Offset_Rest';
    end
 
    ax1 = axes('Parent', hPan, 'Position', [0.08 0.35 0.88 0.58]);
    ax2 = axes('Parent', hPan, 'Position', [0.08 0.04 0.88 0.24]);
 
    % --- Load trials (no min) ---
    try
        scoutStruct = load_fooof_file(subName, scoutName,  'Immagine');
        contraStruct= load_fooof_file(subName, contraName, 'Immagine');
        restStruct  = load_fooof_file(subName, restName,   'Rest');
        dTask  = load_fooof_features(scoutStruct,  roiName);
        dRest  = load_fooof_features(restStruct,   roiName);
        dContra= load_fooof_features(contraStruct, roiName);
    catch ME
        plot_wilcoxon_summary(ax1, allRes, sName, field, featName);
        axis(ax2,'off');
        return
    end
 
    % --- Boxplot + jittered dots ---
    cols = [0.15 0.45 0.8; 0.5 0.5 0.5; 0.85 0.35 0.1];
    lbls = {'Task ipsi','Rest','Task contra'};
    nT   = numel(dTask); nR = numel(dRest); nC = numel(dContra);
    allD = [dTask; dRest; dContra];
    grp  = [ones(nT,1); 2*ones(nR,1); 3*ones(nC,1)];
 
    bx = boxplot(ax1, allD, grp, 'Labels', lbls, ...
                 'Colors', cols, 'Widths', 0.45, 'Symbol', '');
    set(findobj(ax1,'Tag','Box'),      'LineWidth', 1.8);
    set(findobj(ax1,'Tag','Median'),   'LineWidth', 2.0, 'Color', 'k');
    set(findobj(ax1,'Tag','Whisker'),  'LineWidth', 1.2);
    set(findobj(ax1,'Tag','Outliers'), 'Marker', 'none');
    hold(ax1,'on');
 
    % Jittered dots per condition
    jitterW = 0.15;
    dataSets = {dTask, dRest, dContra};
    for g = 1:3
        d    = dataSets{g};
        xPos = g + (rand(numel(d),1)-0.5) * 2 * jitterW;
        scatter(ax1, xPos, d, 22, cols(g,:), 'filled', 'MarkerFaceAlpha', 0.55);
    end
 
    grid(ax1,'on'); box(ax1,'off');
    annotate_wilcoxon(ax1, allRes, sName, field, max([nT nR nC]));
    hold(ax1,'off');
 
    title(ax1, sprintf('%s — %s', strrep(featName,'_',' '), strrep(roiName,'_',' ')), ...
          'Interpreter','none','FontSize',9,'FontWeight','bold');
    ylabel(ax1, srcField);
 
    % --- Stats summary below ---
    plot_wilcoxon_summary(ax2, allRes, sName, field, '');
    axis(ax2,'off');
 
    % --- Mean aperiodic spectre ---
    try
        conds   = {'Immagine','Rest','Immagine'};
        prefixs = {'Scouts','Scouts','Contralateral_Scouts'};
        lbls    = {'Task ipsi','Rest','Task contra'};
        cols    = [0.15 0.45 0.8; 0.5 0.5 0.5; 0.85 0.35 0.1];
        hold(ax2,'on');
        for c = 1:3
            fTag = sprintf('%s_FOOOF_%s', prefixs{c}, conds{c});
            fFiles = bst_process('CallProcess','process_select_files_timefreq',[],[], ...
                'subjectname',subName,'condition',conds{c},'tag',fTag, ...
                'includebad',0,'includeintra',0,'includecommon',0,'outprocesstab','no');
            fFiles(~cellfun(@(cm) strncmp(cm,[fTag ' '],length(fTag)+1),{fFiles.Comment}))=[];
            if isempty(fFiles), continue; end
            % Media tra trial: carica tutti e media
            allTF = [];
            for t = 1:numel(fFiles)
                d = in_bst_data(fFiles(t).FileName);
                rIdx = find(startsWith(d.RowNames, roiName), 1);
                if isempty(rIdx), continue; end
                allTF(end+1,:) = squeeze(d.TF(rIdx,:,:)); %#ok<AGROW>
            end
            if isempty(allTF), continue; end
            meanTF = mean(allTF, 1, 'omitnan');
            freqV  = d.Freqs;
            plot(ax2, freqV, 10*log10(meanTF), 'Color', cols(c,:), ...
                 'LineWidth', 1.8, 'DisplayName', lbls{c});
        end
        hold(ax2,'off');
        xlabel(ax2,'Frequenza (Hz)'); ylabel(ax2,'Potenza aperiodica (dB)');
        legend(ax2,'Location','northeast','FontSize',7);
        grid(ax2,'on'); box(ax2,'off');
        title(ax2,'Spettro aperiodico medio','FontSize',8);
    catch ME
        cla(ax2);
        text(ax2, 0.5, 0.5, sprintf('Spettro non disponibile\n%s', ME.message), ...
             'HorizontalAlignment','center','Units','normalized', ...
             'FontSize',8,'Interpreter','none');
        axis(ax2,'off');
    end
end


%% ========================================================================
%               PLOT ERSD (t-map band × time)
%% ========================================================================

function plot_ersd(featName, roiName, compFileName, subName, hPan)
    compFile = find_presults_bst(subName, compFileName, ...
                   double(~contains(featName,'VsContra')), 'timefreq');
    if isempty(compFile)
        ax = axes('Parent', hPan);
        text(ax, 0.5, 0.5, sprintf('File non trovato:\n%s', compFileName), ...
             'HorizontalAlignment','center','Units','normalized', ...
             'FontSize',9,'Interpreter','none');
        axis(ax,'off'); return
    end

    data = in_bst_data(compFile{1});
    ax   = axes('Parent', hPan, 'Position', [0.09 0.12 0.85 0.80]);

    roiIdx = find(startsWith(data.RowNames, roiName), 1);
    if isempty(roiIdx)
        roiIdx = find(strcmpi(data.RowNames, roiName), 1);
    end
    if isempty(roiIdx)
        text(ax,0.5,0.5,sprintf('ROI "%s" non trovata',roiName), ...
             'HorizontalAlignment','center','Units','normalized','Interpreter','none');
        axis(ax,'off'); return
    end

    tmap = squeeze(data.tmap(roiIdx, :, :));   % [nTime × nBand]
    pmap = squeeze(data.pmap(roiIdx, :, :));
    tmap(pmap > 0.05) = NaN;

    imagesc(ax, data.Time, 1:size(tmap,2), tmap');
    set(ax,'YDir','normal');

    if iscell(data.Freqs)
        bandLabels = data.Freqs(:,1);
    else
        bandLabels = arrayfun(@(f) sprintf('%.0f Hz',f), data.Freqs, 'UniformOutput',false);
    end
    set(ax,'YTick',1:numel(bandLabels),'YTickLabel',bandLabels);

    tTicks = linspace(data.Time(1), data.Time(end), 7);
    set(ax,'XTick',tTicks,'XTickLabel',arrayfun(@(t) sprintf('%.2f',t),tTicks,'UniformOutput',false));

    xlabel(ax,'Tempo (s)'); ylabel(ax,'Banda');
    title(ax, sprintf('%s — %s', strrep(featName,'_',' '), strrep(roiName,'_',' ')), ...
          'Interpreter','none','FontSize',9,'FontWeight','bold');

    cm = colormap(ax,'turbo');
    cm = [0.93 0.93 0.93; cm];   % grey for NaN
    colormap(ax, cm); colorbar(ax);
    climV = max(abs(tmap(:)),[],'omitnan');
    if ~isnan(climV) && climV > 0, clim(ax,[-climV climV]); end
end


%% ========================================================================
%               PLOT ITPC (TF time×freq with bands)
%% ========================================================================

function plot_itpc(featName, roiName, hPan, allRes)
    if ~isfield(allRes,'ITPC') || isempty(allRes.ITPC)
        ax = axes('Parent', hPan);
        text(ax,0.5,0.5,'ITPC unavailable','HorizontalAlignment','center', ...
             'Units','normalized','FontSize',11);
        axis(ax,'off'); return
    end

    itpcData = allRes.ITPC;
    roiIdx   = find(strcmpi(itpcData.roiNames, roiName), 1);
    if isempty(roiIdx)
        roiIdx = find(startsWith(itpcData.roiNames, roiName), 1);
    end

    % Layout: TF map at left, contrasts at right
    ax1 = axes('Parent', hPan, 'Position', [0.06 0.12 0.55 0.78]);
    ax2 = axes('Parent', hPan, 'Position', [0.66 0.12 0.30 0.78]);

    % --- TF map ITPC task ---
    if ~isempty(roiIdx)
        tfData = squeeze(itpcData.ITPC(roiIdx, :, :));   % [nTime × nFreq]  (BST: time×freq)
        % If swapped dimensions, correct
        if size(tfData,1) ~= numel(itpcData.times)
            tfData = tfData';
        end
        imagesc(ax1, itpcData.times, itpcData.freqs, tfData');
        set(ax1,'YDir','normal');
        colormap(ax1,'hot'); cb = colorbar(ax1);
        cb.Label.String = 'ITPC'; clim(ax1,[0 1]);
        xlabel(ax1,'Tempo (s)'); ylabel(ax1,'Frequenza (Hz)');
        title(ax1, sprintf('ITPC Task — %s', strrep(roiName,'_',' ')), ...
              'Interpreter','none','FontSize',9,'FontWeight','bold');

        % Band rectangles
        bandDefs = struct('name',{'Theta','Alpha','Beta'}, ...
                          'lims',{[4,8],[8,13],[13,30]}, ...
                          'color',{[0 0.6 1],[0 0.8 0],[1 0.4 0]});
        hold(ax1,'on');
        for b = 1:numel(bandDefs)
            bN = bandDefs(b).name;
            if ~isfield(itpcData,bN), continue; end
            fl = bandDefs(b).lims(1); fh = bandDefs(b).lims(2);
            clr= bandDefs(b).color;
            isSigB = itpcData.(bN).RelRois(roiIdx);
            val    = itpcData.(bN).Value(roiIdx);
            thr    = itpcData.(bN).Thr;
            alph   = 0.12 + 0.15*double(isSigB);
            tL = itpcData.times(1); tR = itpcData.times(end);
            fill(ax1,[tL tR tR tL],[fl fl fh fh],clr,'FaceAlpha',alph, ...
                 'EdgeColor',clr,'LineWidth',1+double(isSigB));
            sigMk = ''; if isSigB, sigMk = ' ✓'; end
            text(ax1, tL+0.02*(tR-tL), fl+0.3*(fh-fl), ...
                 sprintf('%s %.2f/%.2f%s',bN,val,thr,sigMk), ...
                 'Color',clr,'FontSize',7,'FontWeight','bold', ...
                 'BackgroundColor',[1 1 1 0.5],'Interpreter','none');
        end
        hold(ax1,'off');
    else
        text(ax1,0.5,0.5,sprintf('ROI "%s" not found in ITPC',roiName), ...
             'HorizontalAlignment','center','Units','normalized','Interpreter','none');
        axis(ax1,'off');
    end

    % --- Wilcoxon contrasts (bar chart p-values per band) ---
    sName = matlab.lang.makeValidName(roiName);
    if isfield(allRes,'itpcComp') && isfield(allRes.itpcComp, sName)
        comp = allRes.itpcComp.(sName);
        bNames = {'Theta','Alpha','Beta'};
        compTypes = {'VsRest','VsContra'};
        ctLabels  = {'vs Rest','vs Contra'};
        pVals     = nan(numel(bNames), numel(compTypes));
        isSigMat  = false(numel(bNames), numel(compTypes));
        for b = 1:numel(bNames)
            if ~isfield(comp, bNames{b}), continue; end
            bd = comp.(bNames{b});
            if ~bd.gatePass, continue; end
            for ct = 1:numel(compTypes)
                if ~isfield(bd, compTypes{ct}), continue; end
                pVals(b,ct)    = bd.(compTypes{ct}).p;
                isSigMat(b,ct) = bd.(compTypes{ct}).sig;
            end
        end
        % Bar chart -log10(p)
        logP = -log10(max(pVals, 1e-4));
        logP(isnan(pVals)) = 0;
        bh = bar(ax2, logP, 'grouped');
        bh(1).FaceColor = [0.2 0.5 0.85];
        bh(2).FaceColor = [0.85 0.45 0.1];
        set(ax2,'XTickLabel',bNames,'FontSize',8);
        yline(ax2, -log10(0.05), '--k', 'p=0.05','LabelHorizontalAlignment','left','FontSize',7);
        legend(ax2, ctLabels,'FontSize',7,'Location','northwest');
        ylabel(ax2,'-log_{10}(p)'); title(ax2,'Contrasti ITPC','FontSize',8);
        grid(ax2,'on'); box(ax2,'off');

        % Asterisks
        hold(ax2,'on');
        nGroups = numel(bNames); nBars = numel(compTypes);
        groupW  = 0.8; barW = groupW/nBars;
        for b = 1:nGroups
            for ct = 1:nBars
                if isSigMat(b,ct)
                    xPos = b + (ct - (nBars+1)/2) * barW;
                    text(ax2, xPos, logP(b,ct)+0.05, '*', ...
                         'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
                end
            end
        end
        hold(ax2,'off');
    else
        text(ax2,0.5,0.5,'Unavailable contrasts','HorizontalAlignment','center', ...
             'Units','normalized','FontSize',9,'Interpreter','none');
        axis(ax2,'off');
    end
end


%% ========================================================================
%                     MEAN ERP + ERSD PANEL
%% ========================================================================

function plot_roi_signals(roiName, subName, hPan)

    erpTag  = 'SN_Scout';
    ersdTag = 'Scouts_ERSD_Immagine';
    cond    = 'Immagine';

    delete(hPan.Children);
    axERP = axes('Parent', hPan, 'Position', [0.15 0.60 0.72 0.33]);
    axTF  = axes('Parent', hPan, 'Position', [0.15 0.09 0.72 0.35]);

    % ---------------- ERP ----------------
    try
        plot_roi_erp(axERP, subName, cond, erpTag, roiName);
    catch ME
        cla(axERP);
        text(axERP, 0.5, 0.5, sprintf('ERP non disponibile\n%s', ME.message), ...
             'HorizontalAlignment','center','Units','normalized', ...
             'FontSize',8,'Interpreter','none');
        axis(axERP,'off');
    end

    % ---------------- ERSD ----------------
    try
        plot_roi_ersd(axTF, subName, cond, ersdTag, roiName);
    catch ME
        cla(axTF);
        text(axTF, 0.5, 0.5, sprintf('ERSD non disponibile\n%s', ME.message), ...
             'HorizontalAlignment','center','Units','normalized', ...
             'FontSize',8,'Interpreter','none');
        axis(axTF,'off');
    end
end


function plot_roi_erp(ax, subName, cond, tag, roiName)

    files = select_bst_files_strict(subName, cond, tag, 'matrix');
    if isempty(files)
        error('Nessun file matrix "%s" in %s/%s', tag, subName, cond);
    end

    d0   = in_bst_data(files(1).FileName);
    rIdx = find_roi_row(get_row_names(d0), roiName);
    if isempty(rIdx)
        error('ROI "%s" not present in RowNames di %s', roiName, tag);
    end

    tVec  = d0.Time(:)';
    nSamp = numel(tVec);
    nTr   = numel(files);
    M     = nan(nTr, nSamp);

    for t = 1:nTr
        if t == 1
            d = d0;
        else
            d = in_bst_data(files(t).FileName);
        end
        rI = find_roi_row(get_row_names(d), roiName);
        if isempty(rI) || size(d.Value,2) ~= nSamp, continue; end
        M(t,:) = d.Value(rI, :);
    end
    M(all(isnan(M),2), :) = [];
    if isempty(M), error('No valid trial for ROI'); end

    mu = mean(M, 1, 'omitnan');
    n  = sum(~isnan(M), 1);
    se = std(M, 0, 1, 'omitnan') ./ max(sqrt(n), 1);

    hold(ax,'on');
    fill(ax, [tVec fliplr(tVec)], [mu-se fliplr(mu+se)], [0.15 0.45 0.8], ...
         'FaceAlpha', 0.20, 'EdgeColor','none', 'HandleVisibility','off');
    plot(ax, tVec, mu, 'Color', [0.10 0.35 0.75], 'LineWidth', 1.8);
    yl = ylim(ax);
    if tVec(1) <= 0 && tVec(end) >= 0
        plot(ax, [0 0], yl, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
        ylim(ax, yl);
    end
    yline(ax, 0, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
    hold(ax,'off');

    xlim(ax, [tVec(1) tVec(end)]);
    xlabel(ax, 'Tempo (s)', 'FontSize', 8);
    ylabel(ax, 'Ampiezza (a.u.)', 'FontSize', 8);
    title(ax, sprintf('Mean ERP — %s (n=%d trial)', strrep(roiName,'_',' '), size(M,1)), ...
          'Interpreter','none','FontSize',9,'FontWeight','bold');
    grid(ax,'on'); box(ax,'off'); set(ax,'FontSize',8);
end


function plot_roi_ersd(ax, subName, cond, tag, roiName)

    files = select_bst_files_strict(subName, cond, tag, 'timefreq');
    if isempty(files)
        error('No timefreq file"%s" in %s/%s', tag, subName, cond);
    end

    d0   = in_bst_data(files(1).FileName);
    rIdx = find_roi_row(get_row_names(d0), roiName);
    if isempty(rIdx)
        error('ROI "%s" not present in RowNames di %s', roiName, tag);
    end

    tVec = d0.Time(:)';
    ref  = squeeze(d0.TF(rIdx, :, :));          % [nTime × nFreq]
    if size(ref,1) ~= numel(tVec) && size(ref,2) == numel(tVec)
        ref = ref.';
    end

    acc = zeros(size(ref));
    cnt = 0;
    for t = 1:numel(files)
        if t == 1
            d = d0;
        else
            d = in_bst_data(files(t).FileName);
        end
        rI = find_roi_row(get_row_names(d), roiName);
        if isempty(rI), continue; end
        tf = squeeze(d.TF(rI, :, :));
        if size(tf,1) ~= size(ref,1) && size(tf,2) == size(ref,1)
            tf = tf.';
        end
        if ~isequal(size(tf), size(ref)), continue; end
        acc = acc + real(tf);
        cnt = cnt + 1;
    end
    if cnt == 0, error('No valid trial for ROI'); end
    meanTF = acc / cnt;                          % [nTime × nFreq]

    % Frequencies axes
    if iscell(d0.Freqs)
        bandLabels = d0.Freqs(:,1);
        yVec = 1:numel(bandLabels);
        isBand = true;
    else
        yVec = d0.Freqs(:)';
        isBand = false;
    end

    imagesc(ax, tVec, yVec, meanTF.');
    set(ax, 'YDir', 'normal');
    if isBand
        set(ax, 'YTick', yVec, 'YTickLabel', bandLabels, 'TickLabelInterpreter','none');
        ylabel(ax, 'Band', 'FontSize', 8);
    else
        ylabel(ax, 'Frequency (Hz)', 'FontSize', 8);
    end

    colormap(ax, 'turbo');
    cb = colorbar(ax); cb.FontSize = 7;
    climV = max(abs(meanTF(:)), [], 'omitnan');
    if ~isempty(climV) && ~isnan(climV) && climV > 0
        clim(ax, [-climV climV]);
    end

    if tVec(1) <= 0 && tVec(end) >= 0
        hold(ax,'on');
        plot(ax, [0 0], ylim(ax), 'k--', 'LineWidth', 0.8);
        hold(ax,'off');
    end

    xlabel(ax, 'Tempo (s)', 'FontSize', 8);
    title(ax, sprintf('Mean ERSD — %s (n=%d trial)', strrep(roiName,'_',' '), cnt), ...
          'Interpreter','none','FontSize',9,'FontWeight','bold');
    set(ax,'FontSize',8);
end


function files = select_bst_files_strict(subName, cond, tag, kind)
% Select BST files which Comment starts with tag followed by non
% alphanumerical delimitator (or end of string)
% Escludes Contralateral tag and '<tag>_something'
    switch lower(kind)
        case 'matrix',   procName = 'process_select_files_matrix';
        case 'timefreq', procName = 'process_select_files_timefreq';
        case 'data',     procName = 'process_select_files_data';
        otherwise, error('kind not supported: %s', kind);
    end

    files = bst_process('CallProcess', procName, [], [], ...
        'subjectname',   subName, ...
        'condition',     cond, ...
        'tag',           tag, ...
        'includebad',    0, ...
        'includeintra',  0, ...
        'includecommon', 0, ...
        'outprocesstab', 'no');
    if isempty(files), return; end

    L    = numel(tag);
    keep = false(1, numel(files));
    for k = 1:numel(files)
        c = files(k).Comment;
        if numel(c) < L || ~strncmpi(c, tag, L), continue; end
        if numel(c) == L
            keep(k) = true;
        else
            nxt = c(L+1);
            keep(k) = isempty(regexp(nxt, '[A-Za-z0-9_]', 'once'));
        end
    end
    files = files(keep);
end


function rn = get_row_names(d)
    rn = {};
    if isfield(d, 'RowNames') && ~isempty(d.RowNames)
        rn = d.RowNames;
    elseif isfield(d, 'Description') && ~isempty(d.Description)
        rn = d.Description;
    end
    if ischar(rn), rn = {rn}; end
    rn = rn(:);
end


function idx = find_roi_row(rowNames, roiName)
    idx = [];
    if isempty(rowNames), return; end

    base = regexprep(rowNames, '\s*@.*$', '');
    base = strtrim(base);

    idx = find(strcmpi(base, roiName), 1);
    if ~isempty(idx), return; end
    
    % Fallback: prefix + non alphanumerical delimitator
    pat = ['^' regexptranslate('escape', roiName) '($|[^A-Za-z0-9_])'];
    hit = ~cellfun(@isempty, regexpi(base, pat, 'once'));
    idx = find(hit, 1);
end


%% ========================================================================
%                               UTILITY
%% ========================================================================

function str = format_pval(p)
    if p < 0.001,      str = 'p<0.001';
    elseif p < 0.01,   str = sprintf('p=%.3f', p);
    else,              str = sprintf('p=%.3f', p);
    end
    if p < 0.05, str = ['*' str]; end
end

function fooofFiles = load_fooof_file(subName, fname, condition)
    fooofFiles = bst_process('CallProcess','process_select_files_timefreq',[],[], ...
                'subjectname',subName,...
                'condition',condition,...
                'tag',fname, ...
                'includebad',0,...
                'includeintra',0,...
                'includecommon',0,...
                'outprocesstab','no');
    fooofFiles(~startsWith({fooofFiles.Comment}, fname)) = [];
    if isempty(fooofFiles)
        error('FOOOF files not found')
    end
end

function values = load_fooof_features(fooofStruct, roiName)
    values = NaN(numel(fooofStruct),1);
    for foo = 1:numel(fooofStruct)
        tmp = in_bst_data(fooofStruct(foo).FileName);
        rIdx = find(startsWith(tmp.RowNames, roiName));
        values(foo) = mean(tmp.TF(rIdx,1,:), 'all', 'omitnan');
    end
end