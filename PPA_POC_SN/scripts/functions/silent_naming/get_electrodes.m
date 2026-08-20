% SCRIPT PER ASSEGNARE GLI ELETTRODI BRAINSTORM LIKE PER IL SILENT NAMING
% (DA RENDERE FUNZIONE)

clear
close all
clc

addpath /gea/home3/dati/Brainstorm_Revopoint/brainstorm_250902_src/brainstorm3
addpath(genpath('/gea/home3/dati/PPA_POC/matscripts'))

if ~brainstorm('status')
    brainstorm
end

%% Acquisici la superficie 
revoStruct = in_bst_data('Giulia/tess_scalp_textured_1006_01_tex.mat');
Vertices = revoStruct.Vertices;
Faces = revoStruct.Faces;

% Plot 
figure
patch('Vertices', revoStruct.Vertices, 'Faces', revoStruct.Faces, 'FaceColor',...
    'interp', 'EdgeColor', 'none', 'FaceVertexCData', revoStruct.Color); axis equal;
hold on
% Acquisici tutti i vertici azzurri
hsv = rgb2hsv(revoStruct.Color);
H = hsv(:,1);
S = hsv(:,2);
V = hsv(:,3);

% Azzurro/ciano: H in [0.5, 0.65], saturazione e valore sufficienti
cyanIdxs = find(H >= 0.50 & H <= 0.65 & S > 0.3 & V > 0.3);
scatter3(revoStruct.Vertices(cyanIdxs,1), revoStruct.Vertices(cyanIdxs,2), revoStruct.Vertices(cyanIdxs,3), ...
    20, 'y', 'filled')
cyanCenter = mean(Vertices(cyanIdxs,:),1);
scatter3(cyanCenter(1), cyanCenter(2), cyanCenter(3), ...
    100, 'b', 'filled')


%% Trova la maschera degli elettrodi verdi
greenIdxs = find(H >= 0.22 & H <= 0.45 & S > 0.2 & V > 0.2);
figure
patch('Vertices', revoStruct.Vertices, 'Faces', revoStruct.Faces, 'FaceColor',...
    'interp', 'EdgeColor', 'none', 'FaceVertexCData', revoStruct.Color); axis equal;
hold on
scatter3(revoStruct.Vertices(greenIdxs,1), revoStruct.Vertices(greenIdxs,2), revoStruct.Vertices(greenIdxs,3), ...
    20, 'y', 'filled')
greenMask = false(size(Vertices));
greenMask(greenIdxs) = true;
greenClusters = dbscan(Vertices(greenIdxs,:), 0.005, 20);
uniqueGreenLbls = unique(greenClusters);
greenCenters = [];
for g = 1:numel(uniqueGreenLbls)
    if uniqueGreenLbls(g) == -1
        continue
    end
    clustIdxs = find(greenClusters == uniqueGreenLbls(g));
    clustCenter = mean(Vertices(greenIdxs(clustIdxs),:), 1);  % <-- greenIdxs(clustIdxs)
    greenCenters = [greenCenters; clustCenter];
end
% Merge dei centri troppo vicini
greenDists = pdist2(greenCenters,greenCenters);
greenDists = triu(greenDists);
th = 0.015;
[i,j] = find(greenDists < th & greenDists > 0);
pairs = [i j];
merged = false(size(greenCenters,1),1);
for k = 1:length(i)
    newMean = mean([greenCenters(i(k),:); greenCenters(j(k),:)],1);
    greenCenters = [greenCenters; newMean];
end
% Rimuovi i centri processati
greenCenters([unique(union(i,j))],:) = [];
for g = 1:length(greenCenters)
    scatter3(greenCenters(g,1),greenCenters(g,2),greenCenters(g,3), 100, 'r','filled')
end

%% Passa alla proiezione 2D
capImg2dSize = 900;
% capRangeinIm = 0.5;
[x, y] = bst_project_2d(Vertices(:,1), Vertices(:,2), Vertices(:,3), 'stereo');
[greenX, greenY] = bst_project_2d(greenCenters(:,1), greenCenters(:,2), greenCenters(:,3), 'stereo');
grayness = revoStruct.Color*[1;1;1]/sqrt(3);
lx=linspace(min(x), max(x), capImg2dSize);
ly=linspace(min(y), max(y), capImg2dSize);
[X,Y]=meshgrid(lx,ly);
capImg2d = 0*X;
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
capImg2d(:) = griddata(x,y,grayness,X(:),Y(:),'linear');
figure
imagesc(lx, ly, capImg2d)
axis xy   % importante: orientamento corretto
hold on
colormap gray

for g = 1:length(greenY)
    scatter(greenX(g), greenY(g), 10, 'red', 'filled')
end

% Centroide azzurro
[cyan2Dx,cyan2Dy] = bst_project_2d(cyanCenter(1), cyanCenter(2), cyanCenter(3), 'stereo');
scatter(cyan2Dx,cyan2Dy, 30, 'blue','filled')

%%
% Raccogli tutti i centroidi rossi
allCenters = [greenX, greenY];
bluCenter = [cyan2Dx, cyan2Dy];
fczPos = bluCenter;

% 1. PCA sui centroidi
[coeff, score] = pca(allCenters);

% 2. Proietta
centered = allCenters - mean(allCenters, 1);
proj = centered * coeff(:,1:2);

% 3. Orienta PC1 usando FcZ come ancora
fczCentered = (fczPos - mean(allCenters,1)) * coeff(:,1);
if fczCentered < 0
    coeff(:,1) = -coeff(:,1);
    proj(:,1)  = -proj(:,1);
end
%% Assegna i centroidi verdi al giusti elettrodi
fczProj = fczCentered * coeff(:,1:2);
workProj = proj;
workCoeff = coeff;
if fczProj(1) < 0
    workCoeff(:,1) = -workCoeff(:,1);
    workProj(1,:) = -workProj(1,:);
end
% Trova Fz

%% Bounding box
% 4 angoli estremi nel piano PCA
[~,margin1maxIdx]= max(proj(:,1));
[~,margin1minIdx]= min(proj(:,1));
[~,margin2maxIdx]= max(proj(:,2));
[~,margin2minzdx]= min(proj(:,2));

bbCorners = allCenters([margin1maxIdx, margin1minIdx, margin2maxIdx, margin2minzdx], :);
xMargin = 0.1;
yMargin = 0.15;
x_min = min(bbCorners(:,1))-xMargin; x_max = max(bbCorners(:,1))+xMargin;
y_min = min(bbCorners(:,2))-yMargin; y_max = max(bbCorners(:,2))+yMargin;

rectangle('Position', [x_min, y_min, x_max-x_min, y_max-y_min], ...
    'EdgeColor', 'c', 'LineWidth', 2)
%% Prendi i canali EEG del file EEG, assegna il loro colore
chanStruct = in_bst_data('Giulia/Immagine/channel.mat');
channels = chanStruct.Channel;
greenChans = {'Fp1', 'Fp2', 'F7','Fc5','F3','T7','Tp9','P7','Cp5','O1', 'Oz', 'O2'...
    'P3', 'C3', 'Cz', 'Fc1','Fz','F4','F8','Fc6','T8','C4','Cp6','P8','Tp10',...
    'P4','Cp2','POz', 'Pz','Cp1', 'C1', 'C2'};
for k = 1:length(channels)
    if any(strcmp(channels(k).Name, greenChans))
        channels(k).Color = 'green';
    else
        channels(k).Color = 'white';
    end
end

locs = [channels.Loc]';
% Converti le posizioni 3D in 2D
[locs2dX,locs2dY] = bst_project_2d(locs(:,1), locs(:,2), locs(:,3), 'stereo');
for l = 1:length(locs2dY)
    scatter(locs2dX(l), locs2dY(l), 20, 'yellow','filled')
end


%% Estrai solo i punti nella bounding box
col_mask = lx >= x_min & lx <= x_max;
row_mask = ly >= y_min & ly <= y_max;
section = capImg2d(row_mask, col_mask);
% minRadius = 0;
% maxRadius = 3;
mask = section < 0.3;
el = strel('disk',5);
% mask = imdilate(mask, el);
mask = imclose(mask, el);
CC = bwconncomp(mask);
sizes = cellfun(@numel, CC.PixelIdxList);
toRemove = sizes > 500 | sizes < 1;
for k = find(toRemove)
    mask(CC.PixelIdxList{k}) = false;
end
props = regionprops(mask, 'Area', 'Perimeter', 'PixelIdxList');
% Circolarità = 4*pi*Area / Perimeter^2  (1 = cerchio perfetto, <1 = meno rotondo)
circularity = arrayfun(@(r) 4*pi*r.Area / r.Perimeter^2, props);

circThresh = [0.7, 3]; % aggiusta empiricamente
toRemove = find(circularity < circThresh(1) | circularity > circThresh(2));
for k = 1:length(toRemove)
    mask(props(toRemove(k)).PixelIdxList) = false;
end

figure
imagesc(mask)
props2 = regionprops(mask, 'Centroid', 'Area');

% Filtra per area minima (rimuovi rumore puntiforme)
minArea = 1;
validProps = props2([props2.Area] > minArea);

centroids = vertcat(validProps.Centroid); % Nx2 in pixel della section

% Converti in coordinate UV
lx_sec = lx(col_mask);
ly_sec = ly(row_mask);

figure
imagesc(lx_sec, ly_sec, section); colormap gray; axis xy; hold on

for k = 1:size(centroids,1)
    cx_px = round(centroids(k,1));
    cy_px = round(centroids(k,2));
    cx_px = max(1, min(cx_px, length(lx_sec)));
    cy_px = max(1, min(cy_px, length(ly_sec)));
    scatter(lx_sec(cx_px), ly_sec(cy_px), 50, 'r', 'filled')
end
[capCenters2d, capRadii2d] = imfindcircles(mask, [2 60]);
figure
imagesc(mask); axis equal;
hold on
for i =1:size(capCenters2d,1)
    scatter(capCenters2d(i,1),capCenters2d(i,2),20,'red', 'filled')
end
figure
imagesc(section)
axis equal
hold on
for i =1:size(capCenters2d,1)
    scatter(capCenters2d(i,1),capCenters2d(i,2),20,'red', 'filled')
end
%% Estrai il file dei canali, assegna i colori agli elettrodi
chanStruct = in_bst_data('')
%% Proietta i vertici sulla superficie 2D
sSurfCap = revoStruct;
    capCenters2d = [];
    capImg2d     = [];
    capRadii2d   = [];
    sSurfCap.u   = [];
    sSurfCap.v   = [];
    % if isempty(sSurfCap.Color)
    %     return
    % end

    % Image size [px]
    capImg2dSize = 900;
    capRangeinIm = 0.5;
    % Hyperparameters for circle detection [px]
    % NOTE: these values can vary for new caps
    minRadius = 10;
    maxRadius = 60;
    
    % Flatten the 3D mesh to 2D space using Stereographic projection
    [sSurfCap.u, sSurfCap.v] = bst_project_2d(sSurfCap.Vertices(:,1), sSurfCap.Vertices(:,2), sSurfCap.Vertices(:,3), 'stereo');
    
    % Perform image processing to detect the electrode locations
    % Convert to grayscale
    grayness = sSurfCap.Color*[1;1;1]/sqrt(3);
    
    % Interpolate and fit flattended mesh image from [-capRangeinIm to capRangeinIm] in a capImg2dSize square grid
    % NOTE: Should work with any flattened cap mesh but needs more testing
    ll=linspace(-capRangeinIm, capRangeinIm, capImg2dSize);
    [X,Y]=meshgrid(ll,ll);
    capImg2d = 0*X;
    warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
    capImg2d(:) = griddata(sSurfCap.u(1:end),sSurfCap.v(1:end),grayness,X(:),Y(:),'linear');
    warning('on','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
    mask = capImg2d < 0.3;
    el = strel('disk',5);
    % mask = imdilate(mask, el);
    mask = imclose(mask, el);
    CC = bwconncomp(mask);
    sizes = cellfun(@numel, CC.PixelIdxList);
    toRemove = sizes > 1500;
    for k = find(toRemove)
        mask(CC.PixelIdxList{k}) = false;
    end
    props = regionprops(mask, 'Area', 'Perimeter', 'PixelIdxList');

    % Circolarità = 4*pi*Area / Perimeter^2  (1 = cerchio perfetto, <1 = meno rotondo)
    circularity = arrayfun(@(r) 4*pi*r.Area / r.Perimeter^2, props);
    
    circThresh = 0.5; % aggiusta empiricamente
    toRemove = find(circularity < circThresh);
    for k = 1:length(toRemove)
        mask(props(toRemove(k)).PixelIdxList) = false;
    end
    [capCenters2d, capRadii2d] = imfindcircles(mask, [minRadius maxRadius]);
    figure
    imagesc(mask); axis equal;
    hold on
    for i =1:length(capCenters2d)
        scatter(capCenters2d(i,1),capCenters2d(i,2),20,'red', 'filled')
    end
    figure
    imagesc(capImg2d)
    hold on
    for i =1:length(capCenters2d)
        scatter(capCenters2d(i,1),capCenters2d(i,2),20,'red', 'filled')
    end
    colormap gray
    keyboard
    % Check if white color cap
    % if IsWhiteCap(sSurfCap.Color)
    %     capImg2d = imcomplement(capImg2d);
    % end
    
    % Detect the centers of the electrodes which appear as circles in the flattened image whose radii are in the range below
    warning('off','images:imfindcircles:warnForSmallRadius');
    warning('off','images:imfindcircles:warnForLargeRadiusRange');
    [capCenters2d, capRadii2d] = imfindcircles(capImg2d, [minRadius maxRadius]);

    figure
    hold on
    imagesc(capImg2d)
    set(gca, 'YDir','reverse')
    for i = 1:length(capRadii2d)
        scatter(capCenters2d(i,1),capCenters2d(i,2),10, 'red', 'filled')
    end
    warning('on','images:imfindcircles:warnForSmallRadius');
    warning('on','images:imfindcircles:warnForLargeRadiusRange');

%% Default channels (eeg cap 65 elettrodi brain products, salvata in brainstorm)
eegDefaults = bst_get('EegDefaults');
iCapFamily = find(strcmp({eegDefaults.name}, 'ICBM152'));
capIdx = find(contains({eegDefaults(iCapFamily).contents.fullpath}, 'BrainProducts_ActiCap_65'));
% Load the cap template
ChannelMat = in_bst_channel(eegDefaults(iCapFamily).contents(capIdx).fullpath);
ChanLoc = [ChannelMat.Channel.Loc]';
[~, iMaxLoc] = max(ChanLoc);
[~, iMinLoc] = min(ChanLoc);
% Find most anterior electrode  ~ FP2
frontElec = ChannelMat.Channel(iMaxLoc(1)).Name;
% Find most left electrode      ~ T7
leftElec  = ChannelMat.Channel(iMaxLoc(2)).Name;
% Find most right Electrode     ~ T8
rightElec = ChannelMat.Channel(iMinLoc(2)).Name;
% Find most posterior electrode ~ Oz
postElec  = ChannelMat.Channel(iMinLoc(1)).Name;
% Find most superior electrode  ~ Cz
topElec   = ChannelMat.Channel(iMaxLoc(3)).Name;
% Final list of landmarks
capLandmarkLabels = unique({frontElec, leftElec, rightElec, postElec, topElec}, 'stable');
knownElectPosDefault = [ChannelMat.Channel(iMaxLoc(1)).Loc, ChannelMat.Channel(iMaxLoc(2)).Loc,...
    ChannelMat.Channel(iMinLoc(2)).Loc, ChannelMat.Channel(iMinLoc(1)).Loc, ChannelMat.Channel(iMaxLoc(3)).Loc]';

%% EEG channels caricati su brainstorm
eegChanStruct = in_bst_data('Giulia/Immagine/channel.mat');
% Estrai le posizioni degli stessi elettrodi di cui sopra
eegLoc = [eegChanStruct.Channel.Loc]';

% Trova le posizioni di Fp2, T7, T8, Oz e Cz
fp2Pos = eegChanStruct.Channel(strcmp({eegChanStruct.Channel.Name}, 'Fp2')).Loc;
t7Pos = eegChanStruct.Channel(strcmp({eegChanStruct.Channel.Name}, 'T7')).Loc;
t8Pos = eegChanStruct.Channel(strcmp({eegChanStruct.Channel.Name}, 'T8')).Loc;
ozPos = eegChanStruct.Channel(strcmp({eegChanStruct.Channel.Name}, 'Oz')).Loc;
czPos = eegChanStruct.Channel(strcmp({eegChanStruct.Channel.Name}, 'Cz')).Loc;
% Final list of landmarks
knownEegElectPos = [fp2Pos, t7Pos, t8Pos, ozPos, czPos]';

%% Rotazione 3D dagli elettrodi sulla revoscan
[R,T] = rot3dfit(knownElectPosDefault, knownEegElectPos);

% Porta TUTTI i vertici della revoscan nello spazio dell'eeg
rotVerts = Vertices*R + ones(size(Vertices,1),1)*T;

%% Figura di prova
figure('WindowState','maximized')
subplot(121)
patch('Vertices', Vertices, 'Faces', Faces, 'FaceColor', 'interp',...
    'EdgeColor', 'none', 'FaceVertexCData', revoStruct.Color); axis equal;

subplot(122)
patch('Vertices', rotVerts, 'Faces', Faces, 'FaceColor', 'interp',...
    'EdgeColor', 'none', 'FaceVertexCData', revoStruct.Color); axis equal;
hold on
for ch = 1:length(eegLoc)
    scatter3(eegLoc(ch,1),eegLoc(ch,2),eegLoc(ch,3),100, 'red', 'filled')
end
hold off
