function [centroid] = get_cyan_eletctrode(revoStruct)

Vertices = revoStruct.Vertices;
Faces    = revoStruct.Faces;

% Estrai vertici azzurri (FcZ)
hsv = rgb2hsv(revoStruct.Color);
H = hsv(:,1); S = hsv(:,2); V = hsv(:,3);
cyanIdxs = find(H >= 0.50 & H <= 0.65 & S > 0.3 & V > 0.3);

clusters = dbscan(Vertices(cyanIdxs,:), 0.05, 5);
nClusters = max(clusters);

if nClusters > 1
    scores = struct('size', [], 'sph',[]);
    
    for c = 1:nClusters
        mask     = clusters == c;
        pts      = Vertices(cyanIdxs(mask), :);
        clustSize = sum(mask);
        
        % Sfericità via convex hull
        sphericity = compute_sphericity(pts);
        
        % Punteggio combinato: normalizza dimensione [0,1] poi pondera
        % (la normalizzazione avviene dopo, qui accumula raw)
        scores(c) = struct('size', clustSize, 'sph', sphericity);
    end
    
    % Raccogli separatamente per normalizzare
    sizes = arrayfun(@(c) sum(clusters == c), 1:nClusters)';
    sphs  = zeros(nClusters, 1);
    for c = 1:nClusters
        pts       = Vertices(cyanIdxs(clusters == c), :);
        sphs(c)   = compute_sphericity(pts);
    end
    
    % Normalizza su [0,1]
    norm_size = (sizes - min(sizes)) / (max(sizes) - min(sizes) + eps);
    norm_sph  = (sphs  - min(sphs))  / (max(sphs)  - min(sphs)  + eps);
    
    % Pesi: puoi regolare (0.5/0.5 = peso uguale)
    w_size = 0;
    w_sph  = 1;
    combined = w_size * norm_size + w_sph * norm_sph;
    
    [~, bestClust] = max(combined);
    cyanIdxs = cyanIdxs(clusters == bestClust);
end

centroid = mean(Vertices(cyanIdxs,:), 1);
end


function sph = compute_sphericity(pts)
% Sfericità = pi^(1/3) * (6*V)^(2/3) / A
% calcolata sul convex hull dei punti
    if size(pts, 1) < 4
        sph = 0;
        return
    end
    try
        [~, V] = convhull(pts, 'Simplify', true);  % volume
        % superficie: somma aree delle facce del convex hull
        [K, ~] = convhull(pts, 'Simplify', true);
        v1 = pts(K(:,2),:) - pts(K(:,1),:);
        v2 = pts(K(:,3),:) - pts(K(:,1),:);
        A  = sum(0.5 * vecnorm(cross(v1, v2, 2), 2, 2));
        sph = (pi^(1/3)) * ((6*V)^(2/3)) / A;
    catch
        sph = 0;
    end
end