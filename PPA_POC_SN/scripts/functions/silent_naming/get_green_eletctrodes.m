function [greenCenters] = get_green_eletctrodes(revoStruct)
Vertices   = revoStruct.Vertices;
Faces      = revoStruct.Faces;
hsv = rgb2hsv(revoStruct.Color);
H = hsv(:,1); S = hsv(:,2); V = hsv(:,3);
greenIdxs = find(H >= 0.22 & H <= 0.45 & S > 0.2 & V > 0.2);

greenClusters    = dbscan(Vertices(greenIdxs,:), 0.005, 20);
uniqueGreenLbls  = unique(greenClusters);
greenCenters     = [];
for g = 1:numel(uniqueGreenLbls)
    if uniqueGreenLbls(g) == -1; continue; end
    clustIdxs  = find(greenClusters == uniqueGreenLbls(g));
    greenCenters = [greenCenters; mean(Vertices(greenIdxs(clustIdxs),:), 1)];
end

% Merge centri troppo vicini
greenDists = triu(pdist2(greenCenters, greenCenters));
th = 0.015;
[ii, jj] = find(greenDists < th & greenDists > 0);
for k = 1:length(ii)
    greenCenters = [greenCenters; mean([greenCenters(ii(k),:); greenCenters(jj(k),:)], 1)];
end
greenCenters(unique(union(ii,jj)),:) = [];

end