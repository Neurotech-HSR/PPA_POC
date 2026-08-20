function fillBsSeeds(surfFile)
surfStruct = in_bst_data(surfFile);
Atlases = surfStruct.Atlas;

Vertices = surfStruct.Vertices;
% For each atlas of the surface, set ROI seed if empty
for atl = 1:numel(Atlases)
    scouts = Atlases(atl).Scouts;
    if ~isempty(scouts)
        for sc = 1:numel(scouts)
            currScout = scouts(sc);
            if isempty(currScout.Seed)
                currScout.Vertices = sort(currScout.Vertices);
                % Get center of the region
                V = Vertices(currScout.Vertices,:);
                center = mean(V, 1);
                % Find the vertex that is closer to the center of the ROI
                [~, imin] = min(sum(bst_bsxfun(@minus, V, center) .^ 2, 2));
                surfStruct.Atlas(atl).Scouts(sc).Seed = currScout.Vertices(imin(1));
            end
        end
    end
end
bst_save(file_fullpath(surfFile), surfStruct, 'v7')
end