function sScouts = bst_setScoutRegion(sScouts,sRef)
     % Get total number of vertices
    Nvert = max([sScouts.Vertices, sRef.Vertices]);
    % Build reference atlas
    map = repmat({'UU'}, 1, Nvert);
    for iRef = 1:length(sRef)
        map(sRef(iRef).Vertices) = {sRef(iRef).Region};
    end
    % Assign region to all scouts with unknown regions
    for iScout = 1:length(sScouts)
        % If region is already defined: skip
        if ~ismember(sScouts(iScout).Region, {'UU','LU','RU','CU'})
            continue;
        end
        % Get the list of reference regions in this scout
        allReg = map(sScouts(iScout).Vertices);
        uniqueReg = unique(allReg);
        % Count the region with the most vertices
        countReg = cellfun(@(c)nnz(strcmpi(allReg, c)), uniqueReg);
        [~, iMax] = max(countReg);
        % Set it as the default region
        sScouts(iScout).Region = uniqueReg{iMax};
    end
end