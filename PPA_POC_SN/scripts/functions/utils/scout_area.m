function area = scout_area(scoutStruct, surf, vertArea)
    verts = scoutStruct.Vertices;
    area = sum(vertArea(verts)) * 100 * 100;
end