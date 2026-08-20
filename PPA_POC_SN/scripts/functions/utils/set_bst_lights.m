function set_bst_lights()
    camlight(gca,   0,  40); camlight(gca, 180, 40)
    camlight(gca,   0, -90); camlight(gca,  90,  0); camlight(gca, -90, 0)
    light('Tag','FrontLight','Color',[0.8 0.8 0.8],'Style','infinite','Parent',gca)
    camlight(gca,'headlight')
    material(gca,'dull')
end