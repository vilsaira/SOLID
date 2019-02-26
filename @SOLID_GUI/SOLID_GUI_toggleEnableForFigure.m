function solidGui = SOLID_GUI_toggleEnableForFigure(solidGui, source, event)
    cList = solidGui.ui.Children;
    for c = 1:length(cList)
        try
            switch solidGui.ui.Children(c).Enable
            case 'on'
                solidGui.ui.Children(c).Enable = 'off';
            otherwise
                solidGui.ui.Children(c).Enable = 'on';
            end
        catch
        end
    end
    drawnow;
end