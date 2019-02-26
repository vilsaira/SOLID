function solidGui = SOLID_GUI_callbackShowMask(solidGui, source, event)
    if source.Value
        source.BackgroundColor = 'green';
    else
        source.BackgroundColor = 0.94.*ones(1,3);
    end
    solidGui.SOLID_GUI_updateAxesEtc();
end