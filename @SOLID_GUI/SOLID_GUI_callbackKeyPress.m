function solidGui = SOLID_GUI_callbackKeyPress(solidGui, source, event)    
    switch event.Key
    case 'space'
        % Toggle outlier on/off on the modified Z-score 2D map
        solidGui.SOLID_GUI_callbackSpacebar;
    case 'rightarrow'
        % Jump to next outlier
        solidGui.SOLID_GUI_callbackRightArrow;
    case 'leftarrow'
        % Jump to previous outlier
        solidGui.SOLID_GUI_callbackLeftArrow;
    case 'jacksparrow'
        disp('Why is the rum always gone?')
    end
end