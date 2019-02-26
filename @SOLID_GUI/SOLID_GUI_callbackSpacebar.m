function solidGui = SOLID_GUI_callbackSpacebar(solidGui, source, event)    
    shell = str2double(solidGui.popup.shell.String(solidGui.popup.shell.Value, :));
    shellInds = solidGui.solid.bval == shell;
    shellNums = zeros(size(shellInds));
    shellNums(shellInds) = 1:sum(shellInds);
    curDWI = solidGui.slider.dwi.Value;
    curDWI = find(shellNums == curDWI);
    curSlice = solidGui.slider.axi.Value;
    
    if ismember('control', get(solidGui.ui, 'currentModifier'))
        z = solidGui.solid.modZ(:, curDWI);
        zorig = solidGui.data.modZorig(:, curDWI);                    
        if any(abs(z - zorig) > eps)
            solidGui.solid.modZ(:, curDWI) = solidGui.data.modZorig(:, curDWI);
        else
            solidGui.solid.modZ(:, curDWI) = 0;
        end
    else
        z = solidGui.solid.modZ(curSlice, curDWI);
        zorig = solidGui.data.modZorig(curSlice, curDWI);
        if abs(z - zorig) > eps
            solidGui.solid.modZ(curSlice, curDWI) = solidGui.data.modZorig(curSlice, curDWI);
        else
            solidGui.solid.modZ(curSlice, curDWI) = 0;
        end
    end
    solidGui.SOLID_GUI_updateAxesEtc;  
end