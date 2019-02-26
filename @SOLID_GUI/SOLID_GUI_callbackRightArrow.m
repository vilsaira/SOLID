function solidGui = SOLID_GUI_callbackRightArrow(solidGui, source, event)    
    shell = str2double(solidGui.popup.shell.String(solidGui.popup.shell.Value, :));
    shellInds = solidGui.solid.bval == shell;
    shellNums = zeros(size(shellInds));
    shellNums(shellInds) = 1:sum(shellInds);
    curDWI = solidGui.slider.dwi.Value;
    curDWI = find(shellNums == curDWI);
    maxDWI = find(shellNums == solidGui.slider.dwi.Max);
    curSlice = solidGui.slider.axi.Value;
    n = 0;
    flag = 0;
    for d = curDWI : 1 : maxDWI
        if n == 0
            cs = curSlice + 1;
            if cs == solidGui.slider.axi.Max
                cs = solidGui.slider.axi.Max;
            end                        
            n = 1;
        else
            cs = 1;                        
        end
        for s = cs : 1 : solidGui.slider.axi.Max
            z = solidGui.solid.modZ(s, d);
            if z > solidGui.solid.thresholdLower
                flag = 1;
                solidGui.slider.dwi.Value = shellNums(d);
                solidGui.txtfield.dwi.String = shellNums(d);
                solidGui.slider.axi.Value = s;                
                solidGui.txtfield.axi.String = s;
                break
            end
        end
        if flag
            break
        end
    end
    solidGui.SOLID_GUI_updateAxesEtc;  
end