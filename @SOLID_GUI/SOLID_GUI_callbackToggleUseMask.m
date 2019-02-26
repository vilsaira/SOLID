function solidGui = SOLID_GUI_callbackToggleUseMask(solidGui, source, event)    
    solidGui.SOLID_GUI_toggleEnableForFigure;
    % Check if user has done manual changes on modified Z-scores and notify them that those will be lost if they proceed.
    if any(solidGui.data.modZorig(~isnan(solidGui.data.modZorig)) ~= solidGui.solid.modZ(~isnan(solidGui.data.modZorig)))
        while(1)    
            choice = menu([{'You are about to enable/disable masking. '}, ...
            {'This will reset modified Z-score map and you will loose all manual changes. '}, ...
            {'Please, press "Reset mZ-map" to continue or "Cancel" to return.'}],'Reset mZ-map','Cancel');
            if choice==2 | choice==0
                solidGui.SOLID_GUI_toggleEnableForFigure;
                return;
            else
                break;
            end
        end
    end

    if solidGui.solid.useMask
        solidGui.solid.useMask = false;
        solidGui.menuSettingMask.Checked = 'off';
    else
        solidGui.solid.useMask = true;
        solidGui.menuSettingMask.Checked = 'on';
    end

    if solidGui.solid.useMask
        solidGui.data.modZorig = solidGui.solid.SOLID_calculateModZ(solidGui.solid.DWI, solidGui.solid.bval, solidGui.solid.metric, solidGui.solid.mask);
    else
        solidGui.data.modZorig = solidGui.solid.SOLID_calculateModZ(solidGui.solid.DWI, solidGui.solid.bval, solidGui.solid.metric);
    end

    solidGui.solid.modZ = solidGui.data.modZorig;
    solidGui.data.qSpace3D_updateFlag = true;
    solidGui.SOLID_GUI_updateAxesEtc();
    solidGui.SOLID_GUI_toggleEnableForFigure;
end