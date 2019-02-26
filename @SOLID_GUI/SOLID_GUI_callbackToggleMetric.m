function solidGui = SOLID_GUI_callbackToggleMetric(solidGui, source, event)
    solidGui.SOLID_GUI_toggleEnableForFigure;
    % Check if user has done manual changes on modified Z-scores and notify them that those will be lost if they proceed.
    if any(solidGui.data.modZorig(~isnan(solidGui.data.modZorig)) ~= solidGui.solid.modZ(~isnan(solidGui.data.modZorig)))
        while(1)    
            choice = menu([{'You are about to change SOLID metric.'}, ...
            {'This will reset modified Z-score map and you will loose all manual changes.'}, ...
            {'Please, press "Reset mZ-map" to continue or "Cancel" to return.'}],'Reset mZ-map','Cancel');
            if choice==2 | choice==0
                solidGui.SOLID_GUI_toggleEnableForFigure;
                return;
            else
                break;
            end
        end
    end
    
    name = source.UserData.Name;
    switch name
        case 'Variance'
            solidGui.menuSettingMetricVar.Checked = 'on';
            solidGui.menuSettingMetricMean.Checked = 'off';
            solidGui.menuSettingMetricIod.Checked = 'off';
        case 'Mean'
            solidGui.menuSettingMetricVar.Checked = 'off';
            solidGui.menuSettingMetricMean.Checked = 'on';
            solidGui.menuSettingMetricIod.Checked = 'off';
        case 'Iod'
            solidGui.menuSettingMetricVar.Checked = 'off';
            solidGui.menuSettingMetricMean.Checked = 'off';
            solidGui.menuSettingMetricIod.Checked = 'on';
    end
    solidGui.solid.metric = name;

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