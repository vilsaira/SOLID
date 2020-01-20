function plugin = SOLID_EDTI_plugin_QA(plugin, source, event)
% Function SOLID_EDTI_plugin_QA passess *.mat file information loaded to ExploreDTI to SOLID GUI.
% If a *.mat file is not loaded to ExploreDTI nothing happens.

    if isempty(plugin.EDTI.UserData.DTI.MatfilePath)
        return;
    end
    plugin.SOLID_GUI = SOLID_GUI;
    plugin.SOLID_GUI.SOLID_GUI_toggleEnableForFigure;

    plugin.SOLID_GUI.data.maskTune = plugin.EDTI.UserData.mask_settings.tune_1;
    plugin.SOLID_GUI.data.maskFilter = plugin.EDTI.UserData.mask_settings.mfs;
    plugin.SOLID_GUI.data.metric = plugin.metric;    
    
    niifPath = plugin.EDTI.UserData.DTI.MatfilePath;
    bvalPath = plugin.EDTI.UserData.DTI.MatfilePath;
    bvecPath = plugin.EDTI.UserData.DTI.MatfilePath;
    maskPath = [];    

    plugin.SOLID_GUI.solid = SOLID('in', niifPath,...
                            'bval', bvalPath,...
                            'bvec', bvecPath,...
                            'mask', maskPath,...
                            'maskTune', plugin.SOLID_GUI.data.maskTune,...
                            'maskFilter', plugin.SOLID_GUI.data.maskFilter,...
                            'metric', plugin.SOLID_GUI.data.metric,...
                            'thrL', plugin.thrL,...
                            'thrU', plugin.thrU,...
                            'save', false);
    
    plugin.SOLID_GUI.data.modZorig = plugin.SOLID_GUI.solid.modZ;    
    plugin.SOLID_GUI.data.nearestVecs = plugin.SOLID_GUI.SOLID_GUI_calculateAngularNHood(plugin.SOLID_GUI.solid.bval, plugin.SOLID_GUI.solid.bvec);    

    plugin.SOLID_GUI.popup.shell.String = plugin.SOLID_GUI.solid.uniqueb;
    plugin.SOLID_GUI.popup.shell.Value = 1;
    plugin.SOLID_GUI.data.qSpace3D_updateFlag = true;

    plugin.SOLID_GUI.SOLID_GUI_initAxesEtc();
    plugin.SOLID_GUI.SOLID_GUI_updateAxesEtc();   
    
    plugin.SOLID_GUI.SOLID_GUI_toggleEnableForFigure;
end