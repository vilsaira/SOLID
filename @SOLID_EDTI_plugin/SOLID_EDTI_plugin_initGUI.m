function plugin = SOLID_EDTI_plugin_initGUI(plugin, source, event)
% Function SOLID_EDTI_plugin_initGUI adds SOLID menus to ExploreDTI Settings and Plugins menus
    
% Find ExploreDTI GUI element
    figHandles = get(groot, 'Children');        
    plugin.EDTI = findobj(figHandles, 'Tag', 'MainExploreDTI');    
    
    if isempty(plugin.EDTI)    
        error('SOLID error: ExploreDTI v. 4.8.6 is not found in the current Matlab session.')
    end

    pluginMenu = findobj(plugin.EDTI.Children, 'Tag', 'plugins');

    oldSMECEPI = findobj(pluginMenu, 'Tag', 'SM_EC_EPI_C');
    oldSMECEPI.ForegroundColor = 0.75*ones(1,3);

    % Add SOLID menu element to Plugins        
    qa = uimenu(pluginMenu);
    qa.Tag = 'SOLID_QA';
    qa.Label = 'SOLID QA GUI';
    try
        qa.MenuSelectedFcn = @plugin.SOLID_EDTI_plugin_QA;
    catch
        qa.Callback = @plugin.SOLID_EDTI_plugin_QA;
    end    

    smecepi = uimenu(pluginMenu);
    smecepi.Tag = 'SOLID_SM_EC_EPI_C';
    smecepi.Label = 'SOLID - Correct for subject motion & EC/EPI distortions';
    smecepi.ForegroundColor = 'red';    
    try
        smecepi.MenuSelectedFcn = @plugin.SOLID_EDTI_plugin_startSMECEPI;
    catch
        smecepi.Callback = @plugin.SOLID_EDTI_plugin_startSMECEPI;
    end

    N = length(pluginMenu.Children);
    pluginMenu.Children = pluginMenu.Children([2:N-1, 1, N]); % Sort menu elements

    % Add SOLID menu element to Settings
    settingsMenu = findobj(plugin.EDTI.Children, 'Tag', 'settings');
    solidSettings = uimenu(settingsMenu);
    solidSettings.Tag = 'SOLID_settings';
    solidSettings.Label = 'SOLID';

    solidSettingsParameters = uimenu(solidSettings);
    solidSettingsParameters.Tag = 'SOLID_settings_set';
    solidSettingsParameters.Label = 'Set parameters';
    try            
        solidSettingsParameters.MenuSelectedFcn = @plugin.SOLID_EDTI_plugin_settings;
    catch
        solidSettingsParameters.Callback = @plugin.SOLID_EDTI_plugin_settings;
    end 
    
    solidSettingsExport = uimenu(solidSettings);
    solidSettingsExport.Tag = 'SOLID_export';
    solidSettingsExport.Label = 'Append SOLID info to Parameters_SM_EC_EPI.txt';
    try            
        solidSettingsExport.MenuSelectedFcn = @plugin.SOLID_EDTI_plugin_settingsExport;
    catch
        solidSettingsExport.Callback = @plugin.SOLID_EDTI_plugin_settingsExport;
    end 

end