function SOLID_EDTI_plugin_startSMECEPI(plugin, source, event)

    plugin.EDTI.UserData.solid.thrL = plugin.thrL;
    plugin.EDTI.UserData.solid.thrU = plugin.thrU;
    plugin.EDTI.UserData.solid.metric = plugin.metric;
    plugin.EDTI.UserData.solid.useMask = plugin.useMask;
    plugin.EDTI.UserData.solid.saveResults = plugin.saveResults;

    plugin.SOLID_EDTI_SMECEPI_Main;

end