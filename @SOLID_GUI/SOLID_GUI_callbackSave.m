function solidGui = SOLID_GUI_callbackSave(solidGui, source, event)
    solidGui.SOLID_GUI_toggleEnableForFigure;

    oname = solidGui.solid.SOLID_save(solidGui.solid);
    solidGui.SOLID_GUI_saveManualLabels(oname, solidGui.solid.modZ, solidGui.data.modZorig);

    solidGui.SOLID_GUI_toggleEnableForFigure;
end