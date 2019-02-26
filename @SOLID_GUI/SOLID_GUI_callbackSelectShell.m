function solidGui = SOLID_GUI_callbackSelectShell(solidGui, source, event)
    solidGui.SOLID_GUI_toggleEnableForFigure;
    solidGui.data.qSpace3D_updateFlag = true;
    solidGui.SOLID_GUI_initAxesEtc();
    solidGui.SOLID_GUI_updateAxesEtc();
    solidGui.SOLID_GUI_toggleEnableForFigure;
end