function solidGui = SOLID_GUI_callbackSlider(solidGui, source, event)
    value = round(source.Value);    
    solidGui.txtfield.(source.UserData.Name).String = value;
    if strcmp(source.UserData.Name, 'axi')
        solidGui.data.qSpace3D_updateFlag = true;
    end
    solidGui.SOLID_GUI_updateAxesEtc();
end