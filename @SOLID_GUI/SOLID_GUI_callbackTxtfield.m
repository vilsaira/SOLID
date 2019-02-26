function solidGui = SOLID_GUI_callbackTxtfield(solidGui, source, event)
    name = source.UserData.Name;
    % Invalid input forces the slider value to txtfield    
    if ismember(name, {'cor', 'sag', 'axi,' 'dwi'}) && ~all(ismember(source.String, '0123456789'))
        disp(['Error: Incorrect input to textfield', name]);
        solidGui.txtfield.(name).String = solidGui.slider.(name).Value;
        return;
    end   
    value = str2double(source.String);    
    if ismember(name, {'cor', 'sag', 'axi,' 'dwi'})
        'test'
        value = round(value);
        if value >= solidGui.slider.(name).Min && value <= solidGui.slider.(name).Max
            solidGui.slider.(name).Value = value;
            solidGui.txtfield.(name).String = value;
        else
            disp(['Error: Input to textfield', name, ' is out of image dimensions.']);
            solidGui.txtfield.(name).String = solidGui.slider.(name).Value;
        end
    elseif strcmp(name, 'thrL')
        if value >= 0 && value <= solidGui.solid.thresholdUpper
            solidGui.txtfield.(name).String = value;
            solidGui.solid.thresholdLower = value;
            solidGui.data.qSpace3D_updateFlag = true;
        end
    elseif strcmp(name, 'thrU')
        if value >= 0 && value > solidGui.solid.thresholdLower
            solidGui.txtfield.(name).String = value;
            solidGui.solid.thresholdUpper = value;
            solidGui.data.qSpace3D_updateFlag = true;
        end
    end
    if strcmp(source.UserData.Name, 'axi')
        solidGui.data.qSpace3D_updateFlag = true;
    end
    solidGui.SOLID_GUI_updateAxesEtc();
end