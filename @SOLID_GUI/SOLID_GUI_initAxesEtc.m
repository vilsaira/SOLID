function solidGui = SOLID_GUI_initAxesEtc(solidGui, event, source)

    shell = str2double(solidGui.popup.shell.String(solidGui.popup.shell.Value, :));
    shellInds = solidGui.solid.bval == shell;

    dims = size(solidGui.solid.DWI(:,:,:,shellInds));
    if size(dims, 2) < 4
        dims(4) = 1;
    end

    solidGui.slider.dwi.Min = 1;
    solidGui.slider.dwi.Max = dims(4);
    solidGui.slider.dwi.SliderStep = [1, 1]./(max([dims(4),2])-1);
    solidGui.slider.dwi.Value = 1;
    
    solidGui.slider.slice.Min = 1;
    solidGui.slider.slice.Max = dims(3);
    solidGui.slider.slice.SliderStep = [1, 1]./(dims(3)-1);
    solidGui.slider.slice.Value = round(dims(3)/2);
    
    solidGui.slider.axi.Min = solidGui.slider.slice.Min;
    solidGui.slider.axi.Max = solidGui.slider.slice.Max;
    solidGui.slider.axi.SliderStep = [1, 1]./(dims(3)-1);
    solidGui.slider.axi.Value = solidGui.slider.slice.Value;
    
    solidGui.slider.cor.Min = 1;
    solidGui.slider.cor.Max = dims(1);
    solidGui.slider.cor.SliderStep = [1, 1]./(dims(1)-1);
    solidGui.slider.cor.Value = round(dims(1)/2);
    
    solidGui.slider.sag.Min = 1;
    solidGui.slider.sag.Max = dims(2);
    solidGui.slider.sag.SliderStep = [1, 1]./(dims(2)-1);
    solidGui.slider.sag.Value = round(dims(2)/2);
    
    solidGui.txtfield.dwi.Min = solidGui.slider.dwi.Min;
    solidGui.txtfield.dwi.Max = solidGui.slider.dwi.Max;           
    solidGui.txtfield.dwi.String = solidGui.slider.dwi.Value;
    
    solidGui.txtfield.axi.Min = solidGui.slider.axi.Min;
    solidGui.txtfield.axi.Max = solidGui.slider.axi.Max;
    solidGui.txtfield.axi.String = solidGui.slider.slice.Value;
    
    solidGui.txtfield.cor.Min = solidGui.slider.cor.Min;
    solidGui.txtfield.cor.Max = solidGui.slider.cor.Max;
    solidGui.txtfield.cor.String = solidGui.slider.cor.Value;
    
    solidGui.txtfield.sag.Min = solidGui.slider.sag.Min;
    solidGui.txtfield.sag.Max = solidGui.slider.sag.Max;
    solidGui.txtfield.sag.String = solidGui.slider.sag.Value;
    
    D = solidGui.solid.DWI(:,:,:,shellInds);
    solidGui.data.clim = [max([1, prctile(D(:), 5)]), ...
        max([1, prctile(D(:), 99)])];
    
    colormap(solidGui.ax.dwiTopLeft, 'gray');
    colormap(solidGui.ax.dwiTopCenter, 'gray');
    colormap(solidGui.ax.dwiTopRight, 'gray');
    colormap(solidGui.ax.dwiMiddleLeft, 'gray');
    colormap(solidGui.ax.dwiMiddleCenter, 'gray');
    colormap(solidGui.ax.dwiMiddleRight, 'gray');
    colormap(solidGui.ax.dwiBottomCenter, 'gray');     
        
    solidGui.txtfield.thresholdLower.String = num2str(solidGui.solid.thresholdLower);
    solidGui.txtfield.thresholdUpper.String = num2str(solidGui.solid.thresholdUpper);
    switch solidGui.solid.metric
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
    
    api = iptgetapi(solidGui.scrollPanel);
    mag = api.findFitMag();
    api.setMagnification(mag);

end