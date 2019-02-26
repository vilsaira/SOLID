function solidGui = SOLID_GUI_callbackMouseClick(solidGui, source, event)
    C = get(solidGui.ax.modZ2D, 'CurrentPoint');
    curDWI = round(C(1,1));
    curAXI = round(C(1,2));
    xmin = solidGui.ax.modZ2D.XLim(1);
    xmax = solidGui.ax.modZ2D.XLim(2);
    ymin = solidGui.ax.modZ2D.YLim(1);
    ymax = solidGui.ax.modZ2D.YLim(2);

    if curDWI >= xmin && curDWI <= xmax && ...
        curAXI >= ymin && curAXI <= ymax
        % This triggers if user clicks modified Z-score 2D map
        solidGui.slider.dwi.Value = curDWI;
        solidGui.txtfield.dwi.String = curDWI;
        solidGui.slider.axi.Value = curAXI;
        solidGui.txtfield.axi.String = curAXI;
        solidGui.SOLID_GUI_updateAxesEtc();    
    else
        % Change window / level
        C = get(solidGui.ui, 'CurrentPoint');
        solidGui.data.window.mousePosI = [C(1,1), C(1,2)];
        clim = solidGui.ax.dwiMiddleCenter.CLim;
        solidGui.data.window.width = clim(2)-clim(1);
        solidGui.data.window.level = solidGui.data.window.width/2 + clim(1);        
        start(solidGui.data.window.timer);
    end
end