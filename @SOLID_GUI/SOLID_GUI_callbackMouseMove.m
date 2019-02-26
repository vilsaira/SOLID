function solidGui = SOLID_GUI_callbackMouseMove(solidGui, source, event)    
    if strcmp(solidGui.data.window.timer.running, 'on')
        clim = solidGui.data.clim;
        level = mean(clim);
        width = clim(2)-clim(1);
        
        C = get(solidGui.ui, 'CurrentPoint');
        x = C(1,1);
        y = C(1,2);
        dx = (x-solidGui.data.window.mousePosI(1))/3+1;
        dy = (y-solidGui.data.window.mousePosI(2))/3+1;        
         
        width = width .* dx;
        level = level .* dy;
        
        if level > solidGui.solid.minmax(2)
            level = solidGui.solid.minmax(2);
        end
        if level < solidGui.solid.minmax(1)
            level = solidGui.solid.minmax(1);
        end
        if width < eps
            width = eps;
        end
        
        clim = [level - width/2, level + width/2];
                        
        solidGui.data.clim = clim;
        solidGui.ax.dwiTopLeft.CLim = clim;
        solidGui.ax.dwiTopCenter.CLim = clim;
        solidGui.ax.dwiTopRight.CLim = clim;
        solidGui.ax.dwiMiddleLeft.CLim = clim;
        solidGui.ax.dwiMiddleCenter.CLim = clim;
        solidGui.ax.dwiMiddleRight.CLim = clim;
        solidGui.ax.dwiBottomCenter.CLim = clim;

    end
end