function solidGui = SOLID_GUI_callbackMouseScroll(solidGui, source, event)
    C = get(solidGui.ui, 'CurrentPoint');
    x = C(1,1);
    y = C(1,2);
    % [x,y]
    Pcor = solidGui.ax.dwiTopLeft.Position;
    Psag = solidGui.ax.dwiTopRight.Position;
    Paxi = solidGui.ax.dwiMiddleCenter.Position;
    PmodZ2D = solidGui.scrollPanel.Position;
    fieldName = [];
    if x > Pcor(1) && x < Pcor(1)+Pcor(3) && y > Pcor(2) && y < Pcor(2)+Pcor(4)
        fieldName = 'cor';            
    elseif x > Psag(1) && x < Psag(1)+Psag(3) && y > Psag(2) && y < Psag(2)+Psag(4)
        fieldName = 'sag';
    elseif x > Paxi(1) && x < Paxi(1)+Paxi(3) && y > Paxi(2) && y < Paxi(2)+Paxi(4)
        fieldName = 'axi';
    elseif x > PmodZ2D(1) && x < PmodZ2D(1)+PmodZ2D(3) && y > PmodZ2D(2) && y < PmodZ2D(2)+PmodZ2D(4)
        % If cursor is above modified Z-score map scroll will change slices
        fieldName = 'axi';
    end
    if ~isempty(fieldName) && ismember('control', get(solidGui.ui, 'currentModifier'))
        % If control is pressed scroll will change only dwi
        fieldName = 'dwi';
    end
    
    if ~isempty(fieldName)
        tmp = solidGui.slider.(fieldName).Value - event.VerticalScrollCount;
        if (tmp <= solidGui.slider.(fieldName).Max) && (tmp >= solidGui.slider.(fieldName).Min)
            solidGui.slider.(fieldName).Value = tmp;
            solidGui.txtfield.(fieldName).String = tmp;
            if strcmp(fieldName, 'axi')
                solidGui.data.qSpace3D_updateFlag = true;
            end
            solidGui.SOLID_GUI_updateAxesEtc;
        end
    end

end