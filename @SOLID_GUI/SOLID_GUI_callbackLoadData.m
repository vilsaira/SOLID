function solidGui = SOLID_GUI_callbackLoadData(solidGui, source, event)
    solidGui.SOLID_GUI_toggleEnableForFigure;

    try
        [dwiName, dwiPath] = uigetfile({'*.nii;*.nii.gz'}, 'Select 4D DWI NIfTI file');    
        if isempty(dwiName)
            return
        end
        niifPath = [dwiPath, filesep, dwiName];
        bvalPath = solidGui.SOLID_GUI_loadBValVec(dwiPath, dwiName, 'bval');
        bvecPath = solidGui.SOLID_GUI_loadBValVec(dwiPath, dwiName, 'bvec');
        maskPath = solidGui.SOLID_GUI_loadMask(dwiPath, dwiName);

        solidGui.solid = SOLID('in', niifPath,...
                                'bval', bvalPath,...
                                'bvec', bvecPath,...
                                'mask', maskPath,...
                                'maskTune', solidGui.data.maskTune,...
                                'maskFilter', solidGui.data.maskFilter,...
                                'metric', solidGui.data.metric,...
                                'thrL', str2double(solidGui.txtfield.thresholdLower.String),...
                                'thrU', str2double(solidGui.txtfield.thresholdUpper.String),...
                                'save', false);

        solidGui.data.modZorig = solidGui.solid.modZ;
        solidGui.data.nearestVecs = solidGui.SOLID_GUI_calculateAngularNHood(solidGui.solid.bval, solidGui.solid.bvec);    

        solidGui.popup.shell.String = solidGui.solid.uniqueb;
        solidGui.popup.shell.Value = 1;        
        solidGui.data.qSpace3D_updateFlag = true;
        solidGui.SOLID_GUI_initAxesEtc();
        solidGui.SOLID_GUI_updateAxesEtc();
    catch
        
    end
    solidGui.SOLID_GUI_toggleEnableForFigure;
end