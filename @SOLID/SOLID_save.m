function oname = SOLID_save(solid)
    i = strfind(solid.fName, '.nii');
    if isempty(i)
        i = strfind(solid.fName, '.mat');
    end
    oname = solid.fName(1:i-1);
    oname = strcat(solid.fDir, filesep, oname, '_L_', ...
        num2str(solid.thresholdLower), '_U_', ...
        num2str(solid.thresholdUpper), '_', ...
        solid.metric, '_masked_', num2str(solid.useMask));
    
    solid.SOLID_saveTxtModZ2D(oname, solid.modZ);    
    solid.SOLID_savePNG(oname, solid.modZ, solid.bval, solid.thresholdLower, solid.thresholdUpper);
    
    if contains(solid.fName, '.mat') % assume ExploreDTI mat file                
        solid.SOLID_saveMatFile([oname, '.nii'], solid.modZ, solid.DWI, solid.dims);
        return
    end
        
    try
        solid.SOLID_save4DNIfTI(oname, solid.modZ, solid.fDir, solid.fName);
    catch
        solid.SOLID_saveUntouchNii(oname, solid.modZ, solid.fDir, solid.fName);
    end
end