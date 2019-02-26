function maskPath = SOLID_GUI_loadMask(solidGui, dwiPath, dwiName)

    maskPath = [];
    if strcmp(solidGui.menuSettingMask.Checked, 'off')
        maskPath = 'no'
        return;
    end

    extInd = strfind(dwiName, '.nii');
    maskPath = [dwiPath, filesep, dwiName(1:extInd-1), ['_mask.nii.gz']];
    if ~exist(maskPath, 'file')
        maskPath = [dwiPath, filesep, dwiName(1:extInd-1), ['_mask.nii']];
    end
    
    if ~exist(maskPath, 'file') || ... % If mask file is not found 
        ismember('contro', get(solidGui.ui, 'currentModifier')) % or user press control
        [fname, fdir] = uigetfile(fullfile(dwiPath, {'*.nii;*.nii.gz'}), 'Select mask NIfTI file / press Cancel for default');
        maskPath = [fdir, filesep, fname];
        if isequal(fname, 0)
            maskPath = [];
        end
    end

end