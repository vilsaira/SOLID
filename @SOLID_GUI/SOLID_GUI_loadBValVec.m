function fpath = SOLID_GUI_loadBValVec(solidGui, dwiPath, dwiName, ftype)
    % Input to SOLID_GUI_loadBValVec function should be the full path to DWI NIfTI file. Function searches for type = bval/bvec file based on respective file suffix.

    fpath = [];   

    extInd = strfind(dwiName, '.nii');
    tmpfpath = [dwiPath, filesep, dwiName(1:extInd-1), ['.', ftype]];
    if ~exist(tmpfpath, 'file') || ...                           % If the file does not exist
        ismember('control', get(solidGui.ui, 'currentModifier')) % or user press control 
        [fname, fdir] = uigetfile(fullfile(dwiPath, ['*.', ftype]), ['Select .', ftype, ' / press Cancel for default']);
        if isequal(fname, 0)
            fpath = [dwiPath, filesep, dwiName(1:extInd-1), ['.', ftype]];
        else
            fpath = [fdir, filesep, fname];
        end
    else
        fpath = tmpfpath;
    end
    
    
end