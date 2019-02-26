function nii = SOLID_loadNifti(fpath)

    if contains(fpath, '.mat')  % Assume ExploreDTI mat-file
        in = load(fpath, 'DWI');
        nii = single( cat( 4, in.DWI{ : } ) ); % Convert EDTI DWI cell to 4D matrix.
%         nii = flip(nii, 1);
        return
    end

    try
        nii = niftiread(fpath);
    catch
        try
            nii = load_untouch_nii(fpath);
            nii = nii.img;
        catch
            error('*.nii(.gz) file not found');
        end
    end            
    if any(size(nii) < 2)
        error('Error in NIfTI image dimensions');
    end    
    nii = permute(nii, [2 1 3 4]);
    nii = flip(nii,1);
end