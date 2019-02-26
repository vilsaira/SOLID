function solid =  SOLID_processData(solid, p)
    disp(' ');
    if isempty(p.in)
        error('DWI file not found');
    end
    tmp = strsplit(p.in, filesep);
    dwiName = tmp{end};
    ext = 4;
    if contains(dwiName, '.nii.gz')
        ext = 7;
    end
    dwiPath = strjoin(tmp(1:end-1), filesep);
    if isempty(dwiPath)
        dwiPath = pwd;
    end
    
    bvalPath = [dwiPath, filesep, dwiName(1:end-ext), '.bval'];
    if ~isempty(p.bval)
        bvalPath = p.bval;
    end
    
    bvecPath = [dwiPath, filesep, dwiName(1:end-ext), '.bval'];
    if ~isempty(p.bvec)
        bvecPath = p.bvec;
    end
    
    if ~isempty(p.mask)
        maskPath = p.mask;
    end

    disp(['Loading ' p.in]);
    solid.DWI = solid.SOLID_loadNifti(p.in);    
    solid.dims = size(solid.DWI);
    solid.minmax = [min(solid.DWI(:)), max(solid.DWI(:))];
    solid.bval = solid.SOLID_loadBVal(bvalPath);
    solid.uniqueb = unique(solid.bval);
    solid.bvec = solid.SOLID_loadBVec(bvecPath);
    
    solid.mask = ones(size(solid.DWI(:,:,:,1)));
    solid.useMask = 0;
    solid.maskTune = p.maskTune;
    solid.maskFilter = p.maskFilter;
    if isempty(p.mask)
        try
            disp(['Calculating brain mask']);
            solid.mask = solid.E_DTI_Create_Mask_From_DWI_enhanced_IND(solid.DWI(:,:,:,1), p.maskTune, p.maskFilter); 
            solid.useMask = 1;
        catch
        end
    elseif ~strcmp(p.mask, 'no')
        solid.mask = solid.SOLID_loadNifti(maskPath);
        solid.useMask = 1;
    end
    
    solid.fName = dwiName;
    solid.fDir = dwiPath;
    solid.fullPath = [dwiPath, filesep, dwiName];    
    solid.metric = p.metric;
    solid.thresholdLower = p.thrL;
    solid.thresholdUpper = p.thrU;
    
    disp('Calculating modified Z-scores');        
    solid.modZ = solid.SOLID_calculateModZ(solid.DWI, solid.bval, solid.metric, solid.mask);    
    if p.save
        disp('Saving results');
        oname = solid.SOLID_save(solid);
    end
    disp('SOLID is done!');
    disp(' ');
end