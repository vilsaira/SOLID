function SOLID_save4DNIfTI(oname, modZ, fpath, fname)
    info = niftiinfo(strcat(fpath, filesep, fname));
    img = nifti1read(info);
    img = single(img);
    for k = 1:size(modZ,1)
        for l = 1:size(modZ,2)
            img(:,:,k,l) = repmat(modZ(k,l), [size(img,1), size(img,2), 1,1]);
        end
    end
    img = single(img);
    info.Description = 'SOLID modified Z-score map';
    info.Datatype = 'single';
    info.BitsPerPixel = 32;
    niftiwrite(img, [oname, '_SOLID_map.nii'], info);
    try
        system(['gzip -f ', oname, '_SOLID_map.nii']);
    catch
    end
end