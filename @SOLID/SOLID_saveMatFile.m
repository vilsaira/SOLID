function SOLID_saveMatFile(oname, modZ, DWI, dims)

img = zeros(size(DWI), 'single');
for k = 1:size(modZ,1)
    for l = 1:size(modZ,2)
        img(:,:,k,l) = repmat(modZ(k,l), [size(img,1), size(img,2), 1,1]);
    end
end

E_DTI_write_nifti_file(img, dims, oname);

end