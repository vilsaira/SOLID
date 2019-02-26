function SOLID_saveUntouchNii(oname, modZ, fpath, fname)
    
    hdr = load_untouch_nii(strcat(fpath, filesep, fname));
    img = zeros(size(hdr.img), 'single');
                
    for k = 1:size(modZ,1)
        for l = 1:size(modZ,2)
            img(:,:,k,l) = repmat(modZ(k,l), [size(img,1), size(img,2), 1,1]);
        end
    end
    
    hdr.img = img;
    hdr.hdr.hist.descrip = 'SOLID modified Z-score map';
    hdr.fileprefix = strcat(oname, '_SOLID_map');

    hdr.img = img;
    s = size(img);
    hdr.hdr.dime.dim(1) = length(s); 
    hdr.hdr.dime.dim(2:length(s)+1) = s; %[size(i,1),size(i,2),size(i,3)];
    if length(s) == 3
        hdr.hdr.dime.dim(5) = 1;
    end
    
    dat = whos('img');
    switch (dat.class)
        case 'logical'
            img = single(img);
            hdr.hdr.dime.datatype = 16;
            hdr.hdr.dime.bitpix = 32;        
        case 'uint8'
            hdr.hdr.dime.datatype = 2;
            hdr.hdr.dime.bitpix = 8;
        case 'int8'
            hdr.hdr.dime.datatype = 256;
            hdr.hdr.dime.bitpix = 8;
        case 'int16'
            hdr.hdr.dime.datatype = 4;
            hdr.hdr.dime.bitpix = 16;
        case 'uint16'
            hdr.hdr.dime.datatype = 512;
            hdr.hdr.dime.bitpix = 16;
        case 'single'
            hdr.hdr.dime.datatype = 16;
            hdr.hdr.dime.bitpix = 32;
        case 'double'
            hdr.hdr.dime.datatype = 64;
            hdr.hdr.dime.bitpix = 64;
        otherwise
            hdr.hdr.dime.datatype = 16;
            hdr.hdr.dime.bitpix = 32;        
    end
    [fpath, filesep, hdr.fileprefix, '.nii']
    save_untouch_nii( hdr, [hdr.fileprefix, '.nii'] );
    try
        system(['gzip -f ', hdr.fileprefix, '.nii']);
    catch
    end
end