function bval = SOLID_loadBVal(fpath)    
    % fpath should be the full path to bval file.
    
    if contains(fpath, '.mat')  % Assume ExploreDTI mat-file        
        in = load(fpath, 'b');        
        bval = round( sum( in.b( :, [ 1, 4, 6 ] ), 2 ), -2 ); % Round b-values to the nearest hundred.
        return
    end
    
    try 
        bval = load(fpath);
        bval = round(bval/100)*100;
    catch
        error('*.bval not found');
    end
end