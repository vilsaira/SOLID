function bvec = SOLID_loadBVec(fpath)
    % fpath should be the full path to bvec file.
    
    if contains(fpath, '.mat') % Assume ExploreDTI mat-file
        in = load(fpath, 'b', 'g');        
        b = round( sum( in.b( :, [ 1, 4, 6 ] ), 2 ), -2 ); % Round b-values to the nearest hundred.
        nonzero = b ~= 0;
        bvec = zeros(length(b), 3);
        bvec(nonzero, :) = in.g;
        return
    end
    
    try 
        bvec = load(fpath);
        if size(bvec,1) < size(bvec,2)
            bvec = bvec';
        end
    catch
        error('*.bvec not found');
    end
end