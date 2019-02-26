function out = SOLID_GUI_sphericalInterpolation(bvec, bval, modZ, shellInds, currentAXI)

    %% Spherical interpolation
    % first rotate all bvecs so that the first direction is on
    % z-axis so interpolation near polar region works better.
    % Everything must be rotated back.
        
    
    
    i = find(bval > 0, 1, 'first');
    tmp = modZ(currentAXI,shellInds); % modZ values for specific SLICE
    
    v0 = [0,0,1];
    v1 = bvec(i,:);
    v = cross(v0,v1);
    s = sqrt(sum(v.^2));
    c = dot(v0,v1);
    V = [0,    -v(3),  v(2);...
            v(3),  0,    -v(1);...
        -v(2),  v(1),  0];
    R = eye(3,3) + V + V^2 * (1-c)/s^2;
    bvec = (R'*bvec')';
    
    tmp = [tmp, tmp];            
    [az, el] = cart2sph([bvec(shellInds,1); -bvec(shellInds,1)],...
        [bvec(shellInds,2); -bvec(shellInds,2)],...
        [bvec(shellInds,3); -bvec(shellInds,3)]);
    [~, inds] = sort(el);
    tmp = tmp(inds);
    az = az(inds);
    el = el(inds);
    % copy polar values to avoid twisting
    m = length(az);
    n = 7;
    azimuth = az;
    elevation = el;
    azimuth(n+1:m+n) = az;
    azimuth(1:n) = 2*pi*(1:-1/n:1/n)-pi;
    azimuth(m+n+1:m+2*n) = 2*pi*(1:-1/n:1/n)-pi;
    elevation(n+1:m+n) = el;
    elevation(1:n) = el(1).*ones(1,n);
    elevation(m+n+1:m+2*n) = el(end).*ones(1,n);
    tmp(n+1:m+n) = tmp;
    tmp(1:n) = tmp(1);
    tmp(m+n+1:m+2*n) = tmp(end);
    % mirror azimuth & elevation to interpolate edges correctly
    azimuth = repmat(azimuth, [3,1]);
    elevation = repmat(elevation, [3,1]);
    tmp = repmat(tmp', [3,1]);
    azimuth(1:m+2*n) = azimuth(1:m+2*n) - 2*pi;
    azimuth(2*m+4*n+1:3*m+6*n) = azimuth(2*m+4*n+1:3*m+6*n) + 2*pi;
    
    azimuth_lin = linspace(0, 2*pi, 3*m+3*n);
    elevation_lin = linspace(-pi/2, pi/2, 3*m+3*n);
    [azimuth_mat, elevation_mat] = meshgrid(azimuth_lin, elevation_lin);
    [Xi,Yi,Zi] = sph2cart(azimuth_mat, elevation_mat, ones(size(azimuth_mat)));                                                
    
    f = scatteredInterpolant(azimuth, elevation, double(tmp), 'natural', 'linear');
    Ci = f(azimuth_mat, elevation_mat);
    Ci(isnan(Ci)) = 0;
    
    % Get non-rotated coordinates
    v = R*[Xi(:)'; Yi(:)'; Zi(:)'];    
    Xi(:) = v(1,:);
    Yi(:) = v(2,:);
    Zi(:) = v(3,:);

    out.Xi = Xi;
    out.Yi = Yi;
    out.Zi = Zi;
    out.Ci = Ci;

    

end