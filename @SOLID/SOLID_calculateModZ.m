function modZ = SOLID_calculateModZ(DWI, b, metric, mask)    
    if nargin < 2
        error('Cannot calculate modified Z-score: incorrect DWI or bval.');
        return;
    end         
    if nargin > 2 && ~ismember(metric, {'Variance', 'Mean', 'Iod'})
        error('Cannot calculate modified Z-score: metric must be Variance/Mean/Iod');
        return;
    end

    DWI(DWI < 0 ) = NaN; % User can opt to use their own mask by setting non-brain voxels to negative numbers.
    dims = size(DWI);
    if nargin == 4
        mask = mask(:,:,:,1) > 0;
        DWI(repmat(~mask, [1,1,1,dims(4)])) = NaN;
    end

    modZ = NaN(dims(3:4), 'single');
    uniqueb = unique(b);

    for i = 1 : length(uniqueb)
        shell = (b == uniqueb(i));
        N = sum(shell);
        % if ( uniqueb(i) == 0 )
        %     continue;
        % end
        if ( N < 2 )
            disp(['SOLID: Not enough data with b-value ', num2str(uniqueb(i))]);
            continue;
        end                
        baseline = double(DWI(:,:,:,shell));
        baseline(baseline < eps) = NaN;
        baseline = reshape(baseline, [dims(1)*dims(2), dims(3), N]);
        switch metric
            case 'Mean'
                y = squeeze(mean(baseline, 'omitnan'));
            case 'Iod'
                y = squeeze(var(baseline, 'omitnan'))./squeeze(mean(baseline, 'omitnan'));
            otherwise % Variance
                y = squeeze(var(baseline, 'omitnan'));
        end
        tmp = repmat(median(y,2,'omitnan'), [1,N]);
        MAD = 1.4826 * median( abs(y-tmp),2,'omitnan');
        modZ(:,shell) = abs(y-tmp)./repmat(MAD,[1,N]);                
    end

end