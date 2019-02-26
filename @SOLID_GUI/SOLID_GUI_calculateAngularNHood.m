function aHood = SOLID_calculateAngularNHood(b, v)
    n = sqrt(sum(v.^2,2));
    v = v./repmat(n, [1,3]);
    v(isnan(v)) = 0;
    ad = inf(length(b), length(b));    
    for i = 1:length(b)
        for j = i:length(b)
            vi = v(i, :);
            vj = v(j, :);
            euc = sqrt(sum( (b(j).^2.*vj-b(i).^2.*vi).^2))./b(i);
            ad(i,j) = euc;
        end
    end
    ad(isnan(ad)) = 0;
    ad = triu(ad)+triu(ad,1)';
    [~, aHood] = sort(ad,2);
end