function X = SOLID_EDTI_WLLS_WW(S0,B,sw) %SOLID

    warning off all
    
    crit = 1;
    iter_max = 10;
    cn=0;
    con = 10^-4;
    
    Log_S0 = log(S0);
    W = ones(length(S0),1);
    X = ones(size(B,2),1);
    
    while crit==1
        
        cn = cn+1;
        W_p = W;
        X_p = X;
        w = diag(W_p.^2 .* sw); %SOLID
        X = (B'*w*B)\(B'*w)*Log_S0;
        W = exp(B*X);
        
        if all(abs(X-X_p) <= con*max(abs(X),abs(X_p))) || cn==iter_max
            crit=0;
        end
        
    end
    
    if any(isnan(X)) || any(isinf(X))
        X = B\S0; 
    end
    
    
    warning on all
    