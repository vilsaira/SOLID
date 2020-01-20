function [X, is_outlier] = SOLID_EDTI_robust_Linear_fit_new(S0,B,kappa,iter_max,con, sw) %SOLID
    % Note: solid weights 'sw' are applied only during the final IWLLS
    % estimate. REKINDLE could be adapted to use SOLID information to flag
    % outliers but I haven't investigated this well enough to actually
    % recommend such thing. -VS
    warning off all
    
    % Initialization
    crit_all = 1;
    cn_all = 0; 
    Max_inter_in = 5;
    
    % Step 1: Initial LLS fit
    w = ones(size(S0));
    LogS0 = log(S0);
    W = diag(w);
    X = (B'*W*B)\(B'*W)*LogS0;
    
    while crit_all==1
        
        % Initialization
        cn_all = cn_all + 1;
        X_p_all = X;
        cn = 0;
        crit = 1;
        
        % Step 2: Compute a robust estimate for homoscedastic regression using IRLS.
        while crit==1
    
            % Initialization
            cn = cn+1;
            X_p = X;
            
            fit = B*X; 
            % a. Calculate the residuals e in the linear domain
            res = (LogS0-fit); 
            if all(res==0)
                break;
            end
            % b. Obtain an estimate of the dispersion of the residuals by
            % calculating the median absolute deviation (MAD).
            C = 1.4826*median(abs(res-median(res))) ;
            % c. Recompute the weights according to Eq. [13].
            w = 1./(1 + (res/C).^2).^2;
            W = diag(w);
            % d. Perform WLLS fit with new weights
            X = (B'*W*B)\(B'*W)*LogS0;
            % e. Check convergence
            if all(abs(X-X_p) <= con*max(abs(X),abs(X_p))) || cn == Max_inter_in
                crit=0;
            end
    
        end
    
        % Step 3: Transform variables for heteroscedasticity
        fit = B*X;
        LogS0_2 = LogS0./(exp(-fit));
        B2 = B./repmat((exp(-fit)),[1,size(B,2)]);
    
        % Initialization
        crit = 1;
        cn = 0;
        
        % Step 4: Initial LLS fit in * domain
        w = ones(size(S0));
        W = diag(w);
        X = (B'*W*B)\(B'*W)*LogS0;  
    
        % Step 5: Compute a robust estimate for homoscedastic regression using IRLS.
        while crit==1
    
            % Initialization
            cn = cn+1;
            X_p = X;
            
            fit = B2*X; 
            % a. Calculate the residuals e* in the linear domain
            res = (LogS0_2-fit); 
            if all(res==0)
                break;
            end
            % b. Obtain an estimate of the dispersion of the residuals by
            % calculating the median absolute deviation (MAD).
            C = 1.4826*median(abs(res-median(res))) ; 
            % c. Recompute the weights according to Eq. [13].
            w = 1./(1 + (res/C).^2).^2;
            W = diag(w);
    
            % d. Perform WLLS fit with new weights
            X = (B2'*W*B2)\(B2'*W)*LogS0_2;
            % e. Check convergence
            if all(abs(X-X_p) <= con*max(abs(X),abs(X_p))) || cn == Max_inter_in
                crit=0;
            end
    
        end
    
        %  Step 6: Check convergence overall loop
        if all(abs(X-X_p_all) <= con*max(abs(X),abs(X_p_all))) || cn_all == iter_max
                crit_all=0;
        end
    end
    
    % Step 7: Identify and exclude outliers
    fit = B2*X; % In the first iteration, this is the first fit
    res = (LogS0_2-fit);
    
    C = 1.4826*median(abs(res-median(res))) ;
    IND = (abs(res)<= kappa*C);
    
    % Step 8: Final fit
    X = SOLID_EDTI_plugin.SOLID_EDTI_WLLS_WW(S0(IND),B(IND,:), sw(IND)); %SOLID
    
    is_outlier = ~IND;
    
    warning on all
 
end    