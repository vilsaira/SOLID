function [DT, DWIB0, outlier, chi_sq, chi_sq_iqr] = SOLID_EDTI_Get_DT_from_DWI_b_mask(DWI,b,mask,par,NrB0, solidWeights) %SOLID

    if par.clean_up_PIS==1
        if NrB0>0
            B0s = E_DTI_Clean_up_B0s_2(DWI, mask, NrB0);
            DWI(1:NrB0)=B0s;
            clear B0s;
        end
    end
    
    if par.TE==1
        disp('SOLID weights are not applied in OLLS estimation!')
        Vm = sum(mask(:));
        DWI_m = repmat(double(0),[length(DWI) Vm]);
        
        for i=1:length(DWI)
            DWI_m(i,:) = double(DWI{i}(mask));
        end
        
        DWI_m(DWI_m<=0)=1;
        b2 = [ones(size(b,1),1) -b];
        X = b2\log(DWI_m);
        
        DWIB0 = nan(size(DWI{1}));
        DWIB0(mask) = exp(X(1,:));
        
        dummy = nan(size(DWI{1}));
        
        for i=1:6
            DT{i} = dummy;
            DT{i}(mask) = X(i+1,:);
        end
        
        dwib0 = DWIB0(mask);
        y = max(dwib0(~isinf(dwib0)));
        dwib0(dwib0>y)=y;
        DWIB0(mask) = dwib0;
    
        [outlier, chi_sq, chi_sq_iqr] = E_DTI_Get_outlier_chi_sq_var(DWI,DT,DWIB0,b,par);
        
    elseif par.TE==2
        
        
        Vm = sum(mask(:));
        DWI_m = repmat(double(0),[length(DWI) Vm]);
        solid_m = zeros(size(DWI_m)); %SOLID
        for i=1:length(DWI)
            DWI_m(i,:) = double(DWI{i}(mask));
            tmp = solidWeights(:,:,:,i); %SOLID
            solid_m(i,:) = tmp(mask); %SOLID
        end       
        
        DWI_m(DWI_m<=0)=1;
        b2 = [ones(size(b,1),1) -b];
        
        X = zeros(7,Vm);
        
        parfor i=1:Vm
            
            dwi = DWI_m(:,i);
            sw = solid_m(:,i); %SOLID
            X(:,i) = SOLID_EDTI_plugin.SOLID_EDTI_WLLS_WW(dwi,b2,sw); %SOLID
            
        end
        
        dummy = nan(size(DWI{1}));
        
        for i=1:6
            DT{i} = dummy;
            DT{i}(mask) = X(i+1,:);
        end
        
        DWIB0 = nan(size(DWI{1}));
        DWIB0(mask) = exp(X(1,:));
        
        dwib0 = DWIB0(mask);
        try
            y = max(dwib0(~isinf(dwib0)));
            dwib0(dwib0>y)=y;
        catch
        end
        DWIB0(mask) = dwib0;
    
        [outlier, chi_sq, chi_sq_iqr] = E_DTI_Get_outlier_chi_sq_var(DWI,DT,DWIB0,b,par);
        
    elseif par.TE==3 || (par.TE==4 && par.ROBUST_option==2)
        
        disp('SOLID weights are not applied in non-linear estimation!')
        Vm = sum(mask(:));
        DWI_m = repmat(double(0),[length(DWI) Vm]);
        
        for i=1:length(DWI)
            DWI_m(i,:) = double(DWI{i}(mask));
        end
        
        DWI_m(DWI_m<=0)=1;
        b2 = [ones(size(b,1),1) -b];
        
        X = zeros(7,Vm);
            
        options = statset('nlinfit');
        if par.TE==4
            options.Robust = 'on';
        end
        
        parfor i=1:Vm
            
            X(:,i) = E_DTI_quick_par(options,DWI_m(:,i),b2);
            
        end
        
        dummy = nan(size(DWI{1}));
        
        for i=1:6
            DT{i} = dummy;
            DT{i}(mask) = X(i+1,:);
        end
        
        DWIB0 = nan(size(DWI{1}));
        DWIB0(mask) = exp(X(1,:));
        
        dwib0 = DWIB0(mask);
        y = max(dwib0(~isinf(dwib0)));
        dwib0(dwib0>y)=y;
        DWIB0(mask) = dwib0;
    
        [outlier, chi_sq, chi_sq_iqr] = E_DTI_Get_outlier_chi_sq_var(DWI,DT,DWIB0,b,par);
            
    elseif par.TE==4 && par.ROBUST_option==1
        
        con = par.RE.rel_convergence;
        iter_max = par.RE.max_iter;
        kappa = par.RE.kappa;
        
        Vm = sum(mask(:));
        DWI_m = repmat(double(0),[length(DWI) Vm]);
        solid_m = zeros(size(DWI_m)); %SOLID        
        for i=1:length(DWI)
            DWI_m(i,:) = double(DWI{i}(mask));
            tmp = solidWeights(:,:,:,i); %SOLID
            solid_m(i,:) = tmp(mask); %SOLID
        end  
        
        DWI_m(DWI_m<=0)=1;
        b2 = [ones(size(b,1),1) -b];
        
        X = zeros(7,Vm);
        outlier_m = false([length(DWI) Vm]);
        outlier = false([size(mask) length(DWI)]);     
        
        parfor i=1:Vm    
            [X(:,i), outlier_m(:,i)] = SOLID_EDTI_plugin.SOLID_EDTI_robust_Linear_fit_new(DWI_m(:,i),b2,kappa,iter_max,con,solid_m(:,i));%SOLID - added argument
        end
        
        dummy = nan(size(DWI{1}));
        
        for i=1:6
            DT{i} = dummy;
            DT{i}(mask) = X(i+1,:);
        end
        
        dummy = false(size(DWI{1}));
        
        for i=1:length(DWI)
            dummy(mask) = outlier_m(i,:);
            outlier(:,:,:,i) = dummy;
        end
        
        DWIB0 = nan(size(DWI{1}));
        DWIB0(mask) = exp(X(1,:));
        
        dwib0 = DWIB0(mask);
        try
            y = max(dwib0(~isinf(dwib0)));
            dwib0(dwib0>y)=y;
        catch
        end
        DWIB0(mask) = dwib0;
    
        [dummy_stuff, chi_sq, chi_sq_iqr] = E_DTI_Get_outlier_chi_sq_var(DWI,DT,DWIB0,b,par);    
        
    end
    
end
    
    
    