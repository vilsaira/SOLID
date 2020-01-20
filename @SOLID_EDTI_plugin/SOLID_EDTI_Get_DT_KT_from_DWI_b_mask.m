function [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = SOLID_EDTI_Get_DT_KT_from_DWI_b_mask(DWI,b,mask,g,par,NrB0,solidWeights) %SOLID

if par.clean_up_PIS==1
    if NrB0>0
        B0s = E_DTI_Clean_up_B0s_2(DWI, mask, NrB0);
        DWI(1:NrB0)=B0s;
        clear B0s;
    end
end

if par.TE==1
    disp('SOLID weights are not applied in OLLS estimation!')
    b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
        4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
        4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
        4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
        6*(g(:,1).^2).*(g(:,2).^2) ...
        6*(g(:,1).^2).*(g(:,3).^2) ...
        6*(g(:,2).^2).*(g(:,3).^2) ...
        12*g(:,2).*g(:,3).*(g(:,1).^2) ...
        12*g(:,1).*g(:,3).*(g(:,2).^2) ...
        12*g(:,1).*g(:,2).*(g(:,3).^2)];
    
    % [W1111 W2222 W3333 W1112 W1113 W2221 W2223 W3331 W3332 W1122 W1133 W2233 W1123 W2213 W3312]
    
    b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
    
    b_final = [ones(length(DWI),1) -b b_kurt];
    
    Vm = sum(mask(:));
    DWI_m = repmat(double(0),[length(DWI) Vm]);
    
    for i=1:length(DWI)
        DWI_m(i,:) = double(DWI{i}(mask));
    end
    
    DWI_m(DWI_m<=0)=1;
    X = b_final\log(DWI_m);
    
    DWIB0 = zeros(size(DWI{1}));
    DWIB0(mask) = exp(X(1,:));
    
    dummy = nan(size(DWI{1}));
    
    for i=1:6
        DT{i} = dummy;
        DT{i}(mask) = X(i+1,:);
    end
    
    for i=7:21;
        KT{i-6} = dummy;
        KT{i-6}(mask) = X(i+1,:)'./((DT{1}(mask)+DT{4}(mask)+DT{6}(mask)).^2);
    end
    
    dwib0 = DWIB0(mask);
    y = max(dwib0(~isinf(dwib0)));
    dwib0(dwib0>y)=y;
    DWIB0(mask) = dwib0;
    
    [outlier, chi_sq, chi_sq_iqr] = E_DTI_Get_outlier_chi_sq_var_KT(DWI,DT,DWIB0,b,par,g,KT);

elseif par.TE==2
    
    
    b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
        4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
        4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
        4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
        6*(g(:,1).^2).*(g(:,2).^2) ...
        6*(g(:,1).^2).*(g(:,3).^2) ...
        6*(g(:,2).^2).*(g(:,3).^2) ...
        12*g(:,2).*g(:,3).*(g(:,1).^2) ...
        12*g(:,1).*g(:,3).*(g(:,2).^2) ...
        12*g(:,1).*g(:,2).*(g(:,3).^2)];
    
    b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
    
    b_final = [ones(length(DWI),1) -b b_kurt];
    
    Vm = sum(mask(:));
    DWI_m = repmat(double(0),[length(DWI) Vm]);
    solid_m = zeros(size(DWI_m)); %SOLID    
    for i=1:length(DWI)
        DWI_m(i,:) = double(DWI{i}(mask));
        tmp = solidWeights(:,:,:,i); %SOLID
        solid_m(i,:) = tmp(mask); %SOLID
    end   
    
    DWI_m(DWI_m<=0)=1;
    X = zeros(22,Vm);
    
    parfor i=1:Vm
        
        X(:,i) = SOLID_EDTI_plugin.SOLID_EDTI_WLLS_WW(DWI_m(:,i),b_final,solid_m(:,i));%SOLID - added argument
        
    end

    DWIB0 = zeros(size(DWI{1}));
    DWIB0(mask) = exp(X(1,:));
    
    dummy = nan(size(DWI{1}));
    
    for i=1:6
        DT{i} = dummy;
        DT{i}(mask) = X(i+1,:);
    end
    
    for i=7:21;
        KT{i-6} = dummy;
        KT{i-6}(mask) = X(i+1,:)'./((DT{1}(mask)+DT{4}(mask)+DT{6}(mask)).^2);
    end
 
    dwib0 = DWIB0(mask);
    y = max(dwib0(~isinf(dwib0)));
    dwib0(dwib0>y)=y;
    DWIB0(mask) = dwib0;
    
    [outlier, chi_sq, chi_sq_iqr] = E_DTI_Get_outlier_chi_sq_var_KT(DWI,DT,DWIB0,b,par,g,KT);
    
elseif par.TE==3 || (par.TE==4 && par.ROBUST_option==2)
    disp('SOLID weights are not applied in non-linear estimation!')
    b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
        4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
        4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
        4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
        6*(g(:,1).^2).*(g(:,2).^2) ...
        6*(g(:,1).^2).*(g(:,3).^2) ...
        6*(g(:,2).^2).*(g(:,3).^2) ...
        12*g(:,2).*g(:,3).*(g(:,1).^2) ...
        12*g(:,1).*g(:,3).*(g(:,2).^2) ...
        12*g(:,1).*g(:,2).*(g(:,3).^2)];
    
    b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
    
    b_final = [ones(length(DWI),1) -b b_kurt];
    
    Vm = sum(mask(:));
    DWI_m = repmat(double(0),[length(DWI) Vm]);
    
    for i=1:length(DWI)
        DWI_m(i,:) = double(DWI{i}(mask));
    end
    
    DWI_m(DWI_m<=0)=1;
    X = zeros(22,Vm);
      
    options = statset('nlinfit');
    if par.TE==4
        options.Robust = 'on';
    end
    
    parfor i=1:Vm
                        
        X(:,i) = E_DTI_quick_par(options,DWI_m(:,i),b_final);
        
    end
    
    DWIB0 = zeros(size(DWI{1}));
    DWIB0(mask) = exp(X(1,:));
    
    
    dummy = nan(size(DWI{1}));
    
    for i=1:6
        DT{i} = dummy;
        DT{i}(mask) = X(i+1,:);
    end
    
    for i=7:21;
        KT{i-6} = dummy;
        KT{i-6}(mask) = X(i+1,:)'./((DT{1}(mask)+DT{4}(mask)+DT{6}(mask)).^2);
    end
    
    dwib0 = DWIB0(mask);
    y = max(dwib0(~isinf(dwib0)));
    dwib0(dwib0>y)=y;
    DWIB0(mask) = dwib0;
    
    [outlier, chi_sq, chi_sq_iqr] = E_DTI_Get_outlier_chi_sq_var_KT(DWI,DT,DWIB0,b,par,g,KT);
    
elseif par.TE==4 && par.ROBUST_option==1
    
    con = par.RE.rel_convergence;
    iter_max = par.RE.max_iter;
    kappa = par.RE.kappa;
        
    b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
        4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
        4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
        4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
        6*(g(:,1).^2).*(g(:,2).^2) ...
        6*(g(:,1).^2).*(g(:,3).^2) ...
        6*(g(:,2).^2).*(g(:,3).^2) ...
        12*g(:,2).*g(:,3).*(g(:,1).^2) ...
        12*g(:,1).*g(:,3).*(g(:,2).^2) ...
        12*g(:,1).*g(:,2).*(g(:,3).^2)];
    
    b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;
    
    b_final = [ones(length(DWI),1) -b b_kurt];
    
    Vm = sum(mask(:));
    DWI_m = repmat(double(0),[length(DWI) Vm]);
    solid_m = zeros(size(DWI_m)); %SOLID    
    for i=1:length(DWI)
        DWI_m(i,:) = double(DWI{i}(mask));
        tmp = solidWeights(:,:,:,i); %SOLID
        solid_m(i,:) = tmp(mask); %SOLID
    end  
    
    DWI_m(DWI_m<=0)=1;
    X = zeros(22,Vm);
    
    outlier_m = false([length(DWI) Vm]);
    outlier = false([size(mask) length(DWI)]);    
    
    parfor i=1:Vm
        
        [X(:,i), outlier_m(:,i)] = SOLID_EDTI_plugin.SOLID_EDTI_robust_Linear_fit_new(DWI_m(:,i),b_final,kappa,iter_max,con,solid_m(:,i));%SOLID - added argument
        
    end
    
    DWIB0 = zeros(size(DWI{1}));
    DWIB0(mask) = exp(X(1,:));
    
    dummy = nan(size(DWI{1}));
    
    for i=1:6
        DT{i} = dummy;
        DT{i}(mask) = X(i+1,:);
    end
    
    for i=7:21;
        KT{i-6} = dummy;
        KT{i-6}(mask) = X(i+1,:)'./((DT{1}(mask)+DT{4}(mask)+DT{6}(mask)).^2);
    end
    
    dummy = false(size(DWI{1}));
    
    for i=1:length(DWI)
        dummy(mask) = outlier_m(i,:);
        outlier(:,:,:,i) = dummy;
    end
    
    dwib0 = DWIB0(mask);
    y = max(dwib0(~isinf(dwib0)));
    dwib0(dwib0>y)=y;
    DWIB0(mask) = dwib0;
    
    [~, chi_sq, chi_sq_iqr] = E_DTI_Get_outlier_chi_sq_var_KT(DWI,DT,DWIB0,b,par,g,KT);
 
end

