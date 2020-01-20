function [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = SOLID_EDTI_Get_DT_KT_from_DWI_b_mask_with_constraints(DWI,b,mask,g,par,NrB0,VDims,solidWeights) %SOLID

    [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = SOLID_EDTI_plugin.SOLID_EDTI_Get_DT_KT_from_DWI_b_mask(DWI,b,mask,g,par,NrB0,solidWeights); %SOLID
    
    if par.DKI_constraints.do_it==1
        
        mask_Fin = mask;
        
        constr1 = par.DKI_constraints.constr1;
        constr2 = par.DKI_constraints.constr2;
        
        FWHM = max(VDims);
        sigma_sm = FWHM/(2*sqrt(2*log(2)));
        
        mask = E_DTI_constrain_KT_mask(KT,constr1,constr2);
        
        crit = sum(single(mask(:)))>0;
        cn = 0;
        
        while crit
            
            cn=cn+1;
            
            for t=1:length(DWI)
                DWI{t} = E_DTI_smooth3_ani_voxel_size(single(DWI{t}), 'gaussian', [3 3 3], sigma_sm, VDims);
            end
            
            [DT_, DWIB0_, KT_, outlier_, chi_sq_, chi_sq_iqr_] = SOLID_EDTI_plugin.SOLID_EDTI_Get_DT_KT_from_DWI_b_mask(DWI,b,mask,g,par,NrB0,solidWeights); %SOLID
            
            for i=1:length(KT)
                KT{i}(mask) = KT_{i}(mask);
            end
    
            for i=1:length(DT)
                DT{i}(mask) = DT_{i}(mask);
            end
            
            DWIB0(mask) = DWIB0_(mask);
            chi_sq(mask) = chi_sq_(mask);
            chi_sq_iqr(mask) = chi_sq_iqr_(mask);
            outlier(repmat(mask,[1 1 1 size(outlier,4)])) = outlier_(repmat(mask,[1 1 1 size(outlier,4)]));
            
            mask = E_DTI_constrain_KT_mask(KT,constr1,constr2);
            
            
            if cn==10 || sum(single(mask(:)))==0
                crit=0;
            end
            
        end
        
    end

end