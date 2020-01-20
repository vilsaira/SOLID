 classdef SOLID_EDTI_plugin < handle
% This class is a fork of ExploreDTI 486 SMECEPI algorithm.
% 1) To be used as a commandline option to SMECEPI EDTI *.mat files
% 2) To be used as a compliment to ExploreDTI GUI

properties (Access = public)
    EDTI = []
    solid = [];
    SOLID_GUI = [];
    thrL = 3.5;
    thrU = 10.0;
    metric = 'Variance';
    useMask = true;
    saveResults = false;
end

methods (Access = private)
    plugin = SOLID_EDTI_plugin_initGUI(obj, source, event);    
    plugin = SOLID_EDTI_plugin_QA(obj, source, event);
    SOLID_EDTI_plugin_startSMECEPI(obj, source, event);        
    SOLID_EDTI_plugin_settings(obj, source, event);
    SOLID_EDTI_plugin_settingsExport(obj, source, event);
end

methods (Static, Access = public)
    SOLID_EDTI_SMECEPI_Single(f_in,par);
    [suc, for_trafo] = SOLID_EDTI_SMECEPI_Single_only_EC(f_in,par);
    SOLID_EDTI_SMECEPI_Single_only_EPI(f_in,for_trafo);   
    [DT, DWIB0, outlier, chi_sq, chi_sq_iqr] = SOLID_EDTI_Get_DT_from_DWI_b_mask(DWI,b,mask,par_temp,NrB0, solidWeights);  
    [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = SOLID_EDTI_Get_DT_KT_from_DWI_b_mask_with_constraints(DWI,b,mask,g,par,NrB0,VDims,solidWeights);        
    [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = SOLID_EDTI_Get_DT_KT_from_DWI_b_mask(DWI,b,mask,g,par,NrB0,solidWeights)
    SOLID_EDTI_SMECEPI_Main(Par_file);  
    X = SOLID_EDTI_WLLS_WW(S0,B,sw);
    [X, is_outlier] = SOLID_EDTI_robust_Linear_fit_new(S0,B,kappa,iter_max,con, sw);
end

methods (Access = public)
    function plugin = SOLID_EDTI_plugin(varargin)
                
        plugin = plugin.SOLID_EDTI_plugin_initGUI();  
        
        if nargout == 0
            warning('on', 'all');
            clear plugin;
        end              
        
    end
end

end