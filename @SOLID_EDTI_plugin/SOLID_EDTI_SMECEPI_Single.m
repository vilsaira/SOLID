function SOLID_EDTI_SMECEPI_Single(f_in,par)

    [FOL_FI,F,ext] = fileparts(f_in);
    
    tic;
    [suc, for_trafo] = SOLID_EDTI_plugin.SOLID_EDTI_SMECEPI_Single_only_EC(f_in,par);
    
    if suc==0
        return;
    end
    
    if par.R2D.type~=0
    
        SOLID_EDTI_plugin.SOLID_EDTI_SMECEPI_Single_only_EPI(f_in,for_trafo);
        
    end
    
    t=toc;
    if t<3600
        m = t/60;
        disp(['Computation time for ''' [F ext] ''' was ' num2str(m) ' minutes.'])
    elseif t>3600 && t<(3600*24)
        h = t/3600;
        disp(['Computation time for ''' [F ext] ''' was ' num2str(h) ' hours.'])
    else
        d = t/(3600*24);
        disp(['Computation time for ''' [F ext] ''' was ' num2str(d) ' days.'])
    end
    
end