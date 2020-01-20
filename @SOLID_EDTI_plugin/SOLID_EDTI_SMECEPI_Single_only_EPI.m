function SOLID_EDTI_SMECEPI_Single_only_EPI(f_in,for_trafo)

    [FOLD,FFO] = fileparts(f_in);

    par = for_trafo.par;
    
    [suc, FN_nii] = E_DTI_SMECEPI_check_data_stuff(f_in,par);
    
    if suc==0
        return;
    end
    
    fn_fixed = [for_trafo.dir_temp filesep 'Trafo_fixed.nii'];
    fn_fixed_mask = [for_trafo.dir_temp filesep 'Trafo_fixed_mask.nii'];
    [Fixed,VDims] = E_DTI_read_nifti_file(FN_nii);
    VDims_E = VDims;
    Fixed = single(Fixed);
    Fixed = 10000*(Fixed/max(Fixed(:)));
    E_DTI_write_nifti_file(Fixed,VDims,fn_fixed);
    
    mask = single(Fixed>0);
    if par.EPI.use_f_mask~=1
        mask(:)=1;
    end
    E_DTI_write_nifti_file(mask,VDims,fn_fixed_mask);
    
    fn_moving = [for_trafo.dir_temp filesep 'Trafo_moving.nii'];
    
    if ~isempty(for_trafo.f_out)
        LoadF = for_trafo.f_out;
    else
        LoadF = f_in;
    end
    
    if par.R2D.contrast==1
        load(LoadF,'DWI','VDims')
        DWI = E_DTI_DWI_mat2cell(DWI);
        Moving = single(DWI{1});
        clear DWI;
    elseif par.R2D.contrast==2
        load(LoadF,'FA','VDims')
        FA(isnan(FA))=0;
        FA = FA/sqrt(3);
        Moving = single(FA);
        clear FA;
    elseif par.R2D.contrast==3
        load(LoadF,'DWI','VDims','NrB0')
        DWI = E_DTI_DWI_mat2cell(DWI);
        Moving = single(E_DTI_mean_DWI(DWI, NrB0));
        clear DWI;
    end
    
    yp = prctile(Moving(:),99);
    Moving(Moving>yp)=yp;
    
    Moving = 10000*(Moving/max(Moving(:)));
    
    E_DTI_write_nifti_file(Moving,VDims,fn_moving);
    
    
    fn_par_trafo_rigid = [for_trafo.dir_temp filesep 'Par_file_Trafo_Rigid.txt'];
    
    par_EPI.Hist_bin = par.EPI.Hist_bin;
    par_EPI.Num_Resol = par.EPI.Num_Resol;
    par_EPI.Num_iter = par.EPI.Num_iter;
    par_EPI.Num_samp = par.EPI.Num_samp;
    par_EPI.Interpol = par.Interpol;
    par_EPI.Grid_Spacing = par.EPI.Grid_Spacing([2 1 3]);
    par_EPI.Par_FN = fn_par_trafo_rigid;
    par_EPI.Deriv_Scales = par.EPI.Deriv_Scales([2 1 3]);
    
    suc = E_DTI_Par_Trafo_DTI_rigid(par_EPI);
    if suc==0
        E_DTI_remove_temp_f(for_trafo.dir_temp);
        return;
    end
    
    if par.R2D.type~=3
        
        if ispc
            
            sys_c = ['"' par.E_path '"' ' -f ' '"' fn_fixed '"'...
                ' -fMask ' '"' fn_fixed_mask '"' ...
                ' -m ' '"' fn_moving '"'...
                ' -out ' '"' for_trafo.dir_temp '"' ...
                ' -p ' '"' fn_par_trafo_rigid '"'];
            [r,st]=system(sys_c);
            
        else
            
            [dir_e,el] = fileparts(par.E_path);
            sys_c = ['./' el ' -f ' fn_fixed ...
                ' -fMask ' fn_fixed_mask ...
                ' -m ' fn_moving ...
                ' -out ' for_trafo.dir_temp ...
                ' -p ' fn_par_trafo_rigid];
            
            tedi = cd(dir_e);
            [r,st]=system(sys_c);
            cd(tedi)
        end
        
    else
        
        fn_par_trafo_non_rigid = [for_trafo.dir_temp filesep 'Par_file_Trafo_Non_Rigid.txt'];
        
        par_EPI.Par_FN = fn_par_trafo_non_rigid;
        suc = E_DTI_Par_Trafo_DTI_non_rigid(par_EPI);
        
        if suc==0
            E_DTI_remove_temp_f(for_trafo.dir_temp);
            return;
        end
        
        if ispc
            
            sys_c = ['"' par.E_path '"' ' -f ' '"' fn_fixed '"'...
                ' -fMask ' '"' fn_fixed_mask '"' ...
                ' -m ' '"' fn_moving '"'...
                ' -out ' '"' for_trafo.dir_temp '"' ...
                ' -p ' '"' fn_par_trafo_rigid '"' ...
                ' -p ' '"' fn_par_trafo_non_rigid '"'];
            [r,st]=system(sys_c);
            
        else
            
            [dir_e,el] = fileparts(par.E_path);
            sys_c = ['./' el ' -f ' fn_fixed ...
                ' -fMask ' fn_fixed_mask ...
                ' -m ' fn_moving ...
                ' -out ' for_trafo.dir_temp ...
                ' -p ' fn_par_trafo_rigid ...
                ' -p ' fn_par_trafo_non_rigid];
            
            tedi = cd(dir_e);
            [r,st]=system(sys_c);
            cd(tedi)
        end
        
        
    end
    
    
    if r~=0
        disp('Error(s) during motion/distortion correction:')
        disp(' ')
        disp(st);
        disp(' ')
        disp('See the forum for a solution...')
        E_DTI_remove_temp_f(for_trafo.dir_temp);
        return;
    end
    
    Trafo_rig_result = [for_trafo.dir_temp filesep 'TransformParameters.0.txt'];
    % Rig = E_DTI_get_Trafo_par_EPI_rigid_rotation_components(Trafo_rig_result);
    
    Q = textread(Trafo_rig_result,'%s');
    Rig{1} = str2num(Q{7});
    Rig{2} = str2num(Q{6});
    Rig{3} = str2num(Q{8});
    
    
    At = textread(Trafo_rig_result,'%s','delimiter','\n','bufsize',2^18);
    
    if par.R2D.type==3
        
        Trafo_nonrig_result = [for_trafo.dir_temp filesep 'TransformParameters.1.txt'];
        At_nr = textread(Trafo_nonrig_result,'%s','delimiter','\n','bufsize',2^18);
        
    end
    
    
    
    for i=1:length(for_trafo.trafo_names)
        
        [Fol,dummy] = fileparts(for_trafo.trafo_names{i});
        
        if par.R2D.type~=3
            FinTrafoN{i} = [Fol filesep 'Final_Trafo.txt'];
            At{4} = ['(InitialTransformParametersFileName "' for_trafo.trafo_names{i} '")'];
            suc = E_DTI_SMECEPI_write_tra_file(At,FinTrafoN{i});
            if suc==0
                E_DTI_remove_temp_f(for_trafo.dir_temp);
                return;
            end
        else
            FinTrafoN_R{i} = [Fol filesep 'Final_Trafo_Rigid_step.txt'];
            FinTrafoN{i} = [Fol filesep 'Final_Trafo.txt'];
            At{4} = ['(InitialTransformParametersFileName "' for_trafo.trafo_names{i} '")'];
            suc = E_DTI_SMECEPI_write_tra_file(At,FinTrafoN_R{i});
            if suc==0
                E_DTI_remove_temp_f(for_trafo.dir_temp);
                return;
            end
            At_nr{4} = ['(InitialTransformParametersFileName "' FinTrafoN_R{i} '")'];
            suc = E_DTI_SMECEPI_write_tra_file(At_nr,FinTrafoN{i});
            if suc==0
                E_DTI_remove_temp_f(for_trafo.dir_temp);
                return;
            end
        end
        
    end
    
    
    load(LoadF,'b','bval','info','g','NrB0')
    
    [b, g] = E_DTI_reorient_grad_and_b_matrix_rigid_rotation(b,g,Rig);
    
    if isnan(bval)
        diff_model=2;
    else
        diff_model=1;
    end
    
    [dummy, g] =  E_DTI_Get_Gradients_and_Bvalue(b, NrB0, diff_model);
    
    parfor i=1:length(FinTrafoN)
        E_DTI_do_the_SMEC_trafo_step(for_trafo.f_moving{i},for_trafo.dir_temp_i{i},FinTrafoN{i},par)
                    
        movefile([for_trafo.dir_temp_i{i} filesep 'result.nii'], [for_trafo.dir_temp_i{i} filesep 'result_DWI.nii']); %SOLID avoid writing over DWI registration result
        E_DTI_do_the_SMEC_trafo_step(for_trafo.solid_moving{i},for_trafo.dir_temp_i{i},FinTrafoN{i},par); %SOLID
        movefile([for_trafo.dir_temp_i{i} filesep 'result.nii'], [for_trafo.dir_temp_i{i} filesep 'SOLID.nii']); %SOLID
        modZ4D(:,:,:,i) = E_DTI_read_nifti_file(for_trafo.solid_names{i}); %SOLID
        movefile([for_trafo.dir_temp_i{i} filesep 'result_DWI.nii'], [for_trafo.dir_temp_i{i} filesep 'result.nii']); %SOLID avoid writing over DWI registration result        
    end     
    result_names = for_trafo.result_names;

    parfor i=1:length(FinTrafoN)
        DWI{i} = E_DTI_read_nifti_file(result_names{i});
    end

    
    for i=1:length(DWI)
        if isempty(DWI{i})
            disp('Errors encountered...')
            E_DTI_remove_temp_f(for_trafo.dir_temp);
            suc = 0;
            return;
        end
    end

    for i=1:length(DWI)
        DWI{i} = single(DWI{i});
        DWI{i}(DWI{i}==-1000)=nan;
    end


    dummy = false(size(DWI{1}));
    for i=1:length(DWI)
        dummy = or(dummy,isnan(DWI{i}));
        DWI{i}(isnan(DWI{i}))=0;
        DWI{i}(DWI{i}<0)=0;
    end

    if isempty(par.cust_mask.TS)
        mask = E_DTI_Create_Mask_From_DWI_enhanced(DWI,NrB0,par.mask_P.TS.NDWI,par.mask_P.TS.DWI,par.mask_P.TS.mfs);
    else
        fn_cm = [FOLD filesep FFO par.cust_mask.TS];
        [mask, VDims, suc] = E_DTI_read_nifti_file(fn_cm);
        if suc==0
            E_DTI_remove_temp_f(for_trafo.dir_temp);
            return;
        end
        mask = mask>0;
    end

    mask(dummy)=0;

    par_temp = par;
    par_temp.TE = par.TE.TS;

    % Scale modZ scores from [3.5 10] to [0 1] to obtain solidWeights
    solidWeights = modZ4D;
    solidWeights(solidWeights < par.SOLID.thrL) = par.SOLID.thrL;
    solidWeights(solidWeights > par.SOLID.thrU) = par.SOLID.thrU;
    solidWeights = (solidWeights - par.SOLID.thrL) ./ (par.SOLID.thrU - par.SOLID.thrL); 
    solidWeights = 1-solidWeights;

    if diff_model==1 
        [DT, DWIB0, outlier, chi_sq, chi_sq_iqr] = SOLID_EDTI_plugin.SOLID_EDTI_Get_DT_from_DWI_b_mask(DWI,b,mask,par_temp,NrB0, solidWeights); 
    elseif diff_model==2
        [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = SOLID_EDTI_plugin.SOLID_EDTI_Get_DT_KT_from_DWI_b_mask_with_constraints(DWI,b,mask,g,par_temp,NrB0,VDims, solidWeights);
    end

    [FEFA, FA, FE, SE, eigval] = E_DTI_eigensystem_analytic(DT);

    g = g(NrB0+1:end,:);
    MDims = size(mask);
    VDims = VDims_E;

    f_out = [par.out_folder filesep FFO par.suff.TS];
    par.Rig = Rig;

    max_DWI = 0;
    for i=1:length(DWI)
        max_DWI = max(max_DWI,max(DWI{i}(:)));
    end

    if max_DWI<=intmax('int16')
        for i=1:length(DWI)
            DWI{i} = round(DWI{i});
            DWI{i} = int16(DWI{i});
        end
    elseif max_DWI<=intmax('uint16')
        for i=1:length(DWI)
            DWI{i} = round(DWI{i});
            DWI{i}(DWI{i}<0)=0;
            DWI{i} = uint16(DWI{i});
        end
    end


    try
        save(f_out,'DWI','VDims','b','bval','g','info','FEFA','NrB0','MDims',...
            'FA','FE','SE','eigval','DT','outlier','DWIB0','chi_sq','chi_sq_iqr','par');
    catch me
        E_DTI_remove_temp_f(for_trafo.dir_temp);
        disp(me.message)
        return;
    end

    if par.SOLID.saveResults
        try
            save(f_out, 'solidWeights', '-append');
        catch
        end
    end

    if diff_model==2
        save(f_out,'KT','-append')
    end

    E_DTI_remove_temp_f(for_trafo.dir_temp);
    
end