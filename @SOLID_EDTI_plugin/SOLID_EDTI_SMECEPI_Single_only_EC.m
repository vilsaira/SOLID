function [suc, for_trafo] = SOLID_EDTI_SMECEPI_Single_only_EC(f_in,par)

    for_trafo = struct;

    [FOL_FI,F] = fileparts(f_in);
    
    dir_temp = [par.temp_folder filesep 'Temp_' F];
    
    for_trafo.dir_temp = dir_temp;
    
    warning off all
    suc = mkdir(dir_temp);
    warning on all
    
    if suc==0
        disp('Error: could not write files/folders!')
        disp('Change permission properties of temp folder:')
        disp(par.temp_folder)
        return;
    end
    
    if ispc
        [suc, mess] = fileattrib(dir_temp,'+w -a -s','','s');
    else
        [suc, mess] = fileattrib(dir_temp,'+w +x','a','s');
    end
    
    if suc==0
        disp(['Problem with folder: ' dir_temp])
        disp(['Error message: ' mess])
        E_DTI_remove_temp_f(dir_temp);
        return;
    end
    
    if par.DOF~=1
        par.Par_SMEC = [dir_temp filesep 'par_file_SMEC.txt'];
        
        suc = E_DTI_Par_MDC_aff_DTI(par);
        if suc==0
            E_DTI_remove_temp_f(dir_temp);
            return;
        end
    end
        
    load(f_in,'VDims','DWI','FA', 'b');    
    runSolid = true;
    niifPath = par.SOLID.modZ_f_in;
    if isempty(niifPath)
        niifPath = f_in;
    else
        [modZ4D, ~] = E_DTI_load_nii(niifPath);
        runSolid = false;
    end
    bvalPath = f_in;
    bvecPath = f_in;
    maskPath = par.SOLID.mask_f_in;
    if isempty(maskPath)
        maskPath = [];
    else
        [solidMask, ~] = E_DTI_load_nii(maskPath);
    end
    if ~par.SOLID.useMask
        maskPath = 'no'; % Masking is not used
    end
    maskTune = par.mask_P.NS.NDWI;
    maskFilter = par.mask_P.NS.mfs;
    metric = par.SOLID.metric;
    thrL = par.SOLID.thrL;
    thrU = par.SOLID.thrU;   
    saveResults = par.SOLID.saveResults;
    
    if runSolid
        solid = SOLID('in', niifPath,...
                    'bval', bvalPath,...
                    'bvec', bvecPath,...
                    'mask', maskPath,...
                    'maskTune', maskTune,...
                    'maskFilter', maskFilter,...
                    'metric', metric,...
                    'thrL', thrL,...
                    'thrU', thrU,...
                    'save', saveResults);
        modZ4D = zeros(size(solid.DWI), 'single');
        for k = 1:size(solid.modZ,1)
            for l = 1:size(solid.modZ,2)
                modZ4D(:,:,k,l) = repmat(solid.modZ(k,l), [size(modZ4D,1), size(modZ4D,2), 1,1]);
            end
        end        
    end
    
    
    if par.DOF~=1
        fixed_SMEC_mask = single(~isnan(FA));
        
        if par.use_f_mask == 1
            a=3;
            se = zeros(a,a,a);
            se((a+1)/2,(a+1)/2,(a+1)/2)=1;
            se = smooth3(se,'gaussian',[a a a],1);
            se = se>=se((a+1)/2,(a+1)/2,end);
            fixed_SMEC_mask = single(imdilate(fixed_SMEC_mask,se)>0);
        else
            fixed_SMEC_mask(:)=1;
        end
    end

    
    DWI = E_DTI_DWI_mat2cell(DWI);

    for i=1:length(DWI)
        DWI{i} = single(DWI{i});
    end

    if par.DOF~=1
        par.fixed_SMEC_fn = [dir_temp filesep 'Fixed_SMEC.nii'];
        par.fixed_SMEC_fn_mask = [dir_temp filesep 'Fixed_SMEC_mask.nii'];
        
        fixed_SMEC = single(DWI{1});
        
        if par.Regul==1
            fixed_SMEC = E_DTI_Anisotr_Gauss_smoothing(fixed_SMEC,VDims,2*max(VDims),[3 3 3]);
        end
        
        E_DTI_write_nifti_file(fixed_SMEC,VDims,par.fixed_SMEC_fn)
        E_DTI_write_nifti_file(fixed_SMEC_mask,VDims,par.fixed_SMEC_fn_mask)
        
        if exist(par.fixed_SMEC_fn,'file')~=2
            suc = 0;
            E_DTI_remove_temp_f(dir_temp);
            return;
        end
    end

    Le = length(DWI);    
    f_solid = cell(Le,1); %SOLID
    result_solid = cell(Le,1); %SOLID    

    
    f_moving = cell(Le,1);
    result = cell(Le,1);
    if par.DOF~=1
        DM_info = cell(Le,1);
    end
    VoxS = cell(Le,1);
    dir_temp_i = cell(Le,1);
    trafo_names = cell(Le,1);

    for i=1:Le
        
        dir_temp_i{i} = [dir_temp filesep 'Temp_' num2str(i)];
        warning off all
        suc = mkdir(dir_temp_i{i});
        warning on all
        
        if suc==0
            disp('Error: could not write folder:')
            disp(dir_temp_i{i})
            E_DTI_remove_temp_f(dir_temp);
            return;
        end
        
        VoxS{i} = VDims;
        trafo_names{i} = [dir_temp_i{i} filesep 'TransformParameters.0.txt'];
        result{i} = [dir_temp_i{i} filesep 'result.nii'];
        f_moving{i} = [dir_temp_i{i} filesep 'moving_' num2str(i) '.nii'];

        f_solid{i} = [dir_temp_i{i} filesep 'SOLID_' num2str(i) '.nii']; %SOLID
        result_solid{i} = [dir_temp_i{i} filesep 'SOLID.nii']; %SOLID
        E_DTI_write_nifti_file(modZ4D(:,:,:,i), VoxS{i}, f_solid{i}); %SOLID
    end

    parfor i=1:Le
        E_DTI_output_TR_files(DWI{i},VoxS{i},f_moving{i},par)
    end

    if par.DOF~=1
    
        suc = E_DTI_SMECEPI_make_unity_transf(trafo_names{1},VoxS{1},size(DWI{1}));
        if suc==0
            return;
        end
        
        parfor i=2:Le
            E_DTI_do_the_SMEC_reg_step(f_moving{i},dir_temp_i{i},par)
        end
        
        if par.Regul==1
            parfor i=2:Le
                E_DTI_write_nifti_file(DWI{i},VoxS{i},f_moving{i});
            end
        end
        
        parfor i=2:Le
            E_DTI_do_the_SMEC_trafo_step(f_moving{i},dir_temp_i{i},trafo_names{i},par)
        end
        
        parfor i=2:Le
            DWI{i} = E_DTI_read_nifti_file(result{i});
        end

        parfor i=2:Le  %SOLID
            movefile([dir_temp_i{i} filesep 'result.nii'], [dir_temp_i{i} filesep 'result_DWI.nii']); %SOLID avoid writing over DWI registration result            
            E_DTI_do_the_SMEC_trafo_step(f_solid{i}, dir_temp_i{i}, trafo_names{i}, par); %SOLID
            movefile([dir_temp_i{i} filesep 'result.nii'], [dir_temp_i{i} filesep 'SOLID.nii']); %SOLID
            modZ4D(:,:,:,i) = E_DTI_read_nifti_file(result_solid{i}); %SOLID
            movefile([dir_temp_i{i} filesep 'result_DWI.nii'], [dir_temp_i{i} filesep 'result.nii']); %SOLID avoid writing over DWI registration result
        end  %SOLID
        
        for i=1:length(DWI)
            if isempty(DWI{i})
                disp('Errors encountered...')
                E_DTI_remove_temp_f(dir_temp);
                suc = 0;
                return;
            end
        end

        for i=1:length(DWI)
            DWI{i} = single(DWI{i});
            DWI{i}(DWI{i}==-1000)=nan;
        end

        for i=1:Le
        
            A = textread(trafo_names{i},'%s');
            
            DM_info{i}{1} = ...
                [str2num(A{7})*(180/pi) str2num(A{6})*(180/pi) str2num(A{8})*(180/pi);
                str2num(A{16}) str2num(A{15}) str2num(A{17}(1:end-1));
                str2num(A{13}) str2num(A{12}) str2num(A{14});
                str2num(A{10})  str2num(A{9}) str2num(A{11})];
            
            DM_info{i}{2} = E_DTI_Tra_Par_2_Tra_Mat(DM_info{i}{1});
            
            
            DWI{i} = DWI{i}*det(DM_info{i}{2}(1:3,1:3));
            
        end
        
        load(f_in,'b','NrB0','bval','info')
        
        b_old = b;
        
        for i=1:length(DWI)
            World_Trafo = E_DTI_Tra_Par_2_Tra_Mat(DM_info{i}{1});
            World_T = World_Trafo(1:3,1:3);
            R = ((World_T*World_T')^(-1/2))*World_T;
            B = [b_old(i,1) b_old(i,2)/2 b_old(i,3)/2;...
                b_old(i,2)/2 b_old(i,4) b_old(i,5)/2;...
                b_old(i,3)/2 b_old(i,5)/2 b_old(i,6)];
            B_rot = R*B*R';
            b(i,:) = [B_rot(1,1) 2*B_rot(1,2) 2*B_rot(1,3)...
                B_rot(2,2) 2*B_rot(2,3) B_rot(3,3)];
        end
        
        if isnan(bval)
            diff_model=2;
        else
            diff_model=1;
        end
        
        [dummy, g] =  E_DTI_Get_Gradients_and_Bvalue(b, NrB0, diff_model);
        
        dummy = false(size(DWI{1}));
        for i=1:length(DWI)
            dummy = or(dummy,isnan(DWI{i}));
            DWI{i}(isnan(DWI{i}))=0;
            DWI{i}(DWI{i}<0)=0;
        end
        
        if isempty(par.cust_mask.NS) 
            mask = E_DTI_Create_Mask_From_DWI_enhanced(DWI,NrB0,par.mask_P.NS.NDWI,par.mask_P.NS.DWI,par.mask_P.NS.mfs);
        else
            fn_cm = [FOL_FI filesep F par.cust_mask.NS];
            [mask, VDimsm, suc] = E_DTI_read_nifti_file(fn_cm);
            if suc==0
                return;
            end
            mask = mask>0;
        end
        mask(dummy)=0;
            
        par_temp = par;
        par_temp.TE = par.TE.NS;
        
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
        
        f_out = [par.out_folder filesep F par.suff.NS];
        
        for_trafo.f_out = f_out;
        
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
                'FA','FE','SE','eigval','DT','outlier','DWIB0','chi_sq','chi_sq_iqr','DM_info','par');
        catch me
            suc = 0;
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
        
        if par.R2D.type==0
            E_DTI_remove_temp_f(dir_temp);
        end
        
    else
        
        for_trafo.f_out=[];
        disp('Unity transform is not supported with SOLID');
        return;
        % for i=1:Le  
        %     suc = E_DTI_SMECEPI_make_unity_transf(trafo_names{i},VoxS{i},size(DWI{i}));
        %     if suc==0
        %         return;
        %     end
        % end
        
    end
    
    for_trafo.f_moving = f_moving;
    for_trafo.dir_temp_i = dir_temp_i;
    for_trafo.trafo_names = trafo_names;
    for_trafo.result_names = result;
    for_trafo.par = par;    
    for_trafo.solid_moving = f_solid; %SOLID
    for_trafo.solid_names = result_solid; % SOLID    

end