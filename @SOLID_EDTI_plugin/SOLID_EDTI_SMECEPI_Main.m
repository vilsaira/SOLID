function SOLID_EDTI_plugin_Main(Par_file)
% SOLID_EDTI_plugin_Main is a fork of E_DTI_SMECEPI_Main (v. 4.8.6)
% This function is accessible via commandline and ExploreDTI plugins menu SOLID SMECEPI button.

    if nargin==0
        
        figHandles = get(groot, 'Children');        
        EDTI = findobj(figHandles, 'Tag', 'MainExploreDTI');

        par = E_DTI_Get_SMECEPI_par;
        if par.suc==0
            return;
        end
        
        % Search for already existing and possibly manually edited modZ Nifti from the same path with *.mat file
        extInd = strfind(par.DTI_f_in, '.mat');
        modZNiiName = [par.DTI_f_in(1:extInd-1), '_L_', ...
                      num2str(EDTI.UserData.solid.thrL), '_U_', ...
                      num2str(EDTI.UserData.solid.thrU), '_',...
                      EDTI.UserData.solid.metric, '_masked_', num2str(EDTI.UserData.solid.useMask), '.nii'];
        if ~exist(modZNiiName, 'file');
            modZNiiName = [modZNiiName, '.gz'];
        end
        if ~exist(modZNiiName, 'file')
            modZNiiName = []; % If file is not found, modZ map will be calculated during SMECEPI
        end            
        
        par.SOLID.modZ_f_in = modZNiiName;
        par.SOLID.mask_f_in = [];
        par.SOLID.thrL = EDTI.UserData.solid.thrL;
        par.SOLID.thrU = EDTI.UserData.solid.thrU;
        par.SOLID.useMask = EDTI.UserData.solid.useMask;
        par.SOLID.metric = EDTI.UserData.solid.metric;
        par.SOLID.saveResults = EDTI.UserData.solid.saveResults;
        
    elseif nargin==1
        
        Af = textread(Par_file,'%s','delimiter','\n');
        for i=1:length(Af)
            eval(Af{i});
        end
        par.no_GUI = 1;
        par.EPI.use_f_mask = 1;
        par.use_f_mask = 1;
        if isfield(par,'MDC_constr_fac')==0
            h_f = findobj('Tag','MainExploreDTI');
            data = get(h_f, 'userdata');
            par.MDC_constr_fac = data.MDC_constr_fac;
        end
        if isfield(par, 'SOLID') == 0
            disp(['No SOLID parameters found in the file ', Par_file, '.']);
            disp('Default SOLID parameters will be used.');
            par.SOLID.modZ_f_in = [];
            par.SOLID.mask_f_in = [];
            par.SOLID.thrL = 3.5;
            par.SOLID.thrU = 10.0;
            par.SOLID.useMask = 1;
            par.SOLID.metric = 'Variance';
            par.SOLID.saveResults = 0;
        end

    else
        
        disp('The input of this function is either:')
        disp('(1) Empty: then, it uses the ExploreDTI GUI parameters.')
        disp('(2) A *.txt file with parameter settings.')
        return;
        
    end
    
    % h_w = my_waitbar(0,'Checking data...');pause(0.5)
    DTI_files = E_DTI_SMECEPI_Get_input_DTI(par);
    % my_waitbar(1);close(h_w);pause(0.01);

    if isempty(DTI_files)
        if par.no_GUI==1
            disp('Could not find relevant *.mat files... ')
        else
            my_msgbox('Could not find relevant *.mat files','Warning...','Modal')
        end
        return;
    end

    % disp('Yet another check of par.E_path:')
    % disp(par.E_path)

    suc = E_DTI_do_initial_check_reg_tra_files(par.E_path,par.T_path);

    if suc==0
        if par.no_GUI==0
            my_msgbox('See command prompt for more info...','Error...','Modal')
        end
        return;    
    end

    TS = tic;
    if par.no_GUI==0
        h_w = my_waitbar(0,'Correcting for subject motion and EC/EPI distortions...');pause(0.01)
    end
    for i=1:length(DTI_files)
        disp('Processing file:')
        disp([DTI_files{i} ' ...'])

        SOLID_EDTI_plugin.SOLID_EDTI_SMECEPI_Single(DTI_files{i},par)

        if par.no_GUI==0
            my_waitbar(i/length(DTI_files));
        end
        disp('Done!')
    end
    if par.no_GUI==0
        my_waitbar(1);close(h_w);pause(0.01);
    end

    disp('Processing finished!')

    t=toc(TS);
    if t<3600
        m = t/60;
        disp(['Total computation time was ' num2str(m) ' minutes.'])
    elseif t>3600 && t<(3600*24)
        h = t/3600;
        disp(['Total computation time was ' num2str(h) ' hours.'])
    else
        d = t/(3600*24);
        disp(['Total computation time was ' num2str(d) ' days.'])
    end

end