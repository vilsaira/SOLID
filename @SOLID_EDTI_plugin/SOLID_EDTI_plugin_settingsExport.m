function SOLID_EDTI_plugin_settingsExport(plugin, source, event)

    [filename, pathname] = uigetfile({'*.txt','Text files (*.txt)'},'Select parameter file for ExploreDTI SM/EC/EPI corrections','Parameters_SM_EC_EPI.txt');
    exist_DTI = isa(filename,'char');
    if exist_DTI
        fname = fullfile(pathname,filename);
    else
        return;
    end

    % Check if SOLID parameters are already present in the file
    Af = textread(fname,'%s','delimiter','\n');
    for i=1:length(Af)
        if any(strfind(Af{i}, 'SOLID'));
            disp(['SOLID parameters already present in the file: ', fname]);
            return;
        end
    end    

    thrL = plugin.thrL;
    thrU = plugin.thrU;
    metric = plugin.metric;
    useMask = plugin.useMask;
    saveResults = plugin.saveResults;

    L = {};
    cn = 0;
    cn = cn+1; L{cn} = '% ';
    cn = cn+1; L{cn} = [' par.SOLID.thrL = ', num2str(thrL), ';'];
    cn = cn+1; L{cn} = '% SOLID - the lower modified Z-score threshold (default 3.5).';
    cn = cn+1; L{cn} = '% ------------------------------------------------------';
    cn = cn+1; L{cn} = '% ';
    cn = cn+1; L{cn} = [' par.SOLID.thrU = ', num2str(thrU), ';'];
    cn = cn+1; L{cn} = '% SOLID - the upper modified Z-score threshold (default 10.0).';
    cn = cn+1; L{cn} = '% ------------------------------------------------------';
    cn = cn+1; L{cn} = '% ';
    cn = cn+1; L{cn} = [' par.SOLID.metric = ''', metric, ''';'];
    cn = cn+1; L{cn} = '% SOLID - metric (Variance/Mean/Iod).';
    cn = cn+1; L{cn} = '% ------------------------------------------------------';
    cn = cn+1; L{cn} = '% ';
    cn = cn+1; L{cn} = [' par.SOLID.useMask = ', num2str(useMask), ';'];
    cn = cn+1; L{cn} = '% SOLID - select if brain masking is used .';
    cn = cn+1; L{cn} = '% ------------------------------------------------------';    
    cn = cn+1; L{cn} = '% ';
    cn = cn+1; L{cn} = [' par.SOLID.modZ_f_in = '''';'];
    cn = cn+1; L{cn} = '% SOLID - path to previously calculated modZ NIfTI.';
    cn = cn+1; L{cn} = '% ------------------------------------------------------';
    cn = cn+1; L{cn} = '% ';
    cn = cn+1; L{cn} = [' par.SOLID.mask_f_in = '''';'];
    cn = cn+1; L{cn} = '% SOLID - path to custom brain mask.';
    cn = cn+1; L{cn} = '% ------------------------------------------------------';
    cn = cn+1; L{cn} = '% ';
    cn = cn+1; L{cn} = [' par.SOLID.saveResults = ', num2str(saveResults), ';'];
    cn = cn+1; L{cn} = '% SOLID - output SOLID results (false/true).';
    cn = cn+1; L{cn} = '% ------------------------------------------------------';

    
    [fid,mess] = fopen(fname,'a+');

    if fid==-1
        disp(['Problem with opening file: ' fname])
        disp(['Error message: ' mess])
        return;
    end

    for i=1:length(L)
        fprintf(fid, '%s\n', L{i});
    end

    try
        fclose(fid);
    catch me
        disp(me.message)
        return;
    end

    if ispc
        [su, mess] = fileattrib(fname,'+w -a -s');
    else
        [su, mess] = fileattrib(fname,'+w +x','a');
    end

    if su==0
        disp(['Problem with file: ' par.Par_FN])
        disp(['Error message: ' mess])
        return;
    end
    

end