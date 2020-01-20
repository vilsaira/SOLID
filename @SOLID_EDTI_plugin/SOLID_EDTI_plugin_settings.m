function plugin = SOLID_EDTI_plugin_settings(plugin, source, event)    

    AddOpts.Resize='on';
    AddOpts.WindowStyle='normal';
    AddOpts.Interpreter='tex';
    
    AddOpts.pos = [443 4 1110 1057];
    
    prompt = {'Metric (Mean / Variance / Iod)',...
        'Lower threshold (default: 3.5)',...
        'Upper threshold (default: 10.0)',...
        'Use masking (1 / 0)',...
        'Output SOLID results (1 / 0)'};
    dlg_title = 'SOLID parameter settings';
    num_lines = 1;
    
    def = {plugin.metric, num2str(plugin.thrL), num2str(plugin.thrU), num2str(1), num2str(0)};
    
    answer = my_inputdlg(prompt,dlg_title,num_lines,def,AddOpts);
    pause(0.05)
    
    if isempty(answer)
        return;
    end

    metric = answer{1};
    thrL = str2double(answer{2});
    thrU = str2double(answer{3});
    useMasking = str2double(answer{4}) > 0;
    saveResults = str2double(answer{5}) > 0;

    if ~any(ismember(metric, {'Variance', 'Mean', 'Iod'}))
        my_msgbox('Incorrect input on line 1...','Modal');
        return;
    end

    if thrL < 0 || thrL > thrU
        my_msgbox('Incorrect input on line 2...','Modal');
        return;
    end

    if thrU < 0 || thrU < thrL
        my_msgbox('Incorrect input on line 3...','Modal');
        return;
    end
    
    if ~islogical(useMasking) 
        my_msgbox('Incorrect input on line 4...','Modal');
        return;
    end

    if ~islogical(saveResults) 
        my_msgbox('Incorrect input on line 5...','Modal');
        return;
    end

    plugin.metric = metric;
    plugin.thrL = thrL;
    plugin.thrU = thrU;
    plugin.useMask = useMasking;
    plugin.saveResults = saveResults;
    
end