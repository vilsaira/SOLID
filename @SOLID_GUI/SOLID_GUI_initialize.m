function solidGui = SOLID_GUI_initialize(solidGui)

    pos = get(0, 'MonitorPositions');
    pos = pos(1,:);
    pos = [pos(3)*0.1 pos(4)*0.1 pos(3)*0.8 pos(4)*0.8];

    solidGui.ui = figure;
    solidGui.ui.Name = 'SOLID GUI';
    solidGui.ui.Tag = 'SOLID GUI';
    solidGui.ui.Color = [0, 0, 0];
    solidGui.ui.Position = pos;
    solidGui.ui.Units = 'normalized';
    solidGui.ui.MenuBar = 'none';
    solidGui.ui.ToolBar = 'figure';
    solidGui.ui.WindowButtonDownFcn = @solidGui.SOLID_GUI_callbackMouseClick;
    solidGui.ui.WindowButtonUpFcn = @solidGui.SOLID_GUI_callbackMouseRelease;
    solidGui.ui.WindowScrollWheelFcn = @solidGui.SOLID_GUI_callbackMouseScroll;
    solidGui.ui.WindowButtonMotionFcn = @solidGui.SOLID_GUI_callbackMouseMove;
    solidGui.ui.WindowKeyPressFcn = @solidGui.SOLID_GUI_callbackKeyPress;
    solidGui.ui.WindowKeyReleaseFcn = @solidGui.SOLID_GUI_callbackKeyRelease;

    % DockFig(solidGui.fig)

    ver = version('-release');
    verYear = str2double(ver(1:4));
    if verYear < 2016
        menuVerLabel = 'Label';
        menuVerCallback = 'Callback';
    else
        menuVerLabel = 'Text';
        menuVerCallback = 'MenuSelectedFcn';
    end

    solidGui.menuFile = uimenu(solidGui.ui);
    solidGui.menuFileOpen = uimenu(solidGui.menuFile);
    solidGui.menuFileSave = uimenu(solidGui.menuFile);

    solidGui.menuSetting = uimenu(solidGui.ui);
    solidGui.menuSettingMetric = uimenu(solidGui.menuSetting);
    solidGui.menuSettingMetricVar = uimenu(solidGui.menuSettingMetric);
    solidGui.menuSettingMetricMean = uimenu(solidGui.menuSettingMetric);
    solidGui.menuSettingMetricIod = uimenu(solidGui.menuSettingMetric);
    solidGui.menuSettingMask = uimenu(solidGui.menuSetting);    

    solidGui.menuHelp = uimenu(solidGui.ui);
    solidGui.menuHelpTheory = uimenu(solidGui.menuHelp);
    solidGui.menuHelpGithub = uimenu(solidGui.menuHelp);

    solidGui.menuFile.(menuVerLabel) = 'File';
    solidGui.menuFileOpen.(menuVerLabel) = 'Open';
    solidGui.menuFileOpen.Accelerator = 'O';
    solidGui.menuFileOpen.(menuVerCallback) = @solidGui.SOLID_GUI_callbackLoadData;
    solidGui.menuFileSave.(menuVerLabel) = 'Save results';
    solidGui.menuFileSave.Accelerator = 'S';
    solidGui.menuFileSave.(menuVerCallback) = @solidGui.SOLID_GUI_callbackSave;

    solidGui.menuSetting.(menuVerLabel) = 'Settings';
    solidGui.menuSettingMetric.(menuVerLabel) = 'Select SOLID metric';
    solidGui.menuSettingMetricVar.(menuVerLabel) = 'Variance';
    solidGui.menuSettingMetricMean.(menuVerLabel) = 'Mean';
    solidGui.menuSettingMetricIod.(menuVerLabel) = 'Iod';
    solidGui.menuSettingMetricVar.(menuVerCallback) = @solidGui.SOLID_GUI_callbackToggleMetric;
    solidGui.menuSettingMetricVar.UserData.Name = 'Variance';
    solidGui.menuSettingMetricMean.(menuVerCallback) = @solidGui.SOLID_GUI_callbackToggleMetric;
    solidGui.menuSettingMetricMean.UserData.Name = 'Mean';
    solidGui.menuSettingMetricIod.(menuVerCallback) = @solidGui.SOLID_GUI_callbackToggleMetric;
    solidGui.menuSettingMetricIod.UserData.Name = 'Iod';
    solidGui.menuSettingMetricVar.Checked = 'on';

    solidGui.menuSettingMask.(menuVerLabel) = 'Use brain masking';
    solidGui.menuSettingMask.(menuVerCallback) = @solidGui.SOLID_GUI_callbackToggleUseMask;
    solidGui.menuSettingMask.Checked = 'on';

    solidGui.menuHelp.(menuVerLabel) = 'Help';
    solidGui.menuHelpTheory.(menuVerLabel) = 'See original publication (www)';
    solidGui.menuHelpTheory.(menuVerCallback) = @theory;
    solidGui.menuHelpGithub.(menuVerLabel) = 'See help on Github (www)';
    solidGui.menuHelpGithub.(menuVerCallback) = @github;
        function github(event, source)
            url = 'https://github.com/vilsaira/SOLID';
            web(url, '-new','-browser');
        end
        function theory(event, source)
            url = 'https://doi.org/10.1016/j.neuroimage.2018.07.003';
            web(url, '-new','-browser');
        end

    % Create axes 3x3 grid for DWI viewing
    W = 1/4;
    H = 1/4;
    X = [0, 1, 2].*W;
    Y = [3*H - 2/21, ...
        2*H - 2/20, ...
        1*H - 2/19];

    solidGui.ax.dwiTopLeft = axes('Parent', solidGui.ui);
    solidGui.ax.dwiTopLeft.Position = [X(1) Y(1) W H];            
    solidGui.img.dwiTopLeft = imagesc('Parent', solidGui.ax.dwiTopLeft, 'CData', NaN(50,50));
    axis(solidGui.ax.dwiTopLeft, 'off', 'tight', 'equal');

    solidGui.ax.dwiTopCenter = axes('Parent', solidGui.ui);
    solidGui.ax.dwiTopCenter.Position = [X(2) Y(1) W H];
    solidGui.img.dwiTopCenter = imagesc('Parent', solidGui.ax.dwiTopCenter, 'CData', NaN(50,50));
    axis(solidGui.ax.dwiTopCenter, 'off', 'tight', 'equal');

    solidGui.ax.dwiTopRight = axes('Parent', solidGui.ui);
    solidGui.ax.dwiTopRight.Position = [X(3) Y(1) W H];
    solidGui.img.dwiTopRight = imagesc('Parent', solidGui.ax.dwiTopRight, 'CData', NaN(50,50));
    axis(solidGui.ax.dwiTopRight, 'off', 'tight', 'equal');

    solidGui.ax.dwiMiddleLeft = axes('Parent', solidGui.ui);
    solidGui.ax.dwiMiddleLeft.Position = [X(1) Y(2) W H];
    solidGui.img.dwiMiddleLeft = imagesc('Parent', solidGui.ax.dwiMiddleLeft, 'CData', NaN(50,50));
    axis(solidGui.ax.dwiMiddleLeft, 'tight', 'equal');
    solidGui.ax.dwiMiddleLeft.XColor = [0, 1, 1];
    solidGui.ax.dwiMiddleLeft.YColor = [0, 1, 1];
    solidGui.ax.dwiMiddleLeft.Box = 'on';
    solidGui.ax.dwiMiddleLeft.LineWidth = 2;
    solidGui.ax.dwiMiddleLeft.XTick = [];
    solidGui.ax.dwiMiddleLeft.YTick = [];

    solidGui.ax.dwiMiddleCenter = axes('Parent', solidGui.ui);
    solidGui.ax.dwiMiddleCenter.Position = [X(2) Y(2) W H];
    solidGui.img.dwiMiddleCenter = imagesc('Parent', solidGui.ax.dwiMiddleCenter, 'CData', NaN(50,50));
    axis(solidGui.ax.dwiMiddleCenter', 'tight', 'equal');
    solidGui.ax.dwiMiddleCenter.XColor = [1, 1, 1];
    solidGui.ax.dwiMiddleCenter.YColor = [1, 1, 1];
    solidGui.ax.dwiMiddleCenter.Box = 'on';
    solidGui.ax.dwiMiddleCenter.LineWidth = 2;
    solidGui.ax.dwiMiddleCenter.XTick = [];
    solidGui.ax.dwiMiddleCenter.YTick = [];

    solidGui.ax.dwiMiddleRight = axes('Parent', solidGui.ui);
    solidGui.ax.dwiMiddleRight.Position = [X(3) Y(2) W H];
    solidGui.img.dwiMiddleRight = imagesc('Parent', solidGui.ax.dwiMiddleRight, 'CData', NaN(50,50));
    axis(solidGui.ax.dwiMiddleRight, 'tight', 'equal');
    solidGui.ax.dwiMiddleRight.XColor = [1, 0, 1];
    solidGui.ax.dwiMiddleRight.YColor = [1, 0, 1];
    solidGui.ax.dwiMiddleRight.Box = 'on';
    solidGui.ax.dwiMiddleRight.LineWidth = 2;
    solidGui.ax.dwiMiddleRight.XTick = [];
    solidGui.ax.dwiMiddleRight.YTick = [];

    solidGui.ax.dwiBottomCenter = axes('Parent', solidGui.ui);
    solidGui.ax.dwiBottomCenter.Position = [X(2) Y(3) W H];
    solidGui.img.dwiBottomCenter = imagesc('Parent', solidGui.ax.dwiBottomCenter, 'CData', NaN(50,50));
    axis(solidGui.ax.dwiBottomCenter, 'off', 'tight', 'equal');

    solidGui.ax.qSpace = axes('Parent', solidGui.ui);
    solidGui.ax.qSpace.Position = [X(3) Y(3) W H];
    solidGui.img.qSpace = imagesc('Parent', solidGui.ax.qSpace, 'CData', NaN(50,50));
    axis(solidGui.ax.qSpace, 'square', 'tight');

    solidGui.ax.modZ2D = axes('Parent', solidGui.ui);    
    solidGui.img.modZ2D  = imagesc('Parent', solidGui.ax.modZ2D, 'CData', rand(100,100));
    solidGui.scrollPanel = imscrollpanel(solidGui.ui, solidGui.img.modZ2D);
    axis(solidGui.ax.modZ2D, 'equal', 'tight');    
    solidGui.ax.modZ2D.Box = 'on';
    solidGui.scrollPanel.Position = [X(3)+W+1/100 2/20 W-2/30 16/20];
    api = iptgetapi(solidGui.scrollPanel);
    mag = api.findFitMag();
    api.setMagnification(mag);    
    solidGui.ax.modZ2D.Units = 'normalized';
    solidGui.ax.modZ2D.Position = [0, 0, 1, 1];

    solidGui.ax.modZ2D_colorbar = axes('Parent', solidGui.ui);
    solidGui.img.modZ2D_colorbar = imagesc('Parent', solidGui.ax.modZ2D_colorbar, 'CData', rand(1000,500))
    axis(solidGui.ax.modZ2D_colorbar, 'tight');
    solidGui.ax.modZ2D_colorbar.Box = 'off';
    ch = colorbar(solidGui.ax.modZ2D_colorbar);
    ylabel(ch, 'Modified Z-score');
    ch.Color = [0.9 0.9 0.9];
    solidGui.ax.modZ2D_colorbar.Position = [X(3)+W+1/100 2/20 W-2/30-1/90 16/20];
    % solidGui.ax.modZ2D_colorbar.Visible = 'off';
    solidGui.img.modZ2D_colorbar.Visible = 'off';
    xlabel(solidGui.ax.modZ2D_colorbar, 'Image volume');
    ylabel(solidGui.ax.modZ2D_colorbar, 'Slice');
    solidGui.ax.modZ2D_colorbar.XColor = [0.9 0.9 0.9];
    solidGui.ax.modZ2D_colorbar.YColor = [0.9 0.9 0.9];
    solidGui.ax.modZ2D_colorbar.XTick = [];
    solidGui.ax.modZ2D_colorbar.YTick = [];
    solidGui.ax.modZ2D_colorbar.Color = [0,0,0];

    h_axes = axes('position', ch.Position, 'ylim', ch.Limits, 'color', 'none', 'visible','off');
    solidGui.lines.modZColorbarLine = line(h_axes.XLim, 0*[1 1], 'color', 'red', 'linewidth', 3, 'parent', h_axes);

    solidGui.ax.modZhistogram = axes('Parent', solidGui.ui);
    solidGui.ax.modZhistogram.Position =  [1/20 1/20 2/6-1/10 2/6];
    solidGui.img.modZhistogram = imagesc('Parent', solidGui.ax.modZhistogram, 'CData', NaN(50,50));
    axis(solidGui.ax.modZhistogram, 'tight', 'equal');
    xlabel(solidGui.ax.modZhistogram, 'Modified Z-score');
    ylabel(solidGui.ax.modZhistogram, 'Counts');
    solidGui.ax.modZhistogram.XColor = [0.9 0.9 0.9];
    solidGui.ax.modZhistogram.YColor = [0.9 0.9 0.9];


    solidGui.slider.cor = uicontrol('Parent', solidGui.ui, 'style', 'slider');
    solidGui.slider.cor.BackgroundColor = 'green';
    solidGui.slider.cor.Units = 'normalized';
    solidGui.slider.cor.Position = [0+1/100, 1-1/20+1/80, 1/4-2/100, 1/40];
    solidGui.slider.cor.Callback = @solidGui.SOLID_GUI_callbackSlider;
    solidGui.slider.cor.UserData.Name = 'cor';

    solidGui.slider.sag = uicontrol('Parent', solidGui.ui, 'style', 'slider');
    solidGui.slider.sag.BackgroundColor = 'red';
    solidGui.slider.sag.Units = 'normalized';
    solidGui.slider.sag.Position = [1/4+1/100, 1-1/20+1/80, 1/4-2/100, 1/40];
    solidGui.slider.sag.Callback = @solidGui.SOLID_GUI_callbackSlider;
    solidGui.slider.sag.UserData.Name = 'sag';

    solidGui.slider.axi = uicontrol('Parent', solidGui.ui, 'style', 'slider');
    solidGui.slider.axi.BackgroundColor = 'blue';
    solidGui.slider.axi.Units = 'normalized';
    solidGui.slider.axi.Position = [2/4+1/100, 1-1/20+1/80, 1/4-2/100, 1/40];
    solidGui.slider.axi.Callback = @solidGui.SOLID_GUI_callbackSlider;
    solidGui.slider.axi.UserData.Name = 'axi';

    solidGui.slider.dwi = uicontrol('Parent', solidGui.ui, 'style', 'slider');
    solidGui.slider.dwi.BackgroundColor = 'white';
    solidGui.slider.dwi.Units = 'normalized';
    solidGui.slider.dwi.Position = [3/4+1/100, 1-1/20+1/80, 1/4-2/30, 1/40];
    solidGui.slider.dwi.Callback = @solidGui.SOLID_GUI_callbackSlider;
    solidGui.slider.dwi.UserData.Name = 'dwi';

    solidGui.slider.slice = uicontrol('Parent', solidGui.ui, 'style', 'slider');
    solidGui.slider.slice.BackgroundColor = 'white';
    solidGui.slider.slice.Units = 'normalized';
    pos = solidGui.ax.modZ2D.Position;
    solidGui.slider.slice.Position = [pos(1)+pos(3)+1/15, pos(2), 1/80, pos(4)];
    solidGui.slider.slice.Callback = @solidGui.SOLID_GUI_callbackSlider;
    solidGui.slider.slice.Visible = 'off';
    solidGui.slider.slice.UserData.Name = 'slice';

    % Text fields
    pos = solidGui.slider.cor.Position;
    solidGui.txtfield.cor = uicontrol('Parent', solidGui.ui, 'style', 'edit');
    solidGui.txtfield.cor.BackgroundColor = 'black';            
    solidGui.txtfield.cor.Units = 'normalized';
    solidGui.txtfield.cor.Position = pos + [0, -1/30, 0, 0];
    solidGui.txtfield.cor.String = 'Coronal';
    solidGui.txtfield.cor.ForegroundColor = [0.9 0.9 0.9];
    solidGui.txtfield.cor.Callback = @solidGui.SOLID_GUI_callbackTxtfield;
    solidGui.txtfield.cor.UserData.Name = 'cor';

    pos = solidGui.slider.sag.Position;
    solidGui.txtfield.sag = uicontrol('Parent', solidGui.ui, 'style', 'edit');
    solidGui.txtfield.sag.BackgroundColor = 'black';            
    solidGui.txtfield.sag.Units = 'normalized';
    solidGui.txtfield.sag.Position = pos + [0, -1/30, 0, 0];
    solidGui.txtfield.sag.String = 'Sagittal';
    solidGui.txtfield.sag.ForegroundColor = [0.9 0.9 0.9];
    solidGui.txtfield.sag.Callback = @solidGui.SOLID_GUI_callbackTxtfield;
    solidGui.txtfield.sag.UserData.Name = 'sag';

    pos = solidGui.slider.axi.Position;
    solidGui.txtfield.axi = uicontrol('Parent', solidGui.ui, 'style', 'edit');
    solidGui.txtfield.axi.BackgroundColor = 'black';            
    solidGui.txtfield.axi.Units = 'normalized';
    solidGui.txtfield.axi.Position = pos + [0, -1/30, 0, 0];
    solidGui.txtfield.axi.String = 'Axial';
    solidGui.txtfield.axi.ForegroundColor = [0.9 0.9 0.9];
    solidGui.txtfield.axi.Callback = @solidGui.SOLID_GUI_callbackTxtfield;
    solidGui.txtfield.axi.UserData.Name = 'axi';

    pos = solidGui.slider.dwi.Position;
    solidGui.txtfield.dwi = uicontrol('Parent', solidGui.ui, 'style', 'edit');
    solidGui.txtfield.dwi.BackgroundColor = 'black';            
    solidGui.txtfield.dwi.Units = 'normalized';
    solidGui.txtfield.dwi.Position = pos + [0, -1/30, 0, 0];
    solidGui.txtfield.dwi.String = 'DWI';
    solidGui.txtfield.dwi.ForegroundColor = [0.9 0.9 0.9];
    solidGui.txtfield.dwi.Callback = @solidGui.SOLID_GUI_callbackTxtfield;
    solidGui.txtfield.dwi.UserData.Name = 'dwi';

    solidGui.txtfield.thresholdLowert = uicontrol('Parent', solidGui.ui, 'style', 'text');
    solidGui.txtfield.thresholdLowert.BackgroundColor = 'black';            
    solidGui.txtfield.thresholdLowert.Units = 'normalized';
    solidGui.txtfield.thresholdLowert.Position = [3/4+1/100, 1/20, 1/30, 1/40];
    solidGui.txtfield.thresholdLowert.String = 'Lower_t';
    solidGui.txtfield.thresholdLowert.FontWeight = 'bold';
    solidGui.txtfield.thresholdLowert.ForegroundColor = [0.9 0.9 0.9];    

    pos = solidGui.txtfield.thresholdLowert.Position;
    solidGui.txtfield.thresholdLower = uicontrol('Parent', solidGui.ui, 'style', 'edit');
    solidGui.txtfield.thresholdLower.BackgroundColor = 'black';            
    solidGui.txtfield.thresholdLower.Units = 'normalized';
    solidGui.txtfield.thresholdLower.Position = [pos(1)+pos(3)+1/100, 1/20, 1/30, 1/40];
    solidGui.txtfield.thresholdLower.String = '3.5';
    solidGui.txtfield.thresholdLower.ForegroundColor = [0.9 0.9 0.9];
    solidGui.txtfield.thresholdLower.Callback = @solidGui.SOLID_GUI_callbackTxtfield;
    solidGui.txtfield.thresholdLower.UserData.Name = 'thrL';

    solidGui.txtfield.thresholdUppert = uicontrol('Parent', solidGui.ui, 'style', 'text');
    solidGui.txtfield.thresholdUppert.BackgroundColor = 'black';            
    solidGui.txtfield.thresholdUppert.Units = 'normalized';
    solidGui.txtfield.thresholdUppert.Position = [3/4+1/100, 1/40, 1/30, 1/40];
    solidGui.txtfield.thresholdUppert.String = 'Upper_t';
    solidGui.txtfield.thresholdUppert.FontWeight = 'bold';
    solidGui.txtfield.thresholdUppert.ForegroundColor = [0.9 0.9 0.9];

    pos = solidGui.txtfield.thresholdUppert.Position;
    solidGui.txtfield.thresholdUpper = uicontrol('Parent', solidGui.ui, 'style', 'edit');
    solidGui.txtfield.thresholdUpper.BackgroundColor = 'black';            
    solidGui.txtfield.thresholdUpper.Units = 'normalized';
    solidGui.txtfield.thresholdUpper.Position = [pos(1)+pos(3)+1/100, 1/40, 1/30, 1/40];
    solidGui.txtfield.thresholdUpper.String = '10.0';
    solidGui.txtfield.thresholdUpper.ForegroundColor = [0.9 0.9 0.9];
    solidGui.txtfield.thresholdUpper.Callback = @solidGui.SOLID_GUI_callbackTxtfield;    
    solidGui.txtfield.thresholdUpper.UserData.Name = 'thrU';

    solidGui.togglebtn.mask = uicontrol('Parent', solidGui.ui, 'style', 'togglebutton');
    solidGui.togglebtn.mask.Units = 'normalized';
    solidGui.togglebtn.mask.Position = [pos(1)+pos(3)+1/20, 1/40, 1/20, 1/40];
    solidGui.togglebtn.mask.String = 'Show mask';                  
    solidGui.togglebtn.mask.Callback = @solidGui.SOLID_GUI_callbackShowMask;


    axes(solidGui.ax.dwiMiddleCenter)
    solidGui.lines.axiSAG = line([1,2], [1,2], 'color', 'red');
    solidGui.lines.axiCOR = line([1,2], [1,2], 'color', 'green');
    axes(solidGui.ax.dwiTopLeft)
    solidGui.lines.corSAG = line([1,2], [1,2], 'color', 'red');
    solidGui.lines.corAXI = line([1,2], [1,2], 'color', 'blue');
    axes(solidGui.ax.dwiTopRight)
    solidGui.lines.sagCOR = line([1,2], [1,2], 'color', 'green');
    solidGui.lines.sagAXI = line([1,2], [1,2], 'color', 'blue');
    axes(solidGui.ax.modZ2D)
    solidGui.lines.modZDWI = line([1,2], [1,2], 'color', 'white');
    solidGui.lines.modZAXI = line([1,2], [1,2], 'color', 'white');
    hold(solidGui.ax.modZ2D, 'on');
    solidGui.lines.modZM = plot(1, 1, 'm.', 'markersize', 20);
    solidGui.lines.modZG = plot(1, 1, 'c.', 'markersize', 20);
    hold(solidGui.ax.modZ2D, 'off');

    solidGui.lines.qSpaceVecs = plot3(solidGui.ax.qSpace, 0,0,0, 'ro');
    hold(solidGui.ax.qSpace, 'on');
    solidGui.lines.qSpaceVecsCur = plot3(solidGui.ax.qSpace, 0,0,0, 'gx', 'markersize', 24);
    solidGui.lines.qSpaceVecsNext1 = plot3(solidGui.ax.qSpace, 0,0,0, 'c.', 'markersize', 20);
    solidGui.lines.qSpaceVecsNext2 = plot3(solidGui.ax.qSpace, 0,0,0, 'm.', 'markersize', 20);    
    
    [X,Y,Z] = peaks(15);
    C(:,:,1) = zeros(15); % red
    C(:,:,2) = ones(15).*linspace(0.5,0.6,15); % green
    C(:,:,3) = ones(15).*linspace(0,1,15); % blue
    solidGui.lines.qSpaceSurf = surf(solidGui.ax.qSpace, X, Y, Z, C, 'EdgeColor', 'none');
    
    hold(solidGui.ax.qSpace, 'off')
    axis(solidGui.ax.qSpace, 'square');
    grid(solidGui.ax.qSpace, 'on');
    xlabel(solidGui.ax.qSpace, 'X');
    ylabel(solidGui.ax.qSpace, 'Y');
    zlabel(solidGui.ax.qSpace, 'Z');
    solidGui.ax.qSpace.XColor = [0.5, 0.5, 0.5];
    solidGui.ax.qSpace.YColor = [0.5, 0.5, 0.5];
    solidGui.ax.qSpace.ZColor = [0.5, 0.5, 0.5];    
    
    axis(solidGui.ax.modZhistogram, 'tight');
    solidGui.lines.modZHist = histogram(solidGui.ax.modZhistogram,0,1);
    hold(solidGui.ax.modZhistogram, 'on');
    solidGui.lines.modZHistCurrent = line(solidGui.ax.modZhistogram, [0,0], [0,0], 'color', 'red', 'linewidth', 2);
    solidGui.lines.modZHistLower = line(solidGui.ax.modZhistogram, [0,0], [0,0], 'color', 'green', 'linestyle', '--', 'linewidth', 2);
    solidGui.lines.modZHistUpper = line(solidGui.ax.modZhistogram, [0,0], [0,0], 'color', 'blue', 'linestyle', '--', 'linewidth', 2);
    hold(solidGui.ax.modZhistogram, 'off');
    solidGui.lines.modZHistLegend = legend(solidGui.ax.modZhistogram.Children([1,2,3]), 'Upper threshold', 'Lower threshold', 'Current slice');
    solidGui.ax.modZhistogram.XColor = [0.9 0.9 0.9];
    solidGui.ax.modZhistogram.YColor = [0.9 0.9 0.9];    
    xlabel(solidGui.ax.modZhistogram, 'Modified Z-score');
    ylabel(solidGui.ax.modZhistogram, 'Number of slices');

    pos = solidGui.slider.dwi.Position;
    solidGui.txtfield.shelltxt = uicontrol('Parent', solidGui.ui, 'style', 'text');
    solidGui.txtfield.shelltxt.BackgroundColor = 'black';            
    solidGui.txtfield.shelltxt.Units = 'normalized';
    solidGui.txtfield.shelltxt.Position = [pos(1)+pos(3)+1/300, pos(2)-1/40, 1/20, 1/20];
    solidGui.txtfield.shelltxt.String = 'B-value';
    solidGui.txtfield.shelltxt.FontWeight = 'bold';
    solidGui.txtfield.shelltxt.ForegroundColor = [0.9 0.9 0.9];

    pos = solidGui.txtfield.dwi.Position;
    solidGui.popup.shell = uicontrol('Parent', solidGui.ui, 'Style', 'popup', 'String', {'Shells'});
    solidGui.popup.shell.Units = 'normalized';
    % solidGui.popup.shell.Position = [4/13, 1/8, 1/8, 1/100];
    solidGui.popup.shell.Position = [pos(1)+pos(3)+1/300, pos(2)-1/40, 1/20, 1/20];
    solidGui.popup.shell.Callback = @solidGui.SOLID_GUI_callbackSelectShell;

    solidGui.data.window.changing = false;
    solidGui.data.window.width = 0.5;
    solidGui.data.window.level = 0.5;
    solidGui.data.window.mousePosI = [0,0];
    solidGui.data.window.mousePosF = [0,0];
    solidGui.data.window.timer = timer('Name', 'WindowTimer', ...
        'Period', 0.001, ...
        'StartDelay', 0.001, ...
        'TasksToExecute', inf, ...
        'ExecutionMode', 'fixedSpacing', ...
        'TimerFcn', @solidGui.SOLID_GUI_callbackWindowTimer);
    
end