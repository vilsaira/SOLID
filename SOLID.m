classdef SOLID < handle
% SOLID - Slicewise outlier detection for diffusion weighted MRI data.
%
% General Usage: 
% 
% SOLID values are represented as modified z-scores for each slice and
% volume. You can select a range of SOLID values to flag as outliers by 
% setting a minimum / maximum threshold.
% 
% Once the threshold is set, you will see the z-score plots on the right on
% 2D and interpolated on a spherical surface and distribution of scores 
% plot on left bottom update to reflect this range.  
% 
% For the 2D z-score plot, max values will be in yellow, with all lesser
% values going to dark blue. As an example, a subject with no outliers 
% detected will be completely blue. Perfect! 
% 
% Use the left / right arrow keys to jump between outliers in the 2D plot.
% You can also manually select a single value by clicking directly on that 
% pixel. This will display the corresponding slice in the middle of the 
% GUI. 
% 
% If you do not agree with SOLID slice being an outlier, press SPACEBAR, 
% to toggle the modified Z-score value to 0. You can cancel this action 
% by pressing SPACEBAR again. 
% 
% After checking all SOLID outlier values within the specified threshold
% range, you can now save the results. The saved nifti image can be used in
% further analysis.
% 
% Other features: 
% 
% Hold left mouse button down and move mouse to adjust slice window/level.
% Use the middle button to scroll through slices. 
% 
% Trouble Shooting: 
% 
% You can also turn “Masking On” to show the mask used for SOLID. This is
% useful or checking if sections are included from outside the brain that
% can lead to spurious values.
% 
% Use your own mask:
%
% Use any method to set voxels in DWI NIfTI to negative values and load the
% result in SOLID. Remember to untick the "Use brain mask" option from
% SOLID menu.
%
% Requires ExploreDTI (E_DTI_Create_Mask_From_DWI_enhanced_IND.p)
%
% **********************  Author: Viljami Sairanen  ***********************
% *********************  viljami.sairanen@gmail.com  **********************
% *************  Website: https://github.com/vilsaira/SOLID  **************
% *********************  Last edited: 28 August 2018  *********************  
% 

   events
      Event1
   end
   
   properties (Access = private)
        figureUI
        fileMenu
        fileMenu_Open
        fileMenu_Save
        fileMenu_Metric
        fileMenu_Metric_Var
        fileMenu_Metric_Mean
        fileMenu_Metric_Iod
        fileMenu_UseMask
        popup = struct('shell', []);
        ax = struct('sag', [], ...
                      'cor', [], ...
                      'axiLK', [], ...
                      'axiLMK', [], ...
                      'axiLKM', [], ...
                      'axiLPK', [], ...
                      'axiLKP', [], ...
                      'qspace', [], ...
                      'modZ2D', [], ...
                      'modZhistogram', [], ...
                      'DWIhistogram', []);
        images = struct('sag', [],...
                        'cor', [], ...
                        'axiLK', [], ...
                        'axiLMK', [], ...
                        'axiLKM', [], ...
                        'axiLPK', [], ...
                        'axiLKP', [], ...
                        'qspace', [], ...
                        'modZ2D', [], ...
                        'modZhistogram', [], ...
                        'DWIhistogram', []);
        slider = struct('sag', [], ...
                        'cor', [], ...
                        'axi', [], ...
                        'slice', [], ...
                        'dwi', []);
        txtfield = struct('sag', [], ...
                          'cor', [], ...
                          'axi', [], ...
                          'dwi', [], ...
                          'thresholdLower', [], ...
                          'thresholdUpper', []);
       togglebtn = struct('outlier', [], ...
                          'mask', []);
       data = struct('DWI', [],...
                     'modZ', [],...
                     'modZorig', [],...
                     'bval', [],...
                     'uniqueb', [],...
                     'bvec', [],...
                     'mask', [],...
                     'metric', 'Variance',...
                     'useMask', 1,...
                     'thresholdLower', 3.5,...
                     'thresholdUpper', 10.0,...
                     'clim', [0,1],...
                     'window', [], ...
                     'minmax', [0,1],...
                     'nearestPoints', [],...
                     'fname', [],...
                     'fpath', []);
        lines = struct('axiSAG', [],...
                       'axiCOR', [],...
                       'corSAG', [],...
                       'corAXI', [],...
                       'sagCOR', [],...
                       'sagAXI', [],...
                       'modZAXI', [],...
                       'modZDWI', [],...
                       'modZHistLower', [],...
                       'modZHistUpper', [],...
                       'modZHistCurrent', [],...
                       'dwiHistogramVer', [],...
                       'dwiHistogramHor', []);

   end
   
   
    % App initialization and construction
   methods (Access = private)                
        
       function sliderCorCallback(GUI, source, event)
           value = round(GUI.slider.cor.Value);
           GUI.slider.cor.Value = value;
           GUI.txtfield.cor.String = value;
           updateaxes(GUI);           
       end
       
       function sliderSagCallback(GUI, source, event)
           value = round(GUI.slider.sag.Value);
           GUI.slider.sag.Value = value;
           GUI.txtfield.sag.String = value;
           updateaxes(GUI);             
       end
       
       function sliderAxiCallback(GUI, source, event)
           value = round(GUI.slider.axi.Value);
           GUI.slider.axi.Value = value;
           GUI.slider.slice.Value = value;
           GUI.txtfield.axi.String = value;
           updateaxes(GUI);
       end
       
       function sliderDwiCallback(GUI, source, event)
           value = round(GUI.slider.dwi.Value);
           GUI.slider.dwi.Value = value;
           GUI.txtfield.dwi.String = value;
           updateaxes(GUI);             
       end
       
       function sliderSliceCallback(GUI, source, event)
           value = round(GUI.slider.slice.Value);
           GUI.slider.axi.Value = value;
           GUI.slider.slice.Value = value;
           GUI.txtfield.axi.String = value;
           updateaxes(GUI);           
       end
       
       function txtfieldCorCallback(GUI, source, event)
           value = round(str2double(GUI.txtfield.cor.String));
           if value >= GUI.slider.cpr.Min && value <= GUI.slider.cor.Max
               GUI.slider.cor.Value = value;
               GUI.txtfield.cor.String = value;
           else
               GUI.txtfield.cor.String = GUI.slider.cor.Value;
           end
           updateaxes(GUI);
       end
       
       function txtfieldSagCallback(GUI, source, event)
           value = round(str2double(GUI.txtfield.sag.String));
           if value >= GUI.slider.sag.Min && value <= GUI.slider.sag.Max
               GUI.slider.sag.Value = value;
               GUI.txtfield.sag.String = value;
           else
               GUI.txtfield.sag.String = GUI.slider.sag.Value;
           end
           updateaxes(GUI);
       end
       
       function txtfieldAxiCallback(GUI, source, event)
           value = round(str2double(GUI.txtfield.axi.String));
           if value >= GUI.slider.axi.Min && value <= GUI.slider.axi.Max
               GUI.slider.axi.Value = value;
               GUI.slider.slice.Value = value;
               GUI.txtfield.axi.String = value;
               updateaxes(GUI);
           else
               GUI.txtfield.axi.String = GUI.slider.axi.Value;
           end
       end
       
       function txtfieldDwiCallback(GUI, source, event)
           value = round(str2double(GUI.txtfield.dwi.String));
           if value >= GUI.slider.dwi.Min && value <= GUI.slider.dwi.Max
               GUI.slider.dwi.Value = value;
               GUI.txtfield.dwi.String = value;
               updateaxes(GUI);
           else
               GUI.txtfield.dwi.String = GUI.slider.dwi.Value;
           end
       end
       
       function txtfieldthresholdLowerCallback(GUI, source, event)
           value = str2double(GUI.txtfield.thresholdLower.String);
           if value >= 0 && value <= GUI.data.thresholdUpper
               GUI.txtfield.thresholdLower.String = value;
               GUI.data.thresholdLower = value;
           end              
           updateaxes(GUI);
       end
       
       function txtfieldthresholdUpperCallback(GUI, source, event)
           value = str2double(GUI.txtfield.thresholdUpper.String);
           if value >= 0 && value > GUI.data.thresholdLower
               GUI.txtfield.thresholdUpper.String = value;
               GUI.data.thresholdUpper = value;
           end              
           updateaxes(GUI);        
       end
       
        function toggleVarMean(GUI, source, event)
            if strcmp(version('-release'), '2015a')
                metric = source.Label;
            else
                metric = source.Text;
            end
            switch metric
                case 'Variance'
                    GUI.fileMenu_Metric_Var.Checked = 'on';
                    GUI.fileMenu_Metric_Mean.Checked = 'off';
                    GUI.fileMenu_Metric_Iod.Checked = 'off';
                case 'Mean'
                    GUI.fileMenu_Metric_Var.Checked = 'off';
                    GUI.fileMenu_Metric_Mean.Checked = 'on';
                    GUI.fileMenu_Metric_Iod.Checked = 'off';
                case 'Iod'
                    GUI.fileMenu_Metric_Var.Checked = 'off';
                    GUI.fileMenu_Metric_Mean.Checked = 'off';
                    GUI.fileMenu_Metric_Iod.Checked = 'on';
            end
            GUI.data.metric = metric;
            calculateModZ(GUI);
            updateaxes(GUI);  
        end
        
        function toggleUseMask(GUI, source, event)            
            if strcmp(GUI.fileMenu_UseMask.Checked, 'on')
                GUI.data.useMask = false;
                GUI.fileMenu_UseMask.Checked = 'off';
%                 if strcmp(version('-release'), '2015a')
%                     GUI.fileMenu_UseMask.Label = 'Use mask';
%                 else
%                     GUI.fileMenu_UseMask.Text = 'Use mask';
%                 end
            else
                GUI.data.useMask = true;
                GUI.fileMenu_UseMask.Checked = 'on';
%                 if strcmp(version('-release'), '2015a')
%                     GUI.fileMenu_UseMask.Label = 'Use mask';
%                 else
%                     GUI.fileMenu_UseMask.Text = 'Use mask';
%                 end
            end
            calculateModZ(GUI);
            updateaxes(GUI);
        end
        
        function toggleMaskCallback(GUI, source, event)
            if GUI.togglebtn.mask.Value
                GUI.togglebtn.mask.BackgroundColor = 'green';
            else
                GUI.togglebtn.mask.BackgroundColor = [0.94, 0.94, 0.94];
            end
            updateaxes(GUI);              
        end

        function selectShell(GUI, source, event)
            initializeAxesSlider(GUI);
            updateaxes(GUI);
        end
        
        function initializeAxesSlider(GUI, source, event)
           shell = str2double(GUI.popup.shell.String(GUI.popup.shell.Value, :));
           shellInds = GUI.data.bval == shell;
           
           dims = size(GUI.data.DWI(:,:,:,shellInds));
           if size(dims,2) < 4
               dims(4) = 1;
           end
           
           GUI.slider.dwi.Min = 1;
           GUI.slider.dwi.Max = dims(4);
           GUI.slider.dwi.SliderStep = [1, 1]./(max([dims(4),2])-1);
           GUI.slider.dwi.Value = 1;
           
           GUI.slider.slice.Min = 1;
           GUI.slider.slice.Max = dims(3);
           GUI.slider.slice.SliderStep = [1, 1]./(dims(3)-1);
           GUI.slider.slice.Value = round(dims(3)/2);
           
           GUI.slider.axi.Min = GUI.slider.slice.Min;
           GUI.slider.axi.Max = GUI.slider.slice.Max;
           GUI.slider.axi.SliderStep = [1, 1]./(dims(3)-1);
           GUI.slider.axi.Value = GUI.slider.slice.Value;
           
           GUI.slider.cor.Min = 1;
           GUI.slider.cor.Max = dims(1);
           GUI.slider.cor.SliderStep = [1, 1]./(dims(1)-1);
           GUI.slider.cor.Value = round(dims(1)/2);
           
           GUI.slider.sag.Min = 1;
           GUI.slider.sag.Max = dims(2);
           GUI.slider.sag.SliderStep = [1, 1]./(dims(2)-1);
           GUI.slider.sag.Value = round(dims(2)/2);
           
           GUI.txtfield.dwi.Min = GUI.slider.dwi.Min;
           GUI.txtfield.dwi.Max = GUI.slider.dwi.Max;           
           GUI.txtfield.dwi.String = GUI.slider.dwi.Value;
           
           GUI.txtfield.axi.Min = GUI.slider.axi.Min;
           GUI.txtfield.axi.Max = GUI.slider.axi.Max;
           GUI.txtfield.axi.String = GUI.slider.slice.Value;
           
           GUI.txtfield.cor.Min = GUI.slider.cor.Min;
           GUI.txtfield.cor.Max = GUI.slider.cor.Max;
           GUI.txtfield.cor.String = GUI.slider.cor.Value;
           
           GUI.txtfield.sag.Min = GUI.slider.sag.Min;
           GUI.txtfield.sag.Max = GUI.slider.sag.Max;
           GUI.txtfield.sag.String = GUI.slider.sag.Value;
           
           D = GUI.data.DWI(:,:,:,shellInds);
           GUI.data.clim = [max([1, prctile(D(:), 5)]), ...
               max([1, prctile(D(:), 99)])];
           
           colormap(GUI.ax.sag, 'gray');
           colormap(GUI.ax.cor, 'gray');
           colormap(GUI.ax.axiLK, 'gray');
           colormap(GUI.ax.axiLMK, 'gray');
           colormap(GUI.ax.axiLKM, 'gray');
           colormap(GUI.ax.axiLPK, 'gray');
           colormap(GUI.ax.axiLKP, 'gray');       
           
%            histogram(GUI.ax.DWIhistogram, GUI.data.DWI(:));           
%            GUI.ax.DWIhistogram.YScale = 'log';
           
        end       
        
        function updateaxes(GUI, source, event)            
            shell = str2double(GUI.popup.shell.String(GUI.popup.shell.Value, :));
            shellInds = GUI.data.bval == shell;
            shellNums = zeros(size(shellInds));
            shellNums(shellInds) = 1:sum(shellInds);
            
            
            qspaceView = GUI.ax.qspace.View;
            currentDWI = round(GUI.slider.dwi.Value) == shellNums;
            dwiInd = find(currentDWI,1,'first');
            currentAXI = round(GUI.slider.axi.Value);
            currentSAG = round(GUI.slider.sag.Value);
            currentCOR = round(GUI.slider.cor.Value);
            
            limitsDWI = [GUI.slider.dwi.Min, GUI.slider.dwi.Max];
            limitsAXI = [GUI.slider.axi.Min, GUI.slider.axi.Max];
            limitsCOR = [GUI.slider.cor.Min, GUI.slider.cor.Max];
            limitsSAG = [GUI.slider.sag.Min, GUI.slider.sag.Max];
            
            imgSAG = rot90(squeeze(GUI.data.DWI(:, currentSAG, :, currentDWI)));
%             imgCOR = rot90(squeeze(GUI.data.DWI(currentCOR, end:-1:1, :, currentDWI)));
%             imgAXI = squeeze(GUI.data.DWI(:,end:-1:1,currentAXI, currentDWI));
            imgCOR = rot90(squeeze(GUI.data.DWI(currentCOR, :, :, currentDWI)));
            imgAXI = squeeze(GUI.data.DWI(:,:,currentAXI, currentDWI));

            if dwiInd > 1 && shellInds(dwiInd-1)
                dwiIndM = dwiInd -1;
            else
                dwiIndM = [];
            end
            if dwiInd < size(GUI.data.DWI,4) && shellInds(dwiInd+1)
                dwiIndP = dwiInd+1;
            else
                dwiIndP = [];
            end
%             imgAXI_LMK = squeeze(GUI.data.DWI(:,end:-1:1,currentAXI, dwiIndM));
%             imgAXI_LMK = squeeze(GUI.data.DWI(:,end:-1:1,currentAXI, GUI.data.nearestPoints(dwiInd,2)));
            imgAXI_LMK = squeeze(GUI.data.DWI(:,:,currentAXI, GUI.data.nearestPoints(dwiInd,2)));
            if isempty(imgAXI_LMK)
                imgAXI_LMK = NaN(size(imgAXI_LMK,1), size(imgAXI_LMK,2),1);
            end
%             imgAXI_LPK = squeeze(GUI.data.DWI(:,end:-1:1,currentAXI, dwiIndP));            
%             imgAXI_LPK = squeeze(GUI.data.DWI(:,end:-1:1,currentAXI, GUI.data.nearestPoints(dwiInd,3)));
            imgAXI_LPK = squeeze(GUI.data.DWI(:,:,currentAXI, GUI.data.nearestPoints(dwiInd,3)));
            if isempty(imgAXI_LPK)
                imgAXI_LPK = NaN(size(imgAXI_LPK,1), size(imgAXI_LPK,2),1);
            end
%             imgAXI_LKM = squeeze(GUI.data.DWI(:,end:-1:1,max([currentAXI-1, GUI.slider.axi.Min]), currentDWI));
%             imgAXI_LKP = squeeze(GUI.data.DWI(:,end:-1:1,min([currentAXI+1, GUI.slider.axi.Max]), currentDWI));
            imgAXI_LKM = squeeze(GUI.data.DWI(:,:,max([currentAXI-1, GUI.slider.axi.Min]), currentDWI));
            imgAXI_LKP = squeeze(GUI.data.DWI(:,:,min([currentAXI+1, GUI.slider.axi.Max]), currentDWI));            
            if GUI.togglebtn.mask.Value
                imgSAG = imfuse(imgSAG, rot90(squeeze(GUI.data.mask(:, currentSAG, :))));
%                 imgCOR = imfuse(imgCOR, rot90(squeeze(GUI.data.mask(currentCOR, end:-1:1, :))));
%                 imgAXI = imfuse(imgAXI, squeeze(GUI.data.mask(:,end:-1:1,currentAXI)));
%                 imgAXI_LMK = imfuse(imgAXI_LMK, squeeze(GUI.data.mask(:,end:-1:1,currentAXI)));
%                 imgAXI_LPK = imfuse(imgAXI_LPK, squeeze(GUI.data.mask(:,end:-1:1,currentAXI)));
%                 imgAXI_LKM = imfuse(imgAXI_LKM, squeeze(GUI.data.mask(:,end:-1:1,max([currentAXI-1, GUI.slider.axi.Min]))));
%                 imgAXI_LKP = imfuse(imgAXI_LKP, squeeze(GUI.data.mask(:,end:-1:1,min([currentAXI+1, GUI.slider.axi.Max]))));
                imgCOR = imfuse(imgCOR, rot90(squeeze(GUI.data.mask(currentCOR, :, :))));
                imgAXI = imfuse(imgAXI, squeeze(GUI.data.mask(:,:,currentAXI)));
                imgAXI_LMK = imfuse(imgAXI_LMK, squeeze(GUI.data.mask(:,:,currentAXI)));
                imgAXI_LPK = imfuse(imgAXI_LPK, squeeze(GUI.data.mask(:,:,currentAXI)));
                imgAXI_LKM = imfuse(imgAXI_LKM, squeeze(GUI.data.mask(:,:,max([currentAXI-1, GUI.slider.axi.Min]))));
                imgAXI_LKP = imfuse(imgAXI_LKP, squeeze(GUI.data.mask(:,:,min([currentAXI+1, GUI.slider.axi.Max]))));
            end
            
            GUI.images.sag.CData = flipud(imgSAG);
            set(GUI.ax.sag, 'Ydir', 'normal', 'Clim', GUI.data.clim);  
            GUI.lines.sagAXI.XData = limitsCOR;
            GUI.lines.sagAXI.YData = currentAXI.*[1,1]+0.5;
            GUI.lines.sagCOR.XData = currentCOR.*[1,1]+0.5;
            GUI.lines.sagCOR.YData = limitsAXI;

            GUI.images.cor.CData = flipud(imgCOR);
            set(GUI.ax.cor, 'Ydir', 'normal', 'Clim', GUI.data.clim);  
            GUI.lines.corAXI.XData = limitsSAG;
            GUI.lines.corAXI.YData = currentAXI.*[1,1]+0.5;
            GUI.lines.corSAG.XData = currentSAG.*[1,1]+0.5;
            GUI.lines.corSAG.YData = limitsAXI;
            
            GUI.images.axiLK.CData = imgAXI;
            set(GUI.ax.axiLK, 'Clim', GUI.data.clim);  
            GUI.lines.axiCOR.XData = limitsCOR;
            GUI.lines.axiCOR.YData = currentCOR.*[1,1]+0.5;
            GUI.lines.axiSAG.XData = currentSAG.*[1,1]+0.5;
            GUI.lines.axiSAG.YData = limitsSAG;                        
            
            GUI.images.axiLKM.CData = imgAXI_LKM;
            set(GUI.ax.axiLKM, 'Clim', GUI.data.clim);             
            
            GUI.images.axiLKP.CData = imgAXI_LKP;
            set(GUI.ax.axiLKP, 'Clim', GUI.data.clim); 
            
            GUI.images.axiLMK.CData = imgAXI_LMK;
            set(GUI.ax.axiLMK, 'Clim', GUI.data.clim);             
            
            GUI.images.axiLPK.CData = imgAXI_LPK;
            set(GUI.ax.axiLPK, 'Clim', GUI.data.clim);                        
                
%             GUI.lines.dwiHistogramHor.XData
%             GUI.lines.dwiHistogramHor.XData = GUI.data.clim;
%             GUI.lines.dwiHistogramHor.YData = 0.9.*GUI.ax.DWIhistogram.Ylim(2).*[1,1];
%             GUI.lines.dwiHistogramVer.XData = mean(GUI.data.clim).*[1,1];
%             GUI.lines.dwiHistogramVer.YData = GUI.ax.DWIhistogram.Ylim(2).*[1,1];

            %% Spherical interpolation
            % first rotate all bvecs so that the first direction is on
            % z-axis so interpolation near polar region works better.
            % Everything must be rotated back.
            
            bvec = GUI.data.bvec;           
            i = find(GUI.data.bval > 0, 1, 'first');
            tmp = GUI.data.modZ(currentAXI,shellInds); % modZ values for specific SLICE
            
%             bvec = load('test\003.bvec')'
%             i = find(b > 0, 1, 'first')
%             tmp = load('test\modZ.mat');
%             tmp = tmp.modZ(17,shellInds);
            
            v0 = [0,0,1];
            v1 = bvec(i,:);
            v = cross(v0,v1);
            s = sqrt(sum(v.^2));
            c = dot(v0,v1);
            V = [0,    -v(3),  v(2);...
                 v(3),  0,    -v(1);...
                -v(2),  v(1),  0];
            R = eye(3,3) + V + V^2 * (1-c)/s^2;
            bvec = (R'*bvec')';
            
            
            
            %%
            
%             tmp = modZ(40,:);
            tmp = [tmp, tmp];            
            [az, el] = cart2sph([bvec(shellInds,1); -bvec(shellInds,1)],...
                [bvec(shellInds,2); -bvec(shellInds,2)],...
                [bvec(shellInds,3); -bvec(shellInds,3)]);
            [~, inds] = sort(el);
            tmp = tmp(inds);
            az = az(inds);
            el = el(inds);
            % copy polar values to avoid twisting
            m = length(az);
            n = 21;
            azimuth = az;
            elevation = el;
            azimuth(n+1:m+n) = az;
            azimuth(1:n) = 2*pi*(1:-1/n:1/n)-pi;
            azimuth(m+n+1:m+2*n) = 2*pi*(1:-1/n:1/n)-pi;
            elevation(n+1:m+n) = el;
            elevation(1:n) = el(1).*ones(1,n);
            elevation(m+n+1:m+2*n) = el(end).*ones(1,n);
            tmp(n+1:m+n) = tmp;
            tmp(1:n) = tmp(1);
            tmp(m+n+1:m+2*n) = tmp(end);
            % mirror azimuth & elevation to interpolate edges correctly
            azimuth = repmat(azimuth, [3,1]);
            elevation = repmat(elevation, [3,1]);
            tmp = repmat(tmp', [3,1]);
            azimuth(1:m+2*n) = azimuth(1:m+2*n) - 2*pi;
            azimuth(2*m+4*n+1:3*m+6*n) = azimuth(2*m+4*n+1:3*m+6*n) + 2*pi;
            
            azimuth_lin = linspace(0, 2*pi, 3*m+3*n);
            elevation_lin = linspace(-pi/2, pi/2, 3*m+3*n);
            [azimuth_mat, elevation_mat] = meshgrid(azimuth_lin, elevation_lin);
            [Xi,Yi,Zi] = sph2cart(azimuth_mat, elevation_mat, ones(size(azimuth_mat)));                                                
            
            f = scatteredInterpolant(azimuth, elevation, double(tmp), 'natural', 'linear');
            Ci = f(azimuth_mat, elevation_mat);
            Ci(isnan(Ci)) = 0;
            
            % Get non-rotated coordinates
            v = R*[Xi(:)'; Yi(:)'; Zi(:)'];
            bvec = GUI.data.bvec;
            Xi(:) = v(1,:);
            Yi(:) = v(2,:);
            Zi(:) = v(3,:);                                          
            
%             figure
%             surf(azimuth_mat, elevation_mat, Ci, 'edgecolor', 'none', 'facecolor', 'interp')
%             view(2)
            %%

            plot3(GUI.ax.qspace,[bvec(shellInds,1); -bvec(shellInds,1)],...
                [bvec(shellInds,2); -bvec(shellInds,2)],...
                [bvec(shellInds,3); -bvec(shellInds,3)], 'ro');
            hold(GUI.ax.qspace, 'on')
            
            plot3(GUI.ax.qspace,[bvec(dwiInd,1); -bvec(dwiInd,1)], [bvec(dwiInd,2); -bvec(dwiInd,2)], [bvec(dwiInd,3); -bvec(dwiInd,3)], 'gx', 'markersize', 24);            
            dwiInd2 = GUI.data.nearestPoints(dwiInd,2);
            dwiInd3 = GUI.data.nearestPoints(dwiInd,3);
            plot3(GUI.ax.qspace,[bvec(dwiInd2,1); -bvec(dwiInd2,1)], [bvec(dwiInd2,2); -bvec(dwiInd2,2)], [bvec(dwiInd2,3); -bvec(dwiInd2,3)], 'c.', 'markersize', 20); 
            plot3(GUI.ax.qspace,[bvec(dwiInd3,1); -bvec(dwiInd3,1)], [bvec(dwiInd3,2); -bvec(dwiInd3,2)], [bvec(dwiInd3,3); -bvec(dwiInd3,3)], 'm.', 'markersize', 20);
            surf(GUI.ax.qspace,0.98*Xi, 0.98*Yi, 0.98*Zi, Ci, 'EdgeColor', 'none');                        
            hold(GUI.ax.qspace, 'off')
            shading(GUI.ax.qspace, 'interp');
            axis(GUI.ax.qspace, 'square');
%             legend(GUI.ax.qspace, 'Gradient directions', 'Current direction', 'Location', 'North')
            xlabel(GUI.ax.qspace, 'X');
            ylabel(GUI.ax.qspace, 'Y');
            zlabel(GUI.ax.qspace, 'Z');
            grid(GUI.ax.qspace, 'on');
            GUI.ax.qspace.CLim = [GUI.data.thresholdLower, GUI.data.thresholdUpper];
            GUI.ax.qspace.XColor = [0.5, 0.5, 0.5];
            GUI.ax.qspace.YColor = [0.5, 0.5, 0.5];
            GUI.ax.qspace.ZColor = [0.5, 0.5, 0.5];
            
% view -  The azimuth, az, is the horizontal rotation about the z-axis as
% measured in degrees from the negative y-axis. Positive values indicate 
% counterclockwise rotation of the viewpoint. el is the vertical elevation 
% of the viewpoint in degrees. Positive values of elevation correspond to
% moving above the object; negative values correspond to moving below the object.
%
% cart2sph - Azimuth angle, returned as an array. azimuth is the 
% counterclockwise angle in the x-y plane measured in radians from the
% positive x-axis. The value of the angle is in the range [-pi pi].
% Elevation angle, returned as an array. elevation is the elevation angle 
% in radians from the x-y plane.


            [az, el, ~] = cart2sph(bvec(dwiInd, 1), bvec(dwiInd, 2), bvec(dwiInd, 3));
            az = az./pi.*180;
            el = el./pi.*180;
            az = az - 90;
            az(az < 0) = az + 360;
            el = -el;            
%             GUI.ax.qspace.View = qspaceView;
            GUI.ax.qspace.View = [az, el];

            %%            

            modZclim = [GUI.data.thresholdLower, GUI.data.thresholdUpper];
            GUI.images.modZ2D.CData = GUI.data.modZ(:, shellInds);
            set(GUI.ax.modZ2D, 'YDir', 'normal', 'CLim', modZclim);
            axis(GUI.ax.modZ2D, 'tight');
%             GUI.lines.modZAXI.XData = [find(shellInds, 1, 'first'), find(shellInds,1, 'last')]; %limitsDWI;
            GUI.lines.modZAXI.XData = [0.5, sum(shellInds)+0.5];
            GUI.lines.modZAXI.YData = currentAXI.*[1,1];            
%             GUI.lines.modZDWI.XData = round(GUI.slider.dwi.Value).*[1,1]; %dwiInd.*[1,1];%currentDWI.*[1,1];
            GUI.lines.modZDWI.XData = shellNums(dwiInd).*[1,1];
            GUI.lines.modZDWI.YData = limitsAXI+[-0.5, 0.5];
            
%             shellInds
%             shellNums
%             [dwiInd, dwiInd2, dwiInd3]
%             [shellNums(dwiInd), shellNums(dwiInd2), shellNums(dwiInd3)]
            GUI.lines.modZG.XData = shellNums(dwiInd2);
            GUI.lines.modZG.YData = currentAXI;
            GUI.lines.modZM.XData = shellNums(dwiInd3);
            GUI.lines.modZM.YData = currentAXI;

            
            bins = 0.1:0.2:modZclim(2);
            tmp = GUI.data.modZ(:, shellInds);
            hh = histogram(GUI.ax.modZhistogram, tmp(:), bins);
            axes(GUI.ax.modZhistogram)
            axis(GUI.ax.modZhistogram, 'tight');          
            line(GUI.data.modZ(currentAXI, currentDWI).*[1,1], [1,max(hh.Values)], 'color', 'black', 'linewidth', 2);
            line(modZclim(1).*[1,1], [1,max(hh.Values)], 'color', 'green', 'linestyle', '--');
            line(modZclim(2).*[1,1], [1,max(hh.Values)], 'color', 'red', 'linestyle', '--');                          
%             legend(GUI.ax.modZhistogram.Children([1,2,3]), 'Upper threshold', 'Lower threshold', 'Current slice');
            legend(GUI.ax.modZhistogram.Children([1,2,3]), 'Upper threshold', 'Lower threshold', ['Current ', num2str(GUI.data.modZ(currentAXI, currentDWI),2)]);
            GUI.ax.modZhistogram.XColor = [0.9 0.9 0.9];
            GUI.ax.modZhistogram.YColor = [0.9 0.9 0.9];
            GUI.ax.modZhistogram.XLim = [-0.5, max([modZclim(2), GUI.data.modZ(currentAXI, currentDWI)])+1];
        end
        
        function calculateModZ(GUI, source, event)
            DWI = GUI.data.DWI;
            DWI(DWI < 0) = NaN; % User can opt to use their own mask by setting non-brain voxels to negative numbers.
            b = GUI.data.bval;
            metric = GUI.data.metric;
            
            if GUI.data.useMask
                mask = GUI.data.mask;
                DWI(repmat(~mask, [1,1,1,size(DWI,4)])) = NaN;               
            end
            
            dims = size(DWI);
            modZ = NaN(dims(3:4), 'double');
            uniqueb = unique(b);

            for i = 1 : length(uniqueb)
                shell = (b == uniqueb(i));
                N = sum(shell);
%                 if ( uniqueb(i) == 0 )
%                     continue;
%                 end
                if ( N < 2 )
                    disp(['SOLID: Not enough data with b-value ', num2str(uniqueb(i))]);
                    continue;
                end                
                baseline = double(DWI(:,:,:,shell));
                baseline(baseline < eps) = NaN;
                baseline = reshape(baseline, [dims(1)*dims(2), dims(3), N]);
                switch metric
                    case 'Variance'
                        y = squeeze(var(baseline, 'omitnan'));
                    case 'Mean'
                        y = squeeze(mean(baseline, 'omitnan'));
                    case 'Iod'
                        y = squeeze(var(baseline, 'omitnan'))./squeeze(mean(baseline, 'omitnan'));
                end
                tmp = repmat(median(y,2,'omitnan'), [1,N]);
                MAD = 1.4826 * median( abs(y-tmp),2,'omitnan');
                modZ(:,shell) = abs(y-tmp)./repmat(MAD,[1,N]);                
            end
            
            GUI.data.uniqueb = uniqueb;
            GUI.data.modZ = modZ;
            GUI.data.modZorig = modZ;
            
        end
        
        function loadData(GUI, source, event)           
            
            [fname, fdir] = uigetfile({'*.nii;*.nii.gz'}, 'Select 4D DWI NIfTI file');
            
            try
                DWI = niftiread([fdir, filesep, fname]);
            catch
                try
                    DWI = load_untouch_nii([fdir, filesep, fname]);
                    DWI = DWI.img;
                catch
                    e = errordlg('*.nii file not found', 'File Error');
                    return
                end
            end
            
            if any(size(DWI) < 2)
                e = errordlg('Error in 4D DWI NIfTI image dimensions', 'File error');
                return
            end
            
            try
                bval = load([fdir, filesep, fname(1:end-4), '.bval']);
            catch
                try 
                    [bname, bdir] = uigetfile(fullfile(fdir,'*.bval'), 'Select bval file');
                    bval = load([bdir, filesep, bname]);
                catch
                    e = errordlg('*.bval not found', 'Bval Error');
                    return
                end
            end
            
            try
                bvec = load([fdir, filesep, fname(1:end-4), '.bvec']);
            catch
                try 
                    [vname, vdir] = uigetfile(fullfile(fdir,'*.bvec'), 'Select bvec file');
                    bvec = load([vdir, filesep, vname]);
                catch
                    e = errordlg('*.bvec not found', 'Bvec Error');
                    return
                end
            end
            
            if size(bvec,1) < size(bvec,2)
                bvec = bvec';
            end
            
            img = DWI;
            img = permute(img, [2 1 3 4]);
%             img = flipdim(img,2);
%             img = flipdim(img,1);            
            
            try
                mask = E_DTI_Create_Mask_From_DWI_enhanced_IND(img(:,:,:,1), 0.5, 7); 
            catch
                mask = ones(size(img(:,:,:,1)));
            end
                                 
            bval = round(bval/100)*100;
  
            GUI.data.DWI = img;
            GUI.data.bval = bval;
            GUI.data.bvec = bvec;
            GUI.data.mask = mask;
            GUI.data.minmax = [min(DWI(:)), max(DWI(:))];            
            GUI.data.fname = fname;
            GUI.data.fpath = fdir;
            
            calculateAngularNHood(GUI);
            calculateModZ(GUI);
            GUI.popup.shell.String = GUI.data.uniqueb;
            GUI.popup.shell.Value = 1;
            
            initializeAxesSlider(GUI);
            updateaxes(GUI);
        end        
        
        function calculateAngularNHood(GUI, source, event)
            b = GUI.data.bval;
            v = GUI.data.bvec;
            % normalize v
            n = sqrt(sum(v.^2,2));
            v = v./repmat(n, [1,3]);
            v(isnan(v)) = 0;
            ad = inf(length(b), length(b));
            
            for i = 1:length(b)
                for j = i:length(b)
                    vi = v(i, :);
                    vj = v(j, :);
                    c = cross(vj, vi);
                    d = abs(dot(vj, vi));
                    nc = sqrt(sum(abs(c).^2));
                    euc = sqrt(sum( (b(j).^2.*vj-b(i).^2.*vi).^2))./b(i);
                    ang = atan2(nc, d)./pi;
                    ad(i,j) = euc;
                    % ang = atan2(nc, d)./pi.*180 + e;
%                     ad(i,j) = sqrt(euc.^2 + ang.^2);
%                     [i, j, euc, ang, ad(i,j)]
%                     pause
                end
            end
            ad(isnan(ad)) = 0;
            ad = triu(ad)+triu(ad,1)';
            [ad2, tmp] = sort(ad,2);
            [~, GUI.data.nearestPoints] = sort(ad,2);
            
        end
        
        function saveOutliers(GUI, source, event)
            i = strfind(GUI.data.fname, '.nii');
            oname = GUI.data.fname(1:i-1);
            oname = strcat(GUI.data.fpath, filesep, oname, '_L_', ...
                num2str(GUI.data.thresholdLower), '_U_', ...
                num2str(GUI.data.thresholdUpper), '_', ...
                GUI.data.metric, '_masked', num2str(GUI.data.useMask));
            %% save txt file
            fid = fopen(strcat(oname, '_modZ2D.txt'), 'w');
            for i = 1:size(GUI.data.modZ, 1)
                fprintf(fid, repmat('%f, ', [1, size(GUI.data.modZ,2)-1]), GUI.data.modZ(i, 1:end-1));
                fprintf(fid, '%f\n', GUI.data.modZ(i, end));
            end
            fclose(fid);
            
            %% save 4D nifti
            flag = 0;
            try
                info = niftiinfo(strcat(GUI.data.fpath, filesep, GUI.data.fname));
                img = niftiread(info);                
                flag = 1;
            catch
                try
                    nii = load_untouch_nii(strcat(GUI.data.fpath, filesep, GUI.data.fname));                    
                    img = zeros(size(nii.img));
                    flag = 2;
                catch
                    e = errordlg('*.nii file not found', 'File Error');
                    return
                end
            end
            
            for k = 1:size(GUI.data.modZ,1)
                for l = 1:size(GUI.data.modZ,2)
                    img(:,:,k,l) = repmat(GUI.data.modZ(k,l), [size(img,1), size(img,2), 1,1]);
                end
            end
            
            if flag == 1
                info.Description = 'SOLID modified Z-score map';
                niftiwrite(img, strcat(oname, '_SOLID_map.nii'), info);
            elseif flag == 2
                nii.img = img;
                nii.hdr.hist.descrip = 'SOLID modified Z-score map';
                nii.fileprefix = strcat(oname, '_SOLID_map');
                save_untouch_nii(nii, strcat(oname, '_SOLID_map.nii'));
            end
            
            %% save PNG per shell
            b = GUI.data.bval;
            uniqueb = unique(b);
            for i = 1 : length(uniqueb)
                shell = (b == uniqueb(i));
                modZ = GUI.data.modZ(:,shell);
                modZ(isnan(modZ)) = 0;
                fig = figure('Name', 'SOLID', 'visible','off');
                ax = axes;
                im = imagesc(modZ);
                xlabel('Image volume #');
                ylabel('Slice #');
                cb = colorbar;
                ylabel(cb, 'Modified Z-score');
                title(['SOLID, b-value' num2str(uniqueb(i))]);
                ax.FontSize = 14;
                ax.CLim = [GUI.data.thresholdLower, GUI.data.thresholdUpper];
                axis equal tight;
                print(fig, strcat(oname, '_modZ2D_b', num2str(uniqueb(i)), '.png'), '-dpng', '-r300');
                close(fig);
            end
        end        
        
        function mouseClick(GUI, source, event)            
            C = get(GUI.ax.modZ2D, 'CurrentPoint');
            curDWI = round(C(1,1));
            curSlice = round(C(1,2));
            
            if curDWI >= GUI.ax.modZ2D.XLim(1) && curDWI <= GUI.ax.modZ2D.XLim(2) && ...
                curSlice >= GUI.ax.modZ2D.YLim(1) && curSlice <= GUI.ax.modZ2D.YLim(2)
                
%                 shell = str2double(GUI.popup.shell.String(GUI.popup.shell.Value, :));
%                 shellInds = GUI.data.bval == shell;
%                 shellNums = zeros(size(shellInds));
%                 shellNums(shellInds) = 1:sum(shellInds);
                
%                 curDWI = shellNums(curDWI);
                if curDWI
                    GUI.slider.dwi.Value = curDWI;
                    GUI.txtfield.dwi.String = curDWI;
                    GUI.slider.axi.Value = curSlice;
                    GUI.slider.slice.Value = curSlice;
                    GUI.txtfield.axi.String = curSlice;
                    updateaxes(GUI);
                end
            else
                % Change window / level
                C = get(GUI.figureUI, 'CurrentPoint');
                GUI.data.window.mousePosI = [C(1,1), C(1,2)];
                clim = GUI.ax.axiLK.CLim;
                GUI.data.window.width = clim(2)-clim(1);
                GUI.data.window.level = GUI.data.window.width/2 + clim(1);                
                start(GUI.data.window.timer);
            end
        end
        
        function mouseRelease(GUI, source, event)
            stop(GUI.data.window.timer);
        end
        
        function mouseScroll(GUI, source, event)            
            C = get(GUI.figureUI, 'CurrentPoint');
            x = C(1,1);
            y = C(1,2);
%             [x,y]            
            P = GUI.ax.cor.Position;
            if x > P(1) && x < P(1)+P(3) && y > P(2) && y < P(2)+P(4)
%                 'add cor slice'
                try
                    GUI.slider.cor.Value = GUI.slider.cor.Value - event.VerticalScrollCount;
                    GUI.txtfield.cor.String = GUI.slider.cor.Value;
                    updateaxes(GUI);
                catch
                end
            end
            P = GUI.ax.sag.Position;
            if x > P(1) && x < P(1)+P(3) && y > P(2) && y < P(2)+P(4)
%                 'add sag slice'
                try
                GUI.slider.sag.Value = GUI.slider.sag.Value - event.VerticalScrollCount;
                GUI.txtfield.sag.String = GUI.slider.sag.Value;
                updateaxes(GUI);
                catch
                end
            end
            P = GUI.ax.axiLK.Position;
            if x > P(1) && x < P(1)+P(3) && y > P(2) && y < P(2)+P(4)
%                 'add axi slice'
                try                   
                GUI.slider.axi.Value = GUI.slider.axi.Value - event.VerticalScrollCount;
                GUI.slider.slice.Value = GUI.slider.axi.Value;
                GUI.txtfield.axi.String = GUI.slider.axi.Value;
                updateaxes(GUI);
                catch
                end
            end
            
            
        end       
        
        function windowTimerCallback(GUI, source, event)
%             C = get(GUI.figureUI, 'CurrentPoint');
%             x = C(1,1);
%             y = C(1,2);
%             
%             dx = x-GUI.data.window.mousePosI(1);
%             dy = y-GUI.data.window.mousePosI(2);
%             
%             [x,y, dx,dy]
%             
%             clim = GUI.data.clim;
%             width = clim(2)-clim(1);
%             level = width/2 + clim(1) + dx;
%             width = width + dy;
%             
%             GUI.data.clim = [level-width/2, level+width/2];
% %             [level-width/2, level+width/2]
%                         
%             GUI.ax.sag.CLim = GUI.data.clim;
% %             GUI.ax.cor.Clim = GUI.data.clim;
%             
        end
        
        function mouseMove(GUI, source, event)
            if strcmp(GUI.data.window.timer.running, 'on')                
                
                clim = GUI.data.clim;
                level = mean(clim);
                width = clim(2)-clim(1);
                
                C = get(GUI.figureUI, 'CurrentPoint');
                x = C(1,1);
                y = C(1,2);
                dx = (x-GUI.data.window.mousePosI(1))/3+1;
                dy = (y-GUI.data.window.mousePosI(2))/3+1;                
                 
                
                level = level .* dx;
                width = width .* dy;
                
                if level > GUI.data.minmax(2)
                    level = GUI.data.minmax(2);
                end
                if level < GUI.data.minmax(1)
                    level = GUI.data.minmax(1);
                end
                if width < eps
                    width = eps;
                end
                
                clim = [level - width/2, level + width/2]
                                
                GUI.data.clim = clim;
                                
                GUI.ax.sag.CLim = clim;
                GUI.ax.cor.CLim = clim;
                GUI.ax.axiLK.CLim = clim;
                GUI.ax.axiLMK.CLim = clim;
                GUI.ax.axiLKM.CLim = clim;
                GUI.ax.axiLPK.CLim = clim;
                GUI.ax.axiLKP.CLim = clim;
                
%                 GUI.lines.dwiHistogramHor.XData = GUI.data.clim;
%                 GUI.lines.dwiHistogramHor.YData = 0.9.*GUI.ax.DWIhistogram.Ylim(2).*[1,1];
%                 GUI.lines.dwiHistogramVer.XData = mean(GUI.data.clim).*[1,1];
%                 GUI.lines.dwiHistogramVer.YData = GUI.ax.DWIhistogram.Ylim(2).*[1,1];
                
            end
        end
        
        function keyPress(GUI, source, event)
            % shell list must be unselected !!
            GUI.popup.shell.Enable = 'off';
            GUI.slider.cor.Enable = 'off';
            GUI.slider.sag.Enable = 'off';
            GUI.slider.axi.Enable = 'off';
            GUI.slider.dwi.Enable = 'off';
            pause(0.01);
            if strcmp(event.Key, 'space')
                shell = str2double(GUI.popup.shell.String(GUI.popup.shell.Value, :));
                shellInds = GUI.data.bval == shell;
                shellNums = zeros(size(shellInds));
                shellNums(shellInds) = 1:sum(shellInds);
                curDWI = GUI.slider.dwi.Value;
                curDWI = find(shellNums == curDWI);
                curSlice = GUI.slider.axi.Value;
                z = GUI.data.modZ(curSlice, curDWI);
                zorig = GUI.data.modZorig(curSlice, curDWI);
                if abs(z - zorig) > eps
                    GUI.data.modZ(curSlice, curDWI) = GUI.data.modZorig(curSlice, curDWI);
                else
                    GUI.data.modZ(curSlice, curDWI) = 0;
                end
                updateaxes(GUI);  
            end
            if strcmp(event.Key, 'rightarrow')
                shell = str2double(GUI.popup.shell.String(GUI.popup.shell.Value, :));
                shellInds = GUI.data.bval == shell;
                shellNums = zeros(size(shellInds));
                shellNums(shellInds) = 1:sum(shellInds);
                curDWI = GUI.slider.dwi.Value;
                curDWI = find(shellNums == curDWI);
                maxDWI = find(shellNums == GUI.slider.dwi.Max);
                curSlice = GUI.slider.axi.Value;
                n = 0;
                flag = 0;
                for d = curDWI : 1 : maxDWI
                    if n == 0
                        cs = curSlice + 1;
                        if cs == GUI.slider.axi.Max
                            cs = GUI.slider.axi.Max;
                        end                        
                        n = 1;
                    else
                        cs = 1;                        
                    end
                    for s = cs : 1 : GUI.slider.axi.Max
                        z = GUI.data.modZ(s, d);
                        if z > GUI.data.thresholdLower
                            flag = 1;
                            GUI.slider.dwi.Value = shellNums(d);
                            GUI.txtfield.dwi.String = shellNums(d);
                            GUI.slider.axi.Value = s;
                            GUI.slider.slice.Value = s;
                            GUI.txtfield.axi.String = s;
                            break
                        end
                    end
                    if flag
                        break
                    end
                end
                updateaxes(GUI);  
            end
            if strcmp(event.Key, 'leftarrow')
                shell = str2double(GUI.popup.shell.String(GUI.popup.shell.Value, :));
                shellInds = GUI.data.bval == shell;
                shellNums = zeros(size(shellInds));
                shellNums(shellInds) = 1:sum(shellInds);
                curDWI = GUI.slider.dwi.Value;
                curDWI = find(shellNums == curDWI);
                minDWI = find(shellNums == GUI.slider.dwi.Min);
                curSlice = GUI.slider.axi.Value;
                n = 0;
                flag = 0;
                for d = curDWI : -1 : minDWI
                    if n == 0
                        cs = curSlice - 1;
                        if cs == 0
                            cs = 1;
                        end
                        n = 1;
                    else
                        cs = GUI.slider.axi.Max;                        
                    end
                    for s = cs : -1 : 1
                        z = GUI.data.modZ(s, d);
                        if z > GUI.data.thresholdLower
                            flag = 1;
                            GUI.slider.dwi.Value = shellNums(d);
                            GUI.txtfield.dwi.String = shellNums(d);
                            GUI.slider.axi.Value = s;
                            GUI.slider.slice.Value = s;
                            GUI.txtfield.axi.String = s;
                            break
                        end
                    end
                    if flag
                        break
                    end
                end
                updateaxes(GUI); 
            end
            GUI.popup.shell.Enable = 'on';
            GUI.slider.cor.Enable = 'on';
            GUI.slider.sag.Enable = 'on';
            GUI.slider.axi.Enable = 'on';
            GUI.slider.dwi.Enable = 'on';
            pause(0.01);
        end
        
        % Create UIFigure and components
        function GUI = initializeGUI(GUI)
        
            % Create SOLIDUIFigure
            GUI.figureUI = figure;
            GUI.figureUI.Color = [0 0 0];
            tmp = get(0, 'MonitorPositions');
            tmp = tmp(1,:);
            tmp = [tmp(3)*0.1 tmp(4)*0.1 tmp(3)*0.8 tmp(4)*0.8];
            GUI.figureUI.Position = tmp;
            GUI.figureUI.Name = 'SOLID';
            GUI.figureUI.Units = 'normalized';
%             GUI.figureUI.MenuBar = 'none';
%             GUI.figureUI.ToolBar = 'none';  
            GUI.figureUI.WindowButtonDownFcn = @GUI.mouseClick;
            GUI.figureUI.WindowButtonUpFcn = @GUI.mouseRelease;
            GUI.figureUI.WindowScrollWheelFcn = @GUI.mouseScroll;
            GUI.figureUI.WindowButtonMotionFcn = @GUI.mouseMove;
            GUI.figureUI.WindowKeyPressFcn = @GUI.keyPress;
            
            % Create FileMenu
            GUI.fileMenu = uimenu(GUI.figureUI);
            if strcmp(version('-release'), '2015a')
                GUI.fileMenu.Label = 'SOLID';
            else
                GUI.fileMenu.Text = 'SOLID';
            end
            
            % Create FileMenu -> Open Nifti
            GUI.fileMenu_Open = uimenu(GUI.fileMenu);
            if strcmp(version('-release'), '2015a')
                GUI.fileMenu_Open.Label = 'Open';
            else                
                GUI.fileMenu_Open.Text = 'Open';
            end
            GUI.fileMenu_Open.Accelerator = 'Q';
            if strcmp(version('-release'), '2015a')
                GUI.fileMenu_Open.Callback = @GUI.loadData;
            else
                GUI.fileMenu_Open.MenuSelectedFcn = @GUI.loadData;
            end
            
            % Create FileMenu -> Save outlier map
            GUI.fileMenu_Save = uimenu(GUI.fileMenu);
            if strcmp(version('-release'), '2015a')
                GUI.fileMenu_Save.Label = 'Save results';
            else                
                GUI.fileMenu_Save.Text = 'Save results';
            end            
            GUI.fileMenu_Save.Accelerator = 's';            
            if strcmp(version('-release'), '2015a')
                GUI.fileMenu_Save.Callback = @GUI.saveOutliers;
            else
                GUI.fileMenu_Save.MenuSelectedFcn = @GUI.saveOutliers;
            end
            GUI.fileMenu_Save.Visible = 'on'; 
            
            % Create FileMenu -> Toggle Var/Mean
            GUI.fileMenu_Metric = uimenu(GUI.fileMenu);
            GUI.fileMenu_Metric_Var = uimenu(GUI.fileMenu_Metric);
            GUI.fileMenu_Metric_Mean = uimenu(GUI.fileMenu_Metric);
            GUI.fileMenu_Metric_Iod = uimenu(GUI.fileMenu_Metric);
            if strcmp(version('-release'), '2015a')
                GUI.fileMenu_Metric.Label = 'Select SOLID metric';
                GUI.fileMenu_Metric_Var.Label = 'Variance';
                GUI.fileMenu_Metric_Mean.Label = 'Mean';
                GUI.fileMenu_Metric_Iod.Label = 'Iod';
            else
                GUI.fileMenu_Metric.Text = 'Select SOLID metric';
                GUI.fileMenu_Metric_Var.Text = 'Variance';
                GUI.fileMenu_Metric_Mean.Text = 'Mean';
                GUI.fileMenu_Metric_Iod.Text = 'Iod';
            end           
            if strcmp(version('-release'), '2015a')                
                GUI.fileMenu_Metric_Var.Callback = @GUI.toggleVarMean;
                GUI.fileMenu_Metric_Mean.Callback = @GUI.toggleVarMean;
                GUI.fileMenu_Metric_Iod.Callback = @GUI.toggleVarMean;
            else
                GUI.fileMenu_Metric_Var.MenuSelectedFcn = @GUI.toggleVarMean;
                GUI.fileMenu_Metric_Mean.MenuSelectedFcn = @GUI.toggleVarMean;
                GUI.fileMenu_Metric_Iod.MenuSelectedFcn = @GUI.toggleVarMean;
            end
            GUI.fileMenu_Metric_Var.Checked = 'on';
            
            GUI.fileMenu_UseMask = uimenu(GUI.fileMenu);
            if strcmp(version('-release'), '2015a')
                GUI.fileMenu_UseMask.Label = 'Use brain mask';
            else                
                GUI.fileMenu_UseMask.Text = 'Use brain mask';
            end                         
            if strcmp(version('-release'), '2015a')                
                GUI.fileMenu_UseMask.Callback = @GUI.toggleUseMask;
            else
                GUI.fileMenu_UseMask.MenuSelectedFcn = @GUI.toggleUseMask;
            end            
            GUI.fileMenu_UseMask.Checked = 'on';
            
            % Create axes 3x3 grid DWI
            % (1,1)
            GUI.ax.cor = axes('Parent', GUI.figureUI);
            GUI.ax.cor.Position = [0, 3/4-2/21, 1/4, 1/4];            
            GUI.images.cor = imagesc('Parent', GUI.ax.cor, 'CData', NaN(50,50));
            axis(GUI.ax.cor, 'off', 'tight', 'equal');
            
            % (1,2)
            GUI.ax.axiLKP = axes('Parent', GUI.figureUI);
            GUI.ax.axiLKP.Position = [1/4, 3/4-2/21, 1/4, 1/4];            
            GUI.images.axiLKP = imagesc('Parent', GUI.ax.axiLKP, 'CData', NaN(50,50));
            axis(GUI.ax.axiLKP, 'off', 'tight', 'equal');
            
            % (1,3)
            GUI.ax.sag = axes('Parent', GUI.figureUI);
            GUI.ax.sag.Position = [2/4, 3/4-2/21, 1/4, 1/4];    
            GUI.ax.sag.XDir = 'reverse';
            GUI.images.sag = imagesc('Parent', GUI.ax.sag, 'CData', NaN(50,50));
            axis(GUI.ax.sag, 'off', 'tight', 'equal');
            
            % (2,1)
            GUI.ax.axiLMK = axes('Parent', GUI.figureUI);
            GUI.ax.axiLMK.Position = [0/4, 2/4-2/20, 1/4, 1/4];            
            GUI.images.axiLMK = imagesc('Parent', GUI.ax.axiLMK, 'CData', NaN(50,50));
            axis(GUI.ax.axiLMK, 'tight', 'equal');
            GUI.ax.axiLMK.XColor = [0,1,1];
            GUI.ax.axiLMK.YColor = [0,1,1];
            GUI.ax.axiLMK.Box = 'on';
            GUI.ax.axiLMK.LineWidth = 2;
            GUI.ax.axiLMK.XTick = [];
            GUI.ax.axiLMK.YTick = [];
            
            % (2,2)
            GUI.ax.axiLK = axes('Parent', GUI.figureUI);
            GUI.ax.axiLK.Position = [1/4, 2/4-2/20, 1/4, 1/4];            
            GUI.images.axiLK = imagesc('Parent', GUI.ax.axiLK, 'CData', NaN(50,50));
            axis(GUI.ax.axiLK, 'tight', 'equal');
            GUI.ax.axiLK.XColor = [1,1,1];
            GUI.ax.axiLK.YColor = [1,1,1];
            GUI.ax.axiLK.Box = 'on';
            GUI.ax.axiLK.LineWidth = 2;
            GUI.ax.axiLK.XTick = [];
            GUI.ax.axiLK.YTick = [];
            
            % (2,3)
            GUI.ax.axiLPK = axes('Parent', GUI.figureUI);
            GUI.ax.axiLPK.Position = [2/4, 2/4-2/20, 1/4, 1/4];            
            GUI.images.axiLPK = imagesc('Parent', GUI.ax.axiLPK, 'CData', NaN(50,50));
            axis(GUI.ax.axiLPK, 'tight', 'equal');
            GUI.ax.axiLPK.XColor = [1,0,1];
            GUI.ax.axiLPK.YColor = [1,0,1];
            GUI.ax.axiLPK.Box = 'on';
            GUI.ax.axiLPK.LineWidth = 2;
            GUI.ax.axiLPK.XTick = [];
            GUI.ax.axiLPK.YTick = [];
            
            % (3,2)
            GUI.ax.axiLKM = axes('Parent', GUI.figureUI);
            GUI.ax.axiLKM.Position = [1/4, 1/4-2/19, 1/4, 1/4];            
            GUI.images.axiLKM = imagesc('Parent', GUI.ax.axiLKM, 'CData', NaN(50,50));
            axis(GUI.ax.axiLKM, 'off', 'tight', 'equal');
            
            % Create axes qspace
            GUI.ax.qspace = axes('Parent', GUI.figureUI);
            GUI.ax.qspace.Position = [2/4 1/20 1/4 1/4];
            GUI.images.qspace = imagesc('Parent', GUI.ax.qspace, 'CData', NaN(50,50));
            axis(GUI.ax.qspace, 'square');
%             rotate3d(GUI.ax.qspace, 'on');
            
            % Create axes modZ2D
            GUI.ax.modZ2D = axes('Parent', GUI.figureUI);
            GUI.ax.modZ2D.Position = [3/4+1/100 2/20 1/4-2/30 16/20];
            GUI.images.modZ2D = imagesc('Parent', GUI.ax.modZ2D, 'CData', NaN(50,50));
            axis(GUI.ax.modZ2D, 'tight', 'equal');
            ch = colorbar(GUI.ax.modZ2D);
            ylabel(ch, 'Modified Z-score');
            ch.Color = [0.9 0.9 0.9];            
            GUI.ax.modZ2D.XColor = [0.9 0.9 0.9];
            GUI.ax.modZ2D.YColor = [0.9 0.9 0.9];
            GUI.ax.modZ2D.Box = 'on';
            xlabel(GUI.ax.modZ2D, 'Image volume');
            ylabel(GUI.ax.modZ2D, 'Slice');
            
            % Create axes modZhistogram
            GUI.ax.modZhistogram = axes('Parent', GUI.figureUI);
            GUI.ax.modZhistogram.Position =  [1/20 1/20 2/6-1/10 2/6];
            GUI.images.modZhistogram = imagesc('Parent', GUI.ax.modZhistogram, 'CData', NaN(50,50));
            axis(GUI.ax.modZhistogram, 'tight', 'equal');
            xlabel(GUI.ax.modZhistogram, 'Modified Z-score');
            ylabel(GUI.ax.modZhistogram, 'Counts');
            GUI.ax.modZhistogram.XColor = [0.9 0.9 0.9];
            GUI.ax.modZhistogram.YColor = [0.9 0.9 0.9];

            % Create DWI histogram
            GUI.ax.DWIhistogram = axes('Parent', GUI.figureUI);
            GUI.ax.DWIhistogram.Position = [9/10, 1/150, 1/11, 1/12];
            GUI.images.DWIhistogram = histogram('Parent', GUI.ax.DWIhistogram, NaN(50,1));
            axis(GUI.ax.DWIhistogram, 'off', 'tight');
%             hold(GUI.ax.DWIhistogram, 'on');                       
%             GUI.lines.dwiHistogramVer = line(GUI.ax.DWIhistogram, [1,2], [1,2], 'color', 'red');
%             GUI.lines.dwiHistogramHor = line(GUI.ax.DWIhistogram, [1,2], [1,2], 'color', 'red');
%             hold(GUI.ax.DWIhistogram, 'off');                       
            
            % Create sliders
            GUI.slider.cor = uicontrol('Parent', GUI.figureUI, 'style', 'slider');
            GUI.slider.cor.BackgroundColor = 'green';
            GUI.slider.cor.Units = 'normalized';
            GUI.slider.cor.Position = [0+1/100, 1-1/20+1/80, 1/4-2/100, 1/40];
            GUI.slider.cor.Callback = @GUI.sliderCorCallback;
            
            GUI.slider.sag = uicontrol('Parent', GUI.figureUI, 'style', 'slider');
            GUI.slider.sag.BackgroundColor = 'red';
            GUI.slider.sag.Units = 'normalized';
            GUI.slider.sag.Position = [1/4+1/100, 1-1/20+1/80, 1/4-2/100, 1/40];
            GUI.slider.sag.Callback = @GUI.sliderSagCallback;
            
            GUI.slider.axi = uicontrol('Parent', GUI.figureUI, 'style', 'slider');
            GUI.slider.axi.BackgroundColor = 'blue';
            GUI.slider.axi.Units = 'normalized';
            GUI.slider.axi.Position = [2/4+1/100, 1-1/20+1/80, 1/4-2/100, 1/40];
            GUI.slider.axi.Callback = @GUI.sliderAxiCallback;
            
            GUI.slider.dwi = uicontrol('Parent', GUI.figureUI, 'style', 'slider');
            GUI.slider.dwi.BackgroundColor = 'white';
            GUI.slider.dwi.Units = 'normalized';
            GUI.slider.dwi.Position = [3/4+1/100, 1-1/20+1/80, 1/4-2/30, 1/40];
            GUI.slider.dwi.Callback = @GUI.sliderDwiCallback;
            
            GUI.slider.slice = uicontrol('Parent', GUI.figureUI, 'style', 'slider');
            GUI.slider.slice.BackgroundColor = 'white';
            GUI.slider.slice.Units = 'normalized';
            pos = GUI.ax.modZ2D.Position;
            GUI.slider.slice.Position = [pos(1)+pos(3)+1/15, pos(2), 1/80, pos(4)];
            GUI.slider.slice.Callback = @GUI.sliderSliceCallback;
            GUI.slider.slice.Visible = 'off';
            
            % Text fields
            pos = GUI.slider.cor.Position;
            GUI.txtfield.cor = uicontrol('Parent', GUI.figureUI, 'style', 'edit');
            GUI.txtfield.cor.BackgroundColor = 'black';            
            GUI.txtfield.cor.Units = 'normalized';
            GUI.txtfield.cor.Position = pos + [0, -1/30, 0, 0];
            GUI.txtfield.cor.String = 'Coronal';
            GUI.txtfield.cor.ForegroundColor = [0.9 0.9 0.9];
            GUI.txtfield.cor.Callback = @GUI.txtfieldCorCallback;
            
            pos = GUI.slider.sag.Position;
            GUI.txtfield.sag = uicontrol('Parent', GUI.figureUI, 'style', 'edit');
            GUI.txtfield.sag.BackgroundColor = 'black';            
            GUI.txtfield.sag.Units = 'normalized';
            GUI.txtfield.sag.Position = pos + [0, -1/30, 0, 0];
            GUI.txtfield.sag.String = 'Sagittal';
            GUI.txtfield.sag.ForegroundColor = [0.9 0.9 0.9];
            GUI.txtfield.sag.Callback = @GUI.txtfieldSagCallback;
            
            pos = GUI.slider.axi.Position;
            GUI.txtfield.axi = uicontrol('Parent', GUI.figureUI, 'style', 'edit');
            GUI.txtfield.axi.BackgroundColor = 'black';            
            GUI.txtfield.axi.Units = 'normalized';
            GUI.txtfield.axi.Position = pos + [0, -1/30, 0, 0];
            GUI.txtfield.axi.String = 'Axial';
            GUI.txtfield.axi.ForegroundColor = [0.9 0.9 0.9];
            GUI.txtfield.axi.Callback = @GUI.txtfieldAxiCallback;
            
            pos = GUI.slider.dwi.Position;
            GUI.txtfield.dwi = uicontrol('Parent', GUI.figureUI, 'style', 'edit');
            GUI.txtfield.dwi.BackgroundColor = 'black';            
            GUI.txtfield.dwi.Units = 'normalized';
            GUI.txtfield.dwi.Position = pos + [0, -1/30, 0, 0];
            GUI.txtfield.dwi.String = 'DWI';
            GUI.txtfield.dwi.ForegroundColor = [0.9 0.9 0.9];
            GUI.txtfield.dwi.Callback = @GUI.txtfieldDwiCallback;
            
            GUI.txtfield.thresholdLowert = uicontrol('Parent', GUI.figureUI, 'style', 'text');
            GUI.txtfield.thresholdLowert.BackgroundColor = 'black';            
            GUI.txtfield.thresholdLowert.Units = 'normalized';
            GUI.txtfield.thresholdLowert.Position = [3/4+1/100, 1/20, 1/30, 1/40];
            GUI.txtfield.thresholdLowert.String = 'Lower_t';
            GUI.txtfield.thresholdLowert.FontWeight = 'bold';
            GUI.txtfield.thresholdLowert.ForegroundColor = [0.9 0.9 0.9];
            
            pos = GUI.txtfield.thresholdLowert.Position;
            GUI.txtfield.thresholdLower = uicontrol('Parent', GUI.figureUI, 'style', 'edit');
            GUI.txtfield.thresholdLower.BackgroundColor = 'black';            
            GUI.txtfield.thresholdLower.Units = 'normalized';
            GUI.txtfield.thresholdLower.Position = [pos(1)+pos(3)+1/100, 1/20, 1/30, 1/40];
            GUI.txtfield.thresholdLower.String = '3.5';
            GUI.txtfield.thresholdLower.ForegroundColor = [0.9 0.9 0.9];
            GUI.txtfield.thresholdLower.Callback = @GUI.txtfieldthresholdLowerCallback;
            
            GUI.txtfield.thresholdUppert = uicontrol('Parent', GUI.figureUI, 'style', 'text');
            GUI.txtfield.thresholdUppert.BackgroundColor = 'black';            
            GUI.txtfield.thresholdUppert.Units = 'normalized';
            GUI.txtfield.thresholdUppert.Position = [3/4+1/100, 1/40, 1/30, 1/40];
            GUI.txtfield.thresholdUppert.String = 'Upper_t';
            GUI.txtfield.thresholdUppert.FontWeight = 'bold';
            GUI.txtfield.thresholdUppert.ForegroundColor = [0.9 0.9 0.9];
            
            pos = GUI.txtfield.thresholdUppert.Position;
            GUI.txtfield.thresholdUpper = uicontrol('Parent', GUI.figureUI, 'style', 'edit');
            GUI.txtfield.thresholdUpper.BackgroundColor = 'black';            
            GUI.txtfield.thresholdUpper.Units = 'normalized';
            GUI.txtfield.thresholdUpper.Position = [pos(1)+pos(3)+1/100, 1/40, 1/30, 1/40];
            GUI.txtfield.thresholdUpper.String = '10.0';
            GUI.txtfield.thresholdUpper.ForegroundColor = [0.9 0.9 0.9];
            GUI.txtfield.thresholdUpper.Callback = @GUI.txtfieldthresholdUpperCallback;            
            
            GUI.togglebtn.mask = uicontrol('Parent', GUI.figureUI, 'style', 'togglebutton');
            GUI.togglebtn.mask.Units = 'normalized';
            GUI.togglebtn.mask.Position = [pos(1)+pos(3)+1/20, 1/40, 1/20, 1/40];
            GUI.togglebtn.mask.String = 'Show mask';                  
            GUI.togglebtn.mask.Callback = @GUI.toggleMaskCallback;
           
           axes(GUI.ax.axiLK)
           GUI.lines.axiSAG = line([1,2], [1,2], 'color', 'red');
           GUI.lines.axiCOR = line([1,2], [1,2], 'color', 'green');
           axes(GUI.ax.cor)
           GUI.lines.corSAG = line([1,2], [1,2], 'color', 'red');
           GUI.lines.corAXI = line([1,2], [1,2], 'color', 'blue');
           axes(GUI.ax.sag)
           GUI.lines.sagCOR = line([1,2], [1,2], 'color', 'green');
           GUI.lines.sagAXI = line([1,2], [1,2], 'color', 'blue');
           axes(GUI.ax.modZ2D)
           GUI.lines.modZDWI = line([1,2], [1,2], 'color', 'white');
           GUI.lines.modZAXI = line([1,2], [1,2], 'color', 'white');
           hold(GUI.ax.modZ2D, 'on');
           GUI.lines.modZM = plot(1, 1, 'm.', 'markersize', 20);
           GUI.lines.modZG = plot(1, 1, 'c.', 'markersize', 20);
           hold(GUI.ax.modZ2D, 'off');
           
           GUI.popup.shell = uicontrol('Parent', GUI.figureUI, 'Style', 'popup', 'String', {'Shells'});
           GUI.popup.shell.Units = 'normalized';
           GUI.popup.shell.Position = [4/13, 1/8, 1/8, 1/100];
           GUI.popup.shell.Callback = @GUI.selectShell;                                                                 
           
           GUI.data.window.changing = false;
           GUI.data.window.width = 0.5;
           GUI.data.window.level = 0.5;
           GUI.data.window.mousePosI = [0,0];
           GUI.data.window.mousePosF = [0,0];
           GUI.data.window.timer = timer('Name', 'WindowTimer', ...
               'Period', 0.001, ...
               'StartDelay', 0.001, ...
               'TasksToExecute', inf, ...
               'ExecutionMode', 'fixedSpacing', ...
               'TimerFcn', @GUI.windowTimerCallback);
        

           
        end
        

    end
   
    methods(Access = public)
       
        % Constructor 
        function GUI = SOLID
            
            initializeGUI(GUI);

            if nargout == 0
                clear GUI
            end        
        end
        
        
      
    end
    
    
end