classdef SOLID < handle
% SOLID - Slicewise outlier detection for diffusion weighted MRI data.
%
%             For details, please check the original paper
%           https://doi.org/10.1016/j.neuroimage.2018.07.003
% Use GUI by typing > SOLID
%
% Use COMMANDLINE version by giving following parameters 
% 'in', 'PATH to DWI file'
% 'bval', 'PATH to optional bval file'
% 'bvec', 'PATH to optional bvec file'
% 'mask', 'PATH to optional mask file'
% 'metric', 'Variance/Mean/Iod' | Default Variance
% 'thrL', 'Lower threshold value' | Default 3.5
% 'thrU', 'Upper threshold value' | Default 10.0
% 
% Example: > 
%   SOLID('in', 'DWI.nii.gz', 'mask', 'DWI_mask.nii.gz', 'metric', 'Mean');
% 
% GUI General Usage: 
%
% Import data from menu SOLID > Open. You will be prompted for the
% mandatory DWI NIfTI file and the three optional files (bval/bvec/mask).
% If you do not have specific optional files, press CANCEL, to use default
% values.
% 
% Save results from menu SOLID > Save. Results are saved to the same folder
% where mandatory DWI NIfTI file is located.
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
% obj. 
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
% Automated masking requires ExploreDTI function
%   E_DTI_Create_Mask_From_DWI_enhanced_IND.p
%
% **********************  Author: Viljami Sairanen  ***********************
% *********************  viljami.sairanen@gmail.com  **********************
% *************  Website: https://github.com/vilsaira/SOLID  **************
% *********************  Last edited: 30 August 2018  *********************  
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
       data = struct('UseGUI', 0, ...
                     'DWI', [],...
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
        
       function sliderCorCallback(obj, source, event)
           value = round(obj.slider.cor.Value);
           obj.slider.cor.Value = value;
           obj.txtfield.cor.String = value;
           updateaxes(obj);           
       end
       
       function sliderSagCallback(obj, source, event)
           value = round(obj.slider.sag.Value);
           obj.slider.sag.Value = value;
           obj.txtfield.sag.String = value;
           updateaxes(obj);             
       end
       
       function sliderAxiCallback(obj, source, event)
           value = round(obj.slider.axi.Value);
           obj.slider.axi.Value = value;
           obj.slider.slice.Value = value;
           obj.txtfield.axi.String = value;
           updateaxes(obj);
       end
       
       function sliderDwiCallback(obj, source, event)
           value = round(obj.slider.dwi.Value);
           obj.slider.dwi.Value = value;
           obj.txtfield.dwi.String = value;
           updateaxes(obj);             
       end
       
       function sliderSliceCallback(obj, source, event)
           value = round(obj.slider.slice.Value);
           obj.slider.axi.Value = value;
           obj.slider.slice.Value = value;
           obj.txtfield.axi.String = value;
           updateaxes(obj);           
       end
       
       function txtfieldCorCallback(obj, source, event)
           value = round(str2double(obj.txtfield.cor.String));
           if value >= obj.slider.cpr.Min && value <= obj.slider.cor.Max
               obj.slider.cor.Value = value;
               obj.txtfield.cor.String = value;
           else
               obj.txtfield.cor.String = obj.slider.cor.Value;
           end
           updateaxes(obj);
       end
       
       function txtfieldSagCallback(obj, source, event)
           value = round(str2double(obj.txtfield.sag.String));
           if value >= obj.slider.sag.Min && value <= obj.slider.sag.Max
               obj.slider.sag.Value = value;
               obj.txtfield.sag.String = value;
           else
               obj.txtfield.sag.String = obj.slider.sag.Value;
           end
           updateaxes(obj);
       end
       
       function txtfieldAxiCallback(obj, source, event)
           value = round(str2double(obj.txtfield.axi.String));
           if value >= obj.slider.axi.Min && value <= obj.slider.axi.Max
               obj.slider.axi.Value = value;
               obj.slider.slice.Value = value;
               obj.txtfield.axi.String = value;
               updateaxes(obj);
           else
               obj.txtfield.axi.String = obj.slider.axi.Value;
           end
       end
       
       function txtfieldDwiCallback(obj, source, event)
           value = round(str2double(obj.txtfield.dwi.String));
           if value >= obj.slider.dwi.Min && value <= obj.slider.dwi.Max
               obj.slider.dwi.Value = value;
               obj.txtfield.dwi.String = value;
               updateaxes(obj);
           else
               obj.txtfield.dwi.String = obj.slider.dwi.Value;
           end
       end
       
       function txtfieldthresholdLowerCallback(obj, source, event)
           value = str2double(obj.txtfield.thresholdLower.String);
           if value >= 0 && value <= obj.data.thresholdUpper
               obj.txtfield.thresholdLower.String = value;
               obj.data.thresholdLower = value;
           end              
           updateaxes(obj);
       end
       
       function txtfieldthresholdUpperCallback(obj, source, event)
           value = str2double(obj.txtfield.thresholdUpper.String);
           if value >= 0 && value > obj.data.thresholdLower
               obj.txtfield.thresholdUpper.String = value;
               obj.data.thresholdUpper = value;
           end              
           updateaxes(obj);        
       end
       
        function toggleVarMean(obj, source, event)
            if strcmp(version('-release'), '2015a')
                metric = source.Label;
            else
                metric = source.Text;
            end
            switch metric
                case 'Variance'
                    obj.fileMenu_Metric_Var.Checked = 'on';
                    obj.fileMenu_Metric_Mean.Checked = 'off';
                    obj.fileMenu_Metric_Iod.Checked = 'off';
                case 'Mean'
                    obj.fileMenu_Metric_Var.Checked = 'off';
                    obj.fileMenu_Metric_Mean.Checked = 'on';
                    obj.fileMenu_Metric_Iod.Checked = 'off';
                case 'Iod'
                    obj.fileMenu_Metric_Var.Checked = 'off';
                    obj.fileMenu_Metric_Mean.Checked = 'off';
                    obj.fileMenu_Metric_Iod.Checked = 'on';
            end
            obj.data.metric = metric;
            calculateModZ(obj);
            updateaxes(obj);
        end
        
        function toggleUseMask(obj, source, event)            
            if strcmp(obj.fileMenu_UseMask.Checked, 'on')
                obj.data.useMask = false;
                obj.fileMenu_UseMask.Checked = 'off';
%                 if strcmp(version('-release'), '2015a')
%                     obj.fileMenu_UseMask.Label = 'Use mask';
%                 else
%                     obj.fileMenu_UseMask.Text = 'Use mask';
%                 end
            else
                obj.data.useMask = true;
                obj.fileMenu_UseMask.Checked = 'on';
%                 if strcmp(version('-release'), '2015a')
%                     obj.fileMenu_UseMask.Label = 'Use mask';
%                 else
%                     obj.fileMenu_UseMask.Text = 'Use mask';
%                 end
            end
            calculateModZ(obj);
            updateaxes(obj);
        end
        
        function toggleMaskCallback(obj, source, event)
            if obj.togglebtn.mask.Value
                obj.togglebtn.mask.BackgroundColor = 'green';
            else
                obj.togglebtn.mask.BackgroundColor = [0.94, 0.94, 0.94];
            end
            updateaxes(obj);              
        end

        function selectShell(obj, source, event)
            initializeAxesSlider(obj);
            updateaxes(obj);
        end
        
        function initializeAxesSlider(obj, source, event)
           shell = str2double(obj.popup.shell.String(obj.popup.shell.Value, :));
           shellInds = obj.data.bval == shell;
           
           dims = size(obj.data.DWI(:,:,:,shellInds));
           if size(dims,2) < 4
               dims(4) = 1;
           end
           
           obj.slider.dwi.Min = 1;
           obj.slider.dwi.Max = dims(4);
           obj.slider.dwi.SliderStep = [1, 1]./(max([dims(4),2])-1);
           obj.slider.dwi.Value = 1;
           
           obj.slider.slice.Min = 1;
           obj.slider.slice.Max = dims(3);
           obj.slider.slice.SliderStep = [1, 1]./(dims(3)-1);
           obj.slider.slice.Value = round(dims(3)/2);
           
           obj.slider.axi.Min = obj.slider.slice.Min;
           obj.slider.axi.Max = obj.slider.slice.Max;
           obj.slider.axi.SliderStep = [1, 1]./(dims(3)-1);
           obj.slider.axi.Value = obj.slider.slice.Value;
           
           obj.slider.cor.Min = 1;
           obj.slider.cor.Max = dims(1);
           obj.slider.cor.SliderStep = [1, 1]./(dims(1)-1);
           obj.slider.cor.Value = round(dims(1)/2);
           
           obj.slider.sag.Min = 1;
           obj.slider.sag.Max = dims(2);
           obj.slider.sag.SliderStep = [1, 1]./(dims(2)-1);
           obj.slider.sag.Value = round(dims(2)/2);
           
           obj.txtfield.dwi.Min = obj.slider.dwi.Min;
           obj.txtfield.dwi.Max = obj.slider.dwi.Max;           
           obj.txtfield.dwi.String = obj.slider.dwi.Value;
           
           obj.txtfield.axi.Min = obj.slider.axi.Min;
           obj.txtfield.axi.Max = obj.slider.axi.Max;
           obj.txtfield.axi.String = obj.slider.slice.Value;
           
           obj.txtfield.cor.Min = obj.slider.cor.Min;
           obj.txtfield.cor.Max = obj.slider.cor.Max;
           obj.txtfield.cor.String = obj.slider.cor.Value;
           
           obj.txtfield.sag.Min = obj.slider.sag.Min;
           obj.txtfield.sag.Max = obj.slider.sag.Max;
           obj.txtfield.sag.String = obj.slider.sag.Value;
           
           D = obj.data.DWI(:,:,:,shellInds);
           obj.data.clim = [max([1, prctile(D(:), 5)]), ...
               max([1, prctile(D(:), 99)])];
           
           colormap(obj.ax.sag, 'gray');
           colormap(obj.ax.cor, 'gray');
           colormap(obj.ax.axiLK, 'gray');
           colormap(obj.ax.axiLMK, 'gray');
           colormap(obj.ax.axiLKM, 'gray');
           colormap(obj.ax.axiLPK, 'gray');
           colormap(obj.ax.axiLKP, 'gray');       
           
%            histogram(obj.ax.DWIhistogram, obj.data.DWI(:));           
%            obj.ax.DWIhistogram.YScale = 'log';
           
        end       
        
        function updateaxes(obj, source, event)            
            shell = str2double(obj.popup.shell.String(obj.popup.shell.Value, :));
            shellInds = obj.data.bval == shell;
            shellNums = zeros(size(shellInds));
            shellNums(shellInds) = 1:sum(shellInds);
            
            currentDWI = round(obj.slider.dwi.Value) == shellNums;
            dwiInd = find(currentDWI,1,'first');
            currentAXI = round(obj.slider.axi.Value);
            currentSAG = round(obj.slider.sag.Value);
            currentCOR = round(obj.slider.cor.Value);
            
            limitsDWI = [obj.slider.dwi.Min, obj.slider.dwi.Max];
            limitsAXI = [obj.slider.axi.Min, obj.slider.axi.Max];
            limitsCOR = [obj.slider.cor.Min, obj.slider.cor.Max];
            limitsSAG = [obj.slider.sag.Min, obj.slider.sag.Max];
            
            imgSAG = rot90(squeeze(obj.data.DWI(:, currentSAG, :, currentDWI)));
            imgCOR = rot90(squeeze(obj.data.DWI(currentCOR, :, :, currentDWI)));
            imgAXI = squeeze(obj.data.DWI(:,:,currentAXI, currentDWI));

            if dwiInd > 1 && shellInds(dwiInd-1)
                dwiIndM = dwiInd -1;
            else
                dwiIndM = [];
            end
            if dwiInd < size(obj.data.DWI,4) && shellInds(dwiInd+1)
                dwiIndP = dwiInd+1;
            else
                dwiIndP = [];
            end
            imgAXI_LMK = squeeze(obj.data.DWI(:,:,currentAXI, obj.data.nearestPoints(dwiInd,2)));
            if isempty(imgAXI_LMK)
                imgAXI_LMK = NaN(size(imgAXI_LMK,1), size(imgAXI_LMK,2),1);
            end
            imgAXI_LPK = squeeze(obj.data.DWI(:,:,currentAXI, obj.data.nearestPoints(dwiInd,3)));
            if isempty(imgAXI_LPK)
                imgAXI_LPK = NaN(size(imgAXI_LPK,1), size(imgAXI_LPK,2),1);
            end
            imgAXI_LKM = squeeze(obj.data.DWI(:,:,max([currentAXI-1, obj.slider.axi.Min]), currentDWI));
            imgAXI_LKP = squeeze(obj.data.DWI(:,:,min([currentAXI+1, obj.slider.axi.Max]), currentDWI));            
            if obj.togglebtn.mask.Value
                imgSAG = imfuse(imgSAG, rot90(squeeze(obj.data.mask(:, currentSAG, :))));
                imgCOR = imfuse(imgCOR, rot90(squeeze(obj.data.mask(currentCOR, :, :))));
                imgAXI = imfuse(imgAXI, squeeze(obj.data.mask(:,:,currentAXI)));
                imgAXI_LMK = imfuse(imgAXI_LMK, squeeze(obj.data.mask(:,:,currentAXI)));
                imgAXI_LPK = imfuse(imgAXI_LPK, squeeze(obj.data.mask(:,:,currentAXI)));
                imgAXI_LKM = imfuse(imgAXI_LKM, squeeze(obj.data.mask(:,:,max([currentAXI-1, obj.slider.axi.Min]))));
                imgAXI_LKP = imfuse(imgAXI_LKP, squeeze(obj.data.mask(:,:,min([currentAXI+1, obj.slider.axi.Max]))));
            end
            
            obj.images.sag.CData = flipud(imgSAG);
            set(obj.ax.sag, 'Ydir', 'normal', 'Clim', obj.data.clim);  
            obj.lines.sagAXI.XData = limitsCOR;
            obj.lines.sagAXI.YData = currentAXI.*[1,1]+0.5;
            obj.lines.sagCOR.XData = currentCOR.*[1,1]+0.5;
            obj.lines.sagCOR.YData = limitsAXI;

            obj.images.cor.CData = flipud(imgCOR);
            set(obj.ax.cor, 'Ydir', 'normal', 'Clim', obj.data.clim);  
            obj.lines.corAXI.XData = limitsSAG;
            obj.lines.corAXI.YData = currentAXI.*[1,1]+0.5;
            obj.lines.corSAG.XData = currentSAG.*[1,1]+0.5;
            obj.lines.corSAG.YData = limitsAXI;
            
            obj.images.axiLK.CData = imgAXI;
            set(obj.ax.axiLK, 'Clim', obj.data.clim);  
            obj.lines.axiCOR.XData = limitsCOR;
            obj.lines.axiCOR.YData = currentCOR.*[1,1]+0.5;
            obj.lines.axiSAG.XData = currentSAG.*[1,1]+0.5;
            obj.lines.axiSAG.YData = limitsSAG;                        
            
            obj.images.axiLKM.CData = imgAXI_LKM;
            set(obj.ax.axiLKM, 'Clim', obj.data.clim);             
            
            obj.images.axiLKP.CData = imgAXI_LKP;
            set(obj.ax.axiLKP, 'Clim', obj.data.clim); 
            
            obj.images.axiLMK.CData = imgAXI_LMK;
            set(obj.ax.axiLMK, 'Clim', obj.data.clim);             
            
            obj.images.axiLPK.CData = imgAXI_LPK;
            set(obj.ax.axiLPK, 'Clim', obj.data.clim);                        
                
%             obj.lines.dwiHistogramHor.XData
%             obj.lines.dwiHistogramHor.XData = obj.data.clim;
%             obj.lines.dwiHistogramHor.YData = 0.9.*obj.ax.DWIhistogram.Ylim(2).*[1,1];
%             obj.lines.dwiHistogramVer.XData = mean(obj.data.clim).*[1,1];
%             obj.lines.dwiHistogramVer.YData = obj.ax.DWIhistogram.Ylim(2).*[1,1];

            %% Spherical interpolation
            % first rotate all bvecs so that the first direction is on
            % z-axis so interpolation near polar region works better.
            % Everything must be rotated back.
            
            bvec = obj.data.bvec;           
            i = find(obj.data.bval > 0, 1, 'first');
            tmp = obj.data.modZ(currentAXI,shellInds); % modZ values for specific SLICE
          
            
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
            bvec = obj.data.bvec;
            Xi(:) = v(1,:);
            Yi(:) = v(2,:);
            Zi(:) = v(3,:);                                          
            
%             figure
%             surf(azimuth_mat, elevation_mat, Ci, 'edgecolor', 'none', 'facecolor', 'interp')
%             view(2)
            %%

            plot3(obj.ax.qspace,[bvec(shellInds,1); -bvec(shellInds,1)],...
                [bvec(shellInds,2); -bvec(shellInds,2)],...
                [bvec(shellInds,3); -bvec(shellInds,3)], 'ro');
            hold(obj.ax.qspace, 'on')
            
            plot3(obj.ax.qspace,[bvec(dwiInd,1); -bvec(dwiInd,1)], [bvec(dwiInd,2); -bvec(dwiInd,2)], [bvec(dwiInd,3); -bvec(dwiInd,3)], 'gx', 'markersize', 24);            
            dwiInd2 = obj.data.nearestPoints(dwiInd,2);
            dwiInd3 = obj.data.nearestPoints(dwiInd,3);
            plot3(obj.ax.qspace,[bvec(dwiInd2,1); -bvec(dwiInd2,1)], [bvec(dwiInd2,2); -bvec(dwiInd2,2)], [bvec(dwiInd2,3); -bvec(dwiInd2,3)], 'c.', 'markersize', 20); 
            plot3(obj.ax.qspace,[bvec(dwiInd3,1); -bvec(dwiInd3,1)], [bvec(dwiInd3,2); -bvec(dwiInd3,2)], [bvec(dwiInd3,3); -bvec(dwiInd3,3)], 'm.', 'markersize', 20);
            surf(obj.ax.qspace,0.98*Xi, 0.98*Yi, 0.98*Zi, Ci, 'EdgeColor', 'none');                        
            hold(obj.ax.qspace, 'off')
            shading(obj.ax.qspace, 'interp');
            axis(obj.ax.qspace, 'square');
            xlabel(obj.ax.qspace, 'X');
            ylabel(obj.ax.qspace, 'Y');
            zlabel(obj.ax.qspace, 'Z');
            grid(obj.ax.qspace, 'on');
            obj.ax.qspace.CLim = [obj.data.thresholdLower, obj.data.thresholdUpper];
            obj.ax.qspace.XColor = [0.5, 0.5, 0.5];
            obj.ax.qspace.YColor = [0.5, 0.5, 0.5];
            obj.ax.qspace.ZColor = [0.5, 0.5, 0.5];
            
            [az, el, ~] = cart2sph(bvec(dwiInd, 1), bvec(dwiInd, 2), bvec(dwiInd, 3));
            az = az./pi.*180;
            el = el./pi.*180;
            az = az - 90;
            az(az < 0) = az + 360;
            el = -el;            
            obj.ax.qspace.View = [az, el];

            %%            

            modZclim = [obj.data.thresholdLower, obj.data.thresholdUpper];
            obj.images.modZ2D.CData = obj.data.modZ(:, shellInds);
            set(obj.ax.modZ2D, 'YDir', 'normal', 'CLim', modZclim);
            axis(obj.ax.modZ2D, 'tight');
            obj.lines.modZAXI.XData = [0.5, sum(shellInds)+0.5];
            obj.lines.modZAXI.YData = currentAXI.*[1,1];            
            obj.lines.modZDWI.XData = shellNums(dwiInd).*[1,1];
            obj.lines.modZDWI.YData = limitsAXI+[-0.5, 0.5];
            
            obj.lines.modZG.XData = shellNums(dwiInd2);
            obj.lines.modZG.YData = currentAXI;
            obj.lines.modZM.XData = shellNums(dwiInd3);
            obj.lines.modZM.YData = currentAXI;

            
            bins = 0.1:0.2:modZclim(2);
            tmp = obj.data.modZ(:, shellInds);
            hh = histogram(obj.ax.modZhistogram, tmp(:), bins);
            axes(obj.ax.modZhistogram)
            axis(obj.ax.modZhistogram, 'tight');          
            line(obj.data.modZ(currentAXI, currentDWI).*[1,1], [1,max(hh.Values)], 'color', 'black', 'linewidth', 2);
            line(modZclim(1).*[1,1], [1,max(hh.Values)], 'color', 'green', 'linestyle', '--');
            line(modZclim(2).*[1,1], [1,max(hh.Values)], 'color', 'red', 'linestyle', '--');                          
            legend(obj.ax.modZhistogram.Children([1,2,3]), 'Upper threshold', 'Lower threshold', ['Current ', num2str(obj.data.modZ(currentAXI, currentDWI),2)]);
            obj.ax.modZhistogram.XColor = [0.9 0.9 0.9];
            obj.ax.modZhistogram.YColor = [0.9 0.9 0.9];
            obj.ax.modZhistogram.XLim = [-0.5, max([modZclim(2), obj.data.modZ(currentAXI, currentDWI)])+1];
        end
        
        function calculateModZ(obj, source, event)
            DWI = obj.data.DWI;
            DWI(DWI < 0) = NaN; % User can opt to use their own mask by setting non-brain voxels to negative numbers.
            b = obj.data.bval;
            metric = obj.data.metric;
            
            if obj.data.useMask
                mask = obj.data.mask(:,:,:,1) > 0;                
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
                    case 'Mean'
                        y = squeeze(mean(baseline, 'omitnan'));
                    case 'Iod'
                        y = squeeze(var(baseline, 'omitnan'))./squeeze(mean(baseline, 'omitnan'));
                    otherwise % Variance
                        y = squeeze(var(baseline, 'omitnan'));
                end
                tmp = repmat(median(y,2,'omitnan'), [1,N]);
                MAD = 1.4826 * median( abs(y-tmp),2,'omitnan');
                modZ(:,shell) = abs(y-tmp)./repmat(MAD,[1,N]);                
            end
            
            obj.data.uniqueb = uniqueb;
            obj.data.modZ = modZ;
            obj.data.modZorig = modZ;
            
        end
        
        function nii = loadNifti(obj, fpath)
            try
                nii = niftiread(fpath);
            catch
                try
                    nii = load_untouch_nii(fpath);
                    nii = nii.img;
                catch
                    error('*.nii(.gz) file not found');
                end
            end            
            if any(size(nii) < 2)
                error('Error in NIfTI image dimensions');
            end
            nii = permute(nii, [2 1 3 4]);
        end
        
        function bval = loadBVal(obj, fpath)            
            try 
                bval = load(fpath);
                bval = round(bval/100)*100;
            catch
                error('*.bval not found');
            end
        end
        
        function bvec = loadBVec(obj, fpath)
            try 
                bvec = load(fpath);
                if size(bvec,1) < size(bvec,2)
                    bvec = bvec';
                end
            catch
                error('*.bvec not found');
            end
        end
        
        function loadData(obj, source, event)           
            
            [fnameDWI, fdirDWI] = uigetfile({'*.nii;*.nii.gz'}, 'Select 4D DWI NIfTI file');
            fpath = [fdirDWI, filesep, fnameDWI];
            DWI = obj.loadNifti(fpath);            
            
            [fname, fdir] = uigetfile(fullfile(fdirDWI, '*.bval'), 'Select .bval / press Cancel for default');
            if isequal(fname, 0)
                fpath = [fdirDWI, filesep, fnameDWI(1:end-4), '.bval'];
            else
                fpath = [fdir, filesep, fname];
            end
            bval = obj.loadBVal(fpath);
            
            [fname, fdir] = uigetfile(fullfile(fdirDWI, '*.bvec'), 'Select .bvec / press Cancel for default');
            if isequal(fname, 0)
                fpath = [fdirDWI, filesep, fnameDWI(1:end-4), '.bvec'];
            else
                fpath = [fdir, filesep, fname];
            end
            bvec = obj.loadBVec(fpath);
            
            [fname, fdir] = uigetfile(fullfile(fdirDWI, {'*.nii;*.nii.gz'}), 'Select mask NIfTI file');
            if isequal(fname, 0)
                try
                    mask = E_DTI_Create_Mask_From_DWI_enhanced_IND(DWI(:,:,:,1), 0.5, 7); 
                catch
                    mask = ones(size(DWI(:,:,:,1)));
                end
            else
                fpath = [fdir, filesep, fname];
                mask = obj.loadNifti(fpath);                
            end
            
            obj.data.DWI = DWI;
            obj.data.bval = bval;
            obj.data.bvec = bvec;
            obj.data.mask = mask;
            obj.data.minmax = [min(DWI(:)), max(DWI(:))];            
            obj.data.fname = fnameDWI;
            obj.data.fpath = fdirDWI;
            
            calculateModZ(obj);
            
            if obj.data.UseGUI
                calculateAngularNHood(obj);                        
                obj.popup.shell.String = obj.data.uniqueb;
                obj.popup.shell.Value = 1;            
                initializeAxesSlider(obj);
                updateaxes(obj);
            end
        end        
        
        function calculateAngularNHood(obj, source, event)
            b = obj.data.bval;
            v = obj.data.bvec;
            % normalize v
            n = sqrt(sum(v.^2,2));
            v = v./repmat(n, [1,3]);
            v(isnan(v)) = 0;
            ad = inf(length(b), length(b));
            
            for i = 1:length(b)
                for j = i:length(b)
                    vi = v(i, :);
                    vj = v(j, :);
%                     c = cross(vj, vi);
%                     d = abs(dot(vj, vi));
%                     nc = sqrt(sum(abs(c).^2));
                    euc = sqrt(sum( (b(j).^2.*vj-b(i).^2.*vi).^2))./b(i);
%                     ang = atan2(nc, d)./pi;
                    ad(i,j) = euc;
                    % ang = atan2(nc, d)./pi.*180 + e;
%                     ad(i,j) = sqrt(euc.^2 + ang.^2);
                end
            end
            ad(isnan(ad)) = 0;
            ad = triu(ad)+triu(ad,1)';
            [ad2, tmp] = sort(ad,2);
            [~, obj.data.nearestPoints] = sort(ad,2);
            
        end
        
        function saveOutliers(obj, source, event)
            i = strfind(obj.data.fname, '.nii');
            oname = obj.data.fname(1:i-1);
            oname = strcat(obj.data.fpath, filesep, oname, '_L_', ...
                num2str(obj.data.thresholdLower), '_U_', ...
                num2str(obj.data.thresholdUpper), '_', ...
                obj.data.metric, '_masked', num2str(obj.data.useMask));
            %% save txt file
            fid = fopen(strcat(oname, '_modZ2D.txt'), 'w');
            for i = 1:size(obj.data.modZ, 1)
                fprintf(fid, repmat('%f, ', [1, size(obj.data.modZ,2)-1]), obj.data.modZ(i, 1:end-1));
                fprintf(fid, '%f\n', obj.data.modZ(i, end));
            end
            fclose(fid);
            
            %% save 4D nifti
            flag = 0;
            try
                info = niftiinfo(strcat(obj.data.fpath, filesep, obj.data.fname));
                img = niftiread(info);                
                flag = 1;
            catch
                try
                    strcat(obj.data.fpath, filesep, obj.data.fname)
                    nii = load_untouch_nii(strcat(obj.data.fpath, filesep, obj.data.fname));                    
                    img = zeros(size(nii.img));
                    flag = 2;
                catch
                    error('*.nii file not found');                    
                end
            end
            
            for k = 1:size(obj.data.modZ,1)
                for l = 1:size(obj.data.modZ,2)
                    img(:,:,k,l) = repmat(obj.data.modZ(k,l), [size(img,1), size(img,2), 1,1]);
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
            b = obj.data.bval;
            uniqueb = unique(b);
            for i = 1 : length(uniqueb)
                shell = (b == uniqueb(i));
                modZ = obj.data.modZ(:,shell);
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
                ax.CLim = [obj.data.thresholdLower, obj.data.thresholdUpper];
                axis equal tight;
                print(fig, strcat(oname, '_modZ2D_b', num2str(uniqueb(i)), '.png'), '-dpng', '-r300');
                close(fig);
            end
        end        
        
        function mouseClick(obj, source, event)            
            C = get(obj.ax.modZ2D, 'CurrentPoint');
            curDWI = round(C(1,1));
            curSlice = round(C(1,2));
            
            if curDWI >= obj.ax.modZ2D.XLim(1) && curDWI <= obj.ax.modZ2D.XLim(2) && ...
                curSlice >= obj.ax.modZ2D.YLim(1) && curSlice <= obj.ax.modZ2D.YLim(2)
                
%                 shell = str2double(obj.popup.shell.String(obj.popup.shell.Value, :));
%                 shellInds = obj.data.bval == shell;
%                 shellNums = zeros(size(shellInds));
%                 shellNums(shellInds) = 1:sum(shellInds);
                
%                 curDWI = shellNums(curDWI);
                if curDWI
                    obj.slider.dwi.Value = curDWI;
                    obj.txtfield.dwi.String = curDWI;
                    obj.slider.axi.Value = curSlice;
                    obj.slider.slice.Value = curSlice;
                    obj.txtfield.axi.String = curSlice;
                    updateaxes(obj);
                end
            else
                % Change window / level
                C = get(obj.figureUI, 'CurrentPoint');
                obj.data.window.mousePosI = [C(1,1), C(1,2)];
                clim = obj.ax.axiLK.CLim;
                obj.data.window.width = clim(2)-clim(1);
                obj.data.window.level = obj.data.window.width/2 + clim(1);                
                start(obj.data.window.timer);
            end
        end
        
        function mouseRelease(obj, source, event)
            stop(obj.data.window.timer);
        end
        
        function mouseScroll(obj, source, event)            
            C = get(obj.figureUI, 'CurrentPoint');
            x = C(1,1);
            y = C(1,2);
%             [x,y]            
            P = obj.ax.cor.Position;
            if x > P(1) && x < P(1)+P(3) && y > P(2) && y < P(2)+P(4)
%                 'add cor slice'
                try
                    obj.slider.cor.Value = obj.slider.cor.Value - event.VerticalScrollCount;
                    obj.txtfield.cor.String = obj.slider.cor.Value;
                    updateaxes(obj);
                catch
                end
            end
            P = obj.ax.sag.Position;
            if x > P(1) && x < P(1)+P(3) && y > P(2) && y < P(2)+P(4)
%                 'add sag slice'
                try
                obj.slider.sag.Value = obj.slider.sag.Value - event.VerticalScrollCount;
                obj.txtfield.sag.String = obj.slider.sag.Value;
                updateaxes(obj);
                catch
                end
            end
            P = obj.ax.axiLK.Position;
            if x > P(1) && x < P(1)+P(3) && y > P(2) && y < P(2)+P(4)
%                 'add axi slice'
                try                   
                obj.slider.axi.Value = obj.slider.axi.Value - event.VerticalScrollCount;
                obj.slider.slice.Value = obj.slider.axi.Value;
                obj.txtfield.axi.String = obj.slider.axi.Value;
                updateaxes(obj);
                catch
                end
            end
            
            
        end       
        
        function windowTimerCallback(obj, source, event)
             
        end
        
        function mouseMove(obj, source, event)
            if strcmp(obj.data.window.timer.running, 'on')                
                
                clim = obj.data.clim;
                level = mean(clim);
                width = clim(2)-clim(1);
                
                C = get(obj.figureUI, 'CurrentPoint');
                x = C(1,1);
                y = C(1,2);
                dx = (x-obj.data.window.mousePosI(1))/3+1;
                dy = (y-obj.data.window.mousePosI(2))/3+1;                
                 
                
                level = level .* dx;
                width = width .* dy;
                
                if level > obj.data.minmax(2)
                    level = obj.data.minmax(2);
                end
                if level < obj.data.minmax(1)
                    level = obj.data.minmax(1);
                end
                if width < eps
                    width = eps;
                end
                
                clim = [level - width/2, level + width/2]
                                
                obj.data.clim = clim;
                                
                obj.ax.sag.CLim = clim;
                obj.ax.cor.CLim = clim;
                obj.ax.axiLK.CLim = clim;
                obj.ax.axiLMK.CLim = clim;
                obj.ax.axiLKM.CLim = clim;
                obj.ax.axiLPK.CLim = clim;
                obj.ax.axiLKP.CLim = clim;
                
%                 obj.lines.dwiHistogramHor.XData = obj.data.clim;
%                 obj.lines.dwiHistogramHor.YData = 0.9.*obj.ax.DWIhistogram.Ylim(2).*[1,1];
%                 obj.lines.dwiHistogramVer.XData = mean(obj.data.clim).*[1,1];
%                 obj.lines.dwiHistogramVer.YData = obj.ax.DWIhistogram.Ylim(2).*[1,1];
                
            end
        end
        
        function keyPress(obj, source, event)
            % shell list must be unselected !!
            obj.popup.shell.Enable = 'off';
            obj.slider.cor.Enable = 'off';
            obj.slider.sag.Enable = 'off';
            obj.slider.axi.Enable = 'off';
            obj.slider.dwi.Enable = 'off';
            pause(0.01);
            if strcmp(event.Key, 'space')
                shell = str2double(obj.popup.shell.String(obj.popup.shell.Value, :));
                shellInds = obj.data.bval == shell;
                shellNums = zeros(size(shellInds));
                shellNums(shellInds) = 1:sum(shellInds);
                curDWI = obj.slider.dwi.Value;
                curDWI = find(shellNums == curDWI);
                curSlice = obj.slider.axi.Value;
                z = obj.data.modZ(curSlice, curDWI);
                zorig = obj.data.modZorig(curSlice, curDWI);
                if abs(z - zorig) > eps
                    obj.data.modZ(curSlice, curDWI) = obj.data.modZorig(curSlice, curDWI);
                else
                    obj.data.modZ(curSlice, curDWI) = 0;
                end
                updateaxes(obj);  
            end
            if strcmp(event.Key, 'rightarrow')
                shell = str2double(obj.popup.shell.String(obj.popup.shell.Value, :));
                shellInds = obj.data.bval == shell;
                shellNums = zeros(size(shellInds));
                shellNums(shellInds) = 1:sum(shellInds);
                curDWI = obj.slider.dwi.Value;
                curDWI = find(shellNums == curDWI);
                maxDWI = find(shellNums == obj.slider.dwi.Max);
                curSlice = obj.slider.axi.Value;
                n = 0;
                flag = 0;
                for d = curDWI : 1 : maxDWI
                    if n == 0
                        cs = curSlice + 1;
                        if cs == obj.slider.axi.Max
                            cs = obj.slider.axi.Max;
                        end                        
                        n = 1;
                    else
                        cs = 1;                        
                    end
                    for s = cs : 1 : obj.slider.axi.Max
                        z = obj.data.modZ(s, d);
                        if z > obj.data.thresholdLower
                            flag = 1;
                            obj.slider.dwi.Value = shellNums(d);
                            obj.txtfield.dwi.String = shellNums(d);
                            obj.slider.axi.Value = s;
                            obj.slider.slice.Value = s;
                            obj.txtfield.axi.String = s;
                            break
                        end
                    end
                    if flag
                        break
                    end
                end
                updateaxes(obj);  
            end
            if strcmp(event.Key, 'leftarrow')
                shell = str2double(obj.popup.shell.String(obj.popup.shell.Value, :));
                shellInds = obj.data.bval == shell;
                shellNums = zeros(size(shellInds));
                shellNums(shellInds) = 1:sum(shellInds);
                curDWI = obj.slider.dwi.Value;
                curDWI = find(shellNums == curDWI);
                minDWI = find(shellNums == obj.slider.dwi.Min);
                curSlice = obj.slider.axi.Value;
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
                        cs = obj.slider.axi.Max;                        
                    end
                    for s = cs : -1 : 1
                        z = obj.data.modZ(s, d);
                        if z > obj.data.thresholdLower
                            flag = 1;
                            obj.slider.dwi.Value = shellNums(d);
                            obj.txtfield.dwi.String = shellNums(d);
                            obj.slider.axi.Value = s;
                            obj.slider.slice.Value = s;
                            obj.txtfield.axi.String = s;
                            break
                        end
                    end
                    if flag
                        break
                    end
                end
                updateaxes(obj); 
            end
            obj.popup.shell.Enable = 'on';
            obj.slider.cor.Enable = 'on';
            obj.slider.sag.Enable = 'on';
            obj.slider.axi.Enable = 'on';
            obj.slider.dwi.Enable = 'on';
            pause(0.01);
        end
        
        % Create UIFigure and components
        function initializeGUI(obj)
            
            obj.data.UseGUI = true;
            
            % Create SOLIDUIFigure
            obj.figureUI = figure;
            obj.figureUI.Color = [0 0 0];
            tmp = get(0, 'MonitorPositions');
            tmp = tmp(1,:);
            tmp = [tmp(3)*0.1 tmp(4)*0.1 tmp(3)*0.8 tmp(4)*0.8];
            obj.figureUI.Position = tmp;
            obj.figureUI.Name = 'SOLID';
            obj.figureUI.Units = 'normalized';
%             obj.figureUI.MenuBar = 'none';
%             obj.figureUI.ToolBar = 'none';  
            obj.figureUI.WindowButtonDownFcn = @obj.mouseClick;
            obj.figureUI.WindowButtonUpFcn = @obj.mouseRelease;
            obj.figureUI.WindowScrollWheelFcn = @obj.mouseScroll;
            obj.figureUI.WindowButtonMotionFcn = @obj.mouseMove;
            obj.figureUI.WindowKeyPressFcn = @obj.keyPress;
            
            % Create FileMenu
            obj.fileMenu = uimenu(obj.figureUI);
            if strcmp(version('-release'), '2015a')
                obj.fileMenu.Label = 'SOLID';
            else
                obj.fileMenu.Text = 'SOLID';
            end
            
            % Create FileMenu -> Open Nifti
            obj.fileMenu_Open = uimenu(obj.fileMenu);
            if strcmp(version('-release'), '2015a')
                obj.fileMenu_Open.Label = 'Open';
            else                
                obj.fileMenu_Open.Text = 'Open';
            end
            obj.fileMenu_Open.Accelerator = 'Q';
            if strcmp(version('-release'), '2015a')
                obj.fileMenu_Open.Callback = @obj.loadData;
            else
                obj.fileMenu_Open.MenuSelectedFcn = @obj.loadData;
            end
            
            % Create FileMenu -> Save outlier map
            obj.fileMenu_Save = uimenu(obj.fileMenu);
            if strcmp(version('-release'), '2015a')
                obj.fileMenu_Save.Label = 'Save results';
            else                
                obj.fileMenu_Save.Text = 'Save results';
            end            
            obj.fileMenu_Save.Accelerator = 's';            
            if strcmp(version('-release'), '2015a')
                obj.fileMenu_Save.Callback = @obj.saveOutliers;
            else
                obj.fileMenu_Save.MenuSelectedFcn = @obj.saveOutliers;
            end
            obj.fileMenu_Save.Visible = 'on'; 
            
            % Create FileMenu -> Toggle Var/Mean
            obj.fileMenu_Metric = uimenu(obj.fileMenu);
            obj.fileMenu_Metric_Var = uimenu(obj.fileMenu_Metric);
            obj.fileMenu_Metric_Mean = uimenu(obj.fileMenu_Metric);
            obj.fileMenu_Metric_Iod = uimenu(obj.fileMenu_Metric);
            if strcmp(version('-release'), '2015a')
                obj.fileMenu_Metric.Label = 'Select SOLID metric';
                obj.fileMenu_Metric_Var.Label = 'Variance';
                obj.fileMenu_Metric_Mean.Label = 'Mean';
                obj.fileMenu_Metric_Iod.Label = 'Iod';
            else
                obj.fileMenu_Metric.Text = 'Select SOLID metric';
                obj.fileMenu_Metric_Var.Text = 'Variance';
                obj.fileMenu_Metric_Mean.Text = 'Mean';
                obj.fileMenu_Metric_Iod.Text = 'Iod';
            end           
            if strcmp(version('-release'), '2015a')                
                obj.fileMenu_Metric_Var.Callback = @obj.toggleVarMean;
                obj.fileMenu_Metric_Mean.Callback = @obj.toggleVarMean;
                obj.fileMenu_Metric_Iod.Callback = @obj.toggleVarMean;
            else
                obj.fileMenu_Metric_Var.MenuSelectedFcn = @obj.toggleVarMean;
                obj.fileMenu_Metric_Mean.MenuSelectedFcn = @obj.toggleVarMean;
                obj.fileMenu_Metric_Iod.MenuSelectedFcn = @obj.toggleVarMean;
            end
            obj.fileMenu_Metric_Var.Checked = 'on';
            
            obj.fileMenu_UseMask = uimenu(obj.fileMenu);
            if strcmp(version('-release'), '2015a')
                obj.fileMenu_UseMask.Label = 'Use brain mask';
            else                
                obj.fileMenu_UseMask.Text = 'Use brain mask';
            end                         
            if strcmp(version('-release'), '2015a')                
                obj.fileMenu_UseMask.Callback = @obj.toggleUseMask;
            else
                obj.fileMenu_UseMask.MenuSelectedFcn = @obj.toggleUseMask;
            end            
            obj.fileMenu_UseMask.Checked = 'on';
            
            % Create axes 3x3 grid DWI
            % (1,1)
            obj.ax.cor = axes('Parent', obj.figureUI);
            obj.ax.cor.Position = [0, 3/4-2/21, 1/4, 1/4];            
            obj.images.cor = imagesc('Parent', obj.ax.cor, 'CData', NaN(50,50));
            axis(obj.ax.cor, 'off', 'tight', 'equal');
            
            % (1,2)
            obj.ax.axiLKP = axes('Parent', obj.figureUI);
            obj.ax.axiLKP.Position = [1/4, 3/4-2/21, 1/4, 1/4];            
            obj.images.axiLKP = imagesc('Parent', obj.ax.axiLKP, 'CData', NaN(50,50));
            axis(obj.ax.axiLKP, 'off', 'tight', 'equal');
            
            % (1,3)
            obj.ax.sag = axes('Parent', obj.figureUI);
            obj.ax.sag.Position = [2/4, 3/4-2/21, 1/4, 1/4];    
            obj.ax.sag.XDir = 'reverse';
            obj.images.sag = imagesc('Parent', obj.ax.sag, 'CData', NaN(50,50));
            axis(obj.ax.sag, 'off', 'tight', 'equal');
            
            % (2,1)
            obj.ax.axiLMK = axes('Parent', obj.figureUI);
            obj.ax.axiLMK.Position = [0/4, 2/4-2/20, 1/4, 1/4];            
            obj.images.axiLMK = imagesc('Parent', obj.ax.axiLMK, 'CData', NaN(50,50));
            axis(obj.ax.axiLMK, 'tight', 'equal');
            obj.ax.axiLMK.XColor = [0,1,1];
            obj.ax.axiLMK.YColor = [0,1,1];
            obj.ax.axiLMK.Box = 'on';
            obj.ax.axiLMK.LineWidth = 2;
            obj.ax.axiLMK.XTick = [];
            obj.ax.axiLMK.YTick = [];
            
            % (2,2)
            obj.ax.axiLK = axes('Parent', obj.figureUI);
            obj.ax.axiLK.Position = [1/4, 2/4-2/20, 1/4, 1/4];            
            obj.images.axiLK = imagesc('Parent', obj.ax.axiLK, 'CData', NaN(50,50));
            axis(obj.ax.axiLK, 'tight', 'equal');
            obj.ax.axiLK.XColor = [1,1,1];
            obj.ax.axiLK.YColor = [1,1,1];
            obj.ax.axiLK.Box = 'on';
            obj.ax.axiLK.LineWidth = 2;
            obj.ax.axiLK.XTick = [];
            obj.ax.axiLK.YTick = [];
            
            % (2,3)
            obj.ax.axiLPK = axes('Parent', obj.figureUI);
            obj.ax.axiLPK.Position = [2/4, 2/4-2/20, 1/4, 1/4];            
            obj.images.axiLPK = imagesc('Parent', obj.ax.axiLPK, 'CData', NaN(50,50));
            axis(obj.ax.axiLPK, 'tight', 'equal');
            obj.ax.axiLPK.XColor = [1,0,1];
            obj.ax.axiLPK.YColor = [1,0,1];
            obj.ax.axiLPK.Box = 'on';
            obj.ax.axiLPK.LineWidth = 2;
            obj.ax.axiLPK.XTick = [];
            obj.ax.axiLPK.YTick = [];
            
            % (3,2)
            obj.ax.axiLKM = axes('Parent', obj.figureUI);
            obj.ax.axiLKM.Position = [1/4, 1/4-2/19, 1/4, 1/4];            
            obj.images.axiLKM = imagesc('Parent', obj.ax.axiLKM, 'CData', NaN(50,50));
            axis(obj.ax.axiLKM, 'off', 'tight', 'equal');
            
            % Create axes qspace
            obj.ax.qspace = axes('Parent', obj.figureUI);
            obj.ax.qspace.Position = [2/4 1/20 1/4 1/4];
            obj.images.qspace = imagesc('Parent', obj.ax.qspace, 'CData', NaN(50,50));
            axis(obj.ax.qspace, 'square');
%             rotate3d(obj.ax.qspace, 'on');
            
            % Create axes modZ2D
            obj.ax.modZ2D = axes('Parent', obj.figureUI);
            obj.ax.modZ2D.Position = [3/4+1/100 2/20 1/4-2/30 16/20];
            obj.images.modZ2D = imagesc('Parent', obj.ax.modZ2D, 'CData', NaN(50,50));
            axis(obj.ax.modZ2D, 'tight', 'equal');
            ch = colorbar(obj.ax.modZ2D);
            ylabel(ch, 'Modified Z-score');
            ch.Color = [0.9 0.9 0.9];            
            obj.ax.modZ2D.XColor = [0.9 0.9 0.9];
            obj.ax.modZ2D.YColor = [0.9 0.9 0.9];
            obj.ax.modZ2D.Box = 'on';
            xlabel(obj.ax.modZ2D, 'Image volume');
            ylabel(obj.ax.modZ2D, 'Slice');
            
            % Create axes modZhistogram
            obj.ax.modZhistogram = axes('Parent', obj.figureUI);
            obj.ax.modZhistogram.Position =  [1/20 1/20 2/6-1/10 2/6];
            obj.images.modZhistogram = imagesc('Parent', obj.ax.modZhistogram, 'CData', NaN(50,50));
            axis(obj.ax.modZhistogram, 'tight', 'equal');
            xlabel(obj.ax.modZhistogram, 'Modified Z-score');
            ylabel(obj.ax.modZhistogram, 'Counts');
            obj.ax.modZhistogram.XColor = [0.9 0.9 0.9];
            obj.ax.modZhistogram.YColor = [0.9 0.9 0.9];

            % Create DWI histogram
            obj.ax.DWIhistogram = axes('Parent', obj.figureUI);
            obj.ax.DWIhistogram.Position = [9/10, 1/150, 1/11, 1/12];
            obj.images.DWIhistogram = histogram('Parent', obj.ax.DWIhistogram, NaN(50,1));
            axis(obj.ax.DWIhistogram, 'off', 'tight');
%             hold(obj.ax.DWIhistogram, 'on');                       
%             obj.lines.dwiHistogramVer = line(obj.ax.DWIhistogram, [1,2], [1,2], 'color', 'red');
%             obj.lines.dwiHistogramHor = line(obj.ax.DWIhistogram, [1,2], [1,2], 'color', 'red');
%             hold(obj.ax.DWIhistogram, 'off');                       
            
            % Create sliders
            obj.slider.cor = uicontrol('Parent', obj.figureUI, 'style', 'slider');
            obj.slider.cor.BackgroundColor = 'green';
            obj.slider.cor.Units = 'normalized';
            obj.slider.cor.Position = [0+1/100, 1-1/20+1/80, 1/4-2/100, 1/40];
            obj.slider.cor.Callback = @obj.sliderCorCallback;
            
            obj.slider.sag = uicontrol('Parent', obj.figureUI, 'style', 'slider');
            obj.slider.sag.BackgroundColor = 'red';
            obj.slider.sag.Units = 'normalized';
            obj.slider.sag.Position = [1/4+1/100, 1-1/20+1/80, 1/4-2/100, 1/40];
            obj.slider.sag.Callback = @obj.sliderSagCallback;
            
            obj.slider.axi = uicontrol('Parent', obj.figureUI, 'style', 'slider');
            obj.slider.axi.BackgroundColor = 'blue';
            obj.slider.axi.Units = 'normalized';
            obj.slider.axi.Position = [2/4+1/100, 1-1/20+1/80, 1/4-2/100, 1/40];
            obj.slider.axi.Callback = @obj.sliderAxiCallback;
            
            obj.slider.dwi = uicontrol('Parent', obj.figureUI, 'style', 'slider');
            obj.slider.dwi.BackgroundColor = 'white';
            obj.slider.dwi.Units = 'normalized';
            obj.slider.dwi.Position = [3/4+1/100, 1-1/20+1/80, 1/4-2/30, 1/40];
            obj.slider.dwi.Callback = @obj.sliderDwiCallback;
            
            obj.slider.slice = uicontrol('Parent', obj.figureUI, 'style', 'slider');
            obj.slider.slice.BackgroundColor = 'white';
            obj.slider.slice.Units = 'normalized';
            pos = obj.ax.modZ2D.Position;
            obj.slider.slice.Position = [pos(1)+pos(3)+1/15, pos(2), 1/80, pos(4)];
            obj.slider.slice.Callback = @obj.sliderSliceCallback;
            obj.slider.slice.Visible = 'off';
            
            % Text fields
            pos = obj.slider.cor.Position;
            obj.txtfield.cor = uicontrol('Parent', obj.figureUI, 'style', 'edit');
            obj.txtfield.cor.BackgroundColor = 'black';            
            obj.txtfield.cor.Units = 'normalized';
            obj.txtfield.cor.Position = pos + [0, -1/30, 0, 0];
            obj.txtfield.cor.String = 'Coronal';
            obj.txtfield.cor.ForegroundColor = [0.9 0.9 0.9];
            obj.txtfield.cor.Callback = @obj.txtfieldCorCallback;
            
            pos = obj.slider.sag.Position;
            obj.txtfield.sag = uicontrol('Parent', obj.figureUI, 'style', 'edit');
            obj.txtfield.sag.BackgroundColor = 'black';            
            obj.txtfield.sag.Units = 'normalized';
            obj.txtfield.sag.Position = pos + [0, -1/30, 0, 0];
            obj.txtfield.sag.String = 'Sagittal';
            obj.txtfield.sag.ForegroundColor = [0.9 0.9 0.9];
            obj.txtfield.sag.Callback = @obj.txtfieldSagCallback;
            
            pos = obj.slider.axi.Position;
            obj.txtfield.axi = uicontrol('Parent', obj.figureUI, 'style', 'edit');
            obj.txtfield.axi.BackgroundColor = 'black';            
            obj.txtfield.axi.Units = 'normalized';
            obj.txtfield.axi.Position = pos + [0, -1/30, 0, 0];
            obj.txtfield.axi.String = 'Axial';
            obj.txtfield.axi.ForegroundColor = [0.9 0.9 0.9];
            obj.txtfield.axi.Callback = @obj.txtfieldAxiCallback;
            
            pos = obj.slider.dwi.Position;
            obj.txtfield.dwi = uicontrol('Parent', obj.figureUI, 'style', 'edit');
            obj.txtfield.dwi.BackgroundColor = 'black';            
            obj.txtfield.dwi.Units = 'normalized';
            obj.txtfield.dwi.Position = pos + [0, -1/30, 0, 0];
            obj.txtfield.dwi.String = 'DWI';
            obj.txtfield.dwi.ForegroundColor = [0.9 0.9 0.9];
            obj.txtfield.dwi.Callback = @obj.txtfieldDwiCallback;
            
            obj.txtfield.thresholdLowert = uicontrol('Parent', obj.figureUI, 'style', 'text');
            obj.txtfield.thresholdLowert.BackgroundColor = 'black';            
            obj.txtfield.thresholdLowert.Units = 'normalized';
            obj.txtfield.thresholdLowert.Position = [3/4+1/100, 1/20, 1/30, 1/40];
            obj.txtfield.thresholdLowert.String = 'Lower_t';
            obj.txtfield.thresholdLowert.FontWeight = 'bold';
            obj.txtfield.thresholdLowert.ForegroundColor = [0.9 0.9 0.9];
            
            pos = obj.txtfield.thresholdLowert.Position;
            obj.txtfield.thresholdLower = uicontrol('Parent', obj.figureUI, 'style', 'edit');
            obj.txtfield.thresholdLower.BackgroundColor = 'black';            
            obj.txtfield.thresholdLower.Units = 'normalized';
            obj.txtfield.thresholdLower.Position = [pos(1)+pos(3)+1/100, 1/20, 1/30, 1/40];
            obj.txtfield.thresholdLower.String = '3.5';
            obj.txtfield.thresholdLower.ForegroundColor = [0.9 0.9 0.9];
            obj.txtfield.thresholdLower.Callback = @obj.txtfieldthresholdLowerCallback;
            
            obj.txtfield.thresholdUppert = uicontrol('Parent', obj.figureUI, 'style', 'text');
            obj.txtfield.thresholdUppert.BackgroundColor = 'black';            
            obj.txtfield.thresholdUppert.Units = 'normalized';
            obj.txtfield.thresholdUppert.Position = [3/4+1/100, 1/40, 1/30, 1/40];
            obj.txtfield.thresholdUppert.String = 'Upper_t';
            obj.txtfield.thresholdUppert.FontWeight = 'bold';
            obj.txtfield.thresholdUppert.ForegroundColor = [0.9 0.9 0.9];
            
            pos = obj.txtfield.thresholdUppert.Position;
            obj.txtfield.thresholdUpper = uicontrol('Parent', obj.figureUI, 'style', 'edit');
            obj.txtfield.thresholdUpper.BackgroundColor = 'black';            
            obj.txtfield.thresholdUpper.Units = 'normalized';
            obj.txtfield.thresholdUpper.Position = [pos(1)+pos(3)+1/100, 1/40, 1/30, 1/40];
            obj.txtfield.thresholdUpper.String = '10.0';
            obj.txtfield.thresholdUpper.ForegroundColor = [0.9 0.9 0.9];
            obj.txtfield.thresholdUpper.Callback = @obj.txtfieldthresholdUpperCallback;            
            
            obj.togglebtn.mask = uicontrol('Parent', obj.figureUI, 'style', 'togglebutton');
            obj.togglebtn.mask.Units = 'normalized';
            obj.togglebtn.mask.Position = [pos(1)+pos(3)+1/20, 1/40, 1/20, 1/40];
            obj.togglebtn.mask.String = 'Show mask';                  
            obj.togglebtn.mask.Callback = @obj.toggleMaskCallback;
           
           axes(obj.ax.axiLK)
           obj.lines.axiSAG = line([1,2], [1,2], 'color', 'red');
           obj.lines.axiCOR = line([1,2], [1,2], 'color', 'green');
           axes(obj.ax.cor)
           obj.lines.corSAG = line([1,2], [1,2], 'color', 'red');
           obj.lines.corAXI = line([1,2], [1,2], 'color', 'blue');
           axes(obj.ax.sag)
           obj.lines.sagCOR = line([1,2], [1,2], 'color', 'green');
           obj.lines.sagAXI = line([1,2], [1,2], 'color', 'blue');
           axes(obj.ax.modZ2D)
           obj.lines.modZDWI = line([1,2], [1,2], 'color', 'white');
           obj.lines.modZAXI = line([1,2], [1,2], 'color', 'white');
           hold(obj.ax.modZ2D, 'on');
           obj.lines.modZM = plot(1, 1, 'm.', 'markersize', 20);
           obj.lines.modZG = plot(1, 1, 'c.', 'markersize', 20);
           hold(obj.ax.modZ2D, 'off');
           
           obj.popup.shell = uicontrol('Parent', obj.figureUI, 'Style', 'popup', 'String', {'Shells'});
           obj.popup.shell.Units = 'normalized';
           obj.popup.shell.Position = [4/13, 1/8, 1/8, 1/100];
           obj.popup.shell.Callback = @obj.selectShell;                                                                 
           
           obj.data.window.changing = false;
           obj.data.window.width = 0.5;
           obj.data.window.level = 0.5;
           obj.data.window.mousePosI = [0,0];
           obj.data.window.mousePosF = [0,0];
           obj.data.window.timer = timer('Name', 'WindowTimer', ...
               'Period', 0.001, ...
               'StartDelay', 0.001, ...
               'TasksToExecute', inf, ...
               'ExecutionMode', 'fixedSpacing', ...
               'TimerFcn', @obj.windowTimerCallback);

        end
        
        function runCMD(obj, p)
            
            if isempty(p.in)
                error('DWI file not found');
            end
            disp(' ');
            disp(['Loading ' p.in]);
            tmp = strsplit(p.in, filesep);
            fnameDWI = tmp{end};
            fpathDWI = strjoin(tmp(1:end-1), filesep);
            
            obj.data.DWI = obj.loadNifti(p.in);
            
            if isempty(p.bval)
                obj.data.bval = obj.loadBVal([fpathDWI, filesep, fnameDWI(1:end-4), '.bval']);
            else
                obj.data.bval = obj.loadBVal(p.bval);
            end
            if isempty(p.bvec)
                obj.data.bvec = obj.loadBVec([fpathDWI, filesep, fnameDWI(1:end-4), '.bvec']);
            else
                obj.data.bvec = obj.loadBVal(p.bvec);
            end            
            if isempty(p.mask)
                try
                    disp(['Calculating brain mask']);
                    obj.data.mask = E_DTI_Create_Mask_From_DWI_enhanced_IND(obj.data.DWI(:,:,:,1), 0.5, 7); 
                catch
                    obj.data.mask = ones(size(obj.DWI(:,:,:,1)));
                end
            else
                obj.data.mask = obj.loadNifti(p.mask);
            end
            obj.data.minmax = [min(obj.data.DWI(:)), max(obj.data.DWI(:))];            
            obj.data.fname = fnameDWI;
            obj.data.fpath = fpathDWI;
            
            obj.data.metric = p.metric;
            obj.data.thresholdLower = p.thrL;
            obj.data.thresholdUpper = p.thrU;
            
            disp('Calculating modified Z-scores');
            calculateModZ(obj)            
            disp('Saving results');
            saveOutliers(obj)
            disp('Done!');
            disp(' ');
            
        end
        
    end
   
    methods(Access = public)
              
        function obj = SOLID(varargin)
            
            if nargin == 0
                obj.initializeGUI;
            end
            
            if nargin ~= 0

                p = inputParser;
                checkDWI = @(x) exist(x, 'file') && contains(x, '.nii');                
                checkMask = @(x) exist(x, 'file') && contains(x, '.nii');
                defaultMetric = 'Variance';                
                expectedMetrices = {'Variance', 'Mean', 'Iod'};                
                p.addParameter('in', [], checkDWI);                
                p.addParameter('bval', [], @(x) exist(x, 'file'));
                p.addParameter('bvec', [], @(x) exist(x, 'file'));
                p.addParameter('mask', [], checkMask);
                p.addParameter('metric', defaultMetric, @(x) any(validatestring(x, expectedMetrices)));
                p.addParameter('thrL', 3.5, @(x) isnumeric(x) && (x >= 0));
                p.addParameter('thrU', 10.0, @(x) isnumeric(x) && (x >= 0));
                p.parse(varargin{:});                
                                
                obj.runCMD(p.Results);
                
                
            end
            
            if nargout == 0
                clear obj
            end
            
        end
              
    end
    
end