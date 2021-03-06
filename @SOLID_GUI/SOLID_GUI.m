classdef SOLID_GUI < handle
%% The 2-Clause BSD License
% Copyright 2018 Viljami Sairanen
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation 
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
% USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    events

    end
    
    properties (Access = public)

        ui
        scrollPanel
        solid
        ax = struct('dwiTopLeft', [], ... 
                    'dwiTopCenter', [], ... 
                    'dwiTopRight', [],... 
                    'dwiMiddleLeft', [],...
                    'dwiMiddleCenter', [],...
                    'dwiMiddleRight', [],...
                    'dwiBottomCenter', [],...
                    'qSpace', [],...                    
                    'modZHistogram', [],...
                    'modZ2D', [], ...
                    'modZ2D_colorbar', []);                
        img = struct('dwiTopLeft', [], ...
                     'dwiTopCenter', [], ...
                     'dwiTopRight', [],...
                     'dwiMiddleLeft', [],...
                     'dwiMiddleCenter', [],...
                     'dwiMiddleRight', [],...
                     'dwiBottomCenter', [],...
                     'qSpace', [],...                    
                     'modZHistogram', [],...
                     'modZ2D', [], ...
                     'modZ2D_colorbar', []);
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
                       'thresholdUpper', [], ...
                       'maskTune', [], ...
                       'maskFilter', []);
        togglebtn = struct('outlier', [], ...
                       'mask', []);
        data = struct('modZorig', [],...
                       'clim', [0,1],...
                       'window', [], ...
                       'nearestVecs', [],...
                       'qSpace3D', [],...
                       'qSpace3D_updateFlag', true,...
                       'maskTune', 0.5,...
                       'maskFilter', 7,...
                       'metric', 'Variance');
        lines = struct('axiSAG', [],...
                         'axiCOR', [],...
                         'corSAG', [],...
                         'corAXI', [],...
                         'sagCOR', [],...
                         'sagAXI', [],...
                         'modZAXI', [],...
                         'modZDWI', [],...
                         'modZM', [],...
                         'modZG', [],...
                         'qSpaceVecs', [],...
                         'qSpaceVecsCur', [],...
                         'qSpaceVecsNext1', [],...
                         'qSpaceVecsNext2', [],...
                         'qSpaceSurf', [],...
                         'modZHist', [],...
                         'modZHistLower', [],...
                         'modZHistUpper', [],...
                         'modZHistCurrent', [],...
                         'modZHistLegend', [],...
                         'dwiHistogramVer', [],...
                         'dwiHistogramHor', [],...
                         'modZColorbarLine',[]);        
        popup = struct('shell', []);
        menuFile;
        menuFileOpen;
        menuFileSave;
        menuSetting;
        menuSettingMetric;
        menuSettingMetricVar;
        menuSettingMetricMean;
        menuSettingMetricIod;
        menuSettingMask;
        menuHelp;
        menuHelpTheory;
        menuHelpGithub;
    end
    
    methods (Access = private)
        solidGui = SOLID_GUI_initialize(obj);
                
        mask = SOLID_GUI_loadMask(obj, dwiPath, dwiName);
        out = SOLID_GUI_loadBValVec(obj, dwiPath, dwiName, ftype);

        solidGui = SOLID_GUI_callbackLoadData(obj, event, source);
        solidGui = SOLID_GUI_callbackTextfield(obj, event, source);
        solidGui = SOLID_GUI_callbackSlider(obj, event, source); 
        solidGui = SOLID_GUI_callbackToggleMetric(obj, event, source);
        solidGui = SOLID_GUI_callbackSelectShell(obj, event, source);
        solidGui = SOLID_GUI_callbackToggleUseMask(obj, event, source);
        solidGui = SOLID_GUI_callbackShowMask(obj, event, source);   
        solidGui = SOLID_GUI_callbackMouseClick(obj, event, source);
        solidGui = SOLID_GUI_callbackMouseRelease(obj, event, source);
        solidGui = SOLID_GUI_callbackMouseMove(obj, event, source);
        solidGui = SOLID_GUI_callbackWindowTimer(obj, event, source);
        solidGui = SOLID_GUI_callbackKeyPress(obj, event, source);
        solidGui = SOLID_GUI_callbackKeyRelease(obj, event, source);
        solidGui = SOLID_GUI_callbackSpacebar(obj, event, source);
        solidGui = SOLID_GUI_callbackLeftArrow(obj, event, source);
        solidGui = SOLID_GUI_callbackRightArrow(obj, event, source);
        solidGui = SOLID_GUI_callbackMouseScroll(obj, event, source);        
        solidGui = SOLID_GUI_callbackSave(obj, event, source);        
    end

    methods (Static, Access = private)
        mask = E_DTI_Create_Mask_From_DWI_enhanced_IND(b0img, val1, val2);
        nii = SOLID_loadNifti(fpath);
        bval = SOLID_loadBVal(fpath);
        bvec = SOLID_loadBVec(fpath);        
        qSpace3D = SOLID_GUI_sphericalInterpolation(bvec, bval, modZ, shellInds, currentAXI);
        SOLID_saveTxtModZ2D(oname, modZ);
        SOLID_GUI_saveManualLabels(oname, modZ, modZorig);
        SOLID_save4DNIfTI(oname, modZ, fpath, fname);
        SOLID_saveUntouchNii(oname, modZ, fpath, fname);
        SOLID_savePNG(oname, modZ, b, thresholdLower, thresholdUpper);        
    end
    
    methods (Static, Access = public)
        modZ = SOLID_calculateModZ(DWI, b, metric, mask);         
        aHood = SOLID_GUI_calculateAngularNHood(bval, bvec);
    end
    
    methods (Access = public)   
        function solidGui = SOLID_GUI(varargin)
            
            warning off;
            solidGui.SOLID_GUI_initialize;
            if nargout == 0
                warning on;
                clear solidGui;
            end
        end
        
        solidGui = SOLID_GUI_toggleEnableForFigure(obj, event, source);
        solidGui = SOLID_GUI_initAxesEtc(obj, event, source);
        solidGui = SOLID_GUI_updateAxesEtc(obj, event, source);
    end
    
end