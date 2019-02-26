classdef SOLID < handle
% SOLID - Slicewise outlier detection for diffusion weighted MRI data.
%
% @SOLID class contains all COMMAND LINE utilities. 
%
% Allowed arguments are:
% 'in', 'the path to DWI file'.
% 'bval', 'the path to bval file' (optional).
% 'bvec', 'the path to bvec file' (optional).
% 'mask', 'the path to mask NIfTI' (optional).
% 'maskTune', 'custom masking parameter value' (optional, default 0.5).
% 'maskFilter', 'custom masking parameter value' (optional, default 7).
% 'metric', 'Variance/Mean/Iod' (optional, default Variance).
% 'thrL', 'Lower threshold value' (optional, default 3.5).
% 'thrU', 'Upper threshold value' (optional, default 10.0).
% 'save', 'true/false' (optional, default true).
%
% Examples
% 
% 1) You want to generate a SOLID object for custom usage. You have
% DWI.nii.gz, DWI.bval, DWI.bvec files in the Matlab path and you want to
% use SOLID to estimate the brain mask. Simply write
%
% > s = SOLID('in', 'DWI.nii.gz', 'save', false);
%
% and observer results in s-object. Nothing is written on the hard drive.
%
% 2) You want to analyze a DWI with custom mask and save results on the
% hard drive. 
%
% > SOLID('in', 'DWI.nii.gz', 'mask', '/pathToMask/Mask.nii.gz');
%
% Results are saved on the hard drive to the same folder where DWI.nii.gz
% file is.
%
%
% Automated masking requires ExploreDTI function
%   E_DTI_Create_Mask_From_DWI_enhanced_IND.p
% **********************  Author: Viljami Sairanen  ***********************
% *********************  viljami.sairanen@gmail.com  **********************
% *************  Website: https://github.com/vilsaira/SOLID  **************
% *********************  Last edited: 23 October 2018 *********************
    properties (Access = public)
        DWI 
        dims 
        modZ
        bval 
        uniqueb 
        bvec 
        mask
        useMask
        maskTune
        maskFilter
        metric  
        thresholdLower
        thresholdUpper
        minmax
        fName 
        fDir 
        fullPath 
    end

    methods (Access = private)
        solid = SOLID_processData(obj, input);
    end

    methods (Static, Access = private)
        mask = E_DTI_Create_Mask_From_DWI_enhanced_IND(b0img, val1, val2);
        nii = SOLID_loadNifti(fpath);
        bval = SOLID_loadBVal(fpath);
        bvec = SOLID_loadBVec(fpath);         
        SOLID_saveTxtModZ2D(oname, modZ);
        SOLID_GUI_saveManualLabels(oname, modZ, modZorig);
        SOLID_save4DNIfTI(oname, modZ, fpath, fname);
        SOLID_saveUntouchNii(oname, modZ, fpath, fname);
        SOLID_savePNG(oname, modZ, b, thresholdLower, thresholdUpper);
        SOLID_saveMatFile(oname, modZ, fDir, fName);
    end
    
    methods (Static, Access = public)
        oname = SOLID_save(obj);
        modZ = SOLID_calculateModZ(DWI, b, metric, mask);
    end

    methods (Access = public)
        function solid = SOLID(varargin)
            if nargin ~= 0
                p = inputParser;
                checkDWI = @(x) exist(x, 'file') && (contains(x, '.nii') || contains(x, '.mat'));
                checkMask = @(x)  isempty(x) || strcmp(x, 'no') || (exist(x, 'file') && contains(x, '.nii'));
                defaultMetric = 'Variance';
                expectedMetrices = {'Variance', 'Mean', 'Iod'};
                p.addParameter('in', [], checkDWI);
                p.addParameter('bval', [], @(x) exist(x, 'file'));
                p.addParameter('bvec', [], @(x) exist(x, 'file'));
                p.addParameter('mask', [], checkMask);
                p.addParameter('maskTune', 0.5, @(x) isnumeric(x) && x > 0 && x <= 1);
                p.addParameter('maskFilter', 7, @(x) bitget(x,1) && x > 0); % bitget(x,1) returns true only for odd numbers
                p.addParameter('metric', defaultMetric, @(x) any(validatestring(x, expectedMetrices)));
                p.addParameter('thrL', 3.5, @(x) isnumeric(x) && (x >= 0));
                p.addParameter('thrU', 10.0, @(x) isnumeric(x) && (x >= 0));
                p.addParameter('save', true, @(x) islogical(x));
                p.parse(varargin{:});
                solid.SOLID_processData(p.Results);
            end

            if nargout == 0
                clear solid;
            end
        end
    end

end