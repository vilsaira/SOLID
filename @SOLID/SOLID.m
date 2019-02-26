classdef SOLID < handle

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