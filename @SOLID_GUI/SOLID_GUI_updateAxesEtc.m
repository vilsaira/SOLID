function solidGui = SOLID_GUI_updateAxesEtc(solidGui, source, event)    
    shell = str2double(solidGui.popup.shell.String(solidGui.popup.shell.Value, :));
    shellInds = solidGui.solid.bval == shell;
    shellNums = zeros(size(shellInds));
    shellNums(shellInds) = 1:sum(shellInds);
    
    currentDWI = round(solidGui.slider.dwi.Value) == shellNums;
    dwiInd = find(currentDWI,1,'first');
    dwiInd2 = solidGui.data.nearestVecs(dwiInd,2);
    dwiInd3 = solidGui.data.nearestVecs(dwiInd,3);
    currentAXI = round(solidGui.slider.axi.Value);
    currentSAG = round(solidGui.slider.sag.Value);
    currentCOR = round(solidGui.slider.cor.Value);
    
    limitsDWI = [solidGui.slider.dwi.Min, solidGui.slider.dwi.Max];
    limitsAXI = [solidGui.slider.axi.Min, solidGui.slider.axi.Max];
    limitsCOR = [solidGui.slider.cor.Min, solidGui.slider.cor.Max];
    limitsSAG = [solidGui.slider.sag.Min, solidGui.slider.sag.Max];
    
    imgSAG = rot90(squeeze(solidGui.solid.DWI(:, currentSAG, :, currentDWI)));
    imgCOR = rot90(squeeze(solidGui.solid.DWI(currentCOR, :, :, currentDWI)));
    imgAXI = squeeze(solidGui.solid.DWI(:,:,currentAXI, currentDWI));

    if dwiInd > 1 && shellInds(dwiInd-1)
        dwiIndM = dwiInd -1;
    else
        dwiIndM = [];
    end
    if dwiInd < size(solidGui.solid.DWI,4) && shellInds(dwiInd+1)
        dwiIndP = dwiInd+1;
    else
        dwiIndP = [];
    end
    imgAXI_LMK = squeeze(solidGui.solid.DWI(:,:,currentAXI, solidGui.data.nearestVecs(dwiInd,2)));
    if isempty(imgAXI_LMK)
        imgAXI_LMK = NaN(size(imgAXI_LMK,1), size(imgAXI_LMK,2),1);
    end
    imgAXI_LPK = squeeze(solidGui.solid.DWI(:,:,currentAXI, solidGui.data.nearestVecs(dwiInd,3)));
    if isempty(imgAXI_LPK)
        imgAXI_LPK = NaN(size(imgAXI_LPK,1), size(imgAXI_LPK,2),1);
    end
    imgAXI_LKM = squeeze(solidGui.solid.DWI(:,:,max([currentAXI-1, solidGui.slider.axi.Min]), currentDWI));
    imgAXI_LKP = squeeze(solidGui.solid.DWI(:,:,min([currentAXI+1, solidGui.slider.axi.Max]), currentDWI));            
    if solidGui.togglebtn.mask.Value
        imgSAG = imfuse(imgSAG, rot90(squeeze(solidGui.solid.mask(:, currentSAG, :))));
        imgCOR = imfuse(imgCOR, rot90(squeeze(solidGui.solid.mask(currentCOR, :, :))));
        imgAXI = imfuse(imgAXI, squeeze(solidGui.solid.mask(:,:,currentAXI)));
        imgAXI_LMK = imfuse(imgAXI_LMK, squeeze(solidGui.solid.mask(:,:,currentAXI)));
        imgAXI_LPK = imfuse(imgAXI_LPK, squeeze(solidGui.solid.mask(:,:,currentAXI)));
        imgAXI_LKM = imfuse(imgAXI_LKM, squeeze(solidGui.solid.mask(:,:,max([currentAXI-1, solidGui.slider.axi.Min]))));
        imgAXI_LKP = imfuse(imgAXI_LKP, squeeze(solidGui.solid.mask(:,:,min([currentAXI+1, solidGui.slider.axi.Max]))));
    end
    
    solidGui.img.dwiTopLeft.CData = flipud(imgCOR);
    set(solidGui.ax.dwiTopLeft, 'Ydir', 'normal', 'Clim', solidGui.data.clim);  
    solidGui.lines.sagAXI.XData = limitsCOR;
    solidGui.lines.sagAXI.YData = currentAXI.*[1,1]+0.5;
    solidGui.lines.sagCOR.XData = currentCOR.*[1,1]+0.5;
    solidGui.lines.sagCOR.YData = limitsAXI;

    solidGui.img.dwiTopRight.CData = flipud(imgSAG);
    set(solidGui.ax.dwiTopRight, 'Ydir', 'normal', 'Clim', solidGui.data.clim);  
    solidGui.lines.corAXI.XData = limitsSAG;
    solidGui.lines.corAXI.YData = currentAXI.*[1,1]+0.5;
    solidGui.lines.corSAG.XData = currentSAG.*[1,1]+0.5;
    solidGui.lines.corSAG.YData = limitsAXI;
    
    solidGui.img.dwiMiddleCenter.CData = flipud(imgAXI);
    set(solidGui.ax.dwiMiddleCenter, 'Clim', solidGui.data.clim);  
    solidGui.lines.axiCOR.XData = limitsCOR;
    solidGui.lines.axiCOR.YData = currentCOR.*[1,1]+0.5;
    solidGui.lines.axiSAG.XData = currentSAG.*[1,1]+0.5;
    solidGui.lines.axiSAG.YData = limitsSAG;                        
    
    solidGui.img.dwiBottomCenter.CData = flipud(imgAXI_LKM);
    set(solidGui.ax.dwiBottomCenter, 'Clim', solidGui.data.clim);             
    
    solidGui.img.dwiTopCenter.CData = flipud(imgAXI_LKP);
    set(solidGui.ax.dwiTopCenter, 'Clim', solidGui.data.clim); 
    
    solidGui.img.dwiMiddleLeft.CData = flipud(imgAXI_LMK);
    set(solidGui.ax.dwiMiddleLeft, 'Clim', solidGui.data.clim);             
    
    solidGui.img.dwiMiddleRight.CData = flipud(imgAXI_LPK);
    set(solidGui.ax.dwiMiddleRight, 'Clim', solidGui.data.clim);

    % qSpace plot & modZ2D plot 
    modZclim = [solidGui.solid.thresholdLower, solidGui.solid.thresholdUpper];
    
    if solidGui.data.qSpace3D_updateFlag
%         solidGui.data.qSpace3D = solidGui.SOLID_GUI_sphericalInterpolation(solidGui.solid.bvec, solidGui.solid.bval, solidGui.solid.modZ, shellInds, currentAXI);        
% 
        bins = 0.1:0.2:modZclim(2);
        tmp = solidGui.solid.modZ(:, shellInds);        
        [N, edges] = histcounts(tmp(:), bins);
        solidGui.lines.modZHist.Data = tmp(:);
        solidGui.lines.modZHist.BinEdges = edges;
        solidGui.lines.modZHist.BinLimits = [min(edges), max(edges)];
        solidGui.data.qSpace3D_updateFlag = false;
    end    
%     solidGui.lines.qSpaceVecs.XData = [solidGui.solid.bvec(shellInds,1); -solidGui.solid.bvec(shellInds,1)];
%     solidGui.lines.qSpaceVecs.YData = [solidGui.solid.bvec(shellInds,2); -solidGui.solid.bvec(shellInds,2)];
%     solidGui.lines.qSpaceVecs.ZData = [solidGui.solid.bvec(shellInds,3); -solidGui.solid.bvec(shellInds,3)];
%     solidGui.lines.qSpaceVecsCur.XData = [solidGui.solid.bvec(dwiInd,1); -solidGui.solid.bvec(dwiInd,1)];
%     solidGui.lines.qSpaceVecsCur.YData = [solidGui.solid.bvec(dwiInd,2); -solidGui.solid.bvec(dwiInd,2)];
%     solidGui.lines.qSpaceVecsCur.ZData = [solidGui.solid.bvec(dwiInd,3); -solidGui.solid.bvec(dwiInd,3)];
%     solidGui.lines.qSpaceVecsNext1.XData = [solidGui.solid.bvec(dwiInd2,1); -solidGui.solid.bvec(dwiInd2,1)];
%     solidGui.lines.qSpaceVecsNext1.YData = [solidGui.solid.bvec(dwiInd2,2); -solidGui.solid.bvec(dwiInd2,2)];
%     solidGui.lines.qSpaceVecsNext1.ZData = [solidGui.solid.bvec(dwiInd2,3); -solidGui.solid.bvec(dwiInd2,3)];
%     solidGui.lines.qSpaceVecsNext2.XData = [solidGui.solid.bvec(dwiInd3,1); -solidGui.solid.bvec(dwiInd3,1)];
%     solidGui.lines.qSpaceVecsNext2.YData = [solidGui.solid.bvec(dwiInd3,2); -solidGui.solid.bvec(dwiInd3,2)];
%     solidGui.lines.qSpaceVecsNext2.ZData = [solidGui.solid.bvec(dwiInd3,3); -solidGui.solid.bvec(dwiInd3,3)];
    
%     solidGui.lines.qSpaceSurf.XData = 0.98*solidGui.data.qSpace3D.Xi;
%     solidGui.lines.qSpaceSurf.YData = 0.98*solidGui.data.qSpace3D.Yi;
%     solidGui.lines.qSpaceSurf.ZData = 0.98*solidGui.data.qSpace3D.Zi;
%     solidGui.lines.qSpaceSurf.CData = solidGui.data.qSpace3D.Ci;
%     shading(solidGui.ax.qSpace, 'interp');
    
%     solidGui.ax.qSpace.CLim = [solidGui.solid.thresholdLower, solidGui.solid.thresholdUpper];
%     [az, el, ~] = cart2sph(solidGui.solid.bvec(dwiInd, 1), solidGui.solid.bvec(dwiInd, 2), solidGui.solid.bvec(dwiInd, 3));
%     az = az./pi.*180;
%     el = el./pi.*180;
%     az = az - 90;
%     az(az < 0) = az + 360;
%     el = -el;
%     solidGui.ax.qSpace.View = [az, el];

    solidGui.img.modZ2D.CData = solidGui.solid.modZ(:, shellInds);
    set(solidGui.ax.modZ2D, 'YDir', 'normal', 'CLim', modZclim);
    set(solidGui.ax.modZ2D_colorbar, 'CLim', modZclim);
    % axis(solidGui.ax.modZ2D, 'equal');
    solidGui.ax.modZ2D.Position = [0, 0, 1, 1];
    solidGui.lines.modZAXI.XData = [0.5, sum(shellInds)+0.5];
    solidGui.lines.modZAXI.YData = currentAXI.*[1,1];            
    solidGui.lines.modZDWI.XData = shellNums(dwiInd).*[1,1];
    solidGui.lines.modZDWI.YData = limitsAXI+[-0.5, 0.5];
    
    solidGui.lines.modZG.XData = shellNums(dwiInd2);
    solidGui.lines.modZG.YData = currentAXI;
    solidGui.lines.modZM.XData = shellNums(dwiInd3);
    solidGui.lines.modZM.YData = currentAXI;

    curModZ = solidGui.solid.modZ(currentAXI, currentDWI);
    solidGui.lines.modZHistCurrent.XData = curModZ.*[1,1];
    solidGui.lines.modZHistCurrent.YData = [1,max(solidGui.lines.modZHist.Values)];
    solidGui.lines.modZHistLower.XData = modZclim(1).*[1,1];
    solidGui.lines.modZHistLower.YData =  solidGui.lines.modZHistCurrent.YData;
    solidGui.lines.modZHistUpper.XData = modZclim(2).*[1,1];
    solidGui.lines.modZHistUpper.YData =  solidGui.lines.modZHistCurrent.YData;
    solidGui.lines.modZHistLegend.String = {'Upper threshold', 'Lower threshold', ['Current ', num2str(curModZ,2)]};
    solidGui.ax.modZhistogram.XLim = [-0.01, max([modZclim(2), curModZ])+0.5];

    val = (curModZ - modZclim(1))./(modZclim(2)-modZclim(1));
    if curModZ > modZclim(2)
        val = 1;
    elseif curModZ < modZclim(1)
        val = 0;
    end
    solidGui.lines.modZColorbarLine.YData = val.*[1,1];
    
end