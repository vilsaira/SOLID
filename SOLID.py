#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
SOLID - Slicewise outlier detection for diffusion weighted MRI data.

Author: Viljami Sairanen
Website: https://github.com/vilsaira/SOLID
"""

import matplotlib
matplotlib.use("Agg")
import argparse
import numpy as np
import nibabel as nib
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shutil import which
import numpy.matlib as matlib

DESCRIPTION =   "Slicewise outlier detection (SOLID) for diffusion weighted MRI data. Viljami Sairanen 2018."

class SOLID(object):

    def __init__(self):
        
        self._pathDWI = None
        self._pathBVal = None
        self._pathBVec = None
        self._pathMask = None
        self._pathEDDY = None
        self._metric = "var"
        self._pathOut = None
        self._nameOut = None
        self._printPNG = False

        self._lowerThreshold = 3.5
        self._upperThreshold = 6.0
        self._useMask = False
        self._dataDWI_hdr = None
        self._dataDWI_img = None
        self._dataMask_hdr = None
        self._dataMask_img = None
        self._dataBVal = None
        self._dataBVec = None

        self._modZ2D = None
        self._modZ4Dfirst = None

        self.argParser()
    
    def printPNG(self):
        data = self._modZ2D
        data[np.isnan(data)] = 0
        plt.rcParams["toolbar"] = "None"
        fig, ax = plt.subplots(nrows=1, ncols=1)
        fig.canvas.set_window_title(f"SOLID - {self._pathDWI}")
        fig.tight_layout
        im = ax.imshow(data, interpolation="nearest", cmap="plasma", origin="lower")
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        fig.colorbar(im, cax=cax)
        im.set_clim(self._lowerThreshold, self._upperThreshold)
        ax.set_xlabel("Image volume")
        ax.set_ylabel("Slice")
        ax.set_title("Modified Z-score")
        ax.axis("scaled")
        plt.savefig(f"{self._nameOut}_modZ2D.png", dpi=600)
        if self._printPNG:
            plt.show()

    def saveResults(self):

        self._nameOut = f"{self._pathOut}_{self._metric}"
        if self._useMask:
            self._nameOut = f"{self._nameOut}_masked"
                    
        np.savetxt(f"{self._nameOut}_modZ2D.txt", self._modZ2D, delimiter=", ")

        dims = self._dataDWI_img.shape
        img = np.zeros((dims), dtype=np.float32)
        for i in range(0, dims[2]):
            for j in range(0, dims[3]):
                val = self._modZ2D[i,j]
                img[:,:,i,j] = matlib.repmat(val, dims[0], dims[1])
        
        # EDDY needs the first volume average to be > 0
        # Replace the first volume with the original first b0 image
        # Since the first volume is not transformed, the slice information is already correct
        # and can be replaced after eddy.
        self._modZ4Dfirst = img[:,:,:,0].copy()
        if (self._pathEDDY is not None):
            img[:,:,:,0] = self._dataDWI_img[:,:,:,0]
        
        # img[self._dataMask_img == 0] = np.nan

        hdr = self._dataDWI_hdr.header.copy()
        hdr.descrip = "SOLID modified z-score maps"
        nii = nib.Nifti1Image(img, self._dataDWI_hdr.affine, hdr)
        nib.save(nii, f"{self._nameOut}_modZ4D.nii.gz")
        
        self.printPNG()

    def calculateModZ(self):        
        bvals = np.unique(self._dataBVal.round(-2))
        shape = self._dataDWI_img.shape
        modZall = np.zeros((shape[2], shape[3]))
        modZall.fill(np.nan)

        for b in bvals:
            inds = np.where(self._dataBVal.round(-2) == b)[0]
            if inds.size < 5:
                continue
            
            shell = self._dataDWI_img[:,:,:,inds].astype(np.float32)
            if self._useMask:
                mask = self._dataMask_img.astype(np.int) == 0
                shell[mask] = np.nan

            tmp = np.argwhere(np.isnan(shell))
            dims = shell.shape
            shell = shell.reshape((dims[0]*dims[1], dims[2], dims[3]))
            if self._metric == "var":
                y = np.nanvar(shell, axis=0)
            if self._metric == "mean":
                y = np.nanmean(shell, axis=0)
            if self._metric =="iod":
                y = np.nanvar(shell, axis=0) / np.nanmean(shell, axis=0)
            tmp = matlib.repmat(np.nanmedian(y, axis=1), dims[3], 1).T
            MAD = 1.4826 * np.nanmedian( np.abs(y - tmp), axis=1)
            modZ = np.abs(y-tmp) / matlib.repmat(MAD,dims[3],1).T
            modZall[:, inds] = modZ
        
        self._modZ2D = modZall

        return modZall

    def importData(self):
        # self._dataBVal, self._dataBVec = read_bvals_bvecs(self._pathBVal, self._pathBVec)
        self._dataBVal = np.loadtxt(self._pathBVal)
        self._dataDWI_hdr = nib.load(self._pathDWI)
        self._dataDWI_img = self._dataDWI_hdr.get_fdata()
        if self._useMask:
            self._dataMask_hdr = nib.load(self._pathMask)
            self._dataMask_img = self._dataMask_hdr.get_fdata()

    def applyEDDY(self):
        if os.path.exists(f"{self._pathEDDY}.eddy_command_txt"):
            with open(f"{self._pathEDDY}.eddy_command_txt", "r") as f:
                cmd = f.read().replace("\n", "").split(" ")
        new_cmd = cmd[0]
        new_cmd = f"{new_cmd} --imain={self._nameOut}_modZ4D.nii.gz"
        new_cmd = f"{new_cmd} --init={self._pathEDDY}.eddy_parameters"
        new_cmd = f"{new_cmd} --niter=0"
        new_cmd = f"{new_cmd} --out={self._pathEDDY}_solid_modZ4D.nii.gz"
        matches = ["--mask=", "--acqp=", "--index=", "--bvecs=", "--bvals", "--topup="]
        for i in range(1,len(cmd)):
            for j in range(0, len(matches)):
                if cmd[i].startswith(matches[j]):
                    new_cmd = f"{new_cmd} {cmd[i]}"
        os.system(new_cmd)
        # Return the original first volume to the EDDY corrected modZ4D
        modz_eddy = nib.load(f"{self._pathEDDY}_solid_modZ4D.nii.gz")
        img_eddy = modz_eddy.get_fdata().astype(np.float32)
        img_eddy[:,:,:,0] = self._modZ4Dfirst.astype(np.float32)
        hdr = modz_eddy.header.copy()
        hdr.descrip = "SOLID EDDY transformed modified z-score maps"
        nii = nib.Nifti1Image(img_eddy, modz_eddy.affine, hdr)
        nib.save(nii, f"{self._pathEDDY}_solid_modZ4D.nii.gz")
        # Convert modz to reliability weights
        self.modz2weights(f"{self._pathEDDY}_solid_modZ4D.nii.gz", f"{self._pathEDDY}_solid_reliability_weights.nii.gz", self._lowerThreshold, self._upperThreshold, interptype=2, k=0.2) 

    def modz2weights(self, in_solid_modz4dnii=None, out_solid_weightsnii=None, thr_low=3.5, thr_up=6, interptype=1, k=0.2):
        # interptype defines interpolation method 1 is linear and 2 is logistic (sigmoid)
        if os.path.isfile(in_solid_modz4dnii) is not None and out_solid_weightsnii is not None:
            modznii = nib.load(in_solid_modz4dnii)
            solid_weights = modznii.get_fdata()            
            solid_weights[solid_weights < thr_low] = thr_low
            solid_weights[solid_weights > thr_up] = thr_up
            if interptype == 1:
                # Linear scaling
                solid_weights = (solid_weights - thr_low) / (thr_up - thr_low)
            if interptype == 2:
                # Sigmoid scaling
                solid_weights = (solid_weights - thr_low) * 2.0 / (thr_up - thr_low) - 1.0
                solid_weights = 1 / (1 + np.exp( -solid_weights / k))
            solid_weights = 1 - solid_weights
            hdr = modznii.header.copy()
            hdr.descrip = "SOLID voxelwise reliability weights"
            nii = nib.Nifti1Image(solid_weights, modznii.affine, hdr)
            nib.save(nii, out_solid_weightsnii)

    def argParser(self):
        parser = argparse.ArgumentParser()
        requiredGroup = parser.add_argument_group("required arguments")
        requiredGroup.add_argument("--in", help="DWI file (*.nii.gz)", type=str, dest="dwi", action="store", required=True)
        parser.add_argument("--bval", help="bval file (*.bval)", type=str, dest="bval", action="store")
        parser.add_argument("--bvec", help="bvec file (*.bvec)", type=str, dest="bvec", action="store")
        parser.add_argument("--mask", help="Mask file (*.nii.gz)", type=str, dest="mask", action="store")
        parser.add_argument("--eddyout", help="Path to precalculated FSL EDDY output", type=str, dest="EDDY", action="store")
        parser.add_argument("--out", help="Output file", type=str, dest="out", action="store")
        parser.add_argument("--metric", help="Used metric (var/mean/iod)", type=str, dest="metric", action="store")
        parser.add_argument("--PNG", help="Show modified Z-score as an image", action="store_true", default=False, dest="PNG")
        parser.add_argument("--GUI", help="Launch SOLID GUI", action="store_true", default=False, dest="GUI")
        parser.add_argument("--thrU", help="Set upper modified Z-score threshold (defaut 10.0)", action="store", dest="thrU", default=6.0, type=float)
        parser.add_argument("--thrL", help="Set lower modified Z-score threshold (defaut 3.5)", action="store", dest="thrL", default=3.5, type=float)        
        args = parser.parse_args()
        
        if (args.dwi is None) and (args.GUI is None):
            sys.exit("Argument error, --in/--GUI, provide DWI or launch GUI")    

        self._pathDWI = os.path.realpath(args.dwi)        
        
        fpath, fext = args.dwi.split(os.extsep, 1)

        if args.bval is None:
            self._pathBVal = os.path.realpath(fpath + ".bval")
        else:
            self._pathBVal = os.path.realpath(args.bval)
        
        if args.bvec is None:
            self._pathBVec = os.path.realpath(fpath + ".bvec")
        else:
            self._pathBVec = os.path.realpath(args.bvec)
        
        if args.mask is not None:
            self._pathMask = os.path.realpath(args.mask)
            self._useMask = True

        if args.EDDY is not None:
            self._pathEDDY = os.path.realpath(args.EDDY)
        
        if args.out is None:
            self._pathOut = os.path.realpath(fpath + "_SOLID")
        else:
            self._pathOut = os.path.realpath(args.out)

        if args.metric is not None:            
            self._metric = args.metric

        if args.PNG:
            self._printPNG = True

        if args.thrL is not None:
            self._lowerThreshold = args.thrL

        if args.thrU is not None:
            self._upperThreshold = args.thrU

    @property
    def _pathDWI(self):
        return self.__pathDWI
    
    @_pathDWI.setter
    def _pathDWI(self, value):
        if (value is not None) and not os.path.isfile(value):
            sys.exit("Argument error, --in, cannot locate DWI file (*.nii.gz)")
        self.__pathDWI = value

    @property
    def _pathBVal(self):
        return self.__pathBVal

    @_pathBVal.setter
    def _pathBVal(self, value):
        if (value is not None) and not (os.path.isfile(value)):
            sys.exit("Argument error, --bval, cannot locate bval file (*.bval)")
        self.__pathBVal = value
    
    @property
    def _pathBVec(self):
        return self.__pathBVec

    @_pathBVec.setter
    def _pathBVec(self, value):
        if (value is not None) and not (os.path.isfile(value)):
            sys.exit("Argument error, --bvec, cannot locate bvec file (*.bvec)")
        self.__pathBVec = value

    @property
    def _pathMask(self):
        return self.__pathMask

    @_pathMask.setter
    def _pathMask(self, value):
        if (value is not None) and not (os.path.isfile(value)):
            sys.exit("Argument error, --mask, cannot locate mask file (*.nii.gz)")
        self.__pathMask = value

    @property
    def _pathEDDY(self):
        return self.__pathEDDY

    @_pathEDDY.setter
    def _pathEDDY(self, value):
        if (value is not None) and not (os.path.isfile(f"{value}.eddy_parameters")):
            sys.exit("Argument error, --eddy, cannot locate precalculated FSL EDDY results")
        self.__pathEDDY = value

    @property
    def _metric(self):
        return self.__metric

    @_metric.setter
    def _metric(self, value):        
        if (value != "var") and (value != "iod") and (value != "mean"):
            sys.exit(f"Argument error, --metric, {value} is not defined. Select mean/var/iod instead.")
        self.__metric = value
        
    @property
    def _pathOut(self):
        return self.__pathOut

    @_pathOut.setter
    def _pathOut(self, value):
        if (value is not None) and (os.path.isfile(value)):            
            sys.exit("Argument error, --out, {value} exist already")
        self.__pathOut = value

    @property
    def _printPNG(self):
        return self.__printPNG

    @_printPNG.setter
    def _printPNG(self, value):
        if value:
            self.__printPNG = True        
        else:
            self.__printPNG = False

    @property
    def _lowerThreshold(self):
        return self.__lowerThreshold
    
    @_lowerThreshold.setter
    def _lowerThreshold(self, value):
        if (value >= 0):# and (value < self._upperThreshold):
            self.__lowerThreshold = value
    
    @property
    def _upperThreshold(self):
        return self.__upperThreshold
    
    @_upperThreshold.setter
    def _upperThreshold(self, value):
        if (value > 0):# and (value > self._lowerThreshold):
            self.__upperThreshold = value

    @property
    def _useMask(self):
        return self.__useMask

    @_useMask.setter
    def _useMask(self, value):
        if value:
            self.__useMask = True
        else:
            self.__useMask = False
    
    @property
    def _dataDWI_hdr(self):
        return self.__dataDWI_hdr

    @_dataDWI_hdr.setter
    def _dataDWI_hdr(self, value):
        self.__dataDWI_hdr = value

    @property
    def _dataDWI_img(self):
        return self.__dataDWI_img
    
    @_dataDWI_img.setter
    def _dataDWI_img(self, value):
        self.__dataDWI_img = value

    @property
    def _dataMask_hdr(self):
        return self.__dataMask_hdr

    @_dataMask_hdr.setter
    def _dataMask_hdr(self, value):
        self.__dataMask_hdr = value

    @property
    def _dataMask_img(self):
        return self.__dataMask_img
    
    @_dataMask_img.setter
    def _dataMask_img(self, value):
        self.__dataMask_img = value

    @property
    def _dataBVal(self):
        return self.__dataBVal

    @_dataBVal.setter
    def _dataBVal(self, value):
        self.__dataBVal = value

    @property
    def _dataBVec(self):
        return self.__dataBVec

    @_dataBVec.setter
    def _dataBVec(self, value):
        self.__dataBVec = value

    @property
    def _modZ2D(self):
        return self.__modZ2D

    @_modZ2D.setter
    def _modZ2D(self, value):
        self.__modZ2D = value

if __name__ == "__main__":
    solid = SOLID()
    solid.importData()
    modZall = solid.calculateModZ()
    solid.saveResults()
    if solid._pathEDDY is not None:
        print("Running EDDY on SOLID results")
        solid.applyEDDY()
    sys.exit("Done")