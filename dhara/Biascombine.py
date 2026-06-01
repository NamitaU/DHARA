import numpy as np
import math
import matplotlib.pyplot as plt
from photutils.centroids import centroid_com, centroid_quadratic
from photutils.centroids import centroid_1dg, centroid_2dg
from astropy.io import fits
import os
import fnmatch
from lmfit import Parameters,minimize, fit_report
from scipy.stats import chisquare
from astropy.io.fits import getheader
from astropy.stats import sigma_clipped_stats
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit
from scipy import interpolate
import natsort
from photutils.utils import calc_total_error
from scipy import ndimage
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from astropy.table import hstack
from astropy.table import Column
from photutils.segmentation import make_source_mask
from photutils.datasets import make_4gaussians_image
from astropy.stats import sigma_clipped_stats
from photutils.aperture import CircularAperture, ApertureStats
import wget




def biasCombine(biaspath, file_prefix='bias'):
    """
    combine raw bias images to generate masterbias image.
    Parameters
    ----------
    biaspath: str
        Path to the directory containing the raw bias frames.
    file_prefix: str
        Prefix used to identify bias files when the FITS header does
        not contain object/type information. Examples include
        'bias', 'Bias', or similar naming conventions.
        Default is 'bias'.
    """   
    img_list = natsort.natsorted(os.listdir(biaspath))
    Bias = []
    for files in img_list:
        filepath = os.path.join(biaspath, files)
        # Skip non-fits files
        if not files.endswith('.fits'):
            continue
        use_file = False
        try:
            header = getheader(filepath)
            # Try common header keywords
            objname = (str(header.get('OBJECT', '')).strip().lower())
            # Check if header says bias
            if 'bias' in objname:
                use_file = True
        except Exception as e:
            print(f"Could not read header for {files}: {e}")
        # Fallback to filename prefix
        if not use_file:
            if fnmatch.fnmatch(files.lower(), file_prefix.lower() + '*.fits'):
                use_file = True
        if use_file:
            Bias.append(files)
    n = len(Bias)
    if n == 0:
        raise ValueError("No bias frames found.")
    print(f"Found {n} bias frames")
    # Read first image to determine dimensions
    img1 = fits.getdata(os.path.join(biaspath, Bias[0]))
    dim = img1.shape
    l_datacube = len(dim)
    # Allocate array
    if l_datacube == 3:
        biasarr = np.zeros((n, dim[1], dim[2]))
    else:
        biasarr = np.zeros((n, dim[0], dim[1]))
    # Fill array
    for i in range(n):
        img = fits.getdata(os.path.join(biaspath, Bias[i]))
        if l_datacube == 3:
            biasarr[i, :, :] = img[0, :, :]
        else:
            biasarr[i, :, :] = img
    # Use header from first bias frame
    header = getheader(os.path.join(biaspath, Bias[0]))
    # Median combine
    cmb = np.median(biasarr, axis=0)
    # Write master bias
    F = os.path.join(biaspath, 'masterbias.fits')
    fits.writeto(F, cmb, header, overwrite=True)
    print(f"Master bias written to: {F}")
    return cmb
    
    
# Here the script search for OBJECT keyword in header if it has bias in it. Please change the setting according to your header keyword and bias label in line 45. 

'''
def biasCombine(biaspath, file_prefix ):
    img_list = os.listdir(biaspath)
    img_list = natsort.natsorted(img_list)
    Bias = []
    for files in img_list:
        if fnmatch.fnmatch(files, file_prefix+'*.fits'):
            Bias.append(files)
    n = len(Bias)
    print(n)
    img1 = fits.getdata(os.path.join(biaspath, Bias[0]))
    dim = img1.shape
    l_datacube = len(img1.shape)
    if l_datacube==3:
        biasarr = np.zeros((n,dim[1] ,dim[2]))#1300, 1340))
    else:
        biasarr = np.zeros((n,dim[0] ,dim[1]))
    for i in range(n):
        img = fits.getdata(os.path.join(biaspath, Bias[i]))
        if l_datacube ==3:
            biasarr[i,:,:] = img[0,:,:]
        else:
            biasarr[i,:,:] = img
    header = getheader(os.path.join(p, Bias[1]))
    cmb=np.median(biasarr,axis=0)      # first median combine then normalize with median
    print(cmb.shape, cmb.dtype)
    F = os.path.join(biaspath, 'masterbias.fits')
    fits.writeto(F, cmb, header)
    return()
'''
'''       
p = '/Volumes/namita1TB/PDF/ISM/open_clusters/Observations/2025NovAIMPOL/15Nov2025/Bias'
masterbias = biasCombine(biaspath, "Bias" )
'''