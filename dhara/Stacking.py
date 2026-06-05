import numpy as np
import math
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import os
import fnmatch
from astropy.io.fits import getheader
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit
from scipy import interpolate
from photutils.centroids import centroid_com, centroid_quadratic
from photutils.centroids import centroid_1dg, centroid_2dg
import natsort
from astropy.wcs import WCS
from photutils.utils import calc_total_error
from astropy.stats import SigmaClip
from scipy import ndimage
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from astropy.table import hstack
from astropy.table import Column
from photutils.segmentation import make_source_mask
import wget
from astroquery.vizier import Vizier
import re
import astroalign as aa

def COM(image,x,y,r,s):						## Function to find the centre of mass of an image
  x = int(round(x))
  y = int(round(y))
  #image = image - np.median(image)
  subimage = image[y-r:y+s,x-r:x+s]				## Prod Subimge containing only the star of interest by providing appx centres
  i=x-r 							## (0,0) location of subimage , later to be added to COM coord
  j=y-r
  p,q = ndimage.maximum_position(subimage)		## Getting values of actual centres of the subimage  measurements.center_of_mass
  #print "COM of subimage (x,y):", p,q
  a= q+i							## Getting information about the actual coordinates of that particular centre
  b= p+j							##added opposite since COM fn gives coordinates as (y,x)
  #print "Actual centres (x,y) :", a,b
  return (a,b)

def COM_iter(image,a,b,p,q):
  while True:
     a1,b1 = COM(image,a,b,p,q)
     if(math.sqrt((a1-a)**2+(b1-b)**2)<0.1):
        break
     else:
        a = a1
        b = b1
  return a1,b1

def cent(Image,x,y,h,k ):
    x = int(round(x))
    y = int(round(y))
    sub = Image[y-k:y+k, x-h:x+h]
    xs, ys = centroid_quadratic(sub)
    xp, yp = x-h+xs, y-k+ys
    return(xp,yp)
    
def cent_2g(Image,x,y,h,k ):
    x = int(round(x))
    y = int(round(y))
    sub = Image[y-k:y+k, x-h:x+h]
    xs, ys = centroid_2dg(sub)
    xp, yp = x-h+xs, y-k+ys
    return(xp,yp)

def single_click(event):
    global X1,Y1
    if event.button == 3:
        X1, Y1 = event.xdata, event.ydata
        plt.close(1)
def press(event):
    global X2,Y2
    if event.key == 'x':
        X2, Y2 = event.xdata, event.ydata
        plt.close(1)








def classify_hwp_images(datapath, band, exposure , hwp=['p0', 'p1', 'p2', 'p3'], ang_map=[67.5, 45, 22.5, 0]):
    """
    This function identifies and groups polarimetric images based on
    the HWP position keywords, filter band, and exposure time.
    Parameters
    ----------
    datapath : str
        Path to the directory containing the raw polarimetric images.
    band : str
        Photometric filter used during the observations
        (e.g., 'B', 'R', 'V', 'I').
    exposure : int or float
        Exposure time (in seconds) used to select the same exposure images.
    hwp : list of str, optional
        List of HWP position identifiers used in the filenames or
        FITS headers. Default is ['p0', 'p1', 'p2', 'p3'].
    ang_map : list of float, optional
        Polarization angles (in degrees) corresponding to each HWP
        position. The mapping order should match the `hwp` list.
        Default is p0 = 67.5, p1 = 45, p2 = 22.5, p3 = 0, so [67.5, 45, 22.5, 0].
    """
    img_list = os.listdir(datapath)
    img_list = natsort.natsorted(img_list)
    # Mapping between filename tags and HWP angles
    hwp_map = {
        hwp[0]: ang_map[0],
        hwp[1]: ang_map[1],
        hwp[2]: ang_map[2],
        hwp[3]: ang_map[3]
    }
    # Required HWP positions
    hwp_positions = [0.0, 22.5, 45.0, 67.5]
    # Dictionary to store classified files
    # Final structure:
    # objects[exptime][hwp] = [files]
    objects = {}
    for files in img_list:
        if not files.endswith('.fits') or files.startswith("._"):
            continue
        filepath = os.path.join(datapath, files)
        hwp_value = None
        filter_value = None
        exptime = None
        cluster = None
        # ---------------------------------------------------
        # Read header
        # ---------------------------------------------------
        try:
            #header = getheader(filepath, ignore_missing_simple=True)
            with fits.open(filepath, memmap=False) as hdul:
                header = hdul[0].header
            # HWP keywords
            for key in ['HWP', 'HWPANG', 'HWP_POS']:
                if key in header:
                    hwp_value = round(float(header[key]), 1)
                    break
            # Exposure time
            for key in ['EXPTIME', 'EXPOSURE']:
                if key in header:
                    EXPtime = round(float(header[key]), 2) #in AIMPOL it is in millisecon - change next step based on units
                    exptime = int(EXPtime/1000) # converted in sec
                    break
            # Filter/Band
            for key in ['FILTER', 'BAND', 'FILTER1']:
                if key in header:
                    filter_value = str(header[key]).strip()
                    break
            for key in ['OBJECT', 'OBJTYP']:
                if key in header:
                    cluster = str(header[key]).strip()
        except Exception as e:
            print(f'Header read failed for {files}: {e}')
        
        fname = re.split(r'[_.\s]+', files)
        # ---------------------------------------------------
        # Fallback HWP from filename
        # ---------------------------------------------------
        if hwp_value is None:
            for tag, angle in hwp_map.items():
                if tag in fname:
                    hwp_value = angle
                    break
        # ---------------------------------------------------
        # Fallback filter from filename
        # ---------------------------------------------------
        if filter_value is not None:
            if filter_value != band:
                continue
        else:
            if band not in fname:
                continue
        band_match = True
        if cluster is None:
            cluster = fname[0]
        if exptime is None:
            for t in fname:
                if t.lower().endswith('s') and t[:-1].isdigit():
                    exptime = int(t[:-1]) 
                    break
        if exptime != exposure:
            continue
        else:
            exp_match = True
        # ---------------------------------------------------
        # Store only valid files
        # ---------------------------------------------------
        if (hwp_value in hwp_positions and
                exp_match and
                band_match):
            # Create exposure group if needed
            if hwp_value not in objects:
                objects[hwp_value] = []
            objects[hwp_value].append(files)
    return objects


def ShiftnStack(datapath, biaspath,  band, exposure ,  objects, method, ang_map, spath): # method for stacking sum,or median
    """
    Align, shift, and stack polarimetric images to improve the
    signal-to-noise ratio for further analysis.
    Parameters
    ----------
    datapath : str
        Path to the directory containing the science images.
    biaspath : str
        Path to the directory containing the masterbias. 
    band : str
        Photometric filter used during the observations
        (e.g., 'B', 'V', 'R', 'I').
    exposure : int or float
        Exposure time (in seconds) used to select the images.
    objects : list
        List of target object names or identifiers to be processed.
        output of classify_hwp_images function
    method : str
        Image stacking method to use (e.g., 'median' or 'sum').
    ang_map : list of float
        List of polarization angles corresponding to the HWP positions.
    spath : str
        Output directory where the aligned and stacked images will be
        saved.
    """
    hwp = ang_map
    masterbias = fits.getdata(os.path.join(biaspath, 'masterbias.fits'))
    # -----------------------------
    # reference image (first valid)
    # -----------------------------
    ref_file = objects[hwp[0]][0]
    IMG, head = fits.getdata(os.path.join(datapath, ref_file), header= True)
    IMG = IMG.astype('float64')
    dim = IMG.shape
    if IMG.ndim == 3:
        IMG = IMG[0, :, :]
        xdim, ydim = dim[1], dim[2]
    else:
        IMG = IMG
        xdim, ydim = dim[0], dim[1]
    IMG -= masterbias   
    #-----------------------------
    # Select a reference location 
    #----------------------------
    mean, med, std = sigma_clipped_stats(IMG)
    plt.imshow(IMG,cmap='gray', vmin = med-10*std, vmax = med+10*std)
    plt.colorbar()
    cmd=plt.connect('key_press_event', press)
    plt.title('press x on the center of an isolated-bright object'+ ref_file )
    plt.show()
    ax0, by0 = cent(IMG,X2,Y2,10, 10 )
    stacked = []
    shifted = []
    for k in range(4):        # loop on half waveplate position
        files = objects[hwp[k]]
        n = len(files)
        imgarr = np.zeros((n, xdim, ydim ))
        j = 0
        shifted_hwp = []
        for f in files:
            img = fits.getdata(os.path.join(datapath, f)).astype('float64')
            header = fits.getheader(os.path.join(datapath, f))
            if img.ndim == 3:
                img = img[0, :, :]
            img -= masterbias
            try:
                transf, _ = aa.find_transform(img, IMG) # (current, reference)
                IMGS,_ = aa.apply_transform(transf, img, IMG)
                dx, dy = transf.translation
                shift = np.sqrt(dx**2 + dy**2)
                if shift > 10:   # choose safe threshold
                    raise RuntimeError("Large shift")
            except:
                plt.imshow(img,cmap='gray', vmin = med-10*std, vmax = med+10*std)
                plt.colorbar()
                plt.title('click on the center of same star '+f)
                cmd=plt.connect('button_release_event', single_click)
                plt.show()
                apx, bpy = cent(img,X1,Y1,10,10)
                shiftX, shiftY = ax0-apx, by0-bpy
                IMGS=ndimage.shift(img, [shiftY, shiftX], order=3, mode='constant', cval=0.0, prefilter=False) # mode = 'nearest'
            imgarr[j,:,:]= IMGS  #ndimage.shift(img, [shiftY, shiftX], order=3, mode='constant', cval=0.0, prefilter=False) # mode = 'nearest'
            #------------------------
            # Saved shifted images 
            #------------------------
            shifted_file = os.path.join(spath, 'shifted_'+f)
            shifted_hwp.append('shifted_'+f)
            try:
                fits.writeto(shifted_file, imgarr[j,:,:], head)
            except:
                new_hdr = fits.Header()
                for i, card in enumerate(head.cards[:15]):
                    new_hdr.append(card)
                fits.writeto(shifted_file, imgarr[j,:,:], new_hdr)
            j = j+1
            #Cx = actual_Cx
            #Cy = actual_Cy 
        #--------------------
        # Stack and save 
        #-------------------    
        if(method == 'sum'):
            clipped = sigma_clip(imgarr, sigma=3, axis=0) # sigma clipped for the cosmic rays/bad pixels
            cmb = np.ma.sum(clipped, axis=0)
            cmb = cmb.filled(np.nan)
        elif(method == 'median'):
            cmb = np.median(imgarr, axis = 0)
        else:
            raise ValueError("method for stacking should be either 'sum' or 'median' ")
        #print(cmb.dtype)
        stacked_file = os.path.join(spath, 'stacked_'+f)
        try:
            fits.writeto(stacked_file,cmb, head)
        except:
            fits.writeto(stacked_file,cmb, new_hdr)
            
        stacked.append(['stacked_'+f])
        shifted.append(shifted_hwp) 
    return(stacked, shifted)

