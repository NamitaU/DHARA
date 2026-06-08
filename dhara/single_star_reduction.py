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
from astropy.modeling import models, fitting
from photutils.centroids import centroid_com, centroid_quadratic
from photutils.centroids import centroid_1dg, centroid_2dg
from scipy.optimize import curve_fit
from scipy import interpolate
import natsort
from photutils.aperture import CircularAperture, ApertureStats
from astropy.wcs import WCS
from photutils.utils import calc_total_error
from astropy.stats import SigmaClip
from scipy import ndimage
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils import CircularAnnulus
from photutils import RectangularAnnulus
from astropy.table import hstack
from astropy.table import Column
from photutils.segmentation import make_source_mask
import wget
from astropy.table import Table
from astroquery.vizier import Vizier

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
    
def Gauss2D(xy, x0, y0, s, A):                                 # 2D- gaussian fiting is better than ndimage_maximum_posiyion and Center_of_mass 
    x,y=xy
    c=1/(2*s**2)
    f = A*np.exp(-c*((x-x0)**2+(y-y0)**2))
    return(f)
    
def center(subimg,r):
    x = np.arange(0, subimg.shape[1],1)
    y = np.arange(0, subimg.shape[0], 1)
    xx, yy = np.meshgrid(x,y)
    x0 = r
    y0 = r
    s = max(*subimg.shape)*0.1
    A=np.max(subimg)
    initial_guess=[x0,y0,s,A]
    param, err = curve_fit(Gauss2D, (xx.ravel(), yy.ravel()), subimg.ravel(), p0=initial_guess, maxfev=10000)
    return(param)
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
def gauss1D(x, A, width, mu, B): 
    return (B+(A*np.exp(-(x-mu)**2/width)))
    

# For single star ##
def fwhm(image, x0, y0):
    box = 10
    arr = np.arange(0, 2 * box, 1)
    arr1 = np.arange(0, 2 * box, 0.01)
    datax = image[int(y0), int(x0)-box:int(x0)+box]
    datay = image[int(y0)-box:int(y0)+box, int(x0)]
    data = (datax + datay) / 2
    initial = [np.max(data), 2.5, box, np.min(data)]
    popt, _ = curve_fit(gauss1D, arr, data, p0=initial)
    sigma = np.sqrt(popt[1] / 2)
    F = 2.355 * sigma
    # Optional diagnostic plot
    plt.figure()
    plt.plot(arr, data, 'bo', label=f'FWHM = {F:.2f}')
    plt.plot(arr1, gauss1D(arr1, *popt))
    plt.legend()
    plt.xlabel("Pixel")
    plt.ylabel("Counts")
    plt.title("Gaussian Fit")
    plt.show()
    return F
def save_phot_csv(phot_table, filename):
    t = Table(phot_table)
    t.write(filename, format='csv', overwrite=True)
'''
def fwhm(img):
    arr = np.arange(0, 30, 1)
    arr1 = np.arange(0, 30, 0.01)
    vmin, vmax = np.percentile(img, (1, 99))
    plt.imshow(img,cmap='gray', vmin = vmin, vmax = vmax)
    plt.title('select 5 random stars to compute average FWHM in the field \n press x on the center of bright star ')
    cmd=plt.connect('key_press_event', press)
    plt.show()
    ap, bp = cent(img,X2,Y2,10,10)
    apx, bpy = int(round(ap)), int(round(bp)) 
    datax = img[bpy, apx-15:apx+15]
    datay = img[bpy-15:bpy+15, apx]
    data = (datax+datay)/2
    initial = [max(data), 2.5, 15, min(data)]
    popt, pcov = curve_fit(gauss1D, arr, data, p0=initial)
    F = 2.355*np.sqrt(popt[1]/2)
    #plt.plot(arr, datay,'bo', label=str(F))
    #plt.plot(arr1, gauss1D(arr1, *popt))
    #plt.legend()
    #plt.show()
    return(F)
'''

def measure_fwhm_multi(img, N=5):
    arr = np.arange(0, 20, 1)
    arr1 = np.arange(0, 20, 0.01)
    vmin, vmax = np.percentile(img, (1, 99))
    FS = []
    x, y = [], []
    for i in range(N):
        plt.imshow(img,cmap='gray', vmin = vmin, vmax = vmax)
        plt.title('select 5 random stars to compute average FWHM in the field \n press x on the center of bright star ')
        cmd=plt.connect('key_press_event', press)
        plt.show()
        ap, bp = cent(img,X2,Y2,10,10)
        x.append(ap)
        y.append(bp)
        print(f"Click star {i+1}/{N}")
        FS.append(fwhm(img, x0 = ap, y0 = bp))
    return np.median(FS), x, y


def bkg(Image, x, y, win_in, win_out):
    sub1 = Image[y-win_out:y-win_in, x-win_out:x+win_out]
    sub2 = Image[y+win_in:y+win_out, x-win_out:x+win_out]
    sub3 = Image[y-win_in:y+win_in, x-win_out:x-win_in]
    sub4 = Image[y-win_in:y+win_in, x+win_in:x+win_out]
    s1 = sub1.flatten()
    s2= sub2.flatten()
    s3 = sub3.flatten()
    s4= sub4.flatten()
    S = np.r_[s1, s2, s3, s4]
    bg_counts = np.median(S)
    return(bg_counts)



def single_photometry(img_list, spath, mode):
    """
    perform aperture photometry for single star eo-pair
    Parameters
    ----------
    img_list: list
        Nested 2D list of FITS images sorted according to HWP angle
        and observational set number.
        Example:
        [[p0_s0, p0_s1, ...],
         [p1_s0, p1_s1, ...],
         [p2_s0, p2_s1, ...],
         [p3_s0, p3_s1, ...]]
        use stacked or shifted image lists produced by the 'stacking.py' script.
    spath: str
        Directory to save catalog.
        (for convenience use the same path in which all the processed data is being saved)
    mode: str
        photometric method
        'fixed aperture' or 'pol curve of growth' 
    """
    #------------------------------------------
    # defining o-ray & e-ray pairs
    #------------------------------------------
    ref = fits.getdata(os.path.join(spath, img_list[0][0])) #  first image of p0
    mean, med, std = sigma_clipped_stats(ref)
    plt.imshow(ref,cmap='gray', vmin = med-10*std, vmax = med+10*std)
    plt.colorbar()
    plt.title('press x on the center of oray star '+ str(img_list[0][0]))
    cmd=plt.connect('key_press_event', press)
    plt.show()
    try:
    	ox,oy = cent_2g(ref,X2,Y2,10,10)
    except:
    	ox,oy = cent(ref,X2,Y2,10,10)
    plt.show()
    plt.imshow(ref,cmap='gray', vmin = med-10*std, vmax = med+10*std)
    plt.colorbar()
    plt.title('click on the center of eray star ')
    cmd=plt.connect('button_release_event', single_click)
    plt.show()
    try:
    	ex, ey = cent_2g(ref,X1,Y1,10,10)
    except:
    	ex, ey = cent(ref,X1,Y1,10, 10)
    if len(img_list[0])==1:
        prefix = 'stacked'
    else:
        prefix = 'shifted'
    for i in range(len(img_list[0])): # on number of sets
        ref_set_img = fits.getdata(os.path.join(spath, img_list[0][i])) # p0 of each set
        try:
        	oox, ooy = cent_2g(ref_set_img, ox, oy, 10, 10)
        	eex, eey = cent_2g(ref_set_img, ex, ey, 10, 10)
        except:
        	oox, ooy = cent(ref_set_img, ox, oy, 10, 10)
        	eex, eey = cent(ref_set_img, ex, ey, 10, 10)
        O = np.c_[oox, ooy]
        E= np.c_[eex, eey]
        FWHM= fwhm(ref_set_img, oox, ooy)
        R = 1.5*FWHM
        if mode == 'fixed aperture':
            radii = np.array([R]) 
        elif mode == 'pol curve of growth':
            radii = np.linspace(0.5*FWHM,1.5*FWHM,15)
        else:
            raise ValueError(
            "mode must be either 'pol curve of growth' or 'fixed aperture'")
        np.savetxt(os.path.join(spath, "aperture_radii.txt"),radii,header="Aperture radii in pixels")
        apertureO =  [CircularAperture(O, r=r) for r in radii]
        apertureE =  [CircularAperture(E, r=r) for r in radii] 
        #aperture = CircularAperture(positions, r=R)
        annulus_aperturesO = CircularAnnulus(O, r_in=R+10, r_out=R+15.)
        annulus_aperturesE = CircularAnnulus(E, r_in=R+10, r_out=R+15.)
        for j in range(len(img_list)): # on number of HWP positions
            #print(objs[i][j])
            name = img_list[j][i]
            Image = fits.getdata(os.path.join(spath,img_list[j][i]))
            phot_tableO = aperture_photometry(Image, apertureO)
            phot_tableE = aperture_photometry(Image, apertureE)
            phot_tableO.add_column(np.nan, name = 'RA')
            phot_tableO.add_column(np.nan, name = 'DEC')
            phot_tableE.add_column(np.nan, name = 'RA')
            phot_tableE.add_column(np.nan, name = 'DEC')
            aperstatsO = ApertureStats(Image, annulus_aperturesO) 
            med_bgO = aperstatsO.median
            aperstatsE = ApertureStats(Image, annulus_aperturesE) 
            med_bgE = aperstatsE.median
            for rad in range(len(radii)):
                bkgO_counts = np.multiply(med_bgO,apertureO[rad].area)
                bkgE_counts = np.multiply(med_bgE, apertureE[rad].area)
                phot_tableO.add_column(bkgO_counts, name='bkg_Ocounts_ap'+str(rad))
                phot_tableE.add_column(bkgE_counts, name='bkg_Ecounts_ap'+str(rad))
                phot_finalO = np.subtract(phot_tableO['aperture_sum_'+str(rad)], bkgO_counts)  
                phot_finalE = np.subtract(phot_tableE['aperture_sum_'+str(rad)], bkgE_counts)
                phot_tableO.add_column(phot_finalO, name='Final_phot_ap'+str(rad))
                phot_tableE.add_column(phot_finalE, name='Final_phot_ap'+str(rad)) 
            outO = os.path.join(spath, f'{prefix}_phot_Oray_P{j}_s{i}.csv')
            outE = os.path.join(spath, f'{prefix}_phot_Eray_P{j}_s{i}.csv') 
            save_phot_csv(phot_tableO, outO)
            save_phot_csv(phot_tableE, outE)
    return([phot_tableO, phot_tableE])
           
