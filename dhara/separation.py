import numpy as np
import re
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
from scipy.optimize import curve_fit
from scipy import interpolate
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
from photutils.centroids import centroid_com, centroid_quadratic
from photutils.centroids import centroid_1dg, centroid_2dg

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
def gauss1D(x, A, width, mu, B): 
    return(B+(A*np.exp(-(x-mu)**2/width)))
    
def get_fov(astrometric_image):
    hdr = fits.getheader(astrometric_image)
    comments = hdr['COMMENT']
    plate_scale = None
    for c in comments:
        if 'scale:' in c and 'arcsec/pix' in c:
            match = re.search(r'scale:\s*([0-9.]+)', c)
            if match:
                plate_scale = float(match.group(1))
                break
    if plate_scale is None:
        raise ValueError('Plate scale not found in FITS header')
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']
    fov_x = nx * plate_scale / 60
    fov_y = ny * plate_scale / 60
    return fov_x, fov_y, plate_scale





def imageCoord_eo(astrometry_image, limiting_mag, outpath, showplot = True, recenter = True):
    """
    Generate image pixel coordinates for e-ray and o-ray using astrometric image.
    Parameters
    ----------
    astrometric_image : str
        FITS image with valid WCS solution along with the path 
        (can be generated from astrometry.net service on one of the stacked/shifted file).
    limiting_mag : float
        Gaia G-band limiting magnitude.
    outpath : str
        Directory to save catalog.
    showplot : bool
        Display detected Gaia sources on image.
    recenter : bool
        recenter the coordinates for subpixel accuracy.
    """
    f = fits.open(astrometry_image)
    data =  f[0].data
    header = f[0].header
    ydim,xdim = data.shape
    fov, _, _ = get_fov(astrometry_image)
    vizier_cat  = 'I/355/gaiadr3'
    mean, median, sigma = sigma_clipped_stats(data)
    wcs = WCS(header,naxis=2)
    alpha, delta = wcs.all_pix2world(xdim/2, ydim/2, 1)#1300/2, 1340/2, 1)
    Vizier.ROW_LIMIT = -1
    result = Vizier.query_region(SkyCoord(ra = alpha, dec = delta, unit = (u.deg, u.deg)),
                                 width=fov*u.arcmin,
                                 catalog=vizier_cat, column_filters={"Gmag": f"<{limiting_mag}"})  # czernik 3 , 3 arcmin radius (empol field is 3 arcmin in diameter) and II/349 is the panstarrs ps1 catalog
    if len(result) == 0:
        print("No Gaia sources found.")
        return None
    RA = result[0]['RA_ICRS']
    DEC =  result[0]['DE_ICRS']
    Gmag = result[0]['Gmag']
    xo, yo = wcs.all_world2pix(RA, DEC, 1) 
    # the origin is (1,1) but python needs (0,0) so, in order to use this solution in ppython subtract x and y by 1
    Xo = np.array(xo)
    Yo =np.array(yo)
    plt.imshow(data, vmin=median-5*sigma, vmax=median+5*sigma, cmap='gray')
    plt.plot(Xo-1,Yo-1, 'bx', label='Gmag < '+str(limiting_mag), markersize = 5)
    plt.title('astrometry solution corresponds to o-ray or e-ray ?')
    plt.legend()
    plt.show()
    astro_sol = input("the astrometric solution is for which component (o or e)? ")   
    np.savetxt(os.path.join(outpath,'Gaia_Gmag'+str(limiting_mag)+'.txt'), np.c_[Xo,Yo, RA, DEC, Gmag ], header='xeray      yeray     RA     DEC')
    #----------------------------------------------------------
    # Finding the image coordinates corresponding to both e-ray and o-ray images
    #-----------------------------------------------------------
    print('SELECT an isolated pair of e-ray and o-ray in the upcoming figures')
    plt.imshow(data, vmin=median-5*sigma, vmax=median+5*sigma, cmap='gray')
    plt.title('click on o-ray star ')
    cmd=plt.connect('button_release_event', single_click)
    plt.show()
    oox, ooy = cent(data,X1,Y1,10,10 )
    plt.imshow(data, vmin=median-5*sigma, vmax=median+5*sigma, cmap='gray')
    plt.colorbar()
    plt.title('press x on e-ray star')
    cmd=plt.connect('key_press_event', press)
    plt.show()
    eex, eey = cent(data, X2, Y2, 10, 10)
    d = np.sqrt((oox-eex)**2+(ooy-eey)**2)
    #m = (oy-ey)/(ox-ex) #
    Th=np.arctan2(ooy-eey, oox-eex)
    #A = np.sqrt(1+m**2)
    if astro_sol == 'o':
    	x_oray, y_oray = Xo-1, Yo-1
    	x_eray = x_oray - d*np.cos(Th)
    	y_eray = y_oray - d*np.sin(Th)
    	#x_eray = x_oray - (d/A)
    	#y_eray = y_oray - (d*m/A)
    else:
        x_eray, y_eray = Xo-1, Yo-1
        #x_oray = x_eray + (d/A)
        #y_oray = y_eray + (d*m/A)
        x_oray = x_eray + d*np.cos(Th)
        y_oray = y_eray + d*np.sin(Th)
    #indx_inside = np.where((x_eray > 0) & (y_eray > 0) & (x_oray < 500) & (y_oray < 500))[0]
    indx_inside = np.where((x_eray > d) & (y_eray > d) & (x_oray < xdim-d) & (y_oray < ydim-d))[0]
    if showplot == True:
        plt.imshow(data, vmin=median-5*sigma, vmax=median+5*sigma, cmap='gray')
        plt.colorbar()
        plt.plot(x_eray[indx_inside], y_eray[indx_inside], 'b.', label = 'e-ray')
        plt.plot(x_oray[indx_inside], y_oray[indx_inside], 'r.', label = 'o-ray')
        plt.legend()
        plt.show()
    catalog = np.c_[x_oray[indx_inside]+1, y_oray[indx_inside]+1, x_eray[indx_inside]+1, y_eray[indx_inside]+1, 
    RA[indx_inside], DEC[indx_inside], Gmag[indx_inside]]
    if recenter == True:
        k = h= 15
        arr = np.arange(0, 30, 1) # 30 = 2*h)
        arr1 = np.arange(0, 30, 0.01)
        O = np.c_[x_oray[indx_inside], y_oray[indx_inside]]
        E = np.c_[x_eray[indx_inside], y_eray[indx_inside]]
        ap, bp = cent(data,X2,Y2,10,10)
        apx, bpy = int(round(ap)), int(round(bp))
        datax = data[bpy, apx-k:apx+k]
        datay = data[bpy-k:bpy+k, apx]  
        data_av = (datax+datay)/2
        initial = [max(data_av), 2.5, 15, min(data_av)]
        popt, pcov = curve_fit(gauss1D, arr, data_av, p0=initial)
        FWHM = 2.355*np.sqrt(popt[1]/2) # radius
        print('FWHM = ', FWHM)
        plt.plot(arr, data_av,'bo', label=str(FWHM))
        plt.plot(arr1, gauss1D(arr1, *popt))
        plt.legend()
        plt.show()
        R = int(FWHM)
        reOx = []
        reOy = []
        reEx = []
        reEy = []
        for i in range(len(O[:,0])):
            ox = O[i,0]
            oy = O[i,1]
            ex = E[i,0]
            ey = E[i,1]
            #mox, moy = COM_iter(img, ox, oy, 7, 8) # maximum pixel
            Cox, Coy = cent(data, ox, oy, 2,2) # centroid_quad
            if(np.isnan(Cox)==True):
                Cox, Coy = cent_2g(data, ox, oy,2,2)
            if((np.abs(Cox-ox) > 2) | (np.abs(Coy-oy) > 2)):
                Cox = ox
                Coy = oy        
            Cex, Cey = cent(data, ex, ey, 2,2)
            if(np.isnan(Cex)==True):
                Cex, Cey = cent_2g(data, ex, ey,2,2)
            if((np.abs(Cex-ex) > 2) | (np.abs(Cey-ey) > 2)):
                Cex = ex
                Cey = ey
            reOx.append(Cox+1)
            reOy.append(Coy+1)
            reEx.append(Cex+1)
            reEy.append(Cey+1)
           
        reO = np.c_[reOx, reOy]   # according to DS9 
        reE = np.c_[reEx, reEy] 
        plt.imshow(data, vmin=median-5*sigma, vmax=median+5*sigma, cmap='gray')
        plt.plot(O[:,0], O[:,1], 'r.', label = 'o-ray')
        plt.plot(E[:,0], E[:,1], 'g.', label = 'e-ray')
        plt.plot(reO[:,0]-1, reO[:,1]-1, 'kx' , label = 'recenterd o-ray')#, label='max pixel')
        plt.plot(reE[:,0]-1, reE[:,1]-1, 'bx', label = 'recenterd e-ray')
        plt.show()
        catalog = np.c_[reO, reE, RA[indx_inside], DEC[indx_inside], Gmag[indx_inside] ]
    F = os.path.join(outpath, 'image_coord.txt')
    np.savetxt(F, catalog, header = ' X_oray     Y_oray    X_eray     Y_eray     RA     DEC     Gmag')
    return(F)
