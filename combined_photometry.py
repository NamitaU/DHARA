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
    print(FS)
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




################ combined single and crowded field ###################################################################################
def photometry(img_list, img_coord_file, spath, mode, field): # mode: 'fixed aperture' or 'pol curve of growth'
    """
    perform aperture photometry on e- and o-component images
    Parameters
    ----------
    img_list: str
        2D List of fits image sorted according to HWP and number of sets within each HWP. 
        e.g., [[p0_s0, p0_s1,..], [p1_s0, p1_s1,...], [p2_s0, p2_s1,..], [p3_s0, p3_s1,...]]
        'stacked' or 'shifted',  use the output of Stacking.py
    img_coord_file: str
        filename (along with path) of file containing image coordinates of o-ray and e-ray pairs.
        use the output from separation.py 
    spath: str
        Directory to save catalog.
        (for convenience use the same path in which all the processed data is being saved)
    mode: str
        photometric method
        'fixed aperture' or 'pol curve of growth' 
    """   
    ref_img = fits.getdata(os.path.join(spath, img_list[0][0]))
    mean, med, std = sigma_clipped_stats(ref_img)
    vmin, vmax = np.percentile(ref_img, (1, 99))     
    #----------------------------
    # Source Coordinates 
    #----------------------------
    if field == 'crowded':
        source = np.loadtxt(img_coord_file, comments = '#')
        sky = np.c_[source[:,4], source[:,5]]
        O = np.c_[source[:,0]-1, source[:,1]-1]
        E = np.c_[source[:,2]-1, source[:,3]-1]
        FWHM, x, y = measure_fwhm_multi(ref_img, N=5)
    elif field == 'single star':
        plt.imshow(ref_img,cmap='gray', vmin = med-10*std, vmax = med+10*std)
        plt.colorbar()
        plt.title('press x on the center of oray star ')
        cmd=plt.connect('key_press_event', press)
        plt.show()
        #ox, oy = COM_iter(Img,X2,Y2,30,31)
        try:
            ox,oy = cent_2g(ref_img,X2,Y2,10,10)
        except:
            ox,oy = cent(ref_img,X2,Y2,10,10)
        plt.imshow(ref_img,cmap='gray', vmin = med-10*std, vmax = med+10*std)
        plt.colorbar()
        plt.title('click on the center of eray star ')
        cmd=plt.connect('button_release_event', single_click)
        plt.show()
        try:
            ex, ey = cent_2g(ref_img,X1,Y1,10,10)
        except:
            ex, ey = cent(ref_img,X1,Y1,10, 10)
    else:
        raise ValueError( "field must be either 'single star' or 'crowded'")
    
    for i in range(len(img_list[0])): # on number of sets
        ref_set_img = fits.getdata(os.path.join(spath, img_list[0][i])) # p0 of each set
        if field == 'single star':
            try:
                oox, ooy = cent_2g(ref_set_img, ox, oy, 10, 10)
                eex, eey = cent_2g(ref_set_img, ex, ey, 10, 10)
            except:
                oox, ooy = cent(ref_set_img, ox, oy, 10, 10)
                eex, eey = cent(ref_set_img, ex, ey, 10, 10)
            O = np.c_[oox, ooy]
            E = np.c_[eex, eey]
            FWHM= fwhm(ref_set_img, oox, ooy)
        elif field == 'crowded':
            FS = []
            for cords in range(len(x)):
                FS.append(fwhm(ref_set_img,x[cords], y[cords]))
            FWHM = np.median(FS)
        else:
            raise ValueError( "field must be either 'single star' or 'crowded'")
        R = 1.5*FWHM
        #print('3FWHM = ', 1.5*FWHM, ' 2FWHM = ', FWHM, ' 1FWHM = ', 0.5*FWHM)
        #radii = [1.5*FWHM, 1*FWHM, 0.5*FWHM]
        if mode=='pol curve of growth':
            radii = np.linspace(0.5*FWHM,1.5*FWHM,15)
        elif mode == 'fixed aperture': # 3*FWHM
            radii = np.array([R]) 
        else:
            raise ValueError(
            "mode must be either 'pol curve of growth' or 'fixed aperture'")
        np.savetxt(os.path.join(spath, "aperture_radii.txt"),radii,header="Aperture radii in pixels")
        #-------------------
        # Define Apertures 
        #------------------
        apertureO =  [CircularAperture(O, r=r) for r in radii]
        apertureE =  [CircularAperture(E, r=r) for r in radii]    
        w_in = int(R)+3 # rectangular aperture inner part
        w_out = int(R)+5 # rectangular aperture outer part
        #apO = RectangularAnnulus(O, w_in= 2*w_in, w_out=2*w_out, h_out=2*w_out, h_in = 2*w_in, theta=0.0) # rectangular annulus for visualization only
        #apE = RectangularAnnulus(E, w_in= 2*w_in, w_out=2*w_out, h_out=2*w_out, h_in = 2*w_in, theta=0.0)
        annulus_aperturesO = CircularAnnulus(O, r_in=R+10, r_out=R+15.)
        annulus_aperturesE= CircularAnnulus(E, r_in=R+10, r_out=R+15.)
        #plt.imshow(img,cmap='gray', vmin = vmiin, vmax = vmax)
        #for m in range(len(radii)):
        #    apertureO[m].plot(color='red')
        #    apertureE[m].plot(color= 'blue')
        #apO.plot(color='cyan')
        #apE.plot(color='brown')
        #plt.show()
        if len(img_list[0])==1:
            prefix = 'stacked'
        else:
            prefix = 'shifted' 
        for j in range(len(img_list)):
            fname = img_list[j][i]
            Img = fits.getdata(os.path.join(spath,fname))
            phot_tableO = aperture_photometry(Img, apertureO)  
            phot_tableE = aperture_photometry(Img, apertureE)
            if field == 'single star':
                aperstatsO = ApertureStats(Img, annulus_aperturesO) 
                med_bgO = aperstatsO.median
                aperstatsE = ApertureStats(Img, annulus_aperturesE) 
                med_bgE = aperstatsE.median
                phot_tableO.add_column(np.nan, name = 'RA')
                phot_tableO.add_column(np.nan, name = 'DEC')
                phot_tableE.add_column(np.nan, name = 'RA')
                phot_tableE.add_column(np.nan, name = 'DEC')
            elif field == 'crowded':
                phot_tableO.add_column(sky[:,0], name = 'RA')
                phot_tableO.add_column(sky[:,1], name = 'DEC')
                phot_tableE.add_column(sky[:,0], name = 'RA')
                phot_tableE.add_column(sky[:,1], name = 'DEC')
                med_bgO = np.array([bkg(Img, int(np.round(x)), int(np.round(y)),w_in, w_out) for x, y in O])
                med_bgE = np.array([bkg(Img, int(np.round(x)), int(np.round(y)),w_in, w_out) for x, y in E])
            else:
                raise ValueError("field must be either 'single star' or 'crowded'")
            for rad in range(len(radii)):
                bkgO_counts = np.multiply(med_bgO,apertureO[rad].area)
                bkgE_counts = np.multiply(med_bgE, apertureE[rad].area)
                phot_tableO.add_column(bkgO_counts, name='bkg_Ocounts_ap'+str(rad))
                phot_tableE.add_column(bkgE_counts, name='bkg_Ecounts_ap'+str(rad))
                phot_finalO = np.subtract(phot_tableO['aperture_sum_'+str(rad)], bkgO_counts)  
                phot_finalE = np.subtract(phot_tableE['aperture_sum_'+str(rad)], bkgE_counts) 
                phot_tableO.add_column(phot_finalO, name='Final_phot_ap'+str(rad))
                phot_tableE.add_column(phot_finalE, name='Final_phot_ap'+str(rad)) 
            outO = os.path.join(spath, f'{field}_{prefix}_phot_Oray_P{j}_s{i}.csv')
            outE = os.path.join(spath, f'{prefix}_phot_Eray_P{j}_s{i}.csv')
            save_phot_csv(phot_tableO, outO)
            save_phot_csv(phot_tableE, outE)
            #np.savetxt(os.path.join(spath,'local_bkg_multiple_phot_oray_'+band+'_P'+str(i)+'_all.txt'), phot_tableO)#, header='id     Xcenter_pix     Ycenter_Pix     aperure_sum_1     RA     DEC     bkg_counts_ap1     Final_phot_ap1')
            #np.savetxt(os.path.join(spath,'local_bkg_multiple_phot_eray_'+band+'_P'+str(i)+'_all.txt'), phot_tableE)#, header='id     Xcenter_pix     Ycenter_Pix      aperture_sum_1     RA     DEC     bkg_Ecounts_ap1     Final_phot_ap1')
    return([phot_tableO, phot_tableE])
    
path = '/Volumes/namita1TB/PDF/ISM/open_clusters/Observations/2025NovAIMPOL/15Nov2025'
datapath = path +'/NGC2099'
biaspath = '/Volumes/namita1TB/PDF/ISM/open_clusters/Observations/2025NovAIMPOL/15Nov2025/Bias'
spath = datapath+'/reducedR_15s'
img_list = os.listdir(spath)
img_list = natsort.natsorted(img_list)
img_coord_file = spath + '/image_coord.txt'
Stacked = []
for files in img_list:
    if fnmatch.fnmatch(files, 'stacked*.fits'):
        Stacked.append([files])   
phot = photometry(Stacked, img_coord_file, spath, mode = 'pol curve of growth', field = 'crowded')
'''