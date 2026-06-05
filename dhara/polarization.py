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
from .stokesconversion import StokesConversion
from astropy.modeling import models, fitting
from astropy.nddata import Cutout2D
import pandas as pd


def select_plateau_aperture(apertures, polarization,
                            deriv_threshold=0.005,
                            min_plateau_length=3):
    """
    Automatically select the aperture where polarization reaches a plateau.

    Parameters
    ----------
    apertures : array-like
        Aperture radii (e.g., 1–5 × FWHM).
    polarization : array-like
        Degree of polarization at each aperture.
    deriv_threshold : float
        Maximum derivative (ΔP) threshold to consider the curve 'flat'.
    min_plateau_length : int
        Minimum number of consecutive points meeting the flat condition.

    Returns
    -------
    float
        Selected aperture (first in plateau), or the largest aperture if none found.
    """
    apertures = np.array(apertures)
    polarization = np.array(polarization)
    # First derivative
    dP = np.abs(np.diff(polarization))
    # Identify indices where slope is small
    flat = dP < deriv_threshold
    # Search for consecutive flat region
    count = 0
    for i, is_flat in enumerate(flat):
        if is_flat:
            count += 1
            if count >= min_plateau_length:
                # Plateau starts min_plateau_length points earlier
                star_idx = i - min_plateau_length + 1
                if star_idx ==0:
                    continue
                return apertures[star_idx]
        else:
            count = 0
    # If no plateau found, return the largest aperture
    return apertures[-1]

    

def modulation(alpha, q,u): # alpha is in degrees
    Y = q*np.cos(4*alpha*np.pi/180) + u*np.sin(4*alpha*np.pi/180)
    return(Y)

def overlapping_stars(Xo, Yo, Xe, Ye, aperture_rad):
    Indx = []
    for i in range(len(Xo)):
        for j in range(len(Xe)):
            if(i!=j):
                doe = np.sqrt(( Xo[i]-Xe[j])**2 + (Yo[i]-Ye[j])**2 )  # diff b/w oray of one star and e-ray of other
                if(doe < aperture_rad[i] + aperture_rad[j]):
                    #plt.scatter(oray[i,0],oray[i,1], color= 'green', marker='o')
                    #plt.scatter(eray[j,0], eray[j,1], color= 'green', marker='o')
                    Indx.append(i)
                    Indx.append(j)
                doo = np.sqrt((Xo[i]-Xo[j])**2 + (Yo[i]-Yo[j])**2 ) # diff b/w oray of one star and o-ray of other
                dee = np.sqrt((Xe[i]-Xe[j])**2 + (Ye[i]-Ye[j])**2 )# diff b/w eray of one star and e-ray of other
                if(doo < aperture_rad[i] + aperture_rad[j]):
                    #plt.scatter(oray[i,0], oray[i,1], color= 'magenta', marker='o')
                    #plt.scatter(oray[j,0], oray[j,1], color= 'magenta', marker='o')
                    Indx.append(i)
                    Indx.append(j)
                if(dee < aperture_rad[i] + aperture_rad[j]):
                    #plt.scatter(eray[i,0], eray[i,1], color= 'cyan', marker='o')
                    #plt.scatter(eray[j,0], eray[j,1], color= 'cyan', marker='o')
                    Indx.append(i)
                    Indx.append(j)
    fix = np.unique(Indx)
    return(fix)  


def pol_calc(photO, photE, spath, q_inst, sigma_q_inst, u_inst, sigma_u_inst, psi_inst, sigma_psi_inst, ang_map = [67.5, 45.0, 22.5, 0.0], mode='pol curve of growth', aperture_correction = True):
    """
    Compute instrumental corrected linear polarization parameters.
    Parameters
    ----------
    photO : list
        List of filenames containing photometric measurements of the o-ray, 
        ordered by HWP angle, obtained from the photometry task. 
    photE : list
        List of filenames containing photometric measurements of the e-ray, 
        ordered by HWP angle, obtained from the photometry task. 
    spath : str
        Directory where the polarization results catalog will be saved.
    q_inst : float
        Instrumental Stokes q parameter derived from unpolarized standard
        star observations.
        (0 for unpolarized standard star reduction).
    sigma_q_inst : float
        uncertainty in the instrumental Stokes q parameter derived from unpolarized
        standard star observations.
        (0 for unpolarized standard star reduction).
    u_inst : float
        Instrumental Stokes u parameter derived from unpolarized standard
        star observations.
        (0 for unpolarized standard star reduction).
    sigma_u_inst : float
        uncertainty in the instrumental Stokes u parameter derived from unpolarized
        standard star observations.
        (0 for unpolarized standard star reduction).
    psi_inst : float
        Instrumental zero-offset polarization angle (degrees),
        derived from polarized standard star observations by comparing
        measured and literature values.
        (0 for unpolarized and polarized standard star reduction).
    sigma_psi_inst : float
        Uncertainty in the instrumental zero-offset polarization angle correction
        (in degrees), estimated from the dispersion between measured and literature values.
        (0 for unpolarized and polarized standard star reduction).
    ang_map : list of float, optional
        Polarization angles (in degrees) corresponding to the HWP 
        positions (defined in earlier tasks). 
        Default is [67.5, 45.0, 22.5, 0.0].
    mode : str, optional
        Photometric reduction mode used for photometry task.
        Options include:
        - 'fixed aperture'
        - 'pol curve of growth'
    aperture_correction : bool, optional
        If True, applies aperture correction before computing polarization
        parameters. Default is True.
    """
    A = np.arange(0,90, 0.1)
    al = np.array(ang_map)
    ref = pd.read_csv(os.path.join(spath, photO[0]))
    n = len(ref['RA'])
    radii = np.loadtxt(os.path.join(spath, "aperture_radii.txt"), comments = '#')
    psi_inst = np.deg2rad(psi_inst) # instrumental shift in angle
    sigma_psi_inst = np.deg2rad(sigma_psi_inst)
    Final_cat = []
    pol_1FWHM = []
    for i in range(n):
        phot_radii = []
        for k in range(len(radii)):
            ratio = []
            e_ratio=[]
            for j in range(4):
                phO = pd.read_csv(os.path.join(spath, photO[j]))
                phE = pd.read_csv(os.path.join(spath, photE[j]))
                r = phE['Final_phot_ap'+str(k)].iloc[i]/phO['Final_phot_ap'+str(k)].iloc[i]
                # sigma_R = sqrt(I0 +Ie +Ibo +Ibe)/(I0 + Ib)
                a = phO['aperture_sum_'+str(k)].iloc[i]+\
                phE['aperture_sum_'+str(k)].iloc[i]+\
                phO['bkg_Ocounts_ap'+str(k)].iloc[i]+\
                phE['bkg_Ecounts_ap'+str(k)].iloc[i]
                b = phO['aperture_sum_'+str(k)].iloc[i]+phE['aperture_sum_'+str(k)].iloc[i]
                e = np.sqrt(a)/b
                e_ratio.append(e)
                ratio.append(r)
            ratio = np.array(ratio)
            e_ratio = np.array(e_ratio)
            f = 1/(np.prod(np.abs(ratio)))**(1/4)
            mod = ((ratio*f) - 1)/((ratio*f)+1) 
            popt, pcov = curve_fit(modulation, al, mod,sigma=e_ratio, bounds=((-np.inf, -np.inf), (np.inf, np.inf)),  absolute_sigma = False)
            #--------------------------------------
            # instrumental polarization correction
            #------------------------------------
            q_out = popt[0] - q_inst
            u_out = popt[1]- u_inst
            e_Q = np.sqrt((pcov[0,0]**0.5)**2 + sigma_q_inst**2)
            e_U = np.sqrt((pcov[1,1]**0.5)**2 + sigma_u_inst**2)
            #------------------------------------
            #zero-offset in angle correction
            #------------------------------------
            C_full = np.zeros((3,3))
            C_full[0,0] = e_Q**2   # q,u covarianc
            C_full[1,1] = e_U**2
            C_full[2,2] = sigma_psi_inst**2 # offset variance
            c = np.cos(2*psi_inst)
            s = np.sin(2*psi_inst)
            J = np.array([
            [c, -s, -2*q_out*s - 2*u_out*c],
            [s,  c,  2*q_out*c - 2*u_out*s]])
            C_corr = J @ C_full @ J.T
            e_qt = np.sqrt(C_corr[0,0])
            e_ut = np.sqrt(C_corr[1,1])
            cov_qu_corr = C_corr[0,1]
            qt = q_out*c - u_out*s
            ut = q_out*s + u_out*c
            #----------------------------------
            # Polarization calculation 
            #----------------------------------
            converter = StokesConversion()
            Pol, Pol_mas, ePol, PA, ePA = converter.convert(qt, ut, stokes_q_err=e_qt, stokes_u_err=e_ut, unit='deg' )  
            phot_radii.append([radii[k],qt, e_qt, ut, e_ut, Pol*100, ePol*100, PA, ePA ])
        phot_radii = np.array(phot_radii)
        plt.plot(phot_radii[:,0], phot_radii[:,5], 'r.-')
        plt.show()
        if(mode == 'pol curve of growth'):
            ap = select_plateau_aperture(apertures = phot_radii[:,0],  polarization= phot_radii[:,5], deriv_threshold=0.005, min_plateau_length=5)
            indx = np.where(radii==ap)[0][0]
        elif(mode == 'fixed aperture'):
            indx = 0
        flag = indx
        #CO = phO['Final_phot_ap'+str(indx)].iloc[i] # bkg subtracted o-ray counts of the selected aperture
        #CE = phE['Final_phot_ap'+str(indx)].iloc[i] # bkg subtracted e-ray counts of the selected aperture
        #SNR_O = CO/np.sqrt(phO['bkg_Ocounts_ap'+str(indx)].iloc[i])  # o-ray  bkg sub counts above x bkg sigma level
        #SNR_E = CE/np.sqrt(phO['bkg_Ecounts_ap'+str(indx)].iloc[i])
        row = [ref['id'].iloc[i], ref['RA'].iloc[i], ref['DEC'].iloc[i], *phot_radii[indx,:]]#, CO, CE, SNR_O, SNR_E]
        Final_cat.append(row)
        pol_1FWHM.append(phot_radii[0,:].copy())
    columns = ['id','RA','DEC','aperture', 'q','eq','u','eu','pol','epol','PA','ePA']
    df = pd.DataFrame(Final_cat, columns=columns)
    if aperture_correction:
        Xo, Yo = phO['xcenter'], phO['ycenter']
        Xe, Ye = phE['xcenter'], phE['ycenter']
        aper_rad = df['aperture'].values
        indx_overlap = overlapping_stars(Xo, Yo, Xe, Ye, aper_rad)
        df['overlapping_flag'] = 0 # no overlap
        #--------------------------------------------------------------------------
        # replace flagged results with the results corresponds to 1*FWHM
        #--------------------------------------------------------------------------
        pol_1FWHM = np.array(pol_1FWHM)
        for idx in indx_overlap:
            # replace by 1FWHM measurements
            df.loc[idx, 'aperture'] = pol_1FWHM[idx,0]
            df.loc[idx, 'q']    = pol_1FWHM[idx,1]
            df.loc[idx, 'eq']   = pol_1FWHM[idx,2]
            df.loc[idx, 'u']    = pol_1FWHM[idx,3]
            df.loc[idx, 'eu']   = pol_1FWHM[idx,4]
            df.loc[idx, 'pol']  = pol_1FWHM[idx,5]
            df.loc[idx, 'epol'] = pol_1FWHM[idx,6]
            df.loc[idx, 'PA']   = pol_1FWHM[idx,7]
            df.loc[idx, 'ePA']  = pol_1FWHM[idx,8]
            # flag = 1 for the overlapping stars
            df.loc[idx, 'overlapping_flag'] = 1
    df.to_csv(os.path.join(spath, 'pol_cat.csv'),index=False)
    return(df)
    
        

