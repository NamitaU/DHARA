import os
import fnmatch
import natsort
import pandas as pd
from .Biascombine import biasCombine
from .Stacking import classify_hwp_images as classHWP
from .Stacking import ShiftnStack
from .separation import imageCoord_eo
from .photometry import photometry
from .polarization import pol_calc
from .single_star_reduction import single_photometry

def reduction(
    datapath,
    biaspath,
    band,
    exposure,
    hwp_identifier=['p0', 'p1', 'p2', 'p3'],
    HWP_angle_map=[67.5, 45, 22.5, 0],
    field = 'crowded',
    limiting_mag=15.5,
    photometry_imglist = 'stacked',
    stacking_method='sum',
    photometry_mode='fixed aperture',
    q_inst=0,
    sigma_q_inst=0,
    u_inst=0,
    sigma_u_inst=0,
    psi_inst=0,
    sigma_psi_inst=0,
    aperture_correction=True):
    # -------------------------------------------------
    # output directory
    # -------------------------------------------------
    spath = os.path.join(datapath,f'reduced{band}_{exposure}s')
    os.makedirs(spath, exist_ok=True)

    print('\n===================================')
    print(' STARTING POLARIZATION PIPELINE ')
    print('===================================\n')

    # -------------------------------------------------
    # master bias
    # -------------------------------------------------

    print('Creating master bias')

    masterbias = biasCombine(biaspath, file_prefix='bias')

    # -------------------------------------------------
    # classify HWP images
    # -------------------------------------------------

    print('Classifying HWP images')

    objects = classHWP(datapath, 
    	band, 
    	exposure, 
    	hwp=hwp_identifier, 
    	ang_map=HWP_angle_map)

    # -------------------------------------------------
    # shifting and stacking
    # -------------------------------------------------

    print('Shifting and stacking images')

    stacked, shifted = ShiftnStack(datapath, 
    	biaspath, 
    	band, 
    	exposure, 
    	objects=objects, 
    	method=stacking_method, 
    	ang_map=HWP_angle_map, 
    	spath=spath)
    	
    	
    if(photometry_imglist == 'stacked'):
        phot_imglist = stacked
    elif(photometry_imglist == 'shifted'):
        phot_imglist = shifted
    else:
        raise ValueError("photometry_imglist must be either 'stacked' or 'shifted'")

    
    if field == 'crowded':
        # -------------------------------------------------
        # source detection / eo matching
        # -------------------------------------------------
        print('Finding e & o source coordinates')
        #first upload the image having astrometry solution
        print(f'\n run Astromeyrty.net webservice on one of the stacked images present in your directory {spath} \n')
        astromerty_image = input("filename of astrometry image : ")
        astromerty_image = os.path.join(spath, astromerty_image)
        img_cat = imageCoord_eo(
        astromerty_image,
        limiting_mag=limiting_mag,
        outpath=spath,
        showplot=True,
        recenter=True)
        # -------------------------------------------------
        # photometry
        # -------------------------------------------------
        # select the set of images (stacked or shifted)
        
        print('Performing aperture photometry')
        phot = photometry(
        img_list=phot_imglist,
        img_coord_file=img_cat,
        spath=spath,
        mode=photometry_mode
        )
    elif field == 'single star':
        #-------------------------------------------------
        # directly do the photometry 
        #-------------------------------------------------
        print('Performing aperture photometry')
        phot = single_photometry(img_list=phot_imglist,  spath=spath, mode = photometry_mode)
    else:
        raise ValueError(
            "field must be either 'crowded' or 'single star'")

    files = natsort.natsorted(os.listdir(spath))

    photO = [
        f for f in files
        if fnmatch.fnmatch(f, photometry_imglist+'_phot_Oray*')
    ]

    photE = [
        f for f in files
        if fnmatch.fnmatch(f, photometry_imglist+'_phot_Eray*')
    ]

    # -------------------------------------------------
    # polarization calculation
    # -------------------------------------------------

    print('Calculating polarization')

    pol_cat = pol_calc(
        photO,
        photE,
        spath,
        q_inst,
        sigma_q_inst,
        u_inst,
        sigma_u_inst,
        psi_inst,
        sigma_psi_inst,
        HWP_angle_map,
        aperture_correction=aperture_correction
    )




    print('\n===================================')
    print(' PIPELINE COMPLETE ')
    print('===================================')

    print(f'\nFiles and polarization catalog are saved in:\n{spath}\n')

    return pol_cat