from config import *

from AIMPOL_pipeline import reduction


POL = reduction(
    datapath=datapath,
    biaspath=biaspath,

    band=band,
    exposure=exposure,

    hwp_identifier=hwp_identifier,
    HWP_angle_map=HWP_angle_map,
    
    field = field,
    photometry_imglist =photometry_imglist,
    limiting_mag = limiting_mag,

    stacking_method=stacking_method,
    photometry_mode=photometry_mode,

    q_inst=q_inst,
    sigma_q_inst=sigma_q_inst,

    u_inst=u_inst,
    sigma_u_inst=sigma_u_inst,

    psi_inst=psi_inst,
    sigma_psi_inst=sigma_psi_inst,

    aperture_correction=aperture_correction
)

print(POL)