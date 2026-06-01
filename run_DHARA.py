import config.config as cfg

from dhara.pipeline import reduction



POL = reduction(
    datapath = cfg.datapath,
    biaspath = cfg.biaspath,

    band = cfg.band,
    exposure = cfg.exposure,

    hwp_identifier = cfg.hwp_identifier,
    HWP_angle_map = cfg.HWP_angle_map,
    
    field = cfg.field,
    photometry_imglist = cfg.photometry_imglist,
    limiting_mag = cfg.limiting_mag,

    stacking_method = cfg.stacking_method,
    photometry_mode = cfg.photometry_mode,

    q_inst = cfg.q_inst,
    sigma_q_inst = cfg.sigma_q_inst,

    u_inst = cfg.u_inst,
    sigma_u_inst = cfg.sigma_u_inst,

    psi_inst = cfg.psi_inst,
    sigma_psi_inst = cfg.sigma_psi_inst,

    aperture_correction = cfg.aperture_correction
)

print(POL)