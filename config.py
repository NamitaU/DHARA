###############################################
#-------------------------------------------
#user input to be edited here
#-------------------------------------------
###############################################

# =========================================================
# POLARIMETRY REDUCTION CONFIGURATION
# =========================================================


# -------------------------
# PATHS
# -------------------------

# Base directory of observation run
path = '/Volumes/namita1TB/PDF/ISM/open_clusters/Observations/2025NovAIMPOL/15Nov2025'

# Raw data directory (target field)
datapath = path + '/NGC2099'

# Bias frames directory
biaspath = path + '/Bias'


# -------------------------
# OBSERVATION SETUP
# -------------------------

# Filter band (e.g., U, B, V, R, I)
band = 'R'

# Exposure time in seconds
exposure = 15


# -------------------------
# HWP CONFIGURATION
# -------------------------

# HWP identifiers in filenames
hwp_identifier = ['p0', 'p1', 'p2', 'p3']

# Corresponding HWP angles (degrees)
HWP_angle_map = [67.5, 45, 22.5, 0]


# -------------------------
# SOURCE DETECTION
# -------------------------
# Select if you want to do single star reduction or a crowded field using parameter 'field' = 'crowded' or 'single star'
field = 'single star'
# Gaia limiting magnitude for detection
limiting_mag = 15.5


# -------------------------
# REDUCTION OPTIONS
# -------------------------
# imglist where you want to perform reduction ('stacked' or 'shifted')
photometry_imglist = 'stacked'
# Image stacking method: 'sum' or 'median'
stacking_method = 'sum'

# Photometry method:
#   - 'fixed aperture'
#   - 'pol curve of growth'
photometry_mode = 'pol curve of growth'

# Correct overlapping aperture contamination
aperture_correction = False


# -------------------------
# INSTRUMENTAL POLARIZATION CALIBRATION
# -------------------------

# Instrumental Stokes parameter (from unpolarized standard)
q_inst = 0.0
u_inst = 0.0
sigma_q_inst = 0.0
sigma_u_inst = 0.0

# Instrumental angle offset (degrees)
psi_inst = 0.0
sigma_psi_inst = 0.0

