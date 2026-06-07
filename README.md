# DHARA
### Data Handling and Automated Reduction pipeline for AIMPOL.

DHARA is an interactive, automated data-reduction and analysis pipeline developed for polarimetric observations obtained with the AIMPOL instrument. It performs full end-to-end processing, including bias correction, image stacking, source extraction, photometry, and polarization computation.


## User input 

All user-defined parameters (data paths, exposure times, HWP settings, etc.) are provided in:

config/config.py

Change the parameters according to your requirements.

## Pre-requisites

DHARA has been tested with:

- Python 3.10.0
- numpy 1.26.4
- scipy 1.15.3
- matplotlib 3.5.0
- astropy 6.1.7
- photutils 1.5.0
- natsort 8.0.0
- astroquery 0.4.6
- pandas 1.5.3
- astroalign 2.6.2
- os
- fnmatch
- re

## How to run

python run_DHARA.py or pthon3 run_DHARA.py


##  Important Usage Guidelines

To ensure reproducibility and maintain code integrity:

- The pipeline should be executed only through:
  - `python run_DHARA.py`
- All user-specific inputs and parameters must be defined in:
  - `config/config.py`

- Core modules inside `dhara/` should not be modified for standard runs.

## 👤 Author

- **Name:** Namita Uppal 
- **Affiliation:** Institute of Astrophysics, FORTH 
- **Pipeline:** DHARA (Data Handling and Automated Reduction pipeline for AIMPOL)

### Note: This repository contains the pipeline presented in the manuscript:
  DHARA: Data Handling and Automated Reduction pipeline for AIMPOL.

The manuscript has been submitted to the Journal of Astrophysics and Astronomy (JoAA) and is currently under peer review process. 

Please note that the code and documentation may be updated as revisions to the manuscript are completed. If you use this pipeline, please cite the corresponding publication once it becomes available. This file will be updated with the citation information once the DOI is created.





