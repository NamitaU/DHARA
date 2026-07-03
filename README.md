# DHARA
### Data Handling and Automated Reduction pipeline for AIMPOL.

DHARA is an interactive, automated data-reduction and analysis pipeline developed for polarimetric observations obtained with the AIMPOL instrument. It performs full end-to-end processing, including bias correction, image stacking, source extraction, photometry, and polarization computation.

## [![DOI](https://zenodo.org/badge/1255880455.svg)](https://doi.org/10.5281/zenodo.21158246)

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

## Citation 
This repository contains the pipeline presented in the manuscript:
### DHARA: Data Handling and Automated Reduction pipeline for AIMPOL.

The manuscript has been accepted for publication in the Journal of Astrophysics and Astronomy (JoAA) and is currently in the production/publication process.

The code in this repository corresponds to the version described in the accepted manuscript. Updates, bug fixes, and new features may be incorporated in future releases.

If you use this pipeline in your research, please cite DHARA Version 1.0.0, Zenodo. DOI: 10.5281/zenodo.21158246

and the associated publication. Full bibliographic information and the article DOI will be added here once they become available.





