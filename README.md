# Intro

Supplementary data and code for "Nonparametric estimation of age-depth models from sedimentological and stratigraphic data".

## Authors

__Niklas Hohmann__  
Utrecht University  
email: n.h.hohmann [at] uu.nl  
Web page: [www.uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)  
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

## Requirements

Base R (version >= 4) and the RStudion IDE.

## Usage

Open the file `nonparam_est_adm_supp.Rproj` in RStudio. This will set all paths correctly. Then you can inspect and run the code in the folders.

## Repository structure

* _code_ : folder for code
  * _analysis_steinbruch_schmidt.R_ : code for Steinbruch Schmidt example
* _data_ : folder for data
  * _raw_ : folder for raw data, for read only
    * _SbS_XRF_forfactor3.csv_ : data for Steinbruch Schmidt example, from da Silva (2024)
* _figs_ : folder for figures
* _.gitignore_ : untracked files
* _LICENSE_ : Apache 2.0 license text
* _nonparam_est_adm_supp.Rproj_ : R project file
* _README_ : Readme file

## References

Data in `data/raw/SbS_XRF_forfactor3.csv` and parts of the code in `code/analysis_steinbruch_schmidt.R` are from

* da Silva, A.-C. (2024). Anchoring the Late Devonian mass extinction in absolute time by integrating climatic controls and radio-isotopic dating: Supplementary code (v1.0.0). Zenodo. [DOI: 10.5281/zenodo.12516430](https://doi.org/10.5281/zenodo.12516430)

## Copyright

Copyright 2023-2024 Netherlands eScience Center and Utrecht University.

## License

Apache 2.0 License, see LICENSE file for license text.

## Funding information

Funded by the European Union (ERC, MindTheGap, StG project no 101041077). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.
![European Union and European Research Council logos](https://erc.europa.eu/sites/default/files/2023-06/LOGO_ERC-FLAG_FP.png)
