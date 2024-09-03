# Intro

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13639816.svg)](https://doi.org/10.5281/zenodo.13639816)

Supplementary data and code for "Nonparametric estimation of age-depth models from sedimentological and stratigraphic information".

## Authors

__Niklas Hohmann__  
Utrecht University  
email: n.h.hohmann [at] uu.nl  
Web page: [www.uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)  
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

## Requirements

Base R (version >= 4) and the RStudio IDE.

## Usage

Open the file `nonparam_est_adm_supp.Rproj` in RStudio. This will set all paths correctly, and intall the `renv` package if necessary. The run

```R
renv::restore()
```

to install all dependencies required for the analyses.

To reproduce the results, run

```R
source("code/analyses.R")
```

in the RStudio console. This will run both analyses (site 690 and Steinbruch Schmidt), generate all figues and save them in the `figs/` folder, and save the data under `data/res/site690_data.Rdata` and `data/res/sbs_data.Rdata`.

## Repository structure

* _code_ : folder for code
  * _analysis_steinbruch_schmidt.R_ : code for Steinbruch Schmidt example
  * _analysis_site690.R_ : code for the PETM example
  * _analyses.R_ : runs both examples
* _data_ : folder for data
  * _res_ : folder for resulst. Initially empty, filled once `code/analyses.R` is run
  * _raw_ : folder for raw data, read only!
    * _SbS_XRF_forfactor3.csv_ : data for Steinbruch Schmidt example, from da Silva (2020, 2024)
    * _Farley_and_Eltgroth_2003_supp_data_1_site690.csv_ : data for the PETM example, from Farley and Eltgroth (2003 a, b)
* _figs_ : folder for figures. Initially empty, filled once `code/analyses.R` is run
* _renv_ : folder for `renv` package
* _.gitignore_ : untracked files
* _.Rprofile_ : R session settings
* _LICENSE_ : Apache 2.0 license text
* _nonparam_est_adm_supp.Rproj_ : R project file
* _README_ : Readme file
* _renv.lock_ : lock file of the `renv` package

## References

Data in `data/raw/SbS_XRF_forfactor3.csv` and parts of the code in `code/analysis_steinbruch_schmidt.R` are from

* da Silva, A.-C. (2024). Anchoring the Late Devonian mass extinction in absolute time by integrating climatic controls and radio-isotopic dating: Supplementary code (v1.0.0). Zenodo. [DOI: 10.5281/zenodo.12516430](https://doi.org/10.5281/zenodo.12516430), In supplement to Da Silva, AC., Sinnesael, M., Claeys, P. et al. Anchoring the Late Devonian mass extinction in absolute time by integrating climatic controls and radio-isotopic dating. Sci Rep 10, 12940 (2020). [DOI: 10.1038/s41598-020-69097-6](https://doi.org/10.1038/s41598-020-69097-6)

Data in `data/raw/Farley_and_Eltgroth_2003_supp_data_1_site690.csv`
is from  

* Farley, Kenneth A; Eltgroth, Selene F (2003): (Appendix 1) Helium isotopic ratios and sedimentation rate model of ODP Hole 113-690B [dataset]. PANGAEA, [DOI: 10.1594/PANGAEA.723907](https://doi.org/10.1594/PANGAEA.723907), In supplement to: Farley, KA; Eltgroth, SF (2003): An alternative age model for the Paleocene-Eocene thermal maximum using extraterrestrial 3He. Earth and Planetary Science Letters, 208(3-4), 135-148, [DOI: 10.1016/S0012-821X(03)00017-7](https://doi.org/10.1016/S0012-821X(03)00017-7)

## Citation

Please cite this code as

* Hohmann, N. (2024). Supplementary data and code for "Nonparametric estimation of age-depth models from sedimentological and stratigraphic data" (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.13639816

## Copyright

Copyright 2023-2024 Netherlands eScience Center and Utrecht University.

## License

Apache 2.0 License, see LICENSE file for full license text.

## Funding information

Funded by the European Union (ERC, MindTheGap, StG project no 101041077). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.
![European Union and European Research Council logos](https://erc.europa.eu/sites/default/files/2023-06/LOGO_ERC-FLAG_FP.png)
