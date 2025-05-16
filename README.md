---
editor_options: 
  markdown: 
    wrap: 72
---

# Intro

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13639816.svg)](https://doi.org/10.5281/zenodo.13639816)

Supplementary data and code for "Nonparametric estimation of age-depth
models from sedimentological and stratigraphic information" (Hohmann et
al. 2024, [DOI:
10.5194/egusphere-2024-2857](https://doi.org/10.5194/egusphere-2024-2857).

## Authors

**Niklas Hohmann**\
Utrecht University\
email: n.h.hohmann [at] uu.nl\
Web page: [www.uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)\
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

**Emilia Jarochowska**\
Utrecht University\
email: e.b.jarochowska [at] uu.nl\
Web page:
[www.uu.nl/staff/EBJarochowska](https://www.uu.nl/staff/EBJarochowska)\
ORCID: [0000-0001-8937-9405](https://orcid.org/0000-0001-8937-9405)

## Requirements

Base R (version \>= 4) and the RStudio IDE.

## Usage

Open the file `nonparam_est_adm_supp.Rproj` in RStudio. This will set
all paths correctly, and install the `renv` package if necessary. The
run

``` r
renv::restore()
```

to install all dependencies required for the analyses.

To reproduce the results, run

``` r
source("code/analyses.R")
```

in the RStudio console. This will run both analyses (site 1266 and
Steinbruch Schmidt), generate all figues and save them in the `figs/`
folder, and save the data under `data/res/site690_data.Rdata` and
`data/res/sbs_data.Rdata`.

## Repository structure

-   *code* : folder for code
    -   *analysis_steinbruch_schmidt.R* : code for Steinbruch Schmidt
        example
    -   *analysis_site1266.R* : code for the PETM example
    -   *plot_steinbruch_schmidt.R* : generate plots and statistics for
        sbs example
    -   *plots_site_1266.R* : generate plots and statistics for sbs
        example
    -   *analyses.R* : runs all examples
    -   Helper functions: *condensation_plot.R* *sed_rate_plot.R*
        *median_sedrate_l.R* *petm_recovery_stats.R*
-   *data* : folder for data
    -   *res* : folder for results
        -   *sbs_data.RData* : estimated age-depth models for steinbruch
            schmidt
        -   *site1266_data.RData* : estimated age-depth models for site
            1266
    -   *raw* : folder for raw data, read only!
        -   *SbS_XRF_forfactor3.csv* : data for Steinbruch Schmidt
            example, from da Silva (2020, 2024)
        -   *murphy_et_al_2010_1-s2.0-S0016703710003108-mmc3.txt* : raw
            data from ODP site 1266, from Murphy et al. (2010)
        -   *murphy_et_al_2010_1-s2.0-S0016703710003108-mmc3.csv* :
            tabular data from ODP site 1266, from Murphy et al. (2010)
-   *figs* : folder for figures. Initially empty, filled after
    `code/analyses.R` or the plotting scripts (code/plot\*.R) are run
-   *renv* : folder for `renv` package
-   *.gitignore* : untracked files
-   *.Rprofile* : R session settings
-   *LICENSE* : Apache 2.0 license text
-   *nonparam_est_adm_supp.Rproj* : R project file
-   *README* : Readme file
-   *renv.lock* : lock file of the `renv` package

## References

Data in `data/raw/SbS_XRF_forfactor3.csv` and parts of the code in
`code/analysis_steinbruch_schmidt.R` are from

-   da Silva, A.-C. (2024). Anchoring the Late Devonian mass extinction
    in absolute time by integrating climatic controls and radio-isotopic
    dating: Supplementary code (v1.0.0). Zenodo. [DOI:
    10.5281/zenodo.12516430](https://doi.org/10.5281/zenodo.12516430),
    In supplement to Da Silva, AC., Sinnesael, M., Claeys, P. et al.
    Anchoring the Late Devonian mass extinction in absolute time by
    integrating climatic controls and radio-isotopic dating. Sci Rep 10,
    12940 (2020). [DOI:
    10.1038/s41598-020-69097-6](https://doi.org/10.1038/s41598-020-69097-6)

The materials are protected by copyright and shared under the [CC BY 4.0
license](https://creativecommons.org/licenses/by/4.0/legalcode).

Data in `data/raw/murphy_et_al_2010_1-s2.0-S0016703710003108-mmc3.csv`
is from

-   Murphy, B. H., Farley, K. A., & Zachos, J. C. (2010). An
    extraterrestrial 3He-based timescale for the Paleocene--Eocene
    thermal maximum (PETM) from Walvis Ridge, IODP Site 1266. Geochimica
    et Cosmochimica Acta, 74(17), 5098--5108.
    <https://doi.org/10.1016/j.gca.2010.03.039>

The dataset is protected by copyright and shared under the [CC BY 3.0
license](https://creativecommons.org/licenses/by/3.0/legalcode.en)

## Citation

Please cite this code as

-   Hohmann, N. & Jarochowska, E. (2025). Supplementary data and code
    for "Nonparametric estimation of age-depth models from
    sedimentological and stratigraphic data" (v1.0.0). Zenodo.
    <https://doi.org/10.5281/zenodo.13639816>

## Copyright

Copyright 2023-2025 Utrecht University.

## License

Apache 2.0 License, see LICENSE file for full license text.

## Funding information

Funded by the European Union (ERC, MindTheGap, StG project no
101041077). Views and opinions expressed are however those of the
author(s) only and do not necessarily reflect those of the European
Union or the European Research Council. Neither the European Union nor
the granting authority can be held responsible for them. ![European
Union and European Research Council
logos](https://erc.europa.eu/sites/default/files/2023-06/LOGO_ERC-FLAG_FP.png)
