# Repository for *Revealing the role of spatial heterogeneity in accelerating phase-to-trigger wave transitions in frog egg extracts*

[![DOI:10.1101/2024.01.18.576267](http://img.shields.io/badge/DOI-10.1101/2024.01.18.576267-000000.svg)](https://doi.org/10.1101/2024.01.18.576267) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10583185.svg)](https://doi.org/10.5281/zenodo.10583185)

This repository provides the codes used in the article *Revealing the role of spatial heterogeneity in accelerating phase-to-trigger wave transitions in frog egg extracts*. The repository includes codes to analyze experimental data of mitotic waves, codes to perform numerical simulations using the model introduced by Yang and Ferrell (1-2) using a pseudospectral method (3) and codes for reproducing the figures of the manuscript.

The preprint is available on [bioRxiv ](https://www.biorxiv.org/content/10.1101/2024.01.18.576267)[![DOI:10.1101/2024.01.18.576267](http://img.shields.io/badge/DOI-10.1101/2024.01.18.576267-000000.svg)](https://doi.org/10.1101/2024.01.18.576267)

### Python Dependencies
The codes provided depend on the Python packages.

```
numpy
scipy
matplotlib
statmodels
pandas
openpyxl
rasterio
f90nml
skimage
```
Fortran codes make use of the Intel Math Kernel Libraries.

### Install
Python
```
git clone https://github.com/DanielRuizReynes/mitotic_waves.git
```
Fortran
```
https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html
```
Installation takes less than an hour and runtime is typically on the order of minutes.

### References
1. Yang Q, Ferrell JE. The Cdk1-APC/C cell cycle oscillator circuit functions as a time-delayed, ultrasensitive switch. Nat Cell Biol. 2013;15(5):519–25. 
2. Chang JB, Ferrell JE. Mitotic trigger waves and the spatial coordination of the Xenopus cell cycle. Nature. 2013;500(7464):603–7.
3. Montagne R, Hernández-García E, Amengual A, San Miguel M. Wound-up phase turbulence in the complex Ginzburg-Landau equation. Phys Rev E. 1997;56(1):151.
