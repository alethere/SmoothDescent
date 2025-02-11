# SmoothDescent
This repository contains R code for a package implementing the Smooth Descent algorithm. Smooth Descent is a map-cleaning algorithm that can calculate identity by descent probabilities, predict genotyping errors and impute corrected genotypes. By applying this method, new genetic maps can be estimated to improve genetic map quality. This software was published here: [Therese Navarro et al. 2023](https://doi.org/10.3389/fgene.2023.1049988)

The algorithm is based on SMOOTH ([van Os et al. 2005](https:://doi.org/10.1007/s00122-005-0124-y)), but adds some functionality to it. 

# Installation
Install using `devtools`:
```
devtools::install_github("https://github.com/alethere/SmoothDescent")
```

# Vignette
You can preview the [vignette here](https://htmlpreview.github.io/?https://github.com/alethere/SmoothDescent/blob/master/doc/SmoothDescent_vignette.html), although some of the formatting does not work through github. You can also install the vignette and view it locally by using the code below (will take a couple minutes).
```
devtools::install_github("https://github.com/alethere/SmoothDescent", build_vignette = T)
browseVignettes("SmoothDescent")
```
