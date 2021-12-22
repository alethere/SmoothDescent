# SmoothDescent
This repository contains R code for a package implementing the Smooth Descent algorithm. Smooth Descent is a map-cleaning algorithm that can calculate identity by descent probabilities, predict genotyping errors and impute corrected genotypes. By applying this method, new genetic maps can be estimated to improve genetic map quality. The algorithm is based on SMOOTH (van Os et al. 2006), but adds some functionality to it. 

* Os, Hans & Stam, Piet & Visser, Richard & Eck, Herman. (2006). SMOOTH: A statistical method for successful removal of genotyping errors from high-density genetic linkage data. TAG. Theoretical and applied genetics. Theoretische und angewandte Genetik. DOI: 10.1007/s00122-005-0124-y. https://link.springer.com/article/10.1007/s00122-005-0124-y

# Installation
Install using `devtools`:
```
devtools::install_github("https://github.com/Alethere/SmoothDescent")
```
Download the package [binary file](https://github.com/Alethere/SmoothDescent/raw/master/SmoothDescent_0.1.0.tar.gz) and install it:
```
install.packages("file/location/SmoothDescent_0.1.0.tar.gz", repos = NULL)
```
