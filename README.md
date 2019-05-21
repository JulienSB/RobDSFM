# Robust estimation of Dynamic Semiparametric Factor Model (DSFM) for fMRI data

The present repository contains the software used to estimate DSFM models with non-robust and robust methods on fMRI data, as it is performed in the paper "Robust sieve M-estimation with an application to dimensionality reduction".


## Languages
The software is designed for use on MATLAB. 
Notice that the function compute_splines.m implicitly make use of R.
Therefore, if you do not have R, you should install it: see https://cran.r-project.org/.


## Usage and Contents

To perform inference with DSFM, it is needed to convert the fMRI 4 dimensional (space and time) images to a matrix representation. Then one can convert the estimates back to 4D objects to visualize it.
The folder contains the following main functions:
- getDATA: convert a 4D (space and time) dimensional fMRI object into a matrix, usable for inference with DSFM.
- compute_splines: compute the splines and the 3D matrix of covariates, which gives the voxel index.
  The function makes use of the packages orthogonalsplinebasis and R.matlab (available on R CRAN).
- OLS: Least-squares (non-robust) estimation of the DSFM.
- ROB: Robust estimation of the DSFM.
- getImage: convert the matrix data to a 4D representation.
- plotY: plot matrices as images (slices of the brain). It relies on getImage and plot3D.
- plot3D: plot 4D object. The input time fixes the fourth dimension.


## Data

An exemplary fMRI data (for one subject), called data_example.nii.gz, 
can be found in this repository.



## Repository authors

- Julien Bodelet, Ph.D. student at the Research Center for Statistics (RCS), University of Geneva
- Davide La Vecchia, Associate Professor in Statistics at the Research Center for Statistics (RCS), University of Geneva


## Reference

- Bodelet, J. and La Vecchia D. (2019) Robust sieve M-estimation with an application to dimensionality reduction





