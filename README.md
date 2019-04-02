# Robust estimation of Dynamic Semiparametric Factor Model (DSFM) for fMRI data

The present repository contains the software used to estimate DSFM models with non-robust and robust methods on fMRI data, as it is performed in the paper "Robust sieve M-estimation with an application to dimensionality reduction".


## Languages
The software is designed for use on MATLAB. 
Notice that the function compute_splines.m implicitly makes use of R.
Therefore, if you do not have R, you should install it: see https://cran.r-project.org/.


## Usage and Contents

To perform inference with DSFM, it is needed to convert the 4D images to a matrix representation. Then one can convert the estimates back to a 4D object.
The folder contains following main functions:
- getDATA.m: convert a 4D (space and time) dimensional fMRI object into a matrix, usable for inference with DSFM.
- compute_splines.m: compute the splines and the 3D matrix of covariates, which gives the voxel index.
  It makes use of the package orthogonalsplinebasis (available on CRAN).
- OLS: Least-suares (non-robust) estimation of the DSFM.
- ROB: Robust estimation of the DSFM.
- getImage: Convert the matrix to the matrix representation.
- plotY: allow to plot matrices as 3D images. It relies on getImage and plot3D.
- plot3D: plot slices of a 4D object. The input time fixes the fourth dimension.


## Data

An exemplary fMRI data (for one subject), called data_example.nii.gz, 
can be found in this repository.



## Repository authors

- Julien Bodelet Ph.D. student, University of Geneva
- Davide La Vecchia Associate Professor in Statistics, University of Geneva


## Reference

- Robust sieve M-estimation with an application to dimensionality reduction





