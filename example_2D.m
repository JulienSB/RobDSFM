% This is an example of how to use dsfm1 and plot results.
% Basically, one has to:
%   1- download the fMRI data 
%   2- select the voxels and convert it to a matrix format
%   3- perform robust or non-robust estimation
% This is an example with 2 dimensional data (64 * 64 voxels), that is one slice.


% add dsfm1 to path:
addpath(genpath('dsfm1'))


%% 1- Download fMRI Data Information:
%=====================================

fileName = 'data_example.nii.gz';  % example of data 
info = load_untouch_nii(fileName); 
sz = size(info.img);  % size of the data available
disp(sz);


%% 2- Loading matrix format:
%============================

% Choose number of voxels in each dimension:
xaxis = 1:64;               % 1D
yaxis = 1:64;               % 2D
zaxis = 15;                 % 3D here it is only the 15th slice ()
% choose the number of time points:
T = 180;                    
% number of splines in the three dimensions:
K_vec = [17, 17, 1];        % a good rule of thumb is:  round(ni / 4 + 1)
% number of voxels in the three dimensions:
ni = [size(xaxis, 2), size(yaxis, 2), size(zaxis, 2)]  
% get the data in the matrix format: Y is a T*J matrix, J = prod(ni)
[Y, X, Image, Data] = getDATA(xaxis, yaxis, zaxis, info, T); 
% Compute splines and the covariates X (this is the only function using R):
[splines, X] = compute_splines(ni, K_vec);  
% You can provide an explicit link to the R program. For example:
% [splines, X] = compute_splines(ni, K_vec, 'path-to-R\R-3.5.3\bin'); 
plotY(Y, ni);        % plot of Y (useful to plot a matrix as a 3D images)



%% 3- Estimation:
%================

maxiter = 10;  % maximum number of iterations
L = 2;         % number of latent factors

% Least Squares estimates:
[Z, m, A] = OLS(Y, splines, L, maxiter);  % Z is a T * L matrix, m = A * splines and is (L + 1) * J matrix
plotY(m, ni);         % plot of m0 
plotY(m, ni, 2);      % plot of m1
%plotY(m, ni, 3);     % plot of m2 if it exists 
plot(Z)               % plot of the factors 


% Robust estimates:
tuning = 2; % tuning constant (degree of robustness)
[Zr, mr, Ar] = ROB(Y, splines, Z, A, tuning, maxiter);    % Z and A are required and us as starting values
plotY(mr, ni);       % plot of m0
plotY(mr, ni, 2);    % plot of m1
%plotY(mr, ni, 3);   % plot of m2
plot(Zr)             % plot of the factors 



