% This is an example with 3 dimensional data (64 * 64 * 15 voxels), that is 15 slices.

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
zaxis = 1:15;               % 3D here it is only the 15th slice ()

% choose the number of time points:
T = 180;                    
% number of splines in the three dimensions:
K_vec = [17, 17, 5];        % a good rule of thumb is:  round(ni / 4 + 1)
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

maxiter = 4;   % maximum number of iterations
L = 1;         % number of latent factors

% Least Squares estimates:
[Z, m, A] = OLS(Y, splines, L, maxiter);  % Z is a T * L matrix, m = A * splines and is (L + 1) * J matrix
plotY(m, ni);       % plot of m0
plotY(m, ni, 2);    % plot of m1
plot(Z)


% Robust estimates:
tuning = 2; % tuning constant (degree of robustness)
[Zr, mr] = ROB(Y, splines, Z, A, tuning, maxiter);  
plotY(mr, ni, 1);
plot(Zr) % Plot of the factors



