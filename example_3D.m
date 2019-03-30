% directory:
addpath(genpath('dsfm1'))


%% 1- fMRI Data Information:
%============================

fileName = 'data_example.nii.gz';  % example of file (this is in the data folder) 
info = load_untouch_nii(fileName); 
sz = size(info.img);  % size of the data available
disp(sz);


%% 2- Loading Data:
%===================
  
xaxis = 1:64;
yaxis = 1:64;
zaxis = 1:15;
T = 180;                 % choose the number of time points
ni = [size(xaxis, 2), size(yaxis, 2), size(zaxis, 2)] % vector of number of voxels in the three dimensions
K_vec = [17, 17, 5];        % vector of number of splines in each dimension, a good rule of thumb is:  round(ni / 4 + 1)
[Y, X, Image, Data] = getDATA(xaxis, yaxis, zaxis, info, T);  
[splines, X] = compute_splines(ni, K_vec);    % Compute splines and the covariates X 
% You can provide an explicit link for R. For example:
% [splines, X] = compute_splines(ni, K_vec, 'R-3.5.3\bin'); 
plotY(Y, ni);              % plot of Y



%% 3- Estimation:
%================

maxiter = 4;
L = 1; 

% Least Squares estimates:
[Z, m, A] = OLS(Y, splines, L, maxiter);  
plotY(m, ni);       % plot of m0
plotY(m, ni, 2);    % plot of m1
plot(Z)


% Robust estimates:
tuning = 2; % tuning constant
[Zr, mr] = ROB(Y, splines, Z, A, tuning, maxiter); 
plotY(mr, ni, 1);
plot(Zr) % Plot of the factors



