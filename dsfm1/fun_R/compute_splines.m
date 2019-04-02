function [splines, X, Rpath] = compute_splines(ni, K_vec, Rpath)
% Computes the splines function and the covariates X (the coordinate of the voxels).
% Args:
%   ni :
%   K_vec :
    if nargin < 3
        Rpath = [];
    end
    old_cd = cd;    
    filepath = fileparts(fileparts(fileparts(which('location_fmri_dsfm1.m'))));
    cd(filepath)
    save('dsfm1\fun_R\input.mat', 'ni', 'K_vec') % save the intput
    
    try
        Rpath = RunRcode('dsfm1\fun_R\splines.R', Rpath);      
        out = load('dsfm1\fun_R\output.mat');        % load the ouput
        delete('dsfm1\fun_R\input.mat');
        delete('dsfm1\fun_R\output.mat');
        delete('dsfm1\fun_R\splines.R.log');
    catch ME
        error('Cannot find R path. Please provide Rpath or install R in the programm files.')
    end
    
    if exist('.RData') == 2
        delete('.RData')
    end
    if exist('.Rhistory') == 2
        delete('.Rhistory')
    end
    cd(old_cd)
    

    X = out.X; splines = out.splines;
    splines (splines < 0.001)=0;  % simplify the splines matrix
    splines = sparse(splines );   % sparse it
end




