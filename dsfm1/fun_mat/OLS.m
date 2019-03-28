function [Z, m, A, mse, Yhat] = OLS(Y, splines, L, maxiter, Zinit)
    % Compute OLS estimates for the DSFM.
    % Args:
    %   Y:        T*J dimensional matrix.
    %   splines:  a K*J dimensional matrix of splines basis functions.
    %   L:        number of factors.
    %   maxiter:  number of iterations for the algorithm.
    %   Zinit:    an optional starting value for the factor. It should be an
    %             T*L matrix.
    % Output:
    %   Z:        estimates of the latent factors.
    %   m:        the estimated functions: m = A * splines.
    %   A:        the estimated coefficients of the splines.
    %   mse:      Mean-squared error: mean((Y - Yhat).^2, 'all')
    %   Yhat:     estimated Y.

    [K, J] = size(splines);
    T = size(Y,1);
    spl2 = (splines * splines');

    % starting values for Z:
    if nargin<5
        Z = zscore(normrnd(0, 1, T, L));
    else
        Z = Zinit;
    end
    Zfull = [ones(T,1), Z];
        
    % Initializing Convergence measure:
    mse_new = 999;
    mse = zeros(maxiter,1);
    stop_crit=10;     % stoping criterias
    i = 0;
    
    while i < maxiter && mse_new > stop_crit
        i=i+1;            
        % Minimizing with respect to alpha
        alpha = minimize_alpha(Zfull, Y, splines, spl2, "OLS");
        A = reshape(alpha, [L+1 K]);

        % Minimize with respect to z:
        Z = minimize_Z(A, Y, splines, spl2, L, "OLS");
        Zfull = [ones(T,1), Z];
        
        % Check convergence:
        m = A * splines;
        Yhat = Zfull * m;
        mse_new = check_convergence(Y, Yhat, stop_crit, maxiter, i);
        mse(i) = mse_new;
    end    
    % Output:
    mse = mse(1:i);
end
