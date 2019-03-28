function [Z, m, A, mse, Yhat] = ROB(Y, splines, Zols, Aols, tuning, maxiter, displ)
    % Compute ROBUST estimates for the DSFM.
    % Args:
    %   Y:        T*J dimensional matrix.
    %   splines:  a K*J dimensional matrix of splines basis functions.
    %   Zols:     estimated factors. It is used as a starting value for the
    %             algorithm.
    %   Aols:     estimated coefficients. It is used as a starting value 
    %             for the algorithm.
    %   tuning:   tuning constant.
    %   maxiter:  number of iterations for the algorithm.
    %   displ:    optional display argument. If provided (it can be can
    %             anything) the details of the algorithm will be displayed.
    % Output:
    %   Z:        estimates of the latent factors.
    %   m:        the estimated functions: m = A * splines.
    %   A:        the estimated coefficients of the splines.
    %   mse:      Mean-squared error: mean((Y - Yhat).^2, 'all')
    %   Yhat:     estimated Y.

    [K, J] = size(splines);
    T = size(Y,1);
    spl2 = (splines * splines');
    L = size(Aols, 1) - 1;
    med = median(abs(Y - [ones(T, 1), Zols] * Aols * splines))'; % Median

    
    % Optimization options:
    if(exist('displ') == 1) % display
        options = optimoptions('fminunc','Algorithm',...
        'trust-region','GradObj','on','Hessian','on',...
        'MaxIter',10,'TolFun',7.0000e-05,'Display','final');
    else % do not display
        options = optimoptions('fminunc','Algorithm',...
        'trust-region','GradObj','on','Hessian','on',...
        'MaxIter',10,'TolFun',7.0000e-05,'Display','off');
    end

    % Initializing Convergence measure:
    mse_new = 999;
    mse = zeros(maxiter,1);
    stop_crit=10;     % stoping criterias
    i = 0;
    Z = Zols;                % starting value
    Zfull = [ones(T,1), Z];
    alpha = reshape(Aols,[K*(L+1),1]);

    while i < maxiter && mse_new > stop_crit        
        i=i+1;  
        % Minimize with respect to alpha
        alpha = minimize_alpha(Zfull, Y, splines, spl2, "ROB", alpha, med, tuning, options);
        A = reshape(alpha, [L+1 K]);

        % Minimize with respect to z:
        Z = minimize_Z(A, Y, splines, spl2, L, "ROB", Z, med, tuning, options);
        Zfull = [ones(T,1), Z];
        
        % Check convergence:
        m = A * splines;
        Yhat = Zfull * m;
        mse_new = check_convergence(Y, Yhat, stop_crit, maxiter, i);
        mse(i) = mse_new;
    end
    mse = mse(1:i);
end
