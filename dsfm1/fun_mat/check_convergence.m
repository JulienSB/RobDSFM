function mse_new = check_convergence(Y, Yhat, stop_crit, maxiter, i)
    mse_new  = mean(mean((Y-Yhat).^2)); % Mean squared error
    formatSpec = 'Mse: %4.6f \n';
    fprintf(formatSpec,mse_new)            
    if i >= maxiter
        fprintf(['MaxIter reached: ',num2str(i),'\n'])
        formatSpec='Mse is %4.6f \n';
        fprintf(formatSpec,mse_new)            
    end

    if mse_new < stop_crit
        fprintf('Criterion reached \n')
        fprintf(['Iterations: ',num2str(i),'\n'])
    end
end