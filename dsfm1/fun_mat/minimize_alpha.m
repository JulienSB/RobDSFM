function alphanew = minimize_alpha(Zfull, Y, splines, spl2, type, alpha_init, med, tuning, options)
    % Minimize the loss functions with respect to alpha.
    % The following arguments: alpha_init, med, tuning, options
    %   are required only if type = "ROB".
    [T, J] = size(Y);
    K = size(splines, 1);
    L = size(Zfull, 2) - 1;
    
    % Compute the matrix of variance of Z:
    Z2 = zeros(L+1, L+1, T);
    for t =1:T
        Z2(:,:,t) = Zfull(t,:)' * Zfull(t,:);
    end

    
    % Solve for alpha:
    if type == "OLS"  
        S = 0;
        for t = 1:T
            S = 2*kron((splines * Y(t,: )'), Zfull(t, :)') + S;
        end
        Faa = R_F20_ols(spl2, Z2);
        alphanew = Faa \ S;
    else % type == "ROB"  
        obj_a = @(x) obj_alpha(Y, Zfull, x, Z2, splines, tuning, med, T, L, K, J);
        alphanew = fminunc(obj_a, alpha_init, options);
    end
end


function S2 = R_F20_ols(psi2,ZHAT2)
    S2 = 2 * kron(psi2, sum(ZHAT2, 3));
end

function [f, g, h] = obj_alpha(Y,Qfull,alphanew,ZHAT2,psimatrix,tuning,med,T,L,K,J)

    AHAT = reshape(alphanew, [L+1 K]);
    eps = Y-Qfull*AHAT*psimatrix;
    [M1, M2, M3] = wheighted(eps,tuning,med,T,J);

    f = sum(sum(M1));

    if nargout>1
        g = R_F10opt(M2,psimatrix,Qfull,L,K,J);
    end

    if nargout>2
        h = R_F20opt(M3,psimatrix,ZHAT2,T,L,K);
    end
end

function S = R_F10opt(M2, psimatrix, Qfull, L, K, J)
    S=zeros(K*(L+1),1);
    for j =1:J
        M2j = M2(:,j);
        S = kron(psimatrix(:,j), sum(repmat(-M2j,1,L+1).*Qfull)')+S;
    end
end

function S = R_F20opt(M3,psimatrix,ZHAT2,T,L,K)
    S = zeros(K*(L+1),K*(L+1));
    for t =1:T
        M3t = M3(t,:); M3t(M3t<0.0001) = 0; M3t=sparse(M3t);
        temp = kron(psimatrix*diag(M3t) * psimatrix',ZHAT2(:,:,t));
        S = temp + S;
    end
end

