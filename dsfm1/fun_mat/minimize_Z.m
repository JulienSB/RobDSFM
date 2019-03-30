function Znew = minimize_Z(A, Y, splines, psi2, L, type, Z_init, med, tuning, options)
    % Minimize the loss functions with respect to Z.
    % The following arguments: Z_init, med, tuning, options
    %   are required only if type = "ROB".
    K = size(splines, 1);
    [T J] = size(Y);    
    A2 = A(2:(L+1), :);

    if type == "OLS"  
        a0 = A(1, :)'; 
        apsiA2 = a0' * psi2 * A2';
        A2psi = A2 * splines;
        A2psiA2 = A2 * psi2 * A2';
        for t = 1:T
            Znew(t,:) = (Y(t,:) * A2psi' - apsiA2) / A2psiA2;
        end
    else % type == "ROB"  
        ApsiAJ = zeros(L,L,J);
        for j = 1:J
            ApsiAJ(:, :, j) = A2*(splines(:, j) * splines(:, j)') * A2';
        end
        Apsi = A * splines;
        A2psi = A2 * splines;
        obj_z = @(x) obj_zhatP(Y, x, Apsi, A2psi, ApsiAJ, tuning, med, T, L, J);
        z = fminunc(obj_z, reshape(Z_init',[T*L 1]), options);
        Znew = reshape(z,[L T])';
    end
end


function [f, g, h]=obj_zhatP(Y,x,Apsi,A2psi,ApsiAJ,tuning,med,T,L,J)

    X=reshape(x,[L T])';
    Xfull=[ones(T,1) X];

    eps = Y - Xfull * Apsi;
    [M1, M2, M3] = wheighted(eps,tuning,med,T,J);

    f = sum(sum(M1));

    if nargout>1
        g = R_F01opt(M2,A2psi,T,L);
    end

    if nargout>2
     h = R_F02opt(M3,ApsiAJ,T,L,J);
    end
end


function S2= R_F01opt(M2,A2psi,T,L)

    S2=zeros(T*L,1);
    for t=1:T
        M2t=M2(t,:);
        S1=sum((repmat(M2t,L,1).*A2psi),2);
        S2((1:L)+(t-1)*L,:)=-S1;
    end

end 

function S2= R_F02opt(M3,ApsiAJ,T,L,J)
    % Cannot be parallelized.
    S2=zeros(T*L,T*L);

    for t=1:T
        M3t=M3(t,:);
       S2(((t-1)*L+1):(t*L),((t-1)*L+1):(t*L))=sum(reshape(kron(M3t,ones(L,L)),[L,L,J]).*ApsiAJ,3);
    end
end