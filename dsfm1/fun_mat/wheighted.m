function [M1, M2, M3]= wheighted(eps,tuning,med,T,J)
    M1=zeros(T,J);
    M2=M1;
    M3=M1;
    for j=1:J
        medj=med(j);
        for t=1:T
            [a,b,c] = M_huber(eps(t,j),tuning,medj);
            M1(t,j) = a;
            M2(t,j) = b;
            M3(t,j) = c;
        end
    end 
end

function [sol1, sol2, sol3] = M_huber(y,b,w)
    % Huber loss function and derivatives.
    x=y/w;
    index= abs(x)<=b;


    sol1=x^2*index+(2*b*abs(x)-b^2)*(1-index);

    sol2=(2*x*index+2*b*sign(x)*(1-index))/w;

    sol3=2*index/w^2;
end
