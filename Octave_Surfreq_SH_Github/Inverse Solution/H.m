% Calculate two-loop recursion
function H = H(g,sk,yk,scaling)


    [N,M] = size(sk);
    ro=zeros(M,1);
    alpha =zeros(M,1);
    beta =zeros(M,1);
     q = g;
    for i = 1:M
        ro(i,1) = 1/(yk(:,i)'*sk(:,i));
    end



    for i = M:-1:1
        alpha(i) = ro(i)*sk(:,i)'*q;
        q = q-alpha(i)*yk(:,i);
    end

    r = scaling.*q;

    for i = 1:M
        beta(i) = ro(i)*yk(:,i)'*r;
        r = r + sk(:,i)*(alpha(i)-beta(i));
    end
    %
    H=r;
end % end


