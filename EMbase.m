function P = EMbase(q,s,L,maxiter)
    [m, n] = size(q);
    P = ones(1,n)/n;
    q = (q./repmat(L,m,1))*max(L);
    
    for i = 1:maxiter
        P_old = P;
        a = (repmat(P,m,1).*q)./(repmat(q*P',1,n));
        P = sum(a.*repmat(s,1,n))/sum(sum(a.*repmat(s,1,n)));
        if max(abs(P-P_old))<1e-6
            break;
        end
        
    end
end
