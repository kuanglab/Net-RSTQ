function P = EM(P,alpha,q,s,L,maxiter)
    [m, n] = size(q);
    q = (q./repmat(L,m,1))*max(L);
    
    for i = 1:maxiter
        P_old = P;
        a = (repmat(P,m,1).*q)./(repmat(q*P',1,n));
        P = (sum(a.*repmat(s,1,n))+alpha'-1)/(sum(sum(a.*repmat(s,1,n)))+sum(alpha)-length(alpha));
        if max(abs(P-P_old))<1e-5
            break;
        end
    end
end
