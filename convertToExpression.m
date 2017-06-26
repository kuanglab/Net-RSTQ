function [Expression,rho] = convertToExpression(P1,TmpData)
    
    load CancerGeneNetworkwithHCMC1e5OverlapPPI_New L CancerGene_IX
    [m,n] = size(P1);
    Expression = zeros(m,n);
    rho = zeros(m,n);
    for i = 1:n
        Data = TmpData{i,1};
        for j = 1:size(Data,1)
            ix = find(CancerGene_IX==j);
            rho(ix,i) = (P1(ix,i)./L(ix))/sum(P1(ix,i)./L(ix));
            Expression(ix,i) = (P1(ix,i)./L(ix))*sum(Data{j,2})*50;
        end
    end
end