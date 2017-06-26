function [P,TranscriptName,GeneName] = runNetRSTQ(TmpData,lambda)
    
    load CancerGeneNetworkwithHCMC1e5OverlapPPI_New
    TranscriptName = CancerNMList;
    GeneName = CancerGeneList(CancerGene_IX);
    S = cell(length(CancerGeneList),1);
    IX_Gene = [];
    IX_Gene1 = [];
    for i = 1:length(CancerGeneList)
        tmp = find(CancerGene_IX==i);
        if length(tmp)>1
            tmp1 = sum(sum(CancerGeneNetwork(tmp,:)));
            IX_Gene1 = cat(1,IX_Gene1,i);
            if tmp1>0
                IX_Gene = cat(1,IX_Gene,i);
            end
        end
    end

    % Prepare network

    for i = 1:length(CancerGeneList)
        tmp = find(CancerGene_IX==i);
        S1 =CancerGeneNetwork(tmp,:);
        S1(:,tmp) = [];
        neigh = sum(S1,2);
        for j = 1:size(S1,1)
            if neigh(j)>0
                S1(j,:) = S1(j,:)/neigh(j);
            end
        end
        L1 =L';
        L1(tmp) = [];
        S1 = (S1./repmat(L1,size(S1,1),1)).*repmat(L(tmp),1,size(S1,2));
        S{i,1} = S1;
    end

    weight = lambda/(1-lambda);
    maxiter = 1000;
    
    %%%%%%%%%%% base EM
    nSample = length(TmpData);
    
    P_initial = ones(length(IDX),nSample);

    for i = 1:nSample
        Data = TmpData{i,1};
        for j = 1:length(IX_Gene1)
            tmp1 = find(CancerGene_IX==IX_Gene1(j));
            P_initial(tmp1,i) = ones(1,length(tmp1))/length(tmp1);
        end  
    
        for j = 1:length(IX_Gene1);
            tmp1 = find(CancerGene_IX==IX_Gene1(j));    
            %EM
            if Data{IX_Gene1(j),1}==0|length(Data{IX_Gene1(j),1})==1
                P_New = P_initial(tmp1,i);
            elseif size(Data{IX_Gene1(j),1},1)==1
                P_New = Data{IX_Gene1(j),1}/sum(Data{IX_Gene1(j),1});
            else
                P_New = EMbase(Data{IX_Gene1(j),1},Data{IX_Gene1(j),2},L(tmp1)',maxiter);
                if ~isempty(find(isnan(P_New)))                
                    P_New = P_initial(tmp1,i);
                end
            end
            P_initial(tmp1,i) = P_New;
        end
        clear Data;
    end
    
    %%%%% Net-RSTQ
    Expression = zeros(length(IDX),nSample);
    L_tmp = zeros(maxiter,nSample);
    it=100;
    for i = 1:size(P_initial,2)
        P_ALL = P_initial(:,i);
        Data = TmpData{i,1};
        % number of reads mapped to each gene
        rho = zeros(length(P_ALL),1);
        for j = 1:length(P_ALL)
            rho(j) = sum(Data{CancerGene_IX(j),2});
        end
        rho = rho*weight;
        
        for iter = 1:maxiter

            TotalLikelihood = 0;
            
            ixi = P_ALL==0;
            P_ALL(ixi) = 1e-7;
            
            for j = 1:length(IX_Gene)
                tmp = find(CancerGene_IX==IX_Gene(j));
                nn = length(tmp);
                P_bar = P_ALL;
                P_bar(tmp) = [];
                rho_bar = rho;
                rho_bar(tmp) = [];
                alpha = (S{IX_Gene(j),1}.*repmat(rho_bar',length(tmp),1))*P_bar+1;
                
                TotalLikelihood = TotalLikelihood + gammaln(sum(alpha))-sum(gammaln(alpha)) + sum((alpha-1).*log(P_ALL(tmp)));
                if Data{IX_Gene(j),1}==0
                else
                    q = (Data{IX_Gene(j),1}./repmat(L(tmp)',size(Data{IX_Gene(j),1},1),1))*max(L(tmp));
                    TotalLikelihood = TotalLikelihood + sum(Data{IX_Gene(j),2}.*log(q*P_ALL(tmp)));
                end
            end
            if isinf(TotalLikelihood)
                disp('TotalLikelihood inf');
            end
            L_tmp(iter, i) = TotalLikelihood;
            
            P_old = P_ALL;
            
            for j = 1:length(IX_Gene)
                tmp = find(CancerGene_IX==IX_Gene(j));
                nn = length(tmp);
                
                %%%%Prepare likelihood
                Likelihood = 0;
                %neighbor genes
                Neigh_T = find(sum(CancerGeneNetwork(tmp,:))>0);
                if ~isempty(Neigh_T)
                    Neigh_G = unique(CancerGene_IX(Neigh_T));
                
                    for kk = 1:length(Neigh_G)
                        ixii = find(CancerGene_IX==Neigh_G(kk));
                        P_bar = P_ALL;
                        P_bar(ixii) = [];
                    
                        rho_bar = rho;
                        rho_bar(ixii) = [];
                        alpha = (S{Neigh_G(kk),1}.*repmat(rho_bar',length(ixii),1))*P_bar+1;
                                
                        Likelihood = Likelihood + gammaln(sum(alpha))-sum(gammaln(alpha)) + sum((alpha-1).*log(P_ALL(ixii)));
         
                    end
                end
                %itself
                P_bar = P_ALL;
                P_bar(tmp) = [];
                rho_bar = rho;
                rho_bar(tmp) = [];
                alpha = (S{IX_Gene(j),1}.*repmat(rho_bar',length(tmp),1))*P_bar+1;
            
                Likelihood = Likelihood + gammaln(sum(alpha))-sum(gammaln(alpha)) + sum((alpha-1).*log(P_ALL(tmp)));
                if Data{IX_Gene(j),1}==0
                else
                    q = (Data{IX_Gene(j),1}./repmat(L(tmp)',size(Data{IX_Gene(j),1},1),1))*max(L(tmp));
                    Likelihood = Likelihood + sum(Data{IX_Gene(j),2}.*log(q*P_ALL(tmp)));
                end
                
                %%%%% EM
                if Data{IX_Gene(j),1}==0|length(Data{IX_Gene(j),1})==1
                    P_New = P_initial(tmp,i);
                elseif size(Data{IX_Gene(j),1},1)==1
                    P_New = (Data{IX_Gene(j),2}+alpha'-1)/(sum(Data{IX_Gene(j),2})+sum(alpha)-length(alpha));
                else
                    P_New = EM(P_ALL(tmp)',alpha,Data{IX_Gene(j),1},Data{IX_Gene(j),2},L(tmp)',it);
                end
                %%%%%% Update likelihood
                
                Likelihood_old = Likelihood;
                P_ALL_old = P_ALL;
                P_ALL(tmp) = P_New;
                
                ixi = find(P_ALL==0);
                P_ALL(ixi) = 1e-7;
                
                
                Likelihood = 0;
                %neighbor genes
                Neigh_T = find(sum(CancerGeneNetwork(tmp,:))>0);
                if ~isempty(Neigh_T)
                    Neigh_G = unique(CancerGene_IX(Neigh_T));
                
                    for kk = 1:length(Neigh_G)
                        ixii = find(CancerGene_IX==Neigh_G(kk));
                        P_bar = P_ALL;
                        P_bar(ixii) = [];
                    
                        rho_bar = rho;
                        rho_bar(ixii) = [];
                        alpha = (S{Neigh_G(kk),1}.*repmat(rho_bar',length(ixii),1))*P_bar+1;
                                
                        Likelihood = Likelihood + gammaln(sum(alpha))-sum(gammaln(alpha)) + sum((alpha-1).*log(P_ALL(ixii)));
         
                    end
                end
                %itself
                P_bar = P_ALL;
                P_bar(tmp) = [];
                rho_bar = rho;
                rho_bar(tmp) = [];
                alpha = (S{IX_Gene(j),1}.*repmat(rho_bar',length(tmp),1))*P_bar+1;
            
                Likelihood = Likelihood + gammaln(sum(alpha))-sum(gammaln(alpha)) + sum((alpha-1).*log(P_ALL(tmp)));
                if Data{IX_Gene(j),1}==0
                else
                    q = (Data{IX_Gene(j),1}./repmat(L(tmp)',size(Data{IX_Gene(j),1},1),1))*max(L(tmp));
                    Likelihood = Likelihood + sum(Data{IX_Gene(j),2}.*log(q*P_ALL(tmp)));
                end
                
                if Likelihood>=Likelihood_old
                else
                    P_ALL = P_ALL_old;
                end
  
            end

            if max(abs(P_ALL-P_old))<1e-2
                it = 100;
            elseif max(abs(P_ALL-P_old))<1e-4
                it = 1000;
            else
                it = 20;
            end
            
            if max(abs(P_ALL-P_old))<1e-5
                break;
            end     
        end
        Expression(:,i) = P_ALL;

    end
    P = Expression;
    clear L_tmp;


end