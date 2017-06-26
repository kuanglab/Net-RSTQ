function TmpData = generateDataPaired(files)

    % load the isoforms names in the sanger cancer gene
    load CancerGeneNetworkwithHCMC1e5OverlapPPI_New CancerNMList
    % load human hg19 annotation, each row is an isoform
    load hg19
    [~,~,IX] = intersect(CancerNMList,TranscriptNames_hg19);
    G = GeneNames_hg19(IX);
    G = unique(G);
    clear IX

    Gene_S = zeros(length(G),1);
    Gene_E = zeros(length(G),1);
    Gene_Chr = cell(length(G),1);
    
    % check the start and end position of each gene
    for i = 1:length(G)
        IX = find(strcmp(GeneNames_hg19,G{i,1}));
        a = TranscriptStart_hg19(IX);
        b = TranscriptEnd_hg19(IX);
        c = unique(Chr_hg19(IX));
    
        a = min(a);
        b = max(b);
        
        Gene_S(i,1) = a;
        Gene_E(i,1) = b;
        Gene_Chr{i,1} = c{1,1};
    end
    
    IDX_G = [];
    IDX_T = [];

    for i = 1:length(G)
        ix = find(strcmp(G{i,1},GeneNames_hg19(1:33575)));
        IDX_G = cat(1,IDX_G,ones(length(ix),1)*i);
        IDX_T = cat(1,IDX_T,ix);
    end
    
    TmpData = cell(length(files),1);
    
    for i = 1:length(files)
        Data = cell(length(G),1);
        for j = 1:length(G)
            tmp = strcat('./samtools view -f 0X0040 Data/',files(i,1).name,{' '},Gene_Chr{j,1},':',num2str(Gene_S(j,1)),'-',num2str(Gene_E(j,1)),{' -o Data/'},G{j,1},'.1sam');
            system(tmp{1,1});
            tmp1 = strcat('./samtools view -f 0X0080 Data/',files(i,1).name,{' '},Gene_Chr{j,1},':',num2str(Gene_S(j,1)),'-',num2str(Gene_E(j,1)),{' -o Data/'},G{j,1},'.2sam');
            system(tmp1{1,1});
        
            s1 = dir(strcat('Data/',G{j,1},'.1sam'));
            s2 = dir(strcat('Data/',G{j,1},'.2sam'));
        
        
            if s1.bytes==0||s2.bytes==0
                Data{j,1} = 0;
            else
                fid=fopen(strcat('Data/',G{j,1},'.1sam'));
                data1=textscan(fid,'%s %*s %*s %f %*s %s %s %f %*[^\n]','Delimiter','\t','EndOfline','\n');
                fclose('all');
            
                fid=fopen(strcat('Data/',G{j,1},'.2sam'));
                data2=textscan(fid,'%s %*s %*s %f %*s %s %s %f %*[^\n]','Delimiter','\t','EndOfline','\n');
                fclose('all');
            
                [~,ix1,ix2] = intersect(data1{1,1},data2{1,1});
            
                if ~isempty(ix1)
                
                    for kk = 1:5
                        data1{1,kk} = data1{1,kk}(ix1);
                        data2{1,kk} = data2{1,kk}(ix2);
                    end
                else
                    Data{j,1} = 0;
                end
            
            
                IDX = find(IDX_G==j);
                nn = length(IDX);
                if nn==1
                    Data{j,1} = length(ix1); 
                else
                    tmp1 = cell(nn,1);
                    tmp2 = cell(nn,1);
                    tmp3 = cell(nn,1);
                
                    a = strfind(data1{1,3},'N');
                    idx1 = find(cellfun(@isempty,a));
                    idx1_v = find(~cellfun(@isempty,a));
                
                
                    for k = 1:nn
                        for mm = 1:ExonNum_hg19(IDX_T(IDX(k)))
                            idx = find(data1{1,2}(idx1)>=ExonStart_hg19(IDX_T(IDX(k)),mm)&data1{1,2}(idx1)<=ExonEnd_hg19(IDX_T(IDX(k)),mm));
                            tmp1{k,1} = cat(1,tmp1{k,1},idx1(idx));
                        end
                    end
                
                    if ~isempty(idx1_v)
                        for mm = 1:length(idx1_v)
                        
                            ixi1 = regexp(data1{1,3}(idx1_v(mm)),'\d');
                            ixi1 = setdiff(1:length(data1{1,3}{idx1_v(mm),1}),ixi1{1,1});
                            aa = max(find(ixi1<a{idx1_v(mm),1}(1)));
                            aa = ixi1(aa);
                            num1 = str2double(data1{1,3}{idx1_v(mm),1}(aa+1:a{idx1_v(mm),1}(1)-1));
                            
                            for k = 1:nn
                                aax = find(data1{1,2}(idx1_v(mm))<=ExonEnd_hg19(IDX_T(IDX(k)),1:ExonNum_hg19(IDX_T(IDX(k))))&(data1{1,2}(idx1_v(mm)))>=ExonStart_hg19(IDX_T(IDX(k)),1:ExonNum_hg19(IDX_T(IDX(k)))));
                                if ~isempty(aax);
                                    if aax ~= ExonNum_hg19(IDX_T(IDX(k)))
                                        % approximate the differences
                                        ap = ExonStart_hg19(IDX_T(IDX(k)),aax(1)+1)-ExonEnd_hg19(IDX_T(IDX(k)),aax(1));
                                        if abs(ap-num1)<=10;
                                            tmp1{k,1} = cat(1,tmp1{k,1},idx1_v(mm));
                                        end
                                    end
                                end
                            end
                        end
                    end
                
                    b = strfind(data2{1,3},'N');
                    idx2 = find(cellfun(@isempty,b));
                    idx2_v = find(~cellfun(@isempty,b));
                
                    for k = 1:nn
                        for mm = 1:ExonNum_hg19(IDX_T(IDX(k)))
                            idx = find(data2{1,2}(idx2)>=ExonStart_hg19(IDX_T(IDX(k)),mm)&data2{1,2}(idx2)<=ExonEnd_hg19(IDX_T(IDX(k)),mm));
                            tmp2{k,1} = cat(1,tmp2{k,1},idx2(idx));
                        end
                    end
                
                    if ~isempty(idx2_v)
                        for mm = 1:length(idx2_v)
                            ixi2 = regexp(data2{1,3}(idx2_v(mm)),'\d');
                            ixi2 = setdiff(1:length(data2{1,3}{idx2_v(mm),1}),ixi2{1,1});
                            bb = max(find(ixi2<b{idx2_v(mm),1}(1)));
                            bb = ixi2(bb);
                            num2 = str2double(data2{1,3}{idx2_v(mm),1}(bb+1:b{idx2_v(mm),1}(1)-1));
                            for k = 1:nn
                                bbx = find(data2{1,2}(idx2_v(mm))<=ExonEnd_hg19(IDX_T(IDX(k)),1:ExonNum_hg19(IDX_T(IDX(k))))&(data2{1,2}(idx2_v(mm)))>=ExonStart_hg19(IDX_T(IDX(k)),1:ExonNum_hg19(IDX_T(IDX(k)))));
                                if ~isempty(bbx);
                                    if bbx ~= ExonNum_hg19(IDX_T(IDX(k)))
                                        % approximate the differences
                                        ap = ExonStart_hg19(IDX_T(IDX(k)),bbx(1)+1)-ExonEnd_hg19(IDX_T(IDX(k)),bbx(1));
                                        if abs(ap-num2)<=10;
                                            tmp2{k,1} = cat(1,tmp2{k,1},idx2_v(mm));
                                        end
                                    end
                                end
                            end
                        end
                    end
                    for k = 1:nn
                        tmp3{k,1} = intersect(unique(tmp1{k,1}),unique(tmp2{k,1}));
                    end
                    Data{j,1} = tmp3;
                end  
            end
            delete(strcat('Data/',G{j,1},'.1sam'));
            delete(strcat('Data/',G{j,1},'.2sam'));
        end
 
        Data1 = cell(length(Data),1);
        Data2 = cell(length(Data),2);
        for j = 1:length(Data)
            if length(Data{j,1})==1
                Data1{j,1} = Data{j,1};
                Data2{j,1} = Data{j,1};
                Data2{j,2} = Data{j,1};
            else
                a = 0;
                nn = length(Data{j,1});
                for k = 1:nn
                    b = max(Data{j,1}{k,1});
                    if b>a
                        a=b;
                    end
                end
                if a == 0
                    Data2{j,1} = 0;
                    Data2{j,2} = 0;
                else
                    Data1{j,1} = zeros(a,nn);
                    D = zeros(a,1);
                    for k = 1:nn
                        Data1{j,1}(Data{j,1}{k,1},k) = 1;
                        D = D + Data1{j,1}(:,k)*(2^(nn-k+1));                
                    end
                    d = unique(D);
                    d = sort(d);
                    if d(1)==0
                        d(1) = [];
                    end
                    if isempty(d)
                        Data2{j,1} = 0;
                        Data2{j,2} = 0;
                    else
                        Data2{j,1} = zeros(length(d),nn);
                        Data2{j,2} = zeros(length(d),1);
                        for k = 1:length(d)
                            ix = find(D==d(k));
                            Data2{j,2}(k) = length(ix);
                            Data2{j,1}(k,:) = Data1{j,1}(ix(1),:);
                        end
                    end
                end
            end
        end
        TmpData{i,1} = Data2;  
    end
end