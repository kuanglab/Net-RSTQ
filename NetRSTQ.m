function [Expression,rho,TranscriptName,GeneName,SampleName] = NetRSTQ(lambda,PariedEnd)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (1)* |Install samtools| *
    % Please download the samtools from http://sourceforge.net/projects/samtools/files/samtools/
    % We suggest to install samtools version 0.1.19. Then put the compiled executable samtools file into the folder NetRSTQ
    % (2)* |Put bam files in Data folder| *
    % Please put the aligned bam file(s) with human hg19 annotation in the Data folder under NetRSTQ. 
    % (3) lambda: tuning parameter between prior knowledge and alignment information, 0<lambda<1, the larger the lambda the more we trust network, suggest lambda<0.2;
    % (4) PariedEnd: PariedEnd =1 if it is paried end data, PairedEnd = 0 if it is single end data; 

    files = dir('Data/*.bam');
    if ~isempty(files)
        SampleName = cell(length(files),1);
        for i = 1:length(files)
            SampleName{i,1} = files(i,1).name(1:end-4);
            system(strcat('./samtools index Data/',files(i,1).name));
        end

        if lambda<=0||lambda>=1
            TranscriptName = 0;
            GeneName = 0;
            Expression = 0;
            rho=0;
            disp('Error: lambda should larger than 0 and smaller than 1');
        else
            if PariedEnd==1
                Data = generateDataPaired(files);
                [P1,TranscriptName,GeneName] = runNetRSTQ(Data,lambda);
                [Expression,rho] = convertToExpression(P1,Data);
            elseif PariedEnd==0
                Data = generateDataSingle(files);
                [P1,TranscriptName,GeneName] = runNetRSTQ(Data,lambda); 
                [Expression,rho] = convertToExpression(P1,Data);
            else
                Expression=0;
                TranscriptName = 0;
                GeneName = 0;
                rho=0;
                disp('Error: Please put a valide information, 1: paired end; 0: single end');
            end
        end

    else
        Expression = 0;
        TranscriptName = 0;
        GeneName = 0;
        SampleName = 0;
        rho=0;
        disp('Error: Please put bam files under folder NetRSTQ/Data/');
    end


end
