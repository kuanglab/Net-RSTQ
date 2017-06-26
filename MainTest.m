% (1)* |Install samtools| *
% Please download the samtools from http://sourceforge.net/projects/samtools/files/samtools/
% We suggest to install version 0.1.19 of samtools. Then put the compiled executable samtools file into the sourcecode folder
% (2)* |Put bam files in Data folder| *
% Please put the aligned bam file(s) by human hg19 annotation in the Data folder under NetRSTQ. 
% (3) lambda: tuning parameter between prior knowledge and alignment information, 0<lambda<1, the larger the lambda the more we trust network, suggest lambda<0.2;
% (4) PairedEnd: PairedEnd =1 if it is paired end data, PairedEnd = 0 if it is single end data; 
% Expression: the normalized expression of each transcript
% rho: the relative proportion of transcript in the gene
clc
clear all


% tuning parameter between prior knowledge and alignment information, 0<lambda<1, the larger the lambda the more we trust network, suggest lambda<0.2;
lambda =0.05;

% PairedEnd =1 if it is paired end data, PairedEnd = 0 if it is single end data; 
PairedEnd = 0;

[Expression,rho,TranscriptName,GeneName,SampleName] = NetRSTQ(lambda,PairedEnd);

save NetRSTQ_Result Expression rho TranscriptName GeneName SampleName;
