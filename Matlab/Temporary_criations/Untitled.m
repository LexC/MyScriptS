close all
clear
clc
tic
%% Load and Variables

load D:\GD_UNICAMP\IC_NeuroFisica\Projetos\Coleta_NIRS_fMRI_2015-2017\Processed_data\fMRI\CorrMat_graphs.mat

%% Making mean matrices
tic
for i=16:25%size(CorrMatMean,2)
    %%
    disp(num2str(i))
    CorrMatMean(i).graphs=graphparameters(CorrMatMean(i).map);
    toc
end

save D:\GD_UNICAMP\IC_NeuroFisica\Projetos\Coleta_NIRS_fMRI_2015-2017\Processed_data\fMRI\CorrMat_graphs.mat
disp('--- saved ---')
%%
tic
for i=6:10%size(CorrMat,2)
    %%
    disp(num2str(i))
    CorrMat(i).graphs=graphparameters(CorrMat(i).map);
    toc
end
	save D:\GD_UNICAMP\IC_NeuroFisica\Projetos\Coleta_NIRS_fMRI_2015-2017\Processed_data\fMRI\CorrMat_graphs.mat
    disp('--- saved ---')