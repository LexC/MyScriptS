close all
clear
clc
%% Variables

loc{1}='D:\GD_UNICAMP\IC_NeuroFisica\Projetos\Coleta_NIRS_fMRI_2015-2017\Processed_data\fMRI\Protocolo_01\Patients\Processed';
loc{2}='D:\GD_UNICAMP\IC_NeuroFisica\Projetos\Coleta_NIRS_fMRI_2015-2017\Processed_data\fMRI\Protocolo_02\Patients\Processed';
loc{3}='D:\GD_UNICAMP\IC_NeuroFisica\Projetos\Coleta_NIRS_fMRI_2015-2017\Processed_data\fMRI\Protocolo_02\Volunteers\Processed';

%%
u=1;
for i=1:size(loc,2)
    %%
    [~,b]=dirffin(3,loc{i},'Matrix-VAR.mat');
    %%
    for j=1:size(b,1)
        %%
        c=strsplit(loc{i},'Protocolo_');
        c=strsplit(c{2},'\');
        c=str2num(c{1});
       
        CorrMat(u).Protocolo=c;
        %%
        c=strsplit(loc{i},'\');
        c=c{end-1};
        c=c(1:(end-1));
        
        CorrMat(u).Type=c;
        %%
        c=strsplit(b{j},'\');
        c=strsplit(c{end-2},'_');
        
        CorrMat(u).Subject=str2num(c{3});
        CorrMat(u).Run=str2num(c{4});
        %%
        load(b{j})
        CorrMat(u).map=map;
        u=u+1;
    end
end
%% Making mean matrices

CorrMatMean=CorrMat;

u=2;
for i=2:size(CorrMat,2)
    CorrMatMean(i).flag=0;
    
    if CorrMatMean(i).Subject==CorrMatMean(i-1).Subject
        CorrMatMean(i-u+1).map(:,:,u)=CorrMatMean(i).map;
        CorrMatMean(i).flag=1;
        u=u+1;
    else
        u=2;
    end
end

u=1;
for i=1:size(CorrMat,2)
    
    if CorrMatMean(u).flag==1
        CorrMatMean(u)=[];
    else
        u=u+1;
    end
end

for i=1:size(CorrMatMean ,2)
    if size(CorrMatMean(i).map,3)>1
        CorrMatMean(i).map=fishermean( CorrMatMean(i).map,3);
    end
end

CorrMatMean = rmfield(CorrMatMean,'Run');
CorrMatMean = rmfield(CorrMatMean,'flag');

%% Making Group Mean of matrices
u=1;w=1;
for i=1:size(CorrMat,2)
    switch CorrMat(i).Type
        case 'Patient'
            pmap(:,:,u)=CorrMat(i).map;
            u=u+1;
        case 'Volunteer'
            vmap(:,:,w)=CorrMat(i).map;
            w=w+1;
    end
end

CorrMatGroupMean(1).Type='Patient';
CorrMatGroupMean(1).map=fishermean(pmap,3);
CorrMatGroupMean(2).Type='Volunteer';
CorrMatGroupMean(2).map=fishermean(vmap,3);

%%

save('D:\GD_UNICAMP\IC_NeuroFisica\Projetos\Coleta_NIRS_fMRI_2015-2017\Processed_data\fMRI\CorrMat.mat',...
    'CorrMat','CorrMatMean','CorrMatGroupMean')
