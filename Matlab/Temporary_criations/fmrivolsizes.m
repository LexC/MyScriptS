clear
clc
%%

[~,b]=dirff3('FiltRegrSW*');
%
for i=1:size(b,2)
    %%
    c=strsplit(b(i).name,'\');
    c=c{2};
    c=strsplit(c,'_');
    c{2}=c{2}(2:end);
    d=[str2num(c{2}) str2num(c{3})];
    
    stru=load_nii(b(i).name);
    f=size(stru.img,4);
    
    tabl(d(1),d(2))=f;
end



















