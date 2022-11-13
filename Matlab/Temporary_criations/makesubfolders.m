function makesubfolders( location, numofdiv,strofdiv )
%makesubfolders( location, numofdiv,strofdiv )
%
%   location -> location of the targert folder
%	numofdiv -> number of cells in the name of the product folder
%
% Exemple: makesubfolders( 'C:\TestFolder', 2,'_')
%     if there is a folder named MRI_02_01
%     then the result would be 
%     'C:\TestFolder\MRI_02\MRI_02_01'
%

    p=pwd;
    cd(location)
    %%

    [a,~]=dirff;

    %
    for i=1:size(a,2)
        %%
        b=strsplit(a(i).name,strofdiv);

        foldername=b{1};
        j=1;
        while j<numofdiv
            j=j+1;
            foldername=[foldername,strofdiv,b{j}];
        end
        if isequal(foldername,a(i).name)==0
            if exist(foldername,'dir')==0
                mkdir(foldername)
            end
            disp(foldername)
            movefile(a(i).name,foldername)
        end
    end

    cd(p)

end

