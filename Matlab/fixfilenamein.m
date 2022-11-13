function  fixfilenamein (folders_in,folderpath,oldstr,newstr,name)
% 
% fixfilenamein (folders_in,folderpath,oldstr,newstr,name)
% 
% this function changes unrecognizable characters of files names into
% recognizable ones and\or a specific string (oldstr) into a new one
% (newstr), of an specific folder and it's subfolders.
% 
% folders_in: Number of intended subfolders.
% folderpath = the path of the folder that contains the wanted files
% oldstr = The string that you want to change in the name
% newstr = The new string
% name = Part or full name of the target file or files
% 
% If the variable folderpath is not specificated, the function will act on
% the current folder. If the variables oldstr and newstr are not
% specificated, the function will only change the unrecognizable
% characters. 


%% valiables
if exist('folderpath','var')==1
    [folders,~]=dirffin(folders_in,folderpath);
    folders=sort([folders;folderpath]);
else
    [folders,~]=dirffin(folders_in);
    folders=sort([folders;pwd]);
end

for i=1:size(folders,1)
    if exist('name','var')==1
        fixfilename (folders{i},oldstr,newstr,name)
    elseif exist('newstr','var')==1
        fixfilename (folders{i},oldstr,newstr)
    else
        fixfilename (folders{i})
    end
end

end

