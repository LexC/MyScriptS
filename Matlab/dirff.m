function [ folders,files] = dirff (fullfilepath,PATH_IN,name)
%
% [ folders,files] = dirff (fullfilepath,PATH_IN,name)
% 
%   Returns all the names of folders and files of one folder.
%
%   Explaning variables
% 
% folders: a cell type with the names of the folders.
%     If there is no folders:
%         folders=folders{1}='None'
% files: is a cell type with the names of the files.
%     If there is no files:
%         files=files{1}='None'
% 
% fullfilepath: if equal 1, files and folders will return with the full path.
% PATH_IN: to change the wanted location
% name: lists files that match with it

%% Head

% Changing folder
if exist ('PATH_IN','var') == 1
    PATH_OUT=pwd;
    cd(PATH_IN)
end

% Variables

A=dir;
u=1;w=1;

folders{1}='None';
files{1}='None';

% Excluding unwanted results
for i=1:2
    if isequal(A(1).name,'.')==1 || isequal(A(1).name,'..')==1
        A(1)=[];
    end
end

% If there are no content, stops the script
if size(A,2)==0
    return
end
%% Body

if exist ('name','var') == 0 
%% Results without the variable 'name'
    for i=1:size(A,1)
        if exist(A(i).name,'dir')==7
            folders{u,1}=A(i).name;
            u=u+1;
        else
            files{w,1}=A(i).name;
            w=w+1;
        end
    end
else
%% Results with the variable 'name'
    for i=1:size(A,1)
        if exist(A(i).name,'dir')==7
            folders{u,1}=A(i).name;
            u=u+1;
        end
    end
	B=dir(name);
    for i=1:size(B,1)
        files{w,1}=B(i).name;
        w=w+1;
    end
end
%%  Foot

% Adding path to the variables
if exist ('fullfilepath','var') == 1 && fullfilepath==1 && u~=1
    for i=1:size(folders,1)
        folders{i}=fullfile(pwd,folders{i});
    end
end
if exist ('fullfilepath','var') == 1 && fullfilepath==1 && w~=1
    for i=1:size(files,1)
        files{i}=fullfile(pwd,files{i});
    end
end

% changing the folders back to the original location
if exist ('PATH_IN','var') == 1
	cd(PATH_OUT)
end

end
