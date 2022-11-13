function fixfilename (folderpath,oldstr,newstr,name)
% fixname (oldstr,newstr,folderpath,name)
%
% this function changes unrecognizable characters of files names into
% recognizable ones and\or a specific string (oldstr) into a new one
% (newstr).
% 
% name = Part or full name of the target file or files
% folderpath = the path of the folder that contains the wanted files
% oldstr = The string that you want to change in the name
% newstr = The new string
%
% If the variable folderpath is not specificated, the function will act on
% the current folder. If the variables oldstr and newstr are not
% specificated, the function will only change the unrecognizable
% characters. 

%% valiables

if exist('folderpath','var')==1
    PATH=pwd;
    cd (folderpath)
end

if exist('name','var')==1
    [~,list]=dirff(0,pwd,name);
else
    [~,list]=dirff;
end

%% changing characters

fullfile=list';

if exist('oldstr','var')==1 && exist('newstr','var')==1
    fullfile=strchange(fullfile,newstr,oldstr);
end

% Change space to nothing
fullfile=strchange(fullfile,'_');

%list of problematic characters
oldchar={'(';')';'!';'$';'&';'§';'#';'@';'£';'¢';'¬';
    '¨';'~';'^';'´';'`';
    '+';'-';'=';'%';'º';'ª';'¹';'²';'³';
    'ç';
    'á';'à';'â';'ã';'Á';'À';'Â';'Ã';
    'é';'è';'ê';'É';'È';'Ê';
    'í';'ì';'î';'Í';'Ì';'Î';
    'ó';'ò';'õ';'ô';'Ó';'Ò';'Ô';'Õ';
    'ú';'ù';'û';'Ú';'Ù';'Û'};
% The change of character will happing one to one
newchar={' ';' ';' ';' ';' ';' ';'_';'_';' ';' ';' ';
    ' ';' ';' ';' ';' ';
    '_';'_';'_';'_';'_';'_';'_';'_';'_';
    'c';
    'a';'a';'a';'a';'A';'A';'A';'A';
    'e';'e';'e';'E';'E';'E';
    'i';'i';'i';'I';'I';'I';
    'o';'o';'o';'o';'O';'O';'O';'O';
    'u';'u';'u';'U';'U';'U'};

%Change of the characters
for i=1:size(oldchar,1)
    fullfile=strchange(fullfile,newchar{i},oldchar{i});
end


%% Rewriting name

for i=1:size(list,1)
    if isequal(list{i},fullfile{i})==0
        movefile(list{i},fullfile{i})
    end
end

if exist('PATH','var')==1
    cd (PATH)
end