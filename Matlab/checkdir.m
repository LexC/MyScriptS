function checkdir(PATH)
% Check if a directory exists. If doesn't, the scrit creates it.
if exist(PATH,'dir') ~= 7
    mkdir (PATH)
    disp('*** Folder Created ***')
end

end