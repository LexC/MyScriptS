function fullfile=strchange (list,newchar,oldchar)
% fullfile=strchange (list,newchar,oldchar)
%
% This function change one specific character into anoter in a list of
% strings 
%
% list -> it is cell with size=(1,x)
% newchar -> the new character
% oldchar -> the old character, if this variable is not declere the function
% will act on the space character.
%

for i=1:size(list,2)
    
    % Split the string in two cells
    if exist('oldchar','var')==0 
        a=strsplit(list{i});
    else
        a=strsplit(list{i},oldchar);
    end
    %
    
    fullfile{i}=a{1,1};
    
    % Glue the two cells if the string were slited
    if size(a,2)>1 
        for j=2:size(a,2)
            fullfile{i}=strcat(fullfile{i},newchar,a{1,j});
        end
    end
    %
    
end

end