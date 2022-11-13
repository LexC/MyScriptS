function [ z ] = fishermean( z,dimension,nanflag)
%
%[ z ] = fishermean( z,dimension,nanflag)
%


z=fisher(z);

check1 = exist ('dimension','var');
check2 = exist ('nanflag');

if check1~=0
	check1=10;

end
if check2~=0 
	check2=100;
end

check=check1+check2;

switch check
    case 0
        z=mean (z);
    case 10
        z=mean (z,dimension);
    case 100
        z=mean (z,nanflag);
    case 110
        z=mean (z,dimension,nanflag);
end

z=fisherinv(z);
end

