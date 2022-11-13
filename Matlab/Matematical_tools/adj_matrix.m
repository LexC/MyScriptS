function [ AdjM, AdjMValues ] = adj_matrix ( CM, r )
%This script can be used to Matrix or Vectors
%
%   AdjM = Adjacence Matrix with 0 and 1 as elements
%   AdjMValues  = Adjacence Matrix with 0 and values pf CM as elements
%   CM = Input Matrix
%   r = threshold

[X,Y]=size(CM);
AdjM=zeros(X,Y);
AdjMValues=zeros(X,Y);
R = find (CM >= r );

AdjM(R)=1;
AdjMValues(R)=CM(R);


end