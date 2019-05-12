function [ simMatrix ] = edgelistFile2adjMatrix( networkfile, Ncount )
% Transfer the network edgelist format  in a text file to the adjacency
% matrix format matlab data.
%
% INPUT: 
%   networkfile: A text file contains edgelist data of a network.
%                The format of this file should be:
%                Column 1: id of element i
%                Column 2: id of element j
%                Column 3: edge weight (optional)
%   Ncount: the number of nodes
%
% OUTPUT: The adjacency matrix

xx=load(networkfile);
%Ncount=max(max(xx));
Ecount=size(xx,1);

if (size(xx, 2)<2)
    error('The input network format do not conform to the requirements')
end
simMatrix = zeros(Ncount);
if (size(xx, 2)==3)
    for i=1:Ecount
        ii=xx(i,1);
        jj=xx(i,2);
        simMatrix(ii,jj)=xx(i,3);
        simMatrix(jj,ii)=xx(i,3);
    end
end
if (size(xx, 2)==2)
    for i=1:Ecount
        ii=xx(i,1);
        jj=xx(i,2);
        simMatrix(ii,jj)=1;
        simMatrix(jj,ii)=1;
    end  
end
 
end