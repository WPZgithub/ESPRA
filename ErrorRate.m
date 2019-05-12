function [ errate ] = ErrorRate( x, y )
% Compute the error rate to measure the accuracy of the clustering result.
%
% INPUT: 
%       x, y: The cluster result and the ground truth. Each is a matrix
%       with two columns. The first column is the index of each node. The
%       second column is the corresponding cluster label.
%
% OUTPUT: 
%       errate: the error rate score
%
% Author: Peizhuo Wang <wangpeizhuo_37@163.com>
% Sep. 2016

Ncount=size(x,1);
NCLUSTER_x=max(x(:,2));
NCLUSTER_y=max(y(:,2));

z=zeros(Ncount, NCLUSTER_x);
g=zeros(Ncount, NCLUSTER_y);
for i = 1:Ncount
    z(i, x(i, 2))=1;
    g(i, y(i, 2))=1;
end

errate=trace((z*z'-g*g')'*(z*z'-g*g'));

end