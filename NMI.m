function [ nmi ] = NMI( x, y )
% Compute the normalized mutual information to measure the accuracy of 
% the clustering result.
%
% INPUT: 
%       x, y: The cluster result and the ground truth. Each is a matrix
%       with two columns. The first column is the index of each node. The
%       second column is the corresponding cluster label.
%
% OUTPUT: 
%       nmi: the NMI score (0~1)
%
% Author: Peizhuo Wang (wangpeizhuo_37@163.com)
% Sep. 2016

Ncount = size(x,1);
NCLUSTER_x = max(x(:,2));
NCLUSTER_y = max(y(:,2));
cluster_x = cell(1, NCLUSTER_x);
cluster_y = cell(1, NCLUSTER_y);
for i = 1:NCLUSTER_x
    cluster_x{i} = x(x(:, 2)==i, 1);
end
count = 0;
for i = 1:NCLUSTER_y
    temp = y(y(:, 2)==i, 1);
    if (~isempty(temp))
        count = count+1;
        cluster_y{count} = temp;
    end
end
NCLUSTER_y = count;
cluster_y = cluster_y(1:NCLUSTER_y);

% Compute the confusion matrix
c = zeros(NCLUSTER_x, NCLUSTER_y);
for i = 1:NCLUSTER_x
    for j = 1:NCLUSTER_y
        c(i, j) = length(intersect(cluster_x{i}, cluster_y{j}));
    end
end

% Compute the normalized mutual information
H_x = 0;
H_y = 0;
MI = 0;
for i = 1:NCLUSTER_x
    ci = sum(c(i,:));
    H_x = H_x + ci*log(ci/Ncount);
    for j = 1:NCLUSTER_y
        cj = sum(c(:,j));
        if (i == 1)
            H_y = H_y + cj*log(cj/Ncount);
        end
        if (c(i,j) == 0)
            tempMI = 0;
        else
            tempMI = c(i,j)*log(Ncount*c(i,j)/(ci*cj));
        end
        MI = MI + tempMI;
    end
end

nmi = (-2)*MI/(H_x+H_y);

end