function [ cluster_result, rho, delta ] = cluster( simMatrix, nodeIDs )
% The implementation extends the clustering algorithm proposed by Alex
% Rodriguez and Alessandro Laio (Science, 2014).
% Reference:
%   Rodriguez A, Laio A. Clustering by fast search and find of density 
%   peaks[J]. Science, 2014, 344(6191): 1492-1496.
%
% INPUT: 
%       simMatrix: Similarity matrix of the network
%       nodeIDs: A vector contains the identity of each node
%
% OUTPUT:
%       cluster_result (N, 2): The cluster result
%       rho (N, 1): The local density of each node
%       delta (N, 1): The distance of each node from nodes of higher local density
%
% Modified by Peizhuo Wang <wangpeizhuo_37@163.com>
% Sep. 2016

distance = ones(size(simMatrix))-simMatrix; % transfer similarity to distance
distance = distance-diag(diag(distance));
Ncount = size(distance, 1);
Ecount = Ncount*(Ncount-1)/2;
xx = sparse(distance);
cluster_result = zeros(Ncount, 2);

percent = 2.0;
%fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);
position = round(Ecount*percent/100);
sda = sort(xx(:));%
dc = full(sda(position+Ncount));

% Compute rho
%fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
gaussian_distance = exp(-distance.^2./(dc*dc));
gaussian_distance = gaussian_distance-diag(diag(gaussian_distance));
rho = sum(gaussian_distance,2);
rho = (rho-ones(Ncount, 1)*min(rho))/(max(rho)-min(rho));

% Compute delta
maxd = max(max(distance));
[rho_sorted,ordrho] = sort(rho,'descend');
delta = zeros(Ncount, 1);
delta(ordrho(1)) = -1.;
nneigh(ordrho(1)) = 0;
for ii = 2:Ncount
   delta(ordrho(ii)) = maxd+0.00001;
   for jj = 1:ii-1
     if(distance(ordrho(ii),ordrho(jj)) < delta(ordrho(ii)))
        delta(ordrho(ii)) = distance(ordrho(ii),ordrho(jj));
        nneigh(ordrho(ii)) = ordrho(jj);
     end
   end
end
delta(ordrho(1)) = max(delta(:));
delta = (delta-ones(Ncount, 1)*min(delta))/(max(delta)-min(delta));

% Compute gamma
gamma = rho.*delta;
gamma = (gamma-ones(Ncount, 1)*min(gamma))/(max(gamma)-min(gamma));

% Select the cluster centers
NCLUST = 0;
clusters = ones(Ncount, 1)*(-1);
[gamma_sorted,ordgamma] = sort(gamma,'descend');         
gamma_calculated = gamma_sorted(21:end);
gammamin = mean(gamma_calculated)+2.33*std(gamma_calculated);% Z-score>=2.33, p-value<0.99
%gammamin=mean(gamma_calculated)+2*std(gamma_calculated);% Z-score >= 2
for i = 1:Ncount
    if (gamma_sorted(i) < gammamin)
        break;
    end
    index = ordgamma(i);
    
    NCLUST = NCLUST+1;
    clusters(index) = NCLUST;
    centers(NCLUST) = index;
end

% Nodes assignment
cluster_assign = cell(NCLUST, 1);
cl_no = 1:NCLUST;
for i = 1:NCLUST
    cluster_assign{i} = centers(i);
end
for i = 1:Ncount
  if (clusters(ordrho(i)) == -1)
    clusters(ordrho(i)) = clusters(nneigh(ordrho(i)));
    cluster_assign{clusters(ordrho(i))} = [cluster_assign{clusters(ordrho(i))},ordrho(i)];
  end
end

% We merge two candidate communities with the highest average similarity
% and iterate until the modularity is optimal.
if (NCLUST > 2) 
    result = [(1:Ncount)', clusters];
    Q = modularity(result, simMatrix);
    Qold = -100;
    clusters_new = clusters;
    centers_new = centers;
    while(Q > Qold)
        clusters_new = clusters;
        centers_new = centers;
        result = [(1:Ncount)', clusters];
        dist_cl = distance_clusters( result, distance );
        [row,col] = find(dist_cl == min(min(dist_cl)),1);
        
        % Merge these two clusters£¬and the center of the cluster with 
        % the relative high density is assigned as the new cluster center.
        if (gamma(centers(row)) > gamma(centers(col)))
            clusters(clusters==col) = row;
            centers(col) = 0;
            cl_no(col) = 0;
        else
            clusters(clusters==row) = col;
            centers(row) = 0;
            cl_no(row) = 0;
        end    
        
        cl_no = cl_no(cl_no>0);
        centers = centers(centers>0);
        NCLUST = length(centers);
        for i = 1:NCLUST
            clusters(clusters==cl_no(i)) = i;
        end
        cl_no = 1:NCLUST;

        % Compute the modularity
        Qold = Q;
        result = [(1:Ncount)', clusters];
        Q = modularity(result, simMatrix);
    end
    clusters = clusters_new;
    centers = centers_new;
    NCLUST = length(centers);
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);

for i = 1:Ncount
    cluster_result(i,:) = [nodeIDs(i),clusters(i)];
end

end

%% Compute the modularity
function [ Q ] = modularity( result, simMatrix )

NCLUST = max(result(:,2));
TS = sum(sum(simMatrix))/2.0; 
Q = 0;
for i = 1:NCLUST
    node_in = result(result(:,2)==i,1); 
    IS = sum(sum(simMatrix(node_in, node_in)))/2.0;
    DS = sum(sum(simMatrix(node_in, :))); 
    Q = Q + IS/TS-(DS/(2*TS))^2;
end
end

%% Compute the distance between each pair of clusters
function [ dist_cl ] = distance_clusters( result, distance )

NCLUST=max(result(:,2));
dist_cl=ones(NCLUST)*100;
for i=1:NCLUST-1
    cluster_i=result(result(:,2)==i,1);
    for j=i+1:NCLUST
        cluster_j=result(result(:,2)==j,1);
        dist_cl(i,j)=sum(sum(distance(cluster_i, cluster_j)))/(length(cluster_i)*length(cluster_j));
    end
end

end