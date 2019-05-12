function [ result ] = ESPRA( simMatrices, alpha, beta )
% Evolutionary clustering based on Structural Perturbation and Resource
% Allocation similarity
%
% INPUT:
%   simMatrices (N,N,K): A series of symmetric similarity matrix for
%                        dynamic networks
%   alpha: Parameter for balancing the current clustering (=1) and historical
%          influence (=0)
%   beta: Parameter for trading off the emphasis between the structural
%         perturbation and topological similarity
%
% OUTPUT:
%   result: A cell that contains clustering results at every time step
%
% Author: Peizhuo Wang <wangpeizhuo_37@163.com>
% Sep. 2016

T = length(simMatrices);
result = cell(T,1);
simMatrix_perb = cell(T,1);

% clustering for the first network
disp(['Timestep ', num2str(1)])

A1 = simMatrices{1};
A2 = simMatrices{2};
[V,D] = eig(A1);
simMatrix_perb{1} = beta.*PerturbationSim((A2-A1),V,D)+(1-beta).*RA(A1);
nodes = find(sum(simMatrices{1})~=0);
result{1} = cluster(simMatrix_perb{1}(nodes,nodes),nodes);

for i = 2:T
    disp(['Timestep ', num2str(i)])
    
    A1 = simMatrices{i-1};
    A2 = simMatrices{i};
    [V,D] = eig(A2);
    B = beta.*PerturbationSim((A1-A2), V, D)+(1-beta).*RA(A2);
    simMatrix_perb{i} = alpha*B+(1-alpha)*simMatrix_perb{i-1};
    
    nodes = find(sum(simMatrices{i})~=0);
    result{i} = cluster(simMatrix_perb{i}(nodes,nodes),nodes);
end

end