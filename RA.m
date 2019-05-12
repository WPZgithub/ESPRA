function [ simMatrix ] = RA( A )
% Compute the similarity matrix based on resource allocation (RA) index
% Reference:
%   Zhou T, L¨¹ L, Zhang Y C. Predicting missing links via local information[J]. 
%   The European Physical Journal B, 2009, 71(4): 623-630.
%
% INPUT:
%       A: The adjacency matrix of a network
%
% OUTPUT:
%       simMatrix: The result similarity matrix based on RA index
%
% Author: Peizhuo Wang (wangpeizhuo_37@163.com)
% Sep. 2016

AA = A ./ repmat(sum(A,2), [1,size(A, 1)]);

AA(isnan(AA)) = 0;
AA(isinf(AA)) = 0;
simMatrix = A * AA;
simMatrix = simMatrix .* A;% Only the original edge is considered
simMatrix_max = max(max(simMatrix));
simMatrix_min = min(min(simMatrix));
simMatrix = (simMatrix - ones(size(A, 1))*simMatrix_min) ./ (simMatrix_max-simMatrix_min);
simMatrix = simMatrix - diag(diag(simMatrix));

end