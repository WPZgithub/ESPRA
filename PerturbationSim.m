function [ Ap ] = PerturbationSim( deltaA, V, D )
% Compute the perturbated matrix based on a perturbation term and the
% eigenvectors and eigenvectors of the original matrix
%
% INPUT:
%       deltaA (N,N): The pertubation item between two matrices
%       V (N,N): Each column is a eigenvector of the original matrix
%       D (N,N): The diagonal elements are eigenvalues corresponding to the
%                eigenvectors
%
% OUTPUT:
%       Ap (N,N): The final perturbated matrix 
%
% Author: Peizhuo Wang <wangpeizhuo_37@163.com>
% Sep. 2016

n = size(deltaA,1);

% For a large matrix, we can just use the large eigenvalues whose percentage
% of cumulative sum is more than 0.9. This can decrease noise and reduce
% computing time.
[Y,I] = sort(abs(diag(D)),'descend');
cs = cumsum(Y)/sum(Y);
Ymin = Y(cs<0.9); r = length(Ymin);
Ddelta = zeros(n);
for i = 1:r
    index = I(i);
    x = V(1:n,index);
    Ddelta(index,index) = x'*deltaA*x/(x'*x);
end
D(I(r+1:end),I(r+1:end)) = 0;
Ap = V*(D+Ddelta)*V';
Ap(Ap<0) = 0;

% The result is transformed using a logistic function, such that for the
% element less than 0.3, we make it close to 0; for the element more than
% 0.6, we make it close to 1.
Ap_max = max(max(Ap));
Ap_min = min(min(Ap));
Ap = (Ap - ones(n)*Ap_min) ./ (Ap_max-Ap_min);
Ap = 1./(1+exp(log(9999)-2*log(9999)*Ap));

end