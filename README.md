# ESPRA

This algorithm gives a similarity measure betweens nodes for dynamic networks based on structural perturbation theory, and detects dynamic community structure based on density-based clustering and evolutionary clustering. As described in the paper:
> Wang P, Gao L, Ma X. Dynamic community detection based on network structural perturbation and topological similarity[J]. Journal of Statistical Mechanics: Theory and Experiment, 2017, 2017(1): 013401.

This algorithm is implemented primarily in Matlab 2015b.

## Usage

### Main functions

- `ESPRA.m`: The implementation of the ESPRA algorithm.
```matlab
function [ result ] = ESPRA( simMatrices, alpha, beta )
% INPUT:
%   simMatrices (N,N,K): A series of symmetric similarity matrix for dynamic networks
%   alpha: Parameter for balancing the current clustering (=1) and historical influence (=0)
%   beta: Parameter for trading off the emphasis between the structural perturbation and topological similarity
%
% OUTPUT:
%   result: A cell that contains clustering results at every time step
```

- `PerturbationSim.m`: Compute the structural perturbation similarity.
```matlab
function [ Ap ] = PerturbationSim( deltaA, V, D )
% INPUT:
%       deltaA (N,N): The pertubation item between two matrices
%       V (N,N): Each column is a eigenvector of the original matrix
%       D (N,N): The diagonal elements are eigenvalues corresponding to the  eigenvectors
%
% OUTPUT:
%       Ap (N,N): The final perturbated matrix 
```

- `RA.m`: Compute the similarity matrix based on resource allocation (RA) index
```matlab
function [ simMatrix ] = RA( A )
% INPUT:
%       A: The adjacency matrix of a network
%
% OUTPUT:
%       simMatrix: The result similarity matrix based on RA index
```

- `cluster.m`: The implementation extends the clustering algorithm proposed by Alex Rodriguez and Alessandro Laio (Science, 2014).
```matlab
function [ cluster_result, rho, delta ] = cluster( simMatrix, nodeIDs )
% INPUT: 
%       simMatrix: Similarity matrix of the network
%       nodeIDs: A vector contains the identity of each node
%
% OUTPUT:
%       cluster_result (N, 2): The cluster result
%       rho (N, 1): The local density of each node
%       delta (N, 1): The distance of each node from nodes of higher local density
```
    Please reffer to the following reference for more details:
    Rodriguez A, Laio A. Clustering by fast search and find of density peaks[J].Science, 2014, 344(6191): 1492-1496.

### Synthetic data

The code is for generating simulated data:

- `gen_syn2.m`: Generating the synthetic networks with 4 clusters. This generator is implemented by Yu-Ru Lin (http://www.yurulin.com/download/code/facetnet.html).

### Functions for evaluation

The codes for evaluation measures:

- `NMI.m`: The normalized mutual information
- `ErrorRate.m`: The error rate

------------------------------------------------------------------

If you have any questions, please contact `wangpeizhuo_37@163.com`.