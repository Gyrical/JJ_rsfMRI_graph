%% February 2021 Jessica Jesser
%  September 2022, edited by Tianlu Wang
%% Graph Parameter calculation
%
% Description: for every normalized matrix and across a sparsity range, 
% graph parameters are calculated. 
%
% Input: normalized 3D correlation matrix, 3rd dimension holds subjects
%
% Output: for every subject (row) and sparsity level (column): 
%         cc=average cluster coefficient, 
%         cp=characteristic path length
%         preCGcc=cluster coefficient of PreCG healthy, 
%         preCGcb=centrality betweenness of PreCG healthy
%
% Note: set PreCG index (contralesional index=1, ipsilesional index=116)
%
% Scripts: prune.m
% Scripts from Brain connectivity toolbox: weight_conversion.m, 
% distance_wei.m, clustering_coef_wu.m, charpath.m, betweenness_wei.m
%%

function [cc, cp, preCGcc, preCGcb] = calculate_graphparams(matrix,sparsity,idxPreCG)

% Set default parameter values
if nargin < 3
    idxPreCG = 1;
end

% Settings
n_subj = size(matrix,3);

% Prepare output matrices
cc = zeros(n_subj, length(sparsity));
preCGcc = zeros(n_subj, length(sparsity));
cp = zeros(n_subj, length(sparsity));
preCGcb = zeros(n_subj, length(sparsity));

for subj = 1:n_subj
    for spars = 1:length(sparsity)
        
        % Prepare pruned and distance matrix
        pruned_matrix = prune(matrix(:,:,subj),sparsity(spars));
        conv_matrix = distance_wei(weight_conversion(pruned_matrix, 'lengths'));
        
        % Calculate average clustering coefficient
        cluster_coef = clustering_coef_wu(pruned_matrix);
        cc(subj,spars) = mean(cluster_coef);
        
        % Calculate PreCG clustering coefficient 
        preCGcc(subj,spars) = cluster_coef(idxPreCG)/mean(cluster_coef);

        % Calculate characteristic path with the distance matrix
        cp(subj,spars) = charpath(conv_matrix,0,0);
        
        % Calculate PreCG betweenness centrality with the distance matrix
        cb = betweenness_wei(conv_matrix);
        preCGcb(subj,spars) = cb(idxPreCG);
    
    end
end
