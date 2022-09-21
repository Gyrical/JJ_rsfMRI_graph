%% February 2021 Jessica Jesser
%% September 2022, edited by Tianlu Wang
%% Graph Parameter calculation
%
% Description: for every normalized matrix and across a sparsity range, graph
% parameters are calculated. 
%
% Input: normalized 3D correlation matrix, 3rd dimension holds subjects
%
% Output: for every subject (row) and sparsity level (column): 
%         cc=average cluster coefficient, 
%         cp=characteristic path length
%         preCGcc=cluster coefficient of PreCG healthy, 
%         preCGcb=centrality betweenness of PreCG healthy, 
%
% Scripts: prune.m
% Scripts from Brain connectivity toolbox: (weight_conversion.m,)
% distance_wei.m
% clustering_coef_wu.m, charpath.m, betweenness_wei.m, modularity.m
%%

function [cc, cp, preCGcc, preCGcb] = calculate_graphparams(matrix)
sparse = 0.1:0.05:0.9;

cc = zeros(size(matrix,3), length(sparse));
preCGcc = zeros(size(matrix,3), length(sparse));
cp = zeros(size(matrix,3), length(sparse));
preCGcb = zeros(size(matrix,3), length(sparse));

for j = 1:length(sparse)
    vector_cc = 1:size(matrix,3); % vector of mean of clustering coefficient for every matrix
    vector_pcc = 1:size(matrix,3); % vector of PreCG clustering coefficient for every matrix
    vector_cp = 1:size(matrix,3); 
    vector_cb = 1:size(matrix,3); 
 
    for i = 1:(size(matrix,3))
        pruned_matrix = prune(matrix(:,:,i),sparse(j));
        conv_matrix = distance_wei(weight_conversion(pruned_matrix, 'lengths'));
       %script adaption 01/2015 %conv_matrix = weight_conversion(pruned_matrix, 'lengths'); % conversion to distance matrix, necessary for some graphparameters

        cluster_coef = clustering_coef_wu(pruned_matrix); % calculation of cluster coefficient for every node

        vector_cc(i) = mean(cluster_coef); % calculation of average cluster coefficient
        cc(i,j) = vector_cc(i);

        %vector_pcc(i) = cluster_coef(1); % calculation of PreCG healthy (1st node) cluster coefficient, change to 116 for PreCG lesioned
        vector_pcc(i) = cluster_coef(1)/vector_cc(i); % calculation of normalized PreCG healthy (1st node) cluster coefficient, change to 116 for PreCG lesioned

        preCGcc(i,j) = vector_pcc(i);

        vector_cp(i) = charpath(conv_matrix); % charpath requires input distance matrix
        cp(i,j) = vector_cp(i);

        cb = betweenness_wei(conv_matrix); % betweenness_wei requires input distance matrix
        vector_cb(i)= cb(1); % calculation of PreCG healthy (1st node) centrality betweenness, change to 116 for PreCG lesioned
        preCGcb (i,j) = vector_cb(i);
    
    end
end