%% March 2014 Jessica Jesser
% September 2022, edited by Tianlu Wang
%% Prune Matrix
%
% Description: thresholding of normalized matrix for different values of
% sparsity, checked for nlog > % knet
% 
% Input: matrix: normalized 3D Matrix (symmetrical in 1st and 2nd dimension, 3rd dimension are subjects), sparse = sparsity
%
% Output: pruned_matrix = thresholded 3D Matrix, t = number of kept edges.
% Every Matrix in 3rd dimension gets thresholded to a sparsity value and
% normalized to the mean of the matrix values. Mean degree of nodes must be bigger than log(nodes).
%%
      
function [pruned_matrix,t] = prune(matrix, sparse)
N = size(matrix,1);
t = round((1-sparse) * (N^2-N) * .5); % >>> Jessica: changed code here to only consider number of nodes in the lower triangle
pruned_matrix = zeros(size(matrix));

for subj = 1:size(matrix,3)
    M = zeros(N,N);
    
    % Prune matrix
    half_matrix = tril(matrix(:,:,subj),-1); 
    [sorted_vector,sorted_index] = sort(half_matrix(:),'descend');
    for k = 1:t          
        M(sorted_index(k)) = sorted_vector(k);
    end
    pruned_matrix(:,:,subj) = M+M'; 
   
    % Check that the mean degree is bigger than log(nodes)
    mean_deg = mean(degrees_und(pruned_matrix(:,:,subj)));
    log_nodes =  log(size(matrix,1));
%     fprintf('mean degree: %f\nlog nodes:%f\n',mean_deg,log_nodes) % Uncomment this line to show degree and log nodes per subject
    if mean_deg < log_nodes
        error('mean degree is not bigger than log(nodes)');
    end
    
    % Check for matrix symmetry
    if not(issymmetric(pruned_matrix(:,:,subj)))
        error('matrix is not symmetrical!')
%     else % Uncomment this and next line to print confirmation of symmetry
%         disp('matrix is symmetrical') 
    end
    
    
end    