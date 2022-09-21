%% March 2014 Jessica Jesser
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
t = round((size(matrix,2).*size(matrix,2))-size(matrix,2)-(sparse.*(size(matrix,2).*(size(matrix,2))-size(matrix,2)))); %calculate number of edges to keep by sparsity
pruned_matrix = zeros(size(matrix,1),size(matrix,2),size(matrix,3)); %output matrix
for i = 1:(size(matrix,3))
    M = zeros(size(matrix,1),size(matrix,2));
    half_matrix = tril(matrix(:,:,i),-1); %lower half of matrix, because matrix is symmetrical
    vector_matrix = half_matrix(:);
    [sorted_vector,sorted_index] = sort(vector_matrix,'descend');
    for k = 1:t          
        M(sorted_index(k)) = sorted_vector(k);
    end;
    pruned_matrix(:,:,i) = M+M'; % fill output matrix
   
    if mean(degrees_und(pruned_matrix(:,:,i))) < log(size(matrix,1)) %check for mean node degree is bigger than log(nodes)
        disp('mean degree:'); 
        disp (mean(degrees_und(pruned_matrix(:,:,i))));
        disp('log(nodes):');
        disp(log(size(matrix,1)));
        error('mean degree is not bigger than log(nodes)');
    end;
       for m = 1:size(pruned_matrix,1) % check for matrix symmetry
           for n = 1:size(pruned_matrix,2)
               for o = 1: size(pruned_matrix,3)
                   if pruned_matrix(m,n,o) ~= pruned_matrix(n,m,o)
                       disp('matrix is not symmetrical');
                       error ('check matrix symmetry');
                   end;
               end;
           end;
       end;
       
       if mean(degrees_und(pruned_matrix(:,:,i))) > log(size(matrix,1)) %check for mean node degree is bigger than log(nodes)
        disp('mean degree:'); 
        disp (mean(degrees_und(pruned_matrix(:,:,i))));
        disp('log(nodes):');
        disp(log(size(matrix,1)));
       end
       for m = 1:size(pruned_matrix,1) % check for matrix symmetry
           for n = 1:size(pruned_matrix,2)
               for o = 1: size(pruned_matrix,3)
                   if pruned_matrix(m,n,o) == pruned_matrix(n,m,o);
                   end
               end
           end
       end; disp('matrix is symmetrical');

end    