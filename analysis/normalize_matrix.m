%% April 2014 Jessica Jesser
%% September 2022, edited by Tianlu Wang
%% Matrix Normalization
%
% Description: of every matrix, negative values and main diagonal is set to zero. 
%
% Input: matrix: 3D matrix, 3rd dimension containing subjects
%
% Output: norm_matrix: normalized matrix 
%
% Scripts: no additional scripts needed
%%

function [norm_matrix] = normalize_matrix(matrix)

norm_matrix = zeros(size(matrix));

for subj = 1:size(matrix,3)
    
   mat = matrix(:,:,subj);
   
   % Set negative values to zero
   mat(mat<0) = 0;
   
   % Set main diagonal to 0
   mat(eye(size(mat,1))==1) = 0;
   
   % Save
   norm_matrix(:,:,subj) = mat;
   
end
