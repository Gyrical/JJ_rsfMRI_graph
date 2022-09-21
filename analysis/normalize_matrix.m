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
% replace negative values and main diagnoal by 0
for i = 1:size(matrix,3)
   mat = matrix(:,:,i);
   ind = mat<0;
   mat(ind) = 0;
   norm_matrix(:,:,i) = mat;
   
   for j = 1:(size(matrix,1))
         norm_matrix(j,j,i) = 0;
   end
   
end
