masks = { 'wsub_003'; 'wsub_004'; 'wsub_010'; 
          'wsub_011'; 'wsub_012'; 'wsub_013';
          'wsub_048'; 'wsub_089'; 'wsub_133';
          'wsub_140'; 'wsub_155'; 'wsub_176';
          'wsub_193'; 'wsub_207'; 'wsub_241';
          'wsub_258'; 'wsub_263'; 'wsub_332';
          'wsub_363'; 'wsub_400'; 'wsub_412';
          'wsub_453'; 'wsub_504'; 'wsub_510';
          'wsub_516'; 'wsub_519'; 'wsub_523';
          'wsub_533'; 'wsub_536'; 'wsub_546';
          'wsub_554'; 'wsub_563'; 'wsub_567';
          'wsub_574'; 'wsub_578'; 'wsub_579';
          'wsub_584'; 'wsub_587'; 'wsub_593';
          'wsub_602'; 'wsub_610'; 'wsub_612';
          'wsub_613'; 'wsub_615'; 
          };

lesionsize = cell(2+length(masks),1); 
lesionsize(1,1) = {'subjects'};
lesionsize(1,2) = {'total voxels'};

for k = 1:length(masks)
    lesionsize(1+k,1) = masks(k); 
end 

for kk = 1:length(masks)
matrix = spm_read_vols(spm_vol([masks{kk} '_lesion.nii']));
nonzero = length(find(matrix>0));
lesionsize{1+kk,2} = nonzero;
end

clear k kk nonzero matrix
