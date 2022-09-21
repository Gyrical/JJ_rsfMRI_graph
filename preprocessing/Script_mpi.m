% List of open inputs
% Slice Timing: Session - cfg_files
% Coreg: Estimate: Source Image - cfg_files
% Coreg: Estimate: Other Images - cfg_files
% Segment: Data - cfg_files
% Segment: Masking image - cfg_files
% Normalise: Write: Images to Write - cfg_files
subdirs = {   
%    'sub_003'; 'sub_004';
%    'sub_010'; 'sub_013'; 
%    'sub_140'; 'sub_176'; 
%    'sub_193';
%    'sub_400'; 
     'sub_332';
%     'sub_587'; 
%     'sub_602'; 
%    'sub_578';
%    'sub_610';
%    'sub_593'; 
%     'sub_610';
%    'sub_612'; 'sub_613';
%    'sub_615'; 
                    };

 basedir  = 'X:\Analyse\Normalisierung\Patients_BCCI_200\Pre2\MPI\';  % common part of the path
 restdir = '\rest\';    % cell array for subdirs of rest data
 anadir  = '\ana\';     % cell array for subdirs of anatomy
 maskdir = '\mask\';    % cell array for subdirs of mask
                
                
nrun = length(subdirs); % enter the number of runs here
jobfile = {'X:\Analyse\Normalisierung\script\Script_mpi_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(6, nrun);
for crun = 1:nrun
     datadir = fullfile (basedir, subdirs {crun});
  
    inputs{1, crun} = cellstr(spm_select ('FPList', fullfile (datadir, restdir), '.*')); % Slice Timing: Session - cfg_files
    inputs{2, crun} = cellstr(spm_select ('FPList', fullfile (datadir, anadir), '.*')); % Coreg: Estimate: Source Image - cfg_files
    inputs{3, crun} = cellstr(spm_select ('FPList', fullfile (datadir, maskdir), '.*')); % Coreg: Estimate: Other Images - cfg_files
    inputs{4, crun} = cellstr(spm_select ('FPList', fullfile (datadir, anadir), '.*')); % Segment: Data - cfg_files
    inputs{5, crun} = cellstr(spm_select ('FPList', fullfile (datadir, maskdir), '^m.*\.nii')); % Segment: Masking image - cfg_files
    inputs{6, crun} = cellstr(spm_select ('FPList', fullfile (datadir, maskdir), '^sub.*\.nii')); % Normalise: Write: Images to Write - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});
