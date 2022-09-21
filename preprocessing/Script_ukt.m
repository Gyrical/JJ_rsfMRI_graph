% List of open inputs
% Slice Timing: Session - cfg_files
% Coreg: Estimate: Source Image - cfg_files
% Coreg: Estimate: Other Images - cfg_files
% Segment: Data - cfg_files
% Segment: Masking image - cfg_files
% Normalise: Write: Images to Write - cfg_files
subdirs = {   
%'sub_551'; 'sub_621'; 'sub_622'; 'sub_623';
%    'sub_089'; 'sub_207';
%    'sub_363'; 'sub_516';
%    'sub_523'; 'sub_533'; 
%    'sub_536'; 'sub_612'; 
%    'sub_613';
%     'sub_593'; 
%     'sub_612'; 
%     'sub_613'; 
%     'sub_615';
      'sub_258';
                    };

% basedir  = 'X:\neue_Hirne_von_Leo\norm\';  % common part of the path          
 basedir  = 'X:\Analyse\Normalisierung\Patients_BCI_200\Post2\UKT\';  % common part of the path
 restdir = '\rest\';    % cell array for subdirs of rest data
 anadir  = '\ana\';     % cell array for subdirs of anatomy
 maskdir = '\mask\';    % cell array for subdirs of mask
                
                
nrun = length(subdirs); % enter the number of runs here
jobfile = {'X:\Analyse\Normalisierung\script\Script_ukt_job.m'};
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
