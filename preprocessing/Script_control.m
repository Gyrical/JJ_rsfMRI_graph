% List of open inputs
% Slice Timing: Session - cfg_files
% Coreg: Estimate: Source Image - cfg_files
% Segment: Data - cfg_files

subdirs = {   
%    'con_01';
%    'con_02';
%    'con_03';
    'con_04';
    'con_05';
    'con_06';
    'con_07';
    'con_08';
    'con_09';
    'con_10';
    'con_11';
    'con_12';
    'con_13';
    'con_14';
    'con_15';
    'con_16';
    'con_17';
    'con_18';
                    };
                
 basedir  = 'X:\Analyse\Normalisierung\controls_100\';  % common part of the path
 restdir = '\rest\';    % cell array for subdirs of rest data
 anadir  = '\ana\';     % cell array for subdirs of anatomy
 
nrun = length(subdirs); % enter the number of runs here
jobfile = {'X:\Analyse\Normalisierung\script\Script_control_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);
for crun = 1:nrun
    datadir = fullfile (basedir, subdirs {crun});
    inputs{1, crun} = cellstr(spm_select ('FPList', fullfile (datadir, restdir), '.*')); % Slice Timing: Session - cfg_files
    inputs{2, crun} = cellstr(spm_select ('FPList', fullfile (datadir, anadir), '^2.*\.*')); % Coreg: Estimate: Source Image - cfg_files
    inputs{3, crun} = cellstr(spm_select ('FPList', fullfile (datadir, anadir), '^2.*\.*')); % Segment: Data - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});
