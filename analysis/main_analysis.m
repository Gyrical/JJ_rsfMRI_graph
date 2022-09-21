%% Jessica Jesser, Tianlu Wang
%  September 2022
%% Calculates the graph measures and saves them to *.mat 

%% Initialization

clear all; close all; clc
addpath('/home/tianlu/OneDrive/Code/BrainConnectivityToolbox/2019_03_03_BCT')
dir_main = fileparts(pwd);

%% Load data

fn_patients = dir([dir_main, '/data/patients/*.mat']);
tc_patients = cell(length(fn_patients),1);
for idx = 1:length(fn_patients)
    tc_patients(idx) = struct2cell(load([fn_patients(idx).folder,'/',fn_patients(idx).name]));
end

fn_controls = dir([dir_main, '/data/controls/*.mat']);
tc_controls = cell(length(fn_controls),1);
for idx = 1:length(fn_controls)
    tc_controls(idx) = struct2cell(load([fn_controls(idx).folder,'/',fn_controls(idx).name]));
end


%% Calculate connectivity matrices

% Patients
corr_pat = zeros(size(tc_patients{1},2),size(tc_patients{1},2),length(tc_patients));
for idx = 1:length(tc_patients)
    corr_pat(:,:,idx) = corrcoef(tc_patients{idx});
end

% Controls
corr_con = zeros(size(tc_controls{1},2),size(tc_controls{1},2),length(tc_controls));
for idx = 1:length(tc_controls)
    corr_con(:,:,idx) = corrcoef(tc_controls{idx});
end

% Subgroups
pat_ll = {'XR_sub_574_fl', 'XR_sub_089_fl','XR_sub_516_fl','XR_sub_533_fl','XR_sub_554_fl','XR_sub_263_fl','XR_sub_523_fl','XR_sub_536_fl','XR_sub_612_fl','XR_sub_613_fl'};
corr_ll = zeros(size(tc_patients{1},2),size(tc_patients{1},2),length(pat_ll));
for idx = 1:length(pat_ll)
    corr_ll(:,:,idx) = corrcoef(tc_patients{cellfun(@(x) strcmp(x,[pat_ll{idx},'.mat']), {fn_patients.name})});
end

pat_rl = {'XR_sub_133','XR_sub_546','XR_sub_563','XR_sub_176','XR_sub_048','XR_sub_004','XR_sub_567','XR_sub_587','XR_sub_510','XR_sub_400'};
corr_rl = zeros(size(tc_patients{1},2),size(tc_patients{1},2),length(pat_rl));
for idx = 1:length(pat_rl)
    corr_rl(:,:,idx) = corrcoef(tc_patients{cellfun(@(x) strcmp(x,[pat_rl{idx},'.mat']), {fn_patients.name})});
end


%% Process connectivity matrices

% Fisher z-transform
corr_pat = log((1+corr_pat)./(1-corr_pat))./2;
corr_con = log((1+corr_con)./(1-corr_con))./2;
corr_ll = log((1+corr_ll)./(1-corr_ll))./2;
corr_rl = log((1+corr_rl)./(1-corr_rl))./2;

% >>> Jessica: in normalize_matrix.m only negative and diagonal values are
% set to 0 and no normalization is performed, is this correct? <<<

% Normalize connectivity matrices
norm_pat = normalize_matrix(corr_pat);
norm_con = normalize_matrix(corr_con);
norm_ll = normalize_matrix(corr_ll);
norm_rl = normalize_matrix(corr_rl);


%% Calculate whole-brain and homotopic functional connectivity

[pat_fcwb,pat_fch] = calculate_fc(norm_pat);
[con_fcwb,con_fch] = calculate_fc(norm_con);
[ll_fcwb,ll_fch] = calculate_fc(norm_ll);
[rl_fcwb,rl_fch] = calculate_fc(norm_rl);

%% Calculate graph parameters

[pat_cc, pat_cp, pat_preCGcc, pat_preCGcb] = calculate_graphparams(norm_pat);
[con_cc, con_cp, con_preCGcc, con_preCGcb] = calculate_graphparams(norm_con);
[ll_cc, ll_cp, ll_preCGcc, ll_preCGcb] = calculate_graphparams(norm_ll);
[rl_cc, rl_cp, rl_preCGcc, rl_preCGcb] = calculate_graphparams(norm_rl);

%% Statistical analyses


%% Save data for visualizations in Python

fn_results = [dir_main, '/results/fc_graph_data.mat'];
save(fn_results,'pat_fcwb','pat_fch','pat_cc','pat_cp','pat_preCGcc','pat_preCGcb',...
                'con_fcwb','con_fch','con_cc','con_cp','con_preCGcc','con_preCGcb',...
                'll_fcwb','ll_fch','ll_cc','ll_cp','ll_preCGcc','ll_preCGcb',...
                'rl_fcwb','rl_fch','rl_cc','rl_cp','rl_preCGcc','rl_preCGcb')
 


