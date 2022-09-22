%% Jessica Jesser, Tianlu Wang
%  September 2022
%% main analysis script
% Calculates the connectivity and graph measures, performs statistical
% analyses and and saves measures to *.mat 
% 
% Ran with BrainConnectivityToolbox version 2019_03_03_BCT in MATLAB R2016b
%
%% Initialization

clearvars; close all; clc
addpath('/home/tianlu/OneDrive/Code/BrainConnectivityToolbox/2019_03_03_BCT')
dir_main = fileparts(pwd);
fn_results = [dir_main, '/results/fc_graph_data.mat'];

%% Load data

fn_patients = dir([dir_main, '/data/patients/*.mat']);
fn_patients = table2struct(sortrows(struct2table(fn_patients)));
tc_patients = cell(length(fn_patients),1);
for idx = 1:length(fn_patients)
    tc_patients(idx) = struct2cell(load([fn_patients(idx).folder,'/',fn_patients(idx).name]));
end

fn_controls = dir([dir_main, '/data/controls/*.mat']);
fn_controls = table2struct(sortrows(struct2table(fn_controls)));
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
pat_ll = {'XR_sub_574_fl', 'XR_sub_089_fl','XR_sub_516_fl','XR_sub_533_fl','XR_sub_554_fl',...
          'XR_sub_263_fl','XR_sub_523_fl','XR_sub_536_fl','XR_sub_612_fl','XR_sub_613_fl'};
corr_ll = zeros(size(tc_patients{1},2),size(tc_patients{1},2),length(pat_ll));
for idx = 1:length(pat_ll)
    corr_ll(:,:,idx) = corrcoef(tc_patients{cellfun(@(x) strcmp(x,[pat_ll{idx},'.mat']), {fn_patients.name})});
end

pat_rl = {'XR_sub_133','XR_sub_546','XR_sub_563','XR_sub_176','XR_sub_048',...
          'XR_sub_004','XR_sub_567','XR_sub_587','XR_sub_510','XR_sub_400'};
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
% set to 0 and NO normalization is performed, is this correct? <<<

% Normalize connectivity matrices
norm_pat = normalize_matrix(corr_pat);
norm_con = normalize_matrix(corr_con);
norm_ll = normalize_matrix(corr_ll);
norm_rl = normalize_matrix(corr_rl);


%% Calculate whole-brain and homotopic functional connectivity

[data.pat_fcwb,data.pat_fch] = calculate_fc(norm_pat);
[data.con_fcwb,data.con_fch] = calculate_fc(norm_con);
[data.ll_fcwb,data.ll_fch] = calculate_fc(norm_ll);
[data.rl_fcwb,data.rl_fch] = calculate_fc(norm_rl);

%% Calculate graph parameters

[data.pat_cc, data.pat_cp, data.pat_preCGcc, data.pat_preCGcb] = calculate_graphparams(norm_pat);
[data.con_cc, data.con_cp, data.con_preCGcc, data.con_preCGcb] = calculate_graphparams(norm_con);
[data.ll_cc, data.ll_cp, data.ll_preCGcc, data.ll_preCGcb] = calculate_graphparams(norm_ll);
[data.rl_cc, data.rl_cp, data.rl_preCGcc, data.rl_preCGcb] = calculate_graphparams(norm_rl);


%% Save data for visualizations in Python

save(fn_results,'-struct','data')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Statistical analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
data = load(fn_results);
measures = {'fcwb','fch','cc','cp','preCGcc','preCGcb'};
groups = {'pat','con','ll','rl'};
keys_all = fieldnames(data);
keys_fc = keys_all(cellfun(@(x) contains(x,{'_fcwb','_fch'}),fieldnames(data)));
keys_gm = keys_all(cellfun(@(x) contains(x,{'_cc','_cp','_preCGcc','_preCGcb'}),fieldnames(data)));

%% Functional connectivity

% Test for normality
result_fc_normality = cell(length(keys_fc),4);
for idx = 1:length(keys_fc)
    [H,P,KSSTAT] = kstest(data.(keys_fc{idx}));
    result_fc_normality(idx,:) = {keys_fc{idx},H,P,KSSTAT};
end
disp('Kolmogorov-Smirnov goodness-of-fit hypothesis test results:')
disp(cell2table(result_fc_normality,'VariableNames',{'measure','H','P','KSSTAT'}))

% Compare all stroke patients to healthy controls
result_fc_patcon = cell(2,2);
for idx = 1:2
    pat = data.(strcat(groups{1},'_',measures{idx}));
    con = data.(strcat(groups{2},'_',measures{idx}));
    [P] = ranksum(pat, con);
    result_fc_patcon(idx,:) = {measures{idx},P};
end
disp('Wilcoxon rank sum test results patients vs. controls:')
disp(cell2table(result_fc_patcon,'VariableNames',{'measure','P'}))

% Compare stroke subgroups and healthy controls
GROUP = cellfun(@(x) repmat({x},1,length(data.(strcat(x,'_fch')))),groups(2:end),'UniformOutput',false);
GROUP = [GROUP{:}];
result_fc_anova = cell(2,4);
for idx = 1:2
    X = cell2mat(cellfun(@(x) data.(strcat(x,'_',measures{idx})),groups(2:end),'UniformOutput',false)');
    [P, TAB, STATS] = anova1(X, GROUP,'off'); 
    result_fc_anova(idx,1:3) = {measures{idx},P,TAB{2,5}};
    
    if P < 0.05 
        % Posthoc comparisons
        MC = multcompare(STATS);
        result_fc_anova(idx,4) = {MC};
    end
end
disp('One-way ANOVA test results patient subgroups vs. controls:')
disp(cell2table(result_fc_anova(:,1:3),'VariableNames',{'measure','P','F'}))
for idx = 1:2
    if not(isempty(result_fc_anova{idx,4}))
        disp(['Pairwise compairsons for: ',result_fc_anova{idx,1}]);
        disp(array2table(result_fc_anova{idx,4},'VariableNames',{'G1','G2','Lowerlimit','Difference','Upperlimit','Pvalue'})); 
    end
end




%% Graph measures

sparsities = 0.1:0.05:0.9;

% Compare graph metrics between patients and controls
result_gm_patcon = zeros(4,length(sparsities));
for idx = 3:length(measures)
    pat = data.(strcat(groups{1},'_',measures{idx}));
    con = data.(strcat(groups{2},'_',measures{idx}));
    plist = zeros(length(sparsities),1);
    for idxp = 1:size(pat,2)
        plist(idxp) = ranksum(pat(:,idxp), con(:,idxp));
    end
    result_gm_patcon(idx-2,:) = plist;
end
disp('Wilcoxon rank sum test results patients vs. controls (per sparsity level):')
disp(array2table(result_gm_patcon','VariableNames',measures(3:end), 'RowNames', cellstr(string(sparsities))))

% Compare stroke subgroups and healthy controls
GROUP = cellfun(@(x) repmat({x},1,length(data.(strcat(x,'_fch')))),groups(2:end),'UniformOutput',false);
GROUP = [GROUP{:}];
result_gm_anova = zeros(4,length(sparsities));
for idx = 3:length(measures)
    plist = zeros(length(sparsities),1);
    for idxp = 1:length(sparsities)
        X = cell2mat(cellfun(@(x) data.(strcat(x,'_',measures{idx}))(:,idxp),groups(2:end),'UniformOutput',false)');
        [P, TAB, STATS] = anova1(X, GROUP,'off'); 
        plist(idxp)=P;
    end
    result_gm_anova(idx-2,:) = plist;
end
disp('One-way ANOVA test results patient subgroups vs. controls (per sparsity level):')
disp(array2table(result_gm_anova','VariableNames',measures(3:end), 'RowNames', cellstr(string(sparsities))))











