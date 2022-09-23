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
dir_main = fileparts(pwd);
fn_results = [dir_main, '/results/fc_graph_data.mat'];

      
%% Load data

fn_patients = dir([dir_main, '/data/patients/*.mat']);
fn_patients = table2struct(sortrows(struct2table(fn_patients))); % make sure filenames are sorted
tc_pat = cell(length(fn_patients),1);
for subj = 1:length(fn_patients)
    tc_pat(subj) = struct2cell(load([fn_patients(subj).folder,'/',fn_patients(subj).name]));
end


fn_controls = dir([dir_main, '/data/controls/*.mat']);
fn_controls = table2struct(sortrows(struct2table(fn_controls)));
tc_con = cell(length(fn_controls),1);
for subj = 1:length(fn_controls)
    tc_con(subj) = struct2cell(load([fn_controls(subj).folder,'/',fn_controls(subj).name]));
end


%% Calculate connectivity matrices

% Patients
corr_pat = zeros(size(tc_pat{1},2),size(tc_pat{1},2),length(tc_pat));
for subj = 1:length(tc_pat)
    corr_pat(:,:,subj) = corrcoef(tc_pat{subj});
end

% Controls
corr_con = zeros(size(tc_con{1},2),size(tc_con{1},2),length(tc_con));
for subj = 1:length(tc_con)
    corr_con(:,:,subj) = corrcoef(tc_con{subj});
end


%% Process connectivity matrices

% Fisher z-transform
corr_pat_fi = log((1+corr_pat)./(1-corr_pat))./2;
corr_con_fi = log((1+corr_con)./(1-corr_con))./2;

% Normalize connectivity matrices
norm_pat = normalize_matrix(corr_pat_fi);
norm_con = normalize_matrix(corr_con_fi);


%% Calculate network measures
data = struct();

% Calculate whole-brain and homotopic functional connectivity
[data.pat_fcwb,data.pat_fch] = calculate_fc(norm_pat);
[data.con_fcwb,data.con_fch] = calculate_fc(norm_con);

% Calculate graph parameters
[data.pat_cc, data.pat_cp, data.pat_preCGcc, data.pat_preCGcb] = calculate_graphparams(norm_pat);
[data.con_cc, data.con_cp, data.con_preCGcc, data.con_preCGcb] = calculate_graphparams(norm_con);


%% Get subgroup data

pat_all = cellfun(@(x) {x(1:end-4)},{fn_patients.name}); 
pat_ll = {'XR_sub_574_fl', 'XR_sub_089_fl','XR_sub_516_fl','XR_sub_533_fl','XR_sub_554_fl',...
          'XR_sub_263_fl','XR_sub_523_fl','XR_sub_536_fl','XR_sub_612_fl','XR_sub_613_fl'};
pat_rl = {'XR_sub_133','XR_sub_546','XR_sub_563','XR_sub_176','XR_sub_048',...
          'XR_sub_004','XR_sub_567','XR_sub_587','XR_sub_510','XR_sub_400'};
ll_idx = cellfun(@(x) find(strcmp(pat_all,x)),pat_ll);
rl_idx = cellfun(@(x) find(strcmp(pat_all,x)),pat_rl);

measures = {'fcwb','fch','cc','cp','preCGcc','preCGcb'};

for idxm = 1:length(measures)
    meas = measures{idxm};
    data.(strcat('ll_',meas)) = data.(strcat('pat_',meas))(ll_idx);
    data.(strcat('rl_',meas)) = data.(strcat('pat_',meas))(rl_idx);
end


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

results = struct();

%% Functional connectivity

% Test functional connectivity measures for normality
result_fc_normality = cell(length(keys_fc),4);
for meas = 1:length(keys_fc)
    [H,P,KSSTAT] = kstest(data.(keys_fc{meas}));
    result_fc_normality(meas,:) = {keys_fc{meas},H,P,KSSTAT};
end
disp('Kolmogorov-Smirnov goodness-of-fit hypothesis test results:')
disp(cell2table(result_fc_normality,'VariableNames',{'measure','H','P','KSSTAT'}))

% Compare all stroke patients to healthy controls
result_fc_patcon = cell(2,2);
for meas = 1:2
    pat = data.(strcat(groups{1},'_',measures{meas}));
    con = data.(strcat(groups{2},'_',measures{meas}));
    [P] = ranksum(pat, con);
    result_fc_patcon(meas,:) = {measures{meas},P};
end
disp('Wilcoxon rank sum test results patients vs. controls:')
disp(cell2table(result_fc_patcon,'VariableNames',{'measure','P'}))

% Compare stroke subgroups and healthy controls
GROUP = cellfun(@(x) repmat({x},1,length(data.(strcat(x,'_fch')))),groups(2:end),'UniformOutput',false);
GROUP = [GROUP{:}];
result_fc_anova = cell(2,4);
for meas = 1:2
    X = cell2mat(cellfun(@(x) data.(strcat(x,'_',measures{meas})),groups(2:end),'UniformOutput',false)');
    [P, TAB, STATS] = anova1(X, GROUP,'off'); 
    result_fc_anova(meas,1:3) = {measures{meas},P,TAB{2,5}};
    
    if P < 0.05 
        % Posthoc comparisons
        MC = multcompare(STATS);
        result_fc_anova(meas,4) = {MC};
    end
end
disp('One-way ANOVA test results patient subgroups vs. controls:')
disp(cell2table(result_fc_anova(:,1:3),'VariableNames',{'measure','P','F'}))
for meas = 1:2
    if not(isempty(result_fc_anova{meas,4}))
        disp(['Pairwise compairsons for: ',result_fc_anova{meas,1}]);
        disp(array2table(result_fc_anova{meas,4},'VariableNames',{'G1','G2','Lowerlimit','Difference','Upperlimit','Pvalue'})); 
    end
end




%% Graph measures

sparsities = 0.1:0.05:0.9;

% Compare graph metrics between patients and controls
result_gm_patcon = zeros(4,length(sparsities));
for meas = 3:length(measures)
    pat = data.(strcat(groups{1},'_',measures{meas}));
    con = data.(strcat(groups{2},'_',measures{meas}));
    plist = zeros(length(sparsities),1);
    for idxp = 1:size(pat,2)
        plist(idxp) = ranksum(pat(:,idxp), con(:,idxp));
    end
    result_gm_patcon(meas-2,:) = plist;
end
disp('Wilcoxon rank sum test results patients vs. controls (per sparsity level):')
disp(array2table(result_gm_patcon','VariableNames',measures(3:end), 'RowNames', cellstr(string(sparsities))))

% Compare stroke subgroups and healthy controls
GROUP = cellfun(@(x) repmat({x},1,length(data.(strcat(x,'_fch')))),groups(2:end),'UniformOutput',false);
GROUP = [GROUP{:}];
result_gm_anova = zeros(4,length(sparsities));
for meas = 3:length(measures)
    plist = zeros(length(sparsities),1);
    for idxp = 1:length(sparsities)
        X = cell2mat(cellfun(@(x) data.(strcat(x,'_',measures{meas}))(:,idxp),groups(2:end),'UniformOutput',false)');
        [P, TAB, STATS] = anova1(X, GROUP,'off'); 
        plist(idxp)=P;
    end
    result_gm_anova(meas-2,:) = plist;
end
disp('One-way ANOVA test results patient subgroups vs. controls (per sparsity level):')
disp(array2table(result_gm_anova','VariableNames',measures(3:end), 'RowNames', cellstr(string(sparsities))))

%% Save results










