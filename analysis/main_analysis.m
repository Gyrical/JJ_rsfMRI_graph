%% September 2022 Jessica Jesser, Tianlu Wang
%  
%% main analysis script
% Calculates the connectivity and graph measures, performs statistical
% analyses and saves measures to *.mat 
% 
% Ran with BrainConnectivityToolbox version 2019_03_03_BCT in MATLAB R2016b
%
%% Initialization

clearvars; close all; clc

dir_main = fileparts(pwd);
fn_data = [dir_main, '/results/fc_graph_data.mat'];
fn_stats = [dir_main, '/results/fc_graph_stats.mat'];

measures = {'fcwb','fch','cc','cp','preCGcc','preCGcb'};
groups = {'pat','con','ll','rl'};

load([dir_main, '/data/patients/lesionsize.mat'],'lesionsize');
      
%% Load timecourse data
disp('Loading data...')

fn_patients = dir([dir_main, '/data/patients/XR_sub_*.mat']);
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
disp('Calculating connectivity...')

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

% >>> Jessica: since it's confusing to use the 'normalize_matrix' function
% when no normalization takes place, I just did the preprocessing here in
% the main script

% Set diagonal to 0
corr_pat(repmat(eye(size(corr_pat,1)),1,1,size(corr_pat,3))==1) = 0;
corr_con(repmat(eye(size(corr_con,1)),1,1,size(corr_con,3))==1) = 0;

% use_fisher = false;
use_fisher = true;

if use_fisher
    % Fisher z-transform
    corr_pat_fi = log((1+corr_pat)./(1-corr_pat))./2;
    corr_con_fi = log((1+corr_con)./(1-corr_con))./2;
    
    % Set negative values to 0
    norm_pat = corr_pat_fi; norm_pat(norm_pat<0)=0;
    norm_con = corr_con_fi; norm_con(norm_con<0)=0;

%     % Other option: rescale Z-transformed values to the interval [0,1] (only if Z-transformed)
%     norm_pat = zeros(size(corr_pat));
%     for subj = 1:size(corr_pat,3)
%         mat = corr_pat_fi(:,:,subj);
%         norm_pat(:,:,subj) = (mat-min(mat(:))) ./ (max(mat(:))-min(mat(:)));
%     end
%     norm_con = zeros(size(corr_con));
%     for subj = 1:size(corr_con,3)
%         mat = corr_con_fi(:,:,subj);
%         norm_con(:,:,subj) = (mat-min(mat(:))) ./ (max(mat(:))-min(mat(:)));
%     end
else
    % Take absolute correlation between time courses without z-scoring
    norm_pat = abs(corr_pat);
    norm_con = abs(corr_con);
end


%% Calculate network measures
disp('Calculating network measures...')

if exist(fn_data,'file')==2 
    data = load(fn_data);
else
    data = struct();

    % Calculate whole-brain and homotopic functional connectivity
    [data.pat_fcwb,data.pat_fch] = calculate_fc(norm_pat);
    [data.con_fcwb,data.con_fch] = calculate_fc(norm_con);

    % Calculate graph parameters
    [data.pat_cc, data.pat_cp, data.pat_preCGcc, data.pat_preCGcb] = calculate_graphparams(norm_pat);
    [data.con_cc, data.con_cp, data.con_preCGcc, data.con_preCGcb] = calculate_graphparams(norm_con);

    % Get subgroup data
    pat_all = cellfun(@(x) {x(1:end-4)},{fn_patients.name}); 
    pat_ll = {'XR_sub_574_fl', 'XR_sub_089_fl','XR_sub_516_fl','XR_sub_533_fl','XR_sub_554_fl',...
              'XR_sub_263_fl','XR_sub_523_fl','XR_sub_536_fl','XR_sub_612_fl','XR_sub_613_fl'};
    pat_rl = {'XR_sub_133','XR_sub_546','XR_sub_563','XR_sub_176','XR_sub_048',...
              'XR_sub_004','XR_sub_567','XR_sub_587','XR_sub_510','XR_sub_400'};
    ll_idx = cellfun(@(x) find(strcmp(pat_all,x)),pat_ll);
    rl_idx = cellfun(@(x) find(strcmp(pat_all,x)),pat_rl);

    for idxm = 1:length(measures)
        meas = measures{idxm};
        data.(strcat('ll_',meas)) = data.(strcat('pat_',meas))(ll_idx,:);
        data.(strcat('rl_',meas)) = data.(strcat('pat_',meas))(rl_idx,:);
    end

    
    % Save network measures
    fprintf('Saving data to: %s...\n',fn_data)
    save(fn_data,'-struct','data')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                      Statistical analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization

keys_all = fieldnames(data);
keys_fc = keys_all(cellfun(@(x) contains(x,{'_fcwb','_fch'}),fieldnames(data)));
keys_gm = keys_all(cellfun(@(x) contains(x,{'_cc','_cp','_preCGcc','_preCGcb'}),fieldnames(data)));

results = struct();

%% Functional connectivity

% Test functional connectivity measures for normality
result_fc_normality = cell(length(keys_fc),4);
for idxm = 1:length(keys_fc)
    [H,P,KSSTAT] = kstest(data.(keys_fc{idxm}));
    result_fc_normality(idxm,:) = {keys_fc{idxm},H,P,KSSTAT};
end
disp('Kolmogorov-Smirnov goodness-of-fit hypothesis test results:')
disp(cell2table(result_fc_normality,'VariableNames',{'measure','H','Pval','KSSTAT'}))

% Compare all stroke patients to healthy controls
result_fc_patcon = cell(2,2);
for idxm = 1:2
    meas = measures{idxm};
    pat = data.(strcat(groups{1},'_',meas));
    con = data.(strcat(groups{2},'_',meas));
    [P] = ranksum(pat, con);
    result_fc_patcon(idxm,:) = {meas,P};
    results.(strcat('stats_all_',meas)) = P;
end
disp('Wilcoxon rank sum test results patients vs. controls:')
disp(cell2table(result_fc_patcon,'VariableNames',{'measure','Pval'}))

% Compare stroke subgroups and healthy controls
GROUP = cellfun(@(x) repmat({x},1,length(data.(strcat(x,'_fch')))),groups(2:end),'UniformOutput',false);
GROUP = [GROUP{:}];
result_fc_anova = cell(2,4);
for idxm = 1:2
    meas = measures{idxm};
    % Perform ANOVA
    X = cell2mat(cellfun(@(x) data.(strcat(x,'_',meas)),groups(2:end),'UniformOutput',false)');
    [P, TAB, STATS] = anova1(X, GROUP,'off'); 
    result_fc_anova(idxm,1:3) = {meas,P,TAB{2,5}};
    MC = multcompare(STATS);
        
    % Save results
    results.(strcat('stats_con_ll_rl_',meas)) = P;
    for comp = 1:size(MC,1)
        results.(sprintf('stats_%s_%s_%s',groups{MC(comp,1)+1},groups{MC(comp,2)+1},meas)) = MC(comp,end);
    end
    if P < 0.05 
        result_fc_anova(idxm,4) = {MC};
    end
end
disp('One-way ANOVA test results patients LL vs. RL vs. controls:')
disp(cell2table(result_fc_anova(:,1:3),'VariableNames',{'measure','Pval','F'}))
for idxm = 1:2
    if not(isempty(result_fc_anova{idxm,4}))
        disp(['Pairwise comparisons for: ',result_fc_anova{idxm,1}]);
        disp(array2table(result_fc_anova{idxm,4},'VariableNames',{'G1','G2','Lowerlimit','Difference','Upperlimit','Pvalue'})); 
    end
end


%% Graph measures

sparsities = 0.1:0.05:0.9;

% Compare graph metrics between patients and controls
result_gm_patcon = zeros(4,length(sparsities));
for idxm = 3:length(measures)
    meas = measures{idxm};
    pat = data.(strcat(groups{1},'_',meas));
    con = data.(strcat(groups{2},'_',meas));
    plist = zeros(length(sparsities),1);
    for idxs = 1:size(pat,2)
        plist(idxs) = ranksum(pat(:,idxs), con(:,idxs));
    end
    result_gm_patcon(idxm-2,:) = plist;
    results.(strcat('stats_all_',meas)) = plist;
end
disp('Wilcoxon rank sum test results patients vs. controls (per sparsity level):')
disp(array2table(result_gm_patcon','VariableNames',measures(3:end), 'RowNames', cellstr(string(sparsities))))

% Compare stroke subgroups and healthy controls
GROUP = cellfun(@(x) repmat({x},1,length(data.(strcat(x,'_fch')))),groups(2:end),'UniformOutput',false);
GROUP = [GROUP{:}];
result_gm_anova = zeros(4,length(sparsities));

for idxm = 3:length(measures)
    meas = measures{idxm};
    plist = zeros(length(sparsities),1);
    mclist = zeros(length(sparsities),3);
    
    for idxs = 1:length(sparsities)
        X = cell2mat(cellfun(@(x) data.(strcat(x,'_',measures{idxm}))(:,idxs), groups(2:end),'UniformOutput',false)');
        [P, TAB, STATS] = anova1(X, GROUP,'off'); 
        MC = multcompare(STATS);
        
        plist(idxs)=P;
        mclist(idxs,:) = MC(:,end);
    end
    
    result_gm_anova(idxm-2,:) = plist;
    results.(strcat('stats_con_ll_rl_',meas)) = plist;
    results.(strcat('stats_con_ll_',meas)) = mclist(:,1);
    results.(strcat('stats_con_rl_',meas)) = mclist(:,2);
    results.(strcat('stats_ll_rl_',meas)) = mclist(:,3);

end
disp('One-way ANOVA test results patients LL vs. RL vs. controls (per sparsity level):')
disp(array2table(result_gm_anova','VariableNames',measures(3:end), 'RowNames', cellstr(string(sparsities))))

%% Correlations with lesion size

disp('Correlation between lesion size and network measures:')
for idxm = 1:length(measures)
    meas = measures{idxm};
    results.(strcat('stats_lesionsize_',meas)) = zeros(size(data.(strcat('pat_',meas)),2),2);
    for idxs = 1:size(data.(strcat('pat_',meas)),2)
        [R,P] = corrcoef(data.(strcat('pat_',meas))(:,idxs),lesionsize);
        results.(strcat('stats_lesionsize_',meas))(idxs,:) = [R(1,2),P(1,2)];
    end
end
disp(array2table(cell2mat(cellfun(@(x) {results.(strcat('stats_lesionsize_',x))},measures(1:2))'),'RowNames',measures(1:2),'VariableNames',{'R','P'}))
disp('R per sparsity level:')
disp(array2table(cell2mat(cellfun(@(x) {results.(strcat('stats_lesionsize_',x))(:,1)}, measures(3:end))),'VariableNames',measures(3:end),'RowNames', cellstr(string(sparsities))))
disp('p per sparsity level:')
disp(array2table(cell2mat(cellfun(@(x) {results.(strcat('stats_lesionsize_',x))(:,2)}, measures(3:end))),'VariableNames',measures(3:end),'RowNames', cellstr(string(sparsities))))

%% Save results

fprintf('Saving results to: %s...\n',fn_stats)
save(fn_stats,'-struct','results')








