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

% Data
dir_main = fileparts(pwd);
load([dir_main, '/data/patients/lesionsize.mat'],'lesionsize');
fn_data = [dir_main, '/results/fc_graph_data.mat'];
fn_stats = [dir_main, '/results/fc_graph_stats.mat'];

% Settings
measures = {'fcwb','fch','cc','cp','preCGcc','preCGcb'};
groups = {'pat','con','ll','rl'};
n_roi = 116;
sparsities = 0.1:0.05:0.9;

      
%% Load timecourse data
disp('Loading data...')

fn_patients = dir([dir_main, '/data/patients/XR_sub_*.mat']);
fn_patients = table2struct(sortrows(struct2table(fn_patients))); % make sure filenames are sorted
tc_pat = cell(length(fn_patients),1);
for subj = 1:length(fn_patients)
%     fprintf('%d\t%s\n',subj,fn_patients(subj).name) % Uncomment to check order of files 
    tc_pat(subj) = struct2cell(load([fn_patients(subj).folder,'/',fn_patients(subj).name]));
end

fn_controls = dir([dir_main, '/data/controls/*.mat']);
fn_controls = table2struct(sortrows(struct2table(fn_controls)));
tc_con = cell(length(fn_controls),1);
for subj = 1:length(fn_controls)
%     fprintf('%d\t%s\n',subj,fn_controls(subj).name) % Uncomment to check order of files 
    tc_con(subj) = struct2cell(load([fn_controls(subj).folder,'/',fn_controls(subj).name]));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate network properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculate connectivity matrices
disp('Calculating connectivity...')

% Patients
corr_pat = zeros(n_roi,n_roi,length(tc_pat));
for subj = 1:length(tc_pat)
    corr_pat(:,:,subj) = corrcoef(tc_pat{subj});
end

% Controls
corr_con = zeros(n_roi,n_roi,length(tc_con));
for subj = 1:length(tc_con)
    corr_con(:,:,subj) = corrcoef(tc_con{subj});
end


%% Process connectivity matrices

% Fisher z-transform
corr_pat_fi = log((1+corr_pat)./(1-corr_pat))./2;
corr_con_fi = log((1+corr_con)./(1-corr_con))./2;

% Set negative values to 0
norm_pat = corr_pat_fi; norm_pat(norm_pat<0)=0;
norm_con = corr_con_fi; norm_con(norm_con<0)=0;

% Set diagonal values to 0
norm_pat(repmat(eye(n_roi),1,1,size(norm_pat,3))==1) = 0;
norm_con(repmat(eye(n_roi),1,1,size(norm_con,3))==1) = 0;



%% Calculate network measures
disp('Calculating network measures...')

if exist(fn_data,'file') == 2 
    data = load(fn_data);
else
    data = struct();

    % Calculate whole-brain and homotopic functional connectivity
    [data.pat_fcwb,data.pat_fch] = calculate_fc(norm_pat);
    [data.con_fcwb,data.con_fch] = calculate_fc(norm_con);

    % Calculate graph parameters
    [data.pat_cc, data.pat_cp, data.pat_preCGcc, data.pat_preCGcb] = calculate_graphparams(norm_pat,sparsities);
    [data.con_cc, data.con_cp, data.con_preCGcc, data.con_preCGcb] = calculate_graphparams(norm_con,sparsities);

    
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

    data.lesionsize = lesionsize';
    
    
    % Save network measures
    fprintf('\nSaving data to: %s...\n\n',fn_data)
    save(fn_data,'-struct','data')

end

fprintf('\nFinished calculating network measures!\n\n')

%% Check data
dfn = fieldnames(data);
clc
for i = 1:length(dfn)
    if size(data.(dfn{i}),2) == 1
        fprintf('%s\t%.2f\t%.2f\n',dfn{i},mean(data.(dfn{i})),std(data.(dfn{i})))
    else
        disp(dfn{i})
        disp([mean(data.(dfn{i}))',std(data.(dfn{i}))'])
    end
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

result_fc_patcon = cell(2,5);
for idxm = 1:2
    meas = measures{idxm};
    pat = data.(strcat(groups{1},'_',meas));
    con = data.(strcat(groups{2},'_',meas));
    
    % Use measure of effect size
    stats = mes(pat,con,'cles','nBoot',10000); 
    P = stats.t.p;
    E = stats.cles;
    CI = stats.clesCi;

%     % Use ranksum
%     [P] = ranksum(pat, con); 
  
    result_fc_patcon(idxm,:) = {meas,P,E,CI(1),CI(2)};
    results.(strcat('stats_all_',meas)) = [P,E,CI(1),CI(2)];
end
disp('CLE test results patients vs. controls:')
disp(cell2table(result_fc_patcon,'VariableNames',{'measure','Pval','CLES','CIlo','CIhi'}))


% Compare stroke subgroups and healthy controls

GROUP = cellfun(@(x) repmat({x},1,length(data.(strcat(x,'_fch')))),groups(2:end),'UniformOutput',false);
GROUP = [GROUP{:}];
result_fc_anova = cell(2,4);
for idxm = 1:2
    meas = measures{idxm};
    X = cell2mat(cellfun(@(x) data.(strcat(x,'_',meas)),groups(2:end),'UniformOutput',false)');
    
    % Perform ANOVA   
    [P, TAB, STATS] = anova1(X, GROUP,'off'); 
    result_fc_anova(idxm,1:3) = {meas,P,TAB{2,5}};
    results.(strcat('stats_con_ll_rl_',meas)) = [P,TAB{2,5}];
    
    % Perform pairwise compairsons
    
%     % Use anova multcompare 
%     MC = multcompare(STATS);
%     for comp = 1:size(MC,1)
%         results.(sprintf('stats_%s_%s_%s',groups{MC(comp,1)+1},groups{MC(comp,2)+1},meas)) = MC(comp,end);
%     end
%     if P < 0.05 
%         result_fc_anova(idxm,4) = {MC};
%     end
    
    % Use measure of effect size
    groupcomp = [2,3;2,4;3,4];
    MC = cell(3,6);
    for idxc = 1:size(groupcomp,1)
        g1 = groups{groupcomp(idxc,1)}; g2 = groups{groupcomp(idxc,2)};
        
        stats = mes(data.(strcat(g1,'_',meas)),data.(strcat(g2,'_',meas)),'cles','nBoot',10000);
        results.(sprintf('stats_%s_%s_%s',g1,g2,meas)) = [stats.t.p,stats.cles,stats.clesCi(1),stats.clesCi(2)];
        MC(idxc,:) = {g1,g2,stats.t.p,stats.cles,stats.clesCi(1),stats.clesCi(2)};
    end
    result_fc_anova(idxm,4) = {MC};
end

disp('One-way ANOVA test results patients LL vs. RL vs. CONTROLS:')
disp(cell2table(result_fc_anova(:,1:3),'VariableNames',{'measure','Pval','F'}))

for idxm = 1:2
    if not(isempty(result_fc_anova{idxm,4}))
        disp(['Pairwise comparisons for: ',result_fc_anova{idxm,1}]);
%         disp(array2table(result_fc_anova{idxm,4},'VariableNames',{'G1','G2','Lowerlimit','Difference','Upperlimit','Pvalue'})); 
    disp(array2table(result_fc_anova{idxm,4},'VariableNames',{'G1','G2','Pvalue','CLES','CIlo','CIhi'})); 
    end
end


%% Graph measures


% Compare graph metrics between patients and controls

result_gm_patcon = cell(4,length(sparsities));
for idxm = 3:length(measures)
    meas = measures{idxm};
    pat = data.(strcat(groups{1},'_',meas));
    con = data.(strcat(groups{2},'_',meas));
    
%     % Use ranksum
%     plist = zeros(length(sparsities),1);
%     for idxs = 1:size(pat,2)
%         plist(idxs) = ranksum(pat(:,idxs), con(:,idxs));
%     end
%     result_gm_patcon(idxm-2,:) = {plist};

    % Use CLES
    plist = zeros(length(sparsities),4);
    for idxs = 1:size(pat,2)
        stats = mes(pat(:,idxs),con(:,idxs),'cles','nBoot',10000); 
        plist(idxs,:) = [stats.t.p,stats.cles,stats.clesCi(1),stats.clesCi(2)];
    end
    
    results.(strcat('stats_all_',meas)) = plist;
    disp([meas,' patients vs. controls (per sparsity level):'])
    disp(array2table(plist,'VariableNames',{'Pvalue','CLES','CIlo','CIhi'}, 'RowNames', cellstr(string(sparsities))))
end

% disp('CLES test results patients vs. controls (per sparsity level):')
% disp(array2table(result_gm_patcon','VariableNames',measures(3:end), 'RowNames', cellstr(string(sparsities))))


%% Compare stroke subgroups and healthy controls

GROUP = cellfun(@(x) repmat({x},1,length(data.(strcat(x,'_fch')))),groups(2:end),'UniformOutput',false);
GROUP = [GROUP{:}];
result_gm_anova = cell(4,length(sparsities),2);

for idxm = 3:length(measures)
    meas = measures{idxm};
    plist = zeros(length(sparsities),2);
    mclist = cell(3,length(sparsities),6);
    
    for idxs = 1:length(sparsities)
        X = cell2mat(cellfun(@(x) data.(strcat(x,'_',meas))(:,idxs), groups(2:end),'UniformOutput',false)');
        [P, TAB, STATS] = anova1(X, GROUP,'off'); 
        plist(idxs,:)=[P,TAB{2,5}];
        
%         % Use anova multcompare
%         MC = multcompare(STATS);
%         mclist(idxs,:) = MC(:,end);
%         result_fc_anova(idxm,4) = {MC};
        
        % Use measure of effect size
        groupcomp = [2,3;2,4;3,4];
        MC = cell(3,6);
        for idxc = 1:size(groupcomp,1)
            g1 = groups{groupcomp(idxc,1)}; g2 = groups{groupcomp(idxc,2)};

            stats = mes(data.(strcat(g1,'_',meas))(:,idxs),data.(strcat(g2,'_',meas))(:,idxs),'cles','nBoot',10000);
            mclist(idxc,idxs,:) = {g1,g2,stats.t.p,stats.cles,stats.clesCi(1),stats.clesCi(2)};
            
        end
    end
    
    result_gm_anova(idxm-2,:,:) = num2cell(plist);
    results.(strcat('stats_con_ll_rl_',meas)) = plist;
    for idxc = 1:3
        results.(sprintf('stats_%s_%s_%s',mclist{idxc,1,1},mclist{idxc,1,2},meas)) = cell2mat(squeeze(mclist(idxc,:,3:end)));
    end
    
    disp([meas,' CLES pairwise comparisons (per sparsity level):'])
    for idxc = 1:3
        disp(array2table(squeeze(mclist(idxc,:,:)),'VariableNames',{'G1','G2','Pvalue','CLES','CIlo','CIhi'}, 'RowNames', cellstr(string(sparsities))))
    end
end

disp('One-way ANOVA test results patients LL vs. RL vs. controls (per sparsity level):')
disp('P-values')
disp(array2table(result_gm_anova(:,:,1)','VariableNames',measures(3:end), 'RowNames', cellstr(string(sparsities))))
disp('F-stat')
disp(array2table(result_gm_anova(:,:,2)','VariableNames',measures(3:end), 'RowNames', cellstr(string(sparsities))))


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
disp(array2table(cell2mat(cellfun(@(x) {results.(strcat('stats_lesionsize_',x))},measures(1:2))'),...
    'RowNames',measures(1:2),'VariableNames',{'R','P'}))
disp('R per sparsity level:')
disp(array2table(cell2mat(cellfun(@(x) {results.(strcat('stats_lesionsize_',x))(:,1)}, measures(3:end))),...
    'VariableNames',measures(3:end),'RowNames', cellstr(string(sparsities))))
disp('p per sparsity level:')
disp(array2table(cell2mat(cellfun(@(x) {results.(strcat('stats_lesionsize_',x))(:,2)}, measures(3:end))),...
    'VariableNames',measures(3:end),'RowNames', cellstr(string(sparsities))))


%% Save results

% Save arrays to mat for visualization script
fprintf('Saving results to: %s...\n',fn_stats)
save(fn_stats,'-struct','results')

% Save Table to xlsx
res_gm = {}; res_fc = {};
fnres = fieldnames(results);

for idxr = 1:length(fnres)
    if contains(fnres{idxr},'lesionsize')
        colnames = {'R','P'};
    elseif size(results.(fnres{idxr}),2)==4
        colnames = {'P','CLES','CIlo','CIhi'};
    elseif size(results.(fnres{idxr}),2)==2
        colnames = {'P','F'};
    end
    colnames = cellfun(@(x) [fnres{idxr},'_',x],colnames,'UniformOutput',false);
    
    if size(results.(fnres{idxr}),1)==17
        res_gm(end+1) = {array2table(results.(fnres{idxr}),'VariableNames',colnames,'RowNames',cellstr(string(sparsities)))};
    else
        res_fc(end+1) = {array2table(results.(fnres{idxr}),'VariableNames',colnames)};
    end
    
end

writetable(horzcat(res_gm{:}), [fn_stats(1:end-3),'xlsx'],'WriteRowNames',true,'Sheet','Sheet1')
writetable(horzcat(res_fc{:}), [fn_stats(1:end-3),'xlsx'],'Sheet','Sheet2')

fprintf('\n***** Results saved to %sxlsx! *****\n\n',fn_stats(1:end-3))


