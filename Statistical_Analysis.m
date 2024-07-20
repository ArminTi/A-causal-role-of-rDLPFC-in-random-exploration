%% ============== Statistical Analaysis ===================== %%
%% load data and directories
clear
fundir      = pwd;%[maindir 'TMS_code/'];
datadir     = pwd;%[maindir 'TMS_code/'];
savedir     = pwd;%[maindir];
addpath(fundir);
cd(fundir);

defaultPlotParameters

sub = load_TMS_v1([datadir '/TMSanalysis_192_data']);

% remove bad subjects
[sub, sub_bad] = removeBadSubjects_TMS_v1(sub);


%% creating Table 

T = table([], [], [], [], [], [], [], 'VariableNames', ...
    {'SubjectID', 'TMS_Condition', 'Order', 'Gender', 'Horizon', 'DirectExploration', 'RandomExploration'});
for i = 1:length(sub)
    for cond = {'vertex', 'dlpfc'}
        % Find indices for this condition
        idx_cond = strcmp(sub(i).expt_name, cond);

        % Loop over the two game lengths: 5 and 10
        for gLength = [5, 10]
            % Find indices for games of this length
            idx_gLength = sub(i).gameLength == gLength;

            % Combine condition and game length indices
            idx_combined = idx_cond & idx_gLength & (sub(i).n2(:,4)~=sub(i).n1(:,4));
            mean_hi = nanmean(sub(i).hi(idx_combined, 5));
            idx_combined = idx_cond & idx_gLength & (sub(i).n2(:,4)==sub(i).n1(:,4));
            mean_low = nanmean(sub(i).lm(idx_combined, 5));

            % Add a new row to the table for this subject, condition, and game length
            newRow = {sub(i).subjectID, cond{1}, sub(i).order, sub(i).iswoman, gLength, mean_hi, mean_low};
            T = [T; newRow];
        end
    end
end


%% creating table

T_random = removevars(T, "DirectExploration");
T_random = unstack(T_random,'RandomExploration','TMS_Condition');
T_random = unstack(T_random,{'dlpfc','vertex'},'Horizon');


T_direct = removevars(T, "RandomExploration");
T_direct = unstack(T_direct,'DirectExploration','TMS_Condition');
T_direct = unstack(T_direct,{'dlpfc','vertex'},'Horizon');
%% performing repeated measure ANOVA


T_random.Order = categorical(T_random.Order);
T_random.Gender = categorical(T_random.Gender);

WithinStructure = table([1 1 2 2]',[1 2 1 2]','VariableNames',{'TMScondition','Horizon'});
WithinStructure.TMScondition = categorical(WithinStructure.TMScondition);
WithinStructure.Horizon = categorical(WithinStructure.Horizon);
rm_rand = fitrm(T_random, 'vertex_x5,vertex_x10,dlpfc_x5,dlpfc_x10 ~ Gender + Order','WithinDesign',WithinStructure);
ranovatablerand = ranova(rm_rand,'WithinModel','TMScondition*Horizon');
disp(ranovatablerand)

rm_dir = fitrm(T_direct, 'dlpfc_x5,dlpfc_x10,vertex_x5,vertex_x10~ Gender + Order','WithinDesign',WithinStructure);
ranovatabledir = ranova(rm_dir,'WithinModel','TMScondition*Horizon');
disp(ranovatabledir)


%% Calculate Partial Eta Squared for Random Exploration ANOVA
% Sum of Squares for each effect
% Within-subject factors
SS_effect_rand = ranovatablerand.SumSq([5, 9, 13]); 
SS_error_rand = ranovatablerand.SumSq([8, 12, 16]);
% Total Sum of Squares
SS_total_rand = sum(SS_effect_rand) + sum(SS_error_rand);
% Calculate partial eta squared
eta_squared_rand = SS_effect_rand / SS_total_rand;

eta_squared_rand_table = array2table(eta_squared_rand, 'VariableNames', {'PartialEtaSquared'}, 'RowNames', {'TMScondition', 'Horizon', 'TMScondition:Horizon'});
disp('Partial Eta Squared for Random Exploration ANOVA:')
disp(eta_squared_rand_table)

%% Calculate Partial Eta Squared for Direct Exploration ANOVA
SS_effect_dir = ranovatabledir.SumSq([5, 9, 13]); % Sum of Squares for each effect
SS_error_dir = ranovatabledir.SumSq([8, 12, 16]);

SS_total_dir = sum(SS_effect_dir) + sum(SS_error_dir);
% Calculate partial eta squared
eta_squared_dir = SS_effect_dir / SS_total_dir;


eta_squared_dir_table = array2table(eta_squared_dir, 'VariableNames', {'PartialEtaSquared'}, 'RowNames', {'TMScondition', 'Horizon', 'TMScondition:Horizon'});
disp('Partial Eta Squared for Direct Exploration ANOVA:')
disp(eta_squared_dir_table)

%% post hoc with bonferroni

%posthocResults_rand = multcompare(rm_rand, 'TMScondition', 'by','Order', 'ComparisonType', 'bonferroni');
%disp(posthocResults_rand);

posthocResults_rand = multcompare(rm_rand, 'TMScondition','by', 'Horizon', 'ComparisonType', 'bonferroni');
disp(posthocResults_rand);

%posthocResults_rand = multcompare(rm_rand, 'Horizon', 'ComparisonType', 'bonferroni');
%disp(posthocResults_rand);



%% paired t test and cohens d for effect size
clear p_highInfo p_lowMean 
GL = [5 10];
for sn = 1:length(sub)
    for h = 1:length(GL)
        i_DLPFC = strcmp(sub(sn).expt_name, 'dlpfc');
        i_vertex = strcmp(sub(sn).expt_name, 'vertex');
        
        
        idx = i_vertex;
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)~=sub(sn).n1(:,4)) & idx;
        p_highInfo(h,sn,1) = nanmean(sub(sn).hi(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)==sub(sn).n1(:,4)) & idx;
        p_lowMean(h,sn,1) = nanmean(sub(sn).lm(ind,5));
        p_right(h,sn,1) = nanmean(sub(sn).a(ind,5)==2);
        
        
        
        idx = i_DLPFC;
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)~=sub(sn).n1(:,4)) & idx;
        p_highInfo(h,sn,2) = nanmean(sub(sn).hi(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)==sub(sn).n1(:,4)) & idx;
        p_lowMean(h,sn,2) = nanmean(sub(sn).lm(ind,5));
        p_right(h,sn,2) = nanmean(sub(sn).a(ind,5)==2);
        
    end
end


%Direct Exploration between horizon
%[~,p] = ttest(squeeze(p_highInfo(2,:,1)-p_highInfo(1,:,1))')
%[~,p] = ttest(squeeze(p_highInfo(2,:,2)-p_highInfo(1,:,2))')
%Random Exploration between horizon
%[~,p] = ttest(squeeze(p_lowMean(2,:,1)-p_lowMean(1,:,1))')
%[~,p] = ttest(squeeze(p_lowMean(2,:,2)-p_lowMean(1,:,2))')



% Random Exploration between TMS condition
diff_lowMean_5 = squeeze(p_lowMean(1,:,1) - p_lowMean(1,:,2));
diff_lowMean_10 = squeeze(p_lowMean(2,:,1) - p_lowMean(2,:,2));

% Paired t-tests
[~, ph1rand] = ttest(diff_lowMean_5);
[~, ph6rand] = ttest(diff_lowMean_10);

% Calculate Cohen's d for Random Exploration
d_lowMean_5 = mean(diff_lowMean_5) / std(diff_lowMean_5);
d_lowMean_10 = mean(diff_lowMean_10) / std(diff_lowMean_10);

fprintf('Cohen''s d for Random Exploration (5 trials): %.4f\n', d_lowMean_5);
fprintf('Cohen''s d for Random Exploration (10 trials): %.4f\n', d_lowMean_10);

% Direct Exploration between TMS condition
diff_highInfo_5 = squeeze(p_highInfo(1,:,1) - p_highInfo(1,:,2));
diff_highInfo_10 = squeeze(p_highInfo(2,:,1) - p_highInfo(2,:,2));

% Paired t-tests
[~, ph1dr] = ttest(diff_highInfo_5);
[~, ph6dr] = ttest(diff_highInfo_10);

% Calculate Cohen's d for Direct Exploration
d_highInfo_5 = mean(diff_highInfo_5) / std(diff_highInfo_5);
d_highInfo_10 = mean(diff_highInfo_10) / std(diff_highInfo_10);

fprintf('Cohen''s d for Direct Exploration (5 trials): %.4f\n', d_highInfo_5);
fprintf('Cohen''s d for Direct Exploration (10 trials): %.4f\n', d_highInfo_10);
%% mean and sd for each condition (information in Table 1)

Dir_h1_m_v = mean(p_highInfo(1,:,1));
Dir_h1_s_v = std(p_highInfo(1,:,1));
Dir_h1_m_d = mean(p_highInfo(1,:,2));
Dir_h1_s_d = std(p_highInfo(1,:,2));
Dir_h6_m_v = mean(p_highInfo(2,:,1));
Dir_h6_s_v = std(p_highInfo(2,:,1));
Dir_h6_m_d = mean(p_highInfo(2,:,2));
Dir_h6_s_d = std(p_highInfo(2,:,2));


rand_h1_m_v = mean(p_lowMean(1,:,1));
rand_h1_s_v = std(p_lowMean(1,:,1));
rand_h1_m_d = mean(p_lowMean(1,:,2));
rand_h1_s_d = std(p_lowMean(1,:,2));
rand_h6_m_v = mean(p_lowMean(2,:,1));
rand_h6_s_v = std(p_lowMean(2,:,1));
rand_h6_m_d = mean(p_lowMean(2,:,2));
rand_h6_s_d = std(p_lowMean(2,:,2));