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

% Indices of rows to be removed
%rowsToRemove = [1,17, 7];

% Remove specified rows from the structure array
%sub(rowsToRemove) = [];
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
%% post hoc

%posthocResults_rand = multcompare(rm_rand, 'TMScondition', 'by','Order', 'ComparisonType', 'bonferroni');
%disp(posthocResults_rand);

posthocResults_rand = multcompare(rm_rand, 'TMScondition','by', 'Horizon', 'ComparisonType', 'bonferroni');
disp(posthocResults_rand);

%posthocResults_rand = multcompare(rm_rand, 'Horizon', 'ComparisonType', 'bonferroni');
%disp(posthocResults_rand);


