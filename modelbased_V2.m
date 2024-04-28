
%% refrence and citation


%%%%%%%%%%%This code is a modified version of the the available code in Zajkowski et al 2017 article
%%%Zajkowski, Wojciech; Kossut, Malgorzata; Wilson, Robert, 2017, "A causal role for right frontopolar cortex in directed, but not random, exploration", https://doi.org/10.7910/DVN/CZT6EE
%%%, Harvard Dataverse, V2, UNF:6:DJXEJvnH9X1WmZfT3L+j7Q== [fileUNF]%%%



%%properly cited in article%%
%%modified by Armin TOGHI
%% directories
clear
fundir      = pwd;%[maindir 'TMS_code/'];
datadir     = pwd;%[maindir 'TMS_code/'];
savedir     = pwd;%[maindir];
addpath(fundir);
cd(fundir);

defaultPlotParameters
%% load data
sub = load_TMS_v1([datadir '/TMSanalysis_192_data.csv']);

% remove bad subjects
[sub, sub_bad] = removeBadSubjects_TMS_v1(sub);

% Indices of rows to be removed
%rowsToRemove = [1,17, 6];

% Remove specified rows from the structure array
%sub(rowsToRemove) = [];

%% prepare data
clear a
L = unique(sub(1).gameLength);
i = 1;
NS = length(sub);   % number of subjects
T = 4;              % number of forced choices
U = 2;              % number of uncertainty conditions

a = zeros(NS, 384, T);
c5 = nan(NS, 384);
r = zeros(NS, 384, T);
UC = nan(NS, 384);
GL = nan(NS, 384);
dR = zeros(NS, 384);
ds = zeros(NS, 384);

for sn = 1:length(sub)
    
    % choices on forced trials
    dum = sub(sn).a(:,1:4);
    a(sn,1:size(dum,1),:) = dum;
    
    % choices on free trial
    % note a slight hacky feel here - a is 1 or 2, c5 is 0 or 1.
    dum = sub(sn).a(:,5) == 2;
    L(sn) = length(dum);
    c5(sn,1:size(dum,1)) = dum;
    
    % rewards
    dum = sub(sn).r(:,1:4);
    r(sn,1:size(dum,1),:) = dum;
    
    % game length
    dum = sub(sn).gameLength;
    GL(sn,1:size(dum,1)) = dum;
    
    G(sn) = length(dum);
    
    % uncertainty condition 
    dum = abs(sub(sn).uc - 2) + 1;
    UC(sn, 1:size(dum,1)) = dum;

    %spatial bias
    numGames = size(sub(sn).a, 1);
    ds_currentSubject = nan(numGames, 1);
    % Find indices where option b is chosen (on the right) in the 4th trial
    rightIndices = sub(sn).a(:, 4) == 2; % Option b is on the right in the 4th trial
    % Find indices where option a is chosen (on the left) in the 4th trial
    leftIndices = sub(sn).a(:, 4) == 1; % Option a is on the left in the 4th trial
    
    % Assign spatial bias based on the position of option b
    ds_currentSubject(rightIndices) = 1;  
    ds_currentSubject(leftIndices) = -1;  
    
    % Assign the calculated values back to the main 'ds' matrix
    ds(sn, 1:numGames) = ds_currentSubject';

    % difference in information
    dum = sub(sn).uc - 2;
    dI(sn, 1:size(dum,1)) = dum;

    %difference in rewards
    dR(sn, 1:size(dum,1)) = (sub(sn).o1(:,4) - sub(sn).o2(:,4));
    

    % TMS flag
    dum = strcmp(sub(sn).expt_name, 'dlpfc');
    TMS(sn,1:size(dum,1)) = dum;
    
end

dum = GL(:); dum(dum==0) = [];
H = length(unique(dum));
dum = UC(:); dum(dum==0) = [];
U = length(unique(dum));
GL(GL==5) = 1;
GL(GL==10) = 2;

C1 = (GL-1)*2+UC
C2 = TMS + 1;
nC1 = 4;
nC2 = 2;

% meaning of condition 1
% gl uc c1
%  1  1  1 - horizon 1, [1 3]
%  1  2  2 - horizon 1, [2 2]
%  2  1  3 - horizon 6, [1 3]
%  2  2  4 - horizon 6, [2 2]


datastruct = struct(...
    'C1', C1, 'C2', C2, 'nC1', nC1, 'nC2', nC2, ...
    'NS', NS, 'G',  G,  'T',   T, ...
    'dI', dI, 'a',  a,  'c5',  c5, 'r', r, 'dR', dR, 'ds', ds);

%% model-based 

alpha = zeros(NS, nC1, nC2);  % Information bonus
B = zeros(NS, nC1, nC2);      % Spatial bias
sigma = zeros(NS, nC1, nC2);  % Decision noise

% Loop over subjects
for s = 1:NS
    % Loop over levels of Condition 1
    for c1 = 1:nC1
        % Loop over levels of Condition 2
        for c2 = 1:nC2
            % Find indices for the current combination of conditions
            condition_indices = datastruct.C1(s,:) == c1 & datastruct.C2(s,:) == c2;
            
            % Continue only if there are entries for this condition combination
            if any(condition_indices)
                % Extract choices and difference in information for the current conditions
                dR_sub = datastruct.dR(s, condition_indices);
                dI_sub = datastruct.dI(s, condition_indices);
                ds_sub = datastruct.ds(s, condition_indices);
                choices_sub = datastruct.c5(s, condition_indices);
                xfit = model1(dR_sub, dI_sub, ds_sub, choices_sub);
                % Store the estimated parameters for this condition combination
                alpha(s, c1, c2) = xfit.data(1);
                B(s, c1, c2) = xfit.data(2);
                sigma(s, c1, c2) = xfit.data(3);
            end
        end
    end
end


%% Figure 6: plotting decision noise

[h1, p1] = ttest(sigma(:, 1, 1), sigma(:, 1, 2))
[h1, p1] = ttest(sigma(:, 2, 1), sigma(:, 2, 2))
[h6, p6] = ttest(sigma(:, 3, 1), sigma(:, 3, 2))

means = [mean(sigma(:, 1, 1), 'omitnan'), mean(sigma(:, 1, 2), 'omitnan'); % Horizon 1 for Vertex and DLPFC
         mean(sigma(:, 2, 1), 'omitnan'), mean(sigma(:, 2, 2), 'omitnan');
         mean(sigma(:, 3, 1), 'omitnan'), mean(sigma(:, 3, 2), 'omitnan');
         mean(sigma(:, 4, 1), 'omitnan'), mean(sigma(:, 4, 2), 'omitnan')];

% Calculate standard errors (or standard deviations)
sems = [std(sigma(:, 1, 1), 'omitnan') / sqrt(size(sigma(:, 1, 1), 1)), std(sigma(:, 1, 2), 'omitnan') / sqrt(size(sigma(:, 1, 2), 1));
        std(sigma(:, 2, 1), 'omitnan') / sqrt(size(sigma(:, 2, 1), 1)), std(sigma(:, 2, 2), 'omitnan') / sqrt(size(sigma(:, 2, 2), 1));
        std(sigma(:, 3, 1), 'omitnan') / sqrt(size(sigma(:, 3, 1), 1)), std(sigma(:, 3, 2), 'omitnan') / sqrt(size(sigma(:, 3, 2), 1));
        std(sigma(:, 4, 1), 'omitnan') / sqrt(size(sigma(:, 4, 1), 1)), std(sigma(:, 3, 2), 'omitnan') / sqrt(size(sigma(:, 4, 2), 1))];


figure(1); clf;

b = bar(means, 'grouped'); hold on;

ngroups = size(means, 1);
nbars = size(means, 2);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the error bars
errorbar(x', means, sems, 'k', 'linestyle', 'none');

% Customizing the plot with colors
b(1).FaceColor = 'flat';
b(1).CData = repmat(arminn, [ngroups, 1]); % Custom color for all Horizon 1 bars
b(2).FaceColor = 'flat';
b(2).CData = repmat(arminn_pink, [ngroups, 1]); % Custom color for all Horizon 6 bars


condition1Names = {
    'horizon 1, [2 2]', 
    'horizon 1, [1 3]',
    'horizon 6, [2 2]', 
    'horizon 6, [1 3]', 
};

% Labels and legends
set(gca, 'XTickLabel', condition1Names);
ylabel('Decision noise (\sigma)');
legend({'vertex', 'dlpfc'}, 'Location', 'Best');
 



%% Supplementory material S2

[h1, p1] = ttest(alpha(:, 2, 1), alpha(:, 2, 2))
[h6, p6] = ttest(alpha(:, 4, 1), alpha(:, 4, 2))

% Calculate means
means = [mean(alpha(:, 2, 1), 'omitnan'), mean(alpha(:, 2, 2), 'omitnan'); % Horizon 1 for Vertex and DLPFC
         mean(alpha(:, 4, 1), 'omitnan'), mean(alpha(:, 4, 2), 'omitnan')]; % Horizon 6 for Vertex and DLPFC

% Calculate standard errors (or standard deviations)
sems = [std(alpha(:, 2, 1), 'omitnan') / sqrt(size(alpha(:, 2, 1), 1)), std(alpha(:, 2, 2), 'omitnan') / sqrt(size(alpha(:, 2, 2), 1));
        std(alpha(:, 4, 1), 'omitnan') / sqrt(size(alpha(:, 4, 1), 1)), std(alpha(:, 4, 2), 'omitnan') / sqrt(size(alpha(:, 4, 2), 1))];


figure(1); clf;
%set(gcf, 'position', [211   137   700   650])
%set(gcf, 'position', [588 576 600 400]);
%ax = easy_gridOfEqualFigures([0.2 0.1], [0.2 0.2 0.1]);
%ax = easy_gridOfEqualFigures([0.2 0.1], [0.2 0.03]);
%axes(ax(1)); hold on;

b = bar(means, 'grouped'); hold on;
% Calculate the number of groups and number of bars in each group
ngroups = size(means, 1);
nbars = size(means, 2);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the error bars
errorbar(x', means, sems, 'k', 'linestyle', 'none');

% Customizing the plot with colors
b(1).FaceColor = 'flat';
b(1).CData = repmat(arminn, [ngroups, 1]); % Custom color for all Horizon 1 bars
b(2).FaceColor = 'flat';
b(2).CData = repmat(arminn_pink, [ngroups, 1]); % Custom color for all Horizon 6 bars

% Labels and legends
set(gca, 'XTickLabel', {'Horizon 1', 'Horizon 6'});
ylabel('Information Bonus (\alpha)');
legend({'vertex', 'dlpfc'}, 'Location', 'Best');
%title('Information Bonus');
maxYValues = max(means + sems, [], 2); % Max Y-value for each group of bars, for placing stars
significanceLevels = [p1, p6]; % Array of p-values


% Assuming maxYValues is an array with max Y-values for each comparison
%for i = 1:2
%    yPos = maxYValues(i) * 1.1; % Adjust Y position for annotation
%    if significanceLevels(i) < 0.05
%        plot(i + [0.5, 1.5], [1, 1] * maxYValues(i), '-k', 'LineWidth', 1.5); % Adjust line position
%        text(i + 1, yPos, '*', 'FontSize', 14, 'HorizontalAlignment', 'center'); % Asterisk for significance
%    else
%end

yPos2 = 12;
text(2, yPos2, 'NS', 'FontSize', 14, 'HorizontalAlignment', 'center');
yPos3 =  -5;
text(1, yPos3, 'NS', 'FontSize', 14, 'HorizontalAlignment', 'center');

%% Supplementory material S3

[h1, p1] = ttest(B(:, 1, 1), B(:, 1, 2))
[h1, p1] = ttest(B(:, 2, 1), B(:, 2, 2))
[h6, p6] = ttest(B(:, 3, 1), B(:, 3, 2))
[h1, p1] = ttest(B(:, 4, 1), B(:, 4, 2))
for c1 = 1:length(condition1Names)
    BVertex = B(:, c1, 1);
    BDLPFC = B(:, c1, 2);
    
    % Perform paired t-test
    [h, p] = ttest(BVertex, BDLPFC);
    pValues(c1) = p; % Store p-value
    
    % Plot with custom colors
    subplot(1, length(condition1Names), c1);
    b = bar([1, 2], [mean(BVertex, 'omitnan'), mean(BDLPFC, 'omitnan')]); % Mean sigma for vertex and dlpfc
    
    % Apply custom colors
    b.FaceColor = 'flat';
    b.CData(1,:) = arminn; % Custom color for vertex
    b.CData(2,:) = arminn_pink; % Custom color for dlpfc
    
    set(gca, 'XTickLabel', {'Vertex', 'DLPFC'});
    title(condition1Names{c1});
    ylabel('Spatial Bias (B)');
    
    % Annotate significant differences
    if p < 0.05
        hold on;
        maxY = max([mean(BVertex, 'omitnan'), mean(BDLPFC, 'omitnan')]) * 1.1;
        plot([1, 2], [1, 1] * maxY, '-k', 'LineWidth', 1.5);
        text(1.5, maxY, '*', 'FontSize', 14, 'HorizontalAlignment', 'center'); % Asterisk for significance
        ylim([0, maxY * 1.2]);
    end
end


%% normality check

% Assuming sigma is a 3D matrix with dimensions: subjects x condition1 x condition2
% nC1 represents the number of levels in Condition 1

% Preallocate array for p-values from the Shapiro-Wilk test
swPValues = zeros(nC1, 1);

for c1 = 1:nC1
    % Calculate differences between sigma values for "vertex" and "dlpfc"
    diffSigma = sigma(:, c1, 1) - sigma(:, c1, 2);
    
    % Remove NaN values
    diffSigmaClean = diffSigma(~isnan(diffSigma));
    
    % Perform the Lilliefors test for normality on the differences
    [h, p] = lillietest(diffSigmaClean);
    
    % Display results
    fprintf('Condition 1: %s, Lilliefors test p-value: %.4f\n', condition1Names{c1}, p);
end
%% Figure 7: correlation between p(high_info) and alpha_ p(low_mean) and noisy

clear p_highInfo p_lowMean p_repeat p_repeat13
GL = [5 10];
for sn = 1:length(sub)
    for h = 1:length(GL)
        i_RFPC = strcmp(sub(sn).expt_name, 'dlpfc');
        i_vertex = strcmp(sub(sn).expt_name, 'vertex');
        
        
        idx = i_vertex;
        ind = (sub(sn).gameLength == GL(h)) & idx;
        p_highInfo(h,sn,1) = nanmean(sub(sn).hi(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)~=sub(sn).n1(:,4)) & idx;
        p_repeat13(h, sn,1) = nanmean(sub(sn).rep(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)==sub(sn).n1(:,4)) & idx;
        p_lowMean(h,sn,1) = nanmean(sub(sn).lm(ind,5));
        p_repeat(h, sn,1) = nanmean(sub(sn).rep(ind,5));
        p_right(h,sn,1) = nanmean(sub(sn).a(ind,5)==2);
        
        
        idx = i_RFPC;
        ind = (sub(sn).gameLength == GL(h)) & idx;
        p_highInfo(h,sn,2) = nanmean(sub(sn).hi(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)~=sub(sn).n1(:,4)) & idx;
        p_repeat13(h, sn,2) = nanmean(sub(sn).rep(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)==sub(sn).n1(:,4)) & idx;
        p_lowMean(h,sn,2) = nanmean(sub(sn).lm(ind,5));
        p_repeat(h, sn,2) = nanmean(sub(sn).rep(ind,5));
        p_right(h,sn,2) = nanmean(sub(sn).a(ind,5)==2);
        
    end
end

%p_lowmean(2*18*2)
%sigma(18*4*2)




% gl uc c1
%  1  1  1 - horizon 1, [2 2]
%  1  2  2 - horizon 1, [1 3]
%  2  1  3 - horizon 6, [2 2]
%  2  2  4 - horizon 6, [1 3]

figure(4);


ax(1) = subplot(2, 2, 1)
x_randomfree_5 = reshape(p_lowMean(1,:,:), [], 1);
y_randommodel_5 = reshape(sigma(:, 1, :), [], 1);
scatter(x_randomfree_5, y_randommodel_5, 100, arminn_pink, '.'); hold on 
p = polyfit(x_randomfree_5, y_randommodel_5,1)
f = polyval(p, x_randomfree_5)
plot(x_randomfree_5, f, 'color', arminn, 'LineWidth',3)
ylabel('noise Horizon 1')
xlabel('p(low mean) Horizon 1')
%ylim([0,0.1])
set(gca,'FontSize',14)
t = title('Random exploration', 'FontSize', 16)
pos = get(t, 'Position');
pos(2) = pos(2);
set(t, 'Position', pos);
[r, p] = corr(x_randomfree_5, y_randommodel_5);
text(.21,20,strcat('r =', {' '}, num2str(round(r*100)/100),' ; p < 0.001'),'FontWeight','Normal', 'FontSize',15)

ax(2) = subplot(2,2,3)

x_randomfree_10 = reshape(p_lowMean(2,:,:), [], 1);
y_randommodel_10 = reshape(sigma(:, 3, :), [], 1);
scatter(x_randomfree_10, y_randommodel_10, 100, arminn_pink, '.'); hold on 
p = polyfit(x_randomfree_10, y_randommodel_10,1)
f = polyval(p, x_randomfree_10)
plot(x_randomfree_10, f, 'color', arminn, 'LineWidth',3)
ylabel('noise Horizon 6')
xlabel('p(low mean) Horizon 6')
% ylim([0,1.1])
set(gca,'FontSize',14)
%t = title('Random exploration', 'FontSize', 16)
pos = get(t, 'Position');
pos(2) = pos(2) + 2;
set(t, 'Position', pos);
[r, p] = corr(x_randomfree_10, y_randommodel_10);
text(.21,23,strcat('r =', {' '}, num2str(round(r*100)/100),' ; p < 0.001'),'FontWeight','Normal', 'FontSize',15)


ax(3) = subplot(2,2,2)

x_directfree_1 = reshape(p_highInfo(1,:,:), [], 1);
y_directmodel_1 = reshape(alpha(:, 2, :), [], 1);
scatter(x_directfree_1, y_directmodel_1, 100, arminn_pink, '.'); hold on 
p = polyfit(x_directfree_1, y_directmodel_1,1)
f = polyval(p, x_directfree_1)
plot(x_directfree_1, f, 'color', arminn, 'LineWidth',3)
ylabel('info bonus Horizon 1')
xlabel('p(High Info) Horizon 1')
% ylim([0,1.1])
set(gca,'FontSize',14)
t = title('Direct exploration', 'FontSize', 16)
pos = get(t, 'Position');
pos(2) = pos(2) + 2;
set(t, 'Position', pos);
[r, p] = corr(x_directfree_1, y_directmodel_1);
text(.21,23,strcat('r =', {' '}, num2str(round(r*100)/100),' ; p < 0.001'),'FontWeight','Normal', 'FontSize',15)

ax(4) = subplot(2,2,4)

x_directfree_6 = reshape(p_highInfo(2,:,:), [], 1);
y_directmodel_6 = reshape(alpha(:, 4, :), [], 1);
scatter(x_directfree_6, y_directmodel_6, 100, arminn_pink, '.'); hold on 
p = polyfit(x_directfree_6, y_directmodel_6,1)
f = polyval(p, x_directfree_6)
plot(x_directfree_6, f, 'color', arminn, 'LineWidth',3)
ylabel('info bonus Horizon 6')
xlabel('p(High Info) Horizon 6')
% ylim([0,1.1])
set(gca,'FontSize',14)
pos = get(t, 'Position');
pos(2) = pos(2) + 2;
set(t, 'Position', pos);
[r, p] = corr(x_directfree_6, y_directmodel_6);
text(.21,23,strcat('r =', {' '}, num2str(round(r*100)/100),' ; p < 0.001'),'FontWeight','Normal', 'FontSize',15)

addABCs(ax, [-0.05 0.09], 32)


%% functions

function fit = model1(dR, dI, ds, c5)
    X0 = [0 0 10];           % information bonus, spatial bias, decision noise
    LB = [-100 -100 0];
    UB = [100 100 100];
    
    func = @(x) likelihood(dR, dI, c5, ds, x);
    data = fmincon(func, X0, [], [], [], [], LB, UB);
    
    fit.data(1) = data(1);
    fit.data(2) = data(2);
    fit.data(3) = data(3);
end   


function L = likelihood(dR, dI, c5, ds, x)
    alpha = x(1);
    B = x(2);
    sigma = x(3);
    
    dQ = dR + alpha*dI + B* ds;
    P = 1 ./ (1+exp(dQ/(sigma*sqrt(2))));
    
    logP1 = log(P(c5==1));
    logP2 = log(1-P(c5==0));
    
    L_A = gaussianPrior(alpha, 0, 20);
    L_noise = exponentialPrior(sigma, 1/20);
    
    L = -sum(logP1) - sum(logP2) - L_A - L_noise;
end


function logP = exponentialPrior(x, lambda)
    logP = -lambda * x;
end

function logP = gaussianPrior(x, mu, sigma)
    logP = -(x-mu).^2/sigma^2/2 - log(sqrt(2*pi)*sigma);
end

