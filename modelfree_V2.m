%% ============== Model-free analysis ===================== %%
%% refrence and citation

%%%%Project name : A casual role of right dorsolateral prefrontal cortex in
%%%%random exploration


%%%%%%%%%%% This code were modified version of the the available code in Zajkowski et al 2017 article
%%%Zajkowski, Wojciech; Kossut, Malgorzata; Wilson, Robert, 2017, "A causal role for right frontopolar cortex in directed, but not random, exploration", https://doi.org/10.7910/DVN/CZT6EE
%%%, Harvard Dataverse, V2, UNF:6:DJXEJvnH9X1WmZfT3L+j7Q== [fileUNF]%%%

%%properly cited in article%%
%%%%modified by Armin T
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
%% Figure 3: Model-Free Analysis of first Free-Choice Trials. 

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





 

clear X l l2 t m s e
figure(1); clf;
set(gcf, 'position', [211   137   600   250])
dw = 0.02;
DW = 0.12;
ax = easy_gridOfEqualFigures([0.2  0.2], [0.14 0.17 0.05]);
i = 0;
i=i+1; X(:,i) = p_highInfo(1,:,1); 
i=i+1; X(:,i) = p_highInfo(1,:,2); 
i=i+1; X(:,i) = p_highInfo(2,:,1); vName{i} = 'p(high info), vertex';
i=i+1; X(:,i) = p_highInfo(2,:,2); vName{i} = 'p(high info), rDLPFC';

i=i+1; X(:,i) = p_lowMean(1,:,1); 
i=i+1; X(:,i) = p_lowMean(1,:,2); vName{i} = 'p(low mean), rDLPFC, horizon 1';
i=i+1; X(:,i) = p_lowMean(2,:,1); vName{i} = 'p(low mean), vertex, horizon 6';
i=i+1; X(:,i) = p_lowMean(2,:,2); vName{i} = 'p(low mean), rDLPFC, horizon 6';


clear vName
vName{1} = 'vertex';
vName{2} = 'rDLPFC';
vName{3} = 'vertex';
vName{4} = 'rDLPFC';
vName{5} = 'vertex';
vName{6} = 'rDLPFC';
axes(ax(1)); hold on;

dx = 0;
xx = 1:size(X,2);
xx(1:2:end) = xx(1:2:end)+dx;
xx(2:2:end) = xx(2:2:end)-dx;

set(ax, 'xlim', [0.5 2.5])
axes(ax(1)); hold on;
clear m s
m(1,1) = nanmean(X(:,1)); s(1,1) = nanstd(X(:,1))/sqrt(size(X,1));
m(2,1) = nanmean(X(:,2)); s(2,1) = nanstd(X(:,2))/sqrt(size(X,1));
m(1,2) = nanmean(X(:,3)); s(1,2) = nanstd(X(:,3))/sqrt(size(X,1));
m(2,2) = nanmean(X(:,4)); s(2,2) = nanstd(X(:,4))/sqrt(size(X,1));
e = errorbar(m, s)
leg = legend(e([2 1]), {'horizon 6' 'horizon 1'}, ...
    'orientation', 'vertical', 'location', 'southeast')

ylabel('p(high info)')
xlabel('stimulation condition')
t = text(0, 0, 'Directed exploration');


tt = text(1.5, 0.565, '')
tt(2) = text(0.65, (m(1,1)+m(1,2))/2+0.01, '*');
tt(3) = text(2.3, (m(2,1)+m(2,2))/2+0.004, '*');
ll = plot([0.84 0.72 0.72 0.84]+0.05, [m(1,1) m(1,1) m(1,2) m(1,2)]);
ll(2) = plot([2.2 2.32 2.32 2.2]-0.1, [m(2,1) m(2,1) m(2,2) m(2,2)]);


axes(ax(2)); hold on;
clear m s
m(1,1) = nanmean(X(:,5)); s(1,1) = nanstd(X(:,5))/sqrt(size(X,1));
m(2,1) = nanmean(X(:,6)); s(2,1) = nanstd(X(:,6))/sqrt(size(X,1));
m(1,2) = nanmean(X(:,7)); s(1,2) = nanstd(X(:,7))/sqrt(size(X,1));
m(2,2) = nanmean(X(:,8)); s(2,2) = nanstd(X(:,8))/sqrt(size(X,1));
e(2,:) = errorbar(m, s)
ll(3) = plot([0.84 0.72 0.72 0.84]+0.05, [m(1,1) m(1,1) m(1,2) m(1,2)]);
ll(4) = plot([2.2 2.32 2.32 2.2]-0.1, [m(2,1) m(2,1) m(2,2) m(2,2)]);
tt(4) = text(0.65, (m(1,1)+m(1,2))/2+0.01, '*');
tt(5) = text(2.4, (m(2,1)+m(2,2))/2+0.01, '**');
tt(6) = text(1.5, 0.20, '**')
tt(7) = text(1.5, 0.29, '*')

ylabel('p(low mean)')
xlabel('stimulation condition')
t(2) = text(0, 0, 'Random exploration');

set(t, 'fontsize', 18, 'units', 'normalized', 'position', [0.5 1.15], ...
    'horizontalAlignment', 'center', 'fontweight', 'bold')
set(ax(1), 'xticklabel', {vName{1:2}} ,'ylim', [0.4 0.6])
set(ax(2), 'xticklabel', {vName{3:4}},'ylim', [0.1 0.301])
set(ax, 'view', [0 90], 'xtick', xx, ...
    'tickdir', 'out', 'ytick', [0:0.1:1]);
set(e, 'marker', '.', 'markersize', 30, 'linestyle', '-')
set(e(:,1), 'color', arminn_pink)
set(e(:,2), 'color', arminn)
f = 0.6;

set(tt, 'horizontalalignment', 'center', 'fontsize', 20, ...
    'verticalalignment', 'top')
set(ll, 'color', 'k', 'linewidth', 1)

addABCs(ax, [-0.09 0.07], 32)

%data analysis only t test

%Direct Exploration between horizon
[~,p] = ttest(squeeze(p_highInfo(2,:,1)-p_highInfo(1,:,1))')
[~,p] = ttest(squeeze(p_highInfo(2,:,2)-p_highInfo(1,:,2))')

%Direct Exploration between TMS condition

[~,p] = ttest(squeeze(p_highInfo(1,:,1)-p_highInfo(1,:,2))')
[~,p] = ttest(squeeze(p_highInfo(2,:,1)-p_highInfo(2,:,2))')

%Random Exploration between horizon

[~,p] = ttest(squeeze(p_lowMean(2,:,1)-p_lowMean(1,:,1))')
[~,p] = ttest(squeeze(p_lowMean(2,:,2)-p_lowMean(1,:,2))')

%Random Exploration between TMS condition

[~,p] = ttest(squeeze(p_lowMean(1,:,1)-p_lowMean(1,:,2))')
[~,p] = ttest(squeeze(p_lowMean(2,:,1)-p_lowMean(2,:,2))')
%% Supplementory materials 1: correlation between reward and information

clear rr1 rr6
for sn = 1:length(sub)
    
    dn = sign(sub(sn).n2 - sub(sn).n1);
    dm = sign(sub(sn).o2 - sub(sn).o1);
    gl = sub(sn).gameLength;
    i1 = gl == 5;
    i6 = gl == 10;
    iv = strcmp(sub(sn).expt_name, 'vertex'); 
    ir = strcmp(sub(sn).expt_name, 'dlpfc');
    
    for i = 4:10
        if sum(iv) == 0
            rr1(i-3,1,sn) = nan;
            rr6(i-3,1,sn) = nan;
        else
            rr1(i-3,1,sn) = corr(dn(i1&iv,i), dm(i1&iv,i));
            rr6(i-3,1,sn) = corr(dn(i6&iv,i), dm(i6&iv,i));
        end
        rr1(i-3,2,sn) = corr(dn(i1&ir,i), dm(i1&ir,i));
        rr6(i-3,2,sn) = corr(dn(i6&ir,i), dm(i6&ir,i));
    end
    
end

rr1(2,:) = nan;
rr6(end,:) = nan;

figure(1); clf;
set(gcf, 'position', [588   576   400   250]);

ax = easy_gridOfEqualFigures([0.2 0.1], [0.2 0.03]);
axes(ax(1)); hold on;
plot([0 7], [0 0], 'k--', 'linewidth', 1)

M1 = nanmean(rr1,3);
S1 = nanstd(rr1,[],3)/sqrt(size(rr1,3));
M6 = nanmean(rr6,3);
S6 = nanstd(rr6,[],3)/sqrt(size(rr6,3));
e = plot(M1);
e(:,2) = plot(M6);
e = e';
set(e(1,:), 'color', arminn);
set(e(2,:), 'color', arminn_pink);
set(e(:,1), 'linestyle', '-','marker', '+')
set(e(:,2), 'linestyle', '--','marker', 'x')
set(e, 'markersize', 10, 'linewidth', 1)

set(ax, 'xtick', [1:6], 'xlim', [0.5 6.5], 'tickdir', 'out')
xlabel('free-trial number')
ylabel({'correlation between' 'reward and information'})
legend(e([1 3 2 4]), ...
    {'vertex, horizon 1' 'rDLPFC, horizon 1', 'vertex, horizon 6' 'rDLPFC, horizon 6'}, ...
    'location', 'southeast')


%% Figure 5: correct responses in whole trials

clear correct1 correct6
for sn = 1:length(sub)
    
    dum = sub(sn).co(:,:);
    gl = sub(sn).gameLength;
    i1 = gl == 5;
    i6 = gl == 10;
    iv = strcmp(sub(sn).expt_name, 'vertex');
    ir = strcmp(sub(sn).expt_name, 'dlpfc');

   
    for i =5:10
        correct1(i-4,1,sn) = nanmean(dum(i1&iv,i));
        correct6(i-4,1,sn) = nanmean(dum(i6&iv,i));
        correct1(i-4,2,sn) = nanmean(dum(i1&ir,i));
        correct6(i-4,2,sn) = nanmean(dum(i6&ir,i));
    end
    
end
correct1(2,:) = nan;

figure(1); clf;

set(gcf, 'position', [588 576 600 300]);
ax = easy_gridOfEqualFigures([0.2 0.1], [0.2 0.2 0.1]);


axes(ax(1)); hold on;
plot([0 0], [0 0], 'k--', 'linewidth', 1)

M1 = nanmean(correct1,3);
S1 = nanstd(correct1,[],3)/sqrt(size(correct1,3));
M6 = nanmean(correct6,3);
S6 = nanstd(correct6,[],3)/sqrt(size(correct6,3));


e = plot(M1);
e(:,2) = plot(M6);

e = e';
set(e(1,:), 'color', arminn_pink);
set(e(:,1), 'linestyle', '-','marker', '+')
%set(e, 'markersize', 10, 'linewidth', 1)
%set(ax, 'xtick', [1:6], 'xlim', [0.5 8], 'tickdir', 'out')
set(e(2,:), 'color', arminn);
set(e(:,2), 'linestyle', '--','marker', 'x')
set(e, 'markersize', 10, 'linewidth', 1)

xlabel('free-trial number')
ylabel({'Fraction of correct'})
l1 = legend(e([1 3 2 4]), ...
    {'vertex, horizon 1', 'rDLPFC, horizon 1', 'vertex, horizon 6', 'rDLPFC, horizon 6'}, ...
    'location', 'northwest'); % This will place the legend outside the plot on the right
set(l1, 'FontSize', 8); % Adjust font size as needed


set(ax(1), 'xtick', [1:6], 'xlim', [0.5 6.5], 'tickdir', 'out')

GL = [5 10];
for sn = 1:length(sub)
    for h = 1:length(GL)
        for t = 1:6
            i_dlpfc = strcmp(sub(sn).expt_name, 'dlpfc');
            i_vertex = strcmp(sub(sn).expt_name, 'vertex');
            
            
            idx = i_vertex;
            ind = (sub(sn).gameLength == GL(h)) & idx;
            corrdiff(h, t, 1, sn) = nanmean(sub(sn).co(ind,t+4));
            
            idx = i_dlpfc;
            ind = (sub(sn).gameLength == GL(h)) & idx;
            corrdiff(h, t, 2, sn) = nanmean(sub(sn).co(ind,t+4));
        end
    end
end




d_corr = squeeze(corrdiff(:,:,2,:)-corrdiff(:,:,1,:));



axes(ax(2)); hold on;
M = nanmean(d_corr,3);
S = nanstd(d_corr,[],3)/sqrt(length(sub));
plot([0 7], [0 0], 'k--', 'linewidth', 1)
e = errorbar(M', S');
legend(e, {'horizon 1' 'horizon 6'}, 'location', 'southeast')
xlabel('trial number')
ylabel({'\Delta(Correct Response)' 'rDLPFC - vertex'})

set(ax, 'xlim', [0.5 6.5], 'xtick', [1:6], 'tickdir', 'out')
set(e(:,1), 'color', arminn_pink)
set(e(:,2), 'color', arminn)
set(e, 'markersize', 30, 'marker', '.')
set(ax(2), 'ylim', [-0.1 0.05], 'xlim', [0.5 6.5])

[~,p1] = ttest(squeeze(corrdiff(2,:,2,:)-corrdiff(2,:,1,:))')
[~,p2] = ttest(squeeze(corrdiff(1,:,2,:)-corrdiff(1,:,1,:))')


addABCs(ax, [-0.09 0.07], 32)

%squeeze(mean(mean(correct1(1,:,:))))
%squeeze(mean(mean(correct6(1,:,:))))

%squeeze(std(std(correct1(1,:,:))))
%squeeze(std(std(correct6(1,:,:))))

%% Model-Free Analysis in first block (first 64), second block (second block 64-128) and first+second block (1-128 trials). 

clear p_highInfo p_lowMean p_highInfo_64 p_highInfo_128 p_lowMean_64 p_lowMean_128 

GL = [5 10];
for sn = 1:length(sub)
    for h = 1:length(GL)
        i_dlpfc = strcmp(sub(sn).expt_name, 'dlpfc');
        i_vertex = strcmp(sub(sn).expt_name, 'vertex');
        
        % First block (trials 1-64)
        idx = i_vertex;
        ind = (sub(sn).gameLength == GL(h)) & idx & (sub(sn).game <= 64);
        p_highInfo_64(h,sn,1) = nanmean(sub(sn).hi(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & idx & (sub(sn).game <= 128 & sub(sn).game > 64);
        p_highInfo_128(h,sn,1) = nanmean(sub(sn).hi(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & idx & (sub(sn).game <= 128);
        p_highInfo_overal(h,sn,1) = nanmean(sub(sn).hi(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)==sub(sn).n1(:,4)) & idx & (sub(sn).game <= 64);
        p_lowMean_64(h,sn,1) = nanmean(sub(sn).lm(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)==sub(sn).n1(:,4)) & idx & (sub(sn).game <= 128 & sub(sn).game > 64);
        p_lowMean_128(h,sn,1) = nanmean(sub(sn).lm(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)==sub(sn).n1(:,4)) & idx & (sub(sn).game <= 128);
        p_lowMean_overal(h,sn,1) = nanmean(sub(sn).lm(ind,5));
        
        idx = i_dlpfc;
        ind = (sub(sn).gameLength == GL(h)) & idx & (sub(sn).game <= 64);
        p_highInfo_64(h,sn,2) = nanmean(sub(sn).hi(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & idx & (sub(sn).game <= 128 & sub(sn).game > 64);
        p_highInfo_128(h,sn,2) = nanmean(sub(sn).hi(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & idx & (sub(sn).game <= 128);
        p_highInfo_overal(h,sn,2) = nanmean(sub(sn).hi(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)==sub(sn).n1(:,4)) & idx & (sub(sn).game <= 64);
        p_lowMean_64(h,sn,2) = nanmean(sub(sn).lm(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)==sub(sn).n1(:,4)) & idx & (sub(sn).game <= 128 & sub(sn).game > 64);
        p_lowMean_128(h,sn,2) = nanmean(sub(sn).lm(ind,5));
        ind = (sub(sn).gameLength == GL(h)) & (sub(sn).n2(:,4)==sub(sn).n1(:,4)) & idx & (sub(sn).game <= 128);
        p_lowMean_overal(h,sn,2) = nanmean(sub(sn).lm(ind,5));
        
    end
end




%% supplementory material 2: Model-free analysis of direct exploration in the first 64 and 128 trials

mean_highInfo_64 = nanmean(p_highInfo_64, 2);  % Dimensions: [GL, subject, condition]
mean_highInfo_128 = nanmean(p_highInfo_128, 2); % Dimensions: [GL, subject, condition]
mean_highInfo_overal = nanmean(p_highInfo_overal, 2); % Dimensions: [GL, subject, condition]

numSubjects_64 = size(p_highInfo_64, 2);
numSubjects_128 = size(p_highInfo_128, 2);
numSubjects_overal = size(p_highInfo_overal, 2);

std_highInfo_64 = std(p_highInfo_64, 0, 2) / sqrt(numSubjects_64);
std_highInfo_128 = std(p_highInfo_128, 0, 2) / sqrt(numSubjects_128);
std_highInfo_overal = std(p_highInfo_overal, 0, 2) / sqrt(numSubjects_overal);

% Create a figure window
figure;

% Plotting constants
barWidth = 0.35;
errorBarWidth = 1;
vertexColor = arminn;  % Vertex color
dlPFColor = arminn_pink;  % dlPFC color

% Subplot 1: 64 Trials, horizon 1
ax(1) = subplot(3,2,1);
barData1 = [mean_highInfo_64(1,:, 2)', mean_highInfo_64(1,:, 1)'];
b1 = bar(barData1);
set(b1, 'FaceColor', 'flat');
b1.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b1.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(high Info)');
x = 1:2;
hold on
err = [std_highInfo_64(1, :, 2)', std_highInfo_64(1, :, 1)'];
errorbar(x, barData1, err, 'k', 'linestyle', 'none');                                 
hold off
title('Trial 1-64, horizon 1');

% Subplot 2: 64 Trials, horizon 6
ax(2) = subplot(3,2,2);
barData2 = [mean_highInfo_64(2,:, 2)', mean_highInfo_64(2,:, 1)'];
b2 = bar(barData2);
set(b2, 'FaceColor', 'flat');
b2.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b2.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(high Info)');
x = 1:2;
hold on
err = [std_highInfo_64(2, :, 2)', std_highInfo_64(2, :, 1)'];
errorbar(x, barData2, err, 'k', 'linestyle', 'none');                                 
hold off
title('Trial 1-64, horizon 6');

% Subplot 3: 128 Trials, horizon 1
ax(3) = subplot(3,2,3);
barData3 = [mean_highInfo_128(1,:, 2)', mean_highInfo_128(1,:, 1)'];
b3 = bar(barData3);
set(b3, 'FaceColor', 'flat');
b3.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b3.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(high Info)');
x = 1:2;
hold on
err = [std_highInfo_128(1, :, 2)', std_highInfo_128(1, :, 1)'];
errorbar(x, barData3, err, 'k', 'linestyle', 'none');                                 
hold off
title('Trial 64-128, horizon 1');

% Subplot 4: 128 Trials, horizon 6
ax(4) = subplot(3,2,4);
barData4 = [mean_highInfo_128(2,:, 2)', mean_highInfo_128(2,:, 1)'];
b4 = bar(barData4);
set(b4, 'FaceColor', 'flat');
b4.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b4.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(high Info)');
x = 1:2;
hold on
err = [std_highInfo_128(2, :, 2)', std_highInfo_128(2, :, 1)'];
errorbar(x, barData4, err, 'k', 'linestyle', 'none');                                 
hold off
title('Trial 64-128, horizon 6');

% Subplot 5: Overall, horizon 1
ax(5) = subplot(3,2,5);
barData5 = [mean_highInfo_overal(1,:, 2)', mean_highInfo_overal(1,:, 1)'];
b5 = bar(barData5);
set(b5, 'FaceColor', 'flat');
b5.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b5.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(high Info)');
x = 1:2;
hold on
err = [std_highInfo_overal(1, :, 2)', std_highInfo_overal(1, :, 1)'];
errorbar(x, barData5, err, 'k', 'linestyle', 'none');                                 
hold off
title('Trial 1-128, horizon 1');

% Subplot 6: Overall, horizon 6
ax(6) = subplot(3,2,6);
barData6 = [mean_highInfo_overal(2,:, 2)', mean_highInfo_overal(2,:, 1)'];
b6 = bar(barData6);
set(b6, 'FaceColor', 'flat');
b6.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b6.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(high Info)');
x = 1:2;
hold on
err = [std_highInfo_overal(2, :, 2)', std_highInfo_overal(2, :, 1)'];
errorbar(x, barData6, err, 'k', 'linestyle', 'none');                                 
hold off
title('Trial 1-128, horizon 6');

% Optionally adjust subplot positions to avoid overlap and set figure size
set(gcf, 'Position', [100, 100, 800, 1000]); % Modify as needed

% Add an overall title to the figure
%sgtitle('High Information Under Different Conditions and Trials');

% Call the addABCs function after finalizing the plot
addABCs(ax, [-0.09 0.07], 32);

%% figure 4. : Model-free analysis of random exploration in the first 64 and 128 trials

mean_lowmean_64 = nanmean(p_lowMean_64, 2);  % Dimensions: [GL, subject, condition]
mean_lowmean_128 = nanmean(p_lowMean_128, 2); % Dimensions: [GL, subject, condition]
mean_lowmean_overal = nanmean(p_lowMean_overal, 2); % Dimensions: [GL, subject, condition]

numSubjects_64 = size(p_lowMean_64, 2);
numSubjects_128 = size(p_lowMean_128, 2);
numSubjects_overal = size(p_lowMean_overal, 2);
std_lowmean_64 = std(p_lowMean_64, 0, 2) / sqrt(numSubjects_64);
std_lowmean_128 = std(p_lowMean_128, 0, 2) / sqrt(numSubjects_128);
std_lowmean_overal = std(p_lowMean_overal, 0, 2) / sqrt(numSubjects_overal);

% Create a figure window
figure;

% Plotting constants
barWidth = 0.35;
errorBarWidth = 1;
vertexColor = arminn;  % Vertex color
dlPFColor = arminn_pink;  % dlPFC color

% Subplot 1: 64 Trials, horizon 1
ax(1) = subplot(3,2,1);
barData1 = [mean_lowmean_64(1,:, 2)', mean_lowmean_64(1,:, 1)'];
b1 = bar(barData1);
set(b1, 'FaceColor', 'flat');
b1.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b1.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(low Mean)');
x = 1:2;
hold on
err = [std_lowmean_64(1, :, 2)', std_lowmean_64(1, :, 1)'];
errorbar(x, barData1, err, 'k', 'linestyle', 'none');                                 
hold off
title('Trial 1-64, horizon 1');

% Subplot 2: 64 Trials, horizon 6
ax(2) = subplot(3,2,2);
barData2 = [mean_lowmean_64(2,:, 2)', mean_lowmean_64(2,:, 1)'];
b2 = bar(barData2);
set(b2, 'FaceColor', 'flat');
b2.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b2.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(low Mean)');
x = 1:2;
hold on
err = [std_lowmean_64(2, :, 2)', std_lowmean_64(2, :, 1)'];
errorbar(x, barData2, err, 'k', 'linestyle', 'none');                                 
hold off
title('Trial 1-64, horizon 6');

% Subplot 3: 128 Trials, horizon 1
ax(3) = subplot(3,2,3);
barData3 = [mean_lowmean_128(1,:, 2)', mean_lowmean_128(1,:, 1)'];
b3 = bar(barData3);
set(b3, 'FaceColor', 'flat');
b3.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b3.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(low Mean)');
x = 1:2;
hold on
err = [std_lowmean_128(1, :, 2)', std_lowmean_128(1, :, 1)'];
errorbar(x, barData3, err, 'k', 'linestyle', 'none');                                 
hold off
H = sigstar({[1,2]}, [0.004]);
title('Trial 64-128, horizon 1');

% Subplot 4: 128 Trials, horizon 6
ax(4) = subplot(3,2,4);
barData4 = [mean_lowmean_128(2,:, 2)', mean_lowmean_128(2,:, 1)'];
b4 = bar(barData4);
set(b4, 'FaceColor', 'flat');
b4.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b4.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(low Mean)');
x = 1:2;
hold on
err = [std_lowmean_128(2, :, 2)', std_lowmean_128(2, :, 1)'];
errorbar(x, barData4, err, 'k', 'linestyle', 'none');                                 
hold off
H = sigstar({[1,2]}, [0.048]);
title('Trial 64-128, horizon 6');

% Subplot 5: Overall, horizon 1
ax(5) = subplot(3,2,5);
barData5 = [mean_lowmean_overal(1,:, 2)', mean_lowmean_overal(1,:, 1)'];
b5 = bar(barData5);
set(b5, 'FaceColor', 'flat');
b5.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b5.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(low Mean)');
x = 1:2;
hold on
err = [std_lowmean_overal(1, :, 2)', std_lowmean_overal(1, :, 1)'];
errorbar(x, barData5, err, 'k', 'linestyle', 'none');                                 
hold off
H = sigstar({[1,2]}, [0.04]);
title('Trial 1-128, horizon 1');

% Subplot 6: Overall, horizon 6
ax(6) = subplot(3,2,6);
barData6 = [mean_lowmean_overal(2,:, 2)', mean_lowmean_overal(2,:, 1)'];
b6 = bar(barData6);
set(b6, 'FaceColor', 'flat');
b6.CData(1,:) = dlPFColor;  % Color for dlPFC condition
b6.CData(2,:) = vertexColor;  % Color for Vertex condition
set(gca, 'xticklabel', {'rDLPFC', 'Vertex'});
ylabel('p(low Mean)');
x = 1:2;
hold on
err = [std_lowmean_overal(2, :, 2)', std_lowmean_overal(2, :, 1)'];
errorbar(x, barData6, err, 'k', 'linestyle', 'none');                                 
hold off
H = sigstar({[1,2]}, [0.02]);
title('Trial 1-128, horizon 6');

% Optionally adjust subplot positions to avoid overlap and set figure size
set(gcf, 'Position', [100, 100, 800, 1000]); % Modify as needed

% Add an overall title to the figure
%sgtitle('Random Exploration Under Different Conditions and Trials');

% Call the addABCs function after finalizing the plot
addABCs(ax, [-0.09 0.07], 32);

%% paired t test for each block

%Direct Exploration for first second and overall block

%[~,p] = ttest(squeeze(p_highInfo_64(1,:,1)-p_highInfo_64(1,:,2))')
%[~,p] = ttest(squeeze(p_highInfo_64(2,:,1)-p_highInfo_64(2,:,2))')
%[~,p] = ttest(squeeze(p_highInfo_128(1,:,1)-p_highInfo_128(1,:,2))')
%[~,p] = ttest(squeeze(p_highInfo_128(2,:,1)-p_highInfo_128(2,:,2))')
%[~,p] = ttest(squeeze(p_highInfo_overal(1,:,1)-p_highInfo_128(1,:,2))')
%[~,p] = ttest(squeeze(p_highInfo_overal(2,:,1)-p_highInfo_128(2,:,2))')


%Random Exploration
[~,p] = ttest(squeeze(p_lowMean_64(1,:,1)-p_lowMean_64(1,:,2))')
[~,p] = ttest(squeeze(p_lowMean_64(2,:,1)-p_lowMean_64(2,:,2))')
[~,p] = ttest(squeeze(p_lowMean_128(1,:,1)-p_lowMean_128(1,:,2))')
[~,p] = ttest(squeeze(p_lowMean_128(2,:,1)-p_lowMean_128(2,:,2))')
[~,p] = ttest(squeeze(p_lowMean_overal(1,:,1)-p_lowMean_128(1,:,2))')
[~,p] = ttest(squeeze(p_lowMean_overal(2,:,1)-p_lowMean_128(2,:,2))')

%% supplementory material S4: reaction times

clear time 
GL = [5 10];
for sn = 1:length(sub)
    for h = 1:length(GL)
        for t = 1:6
            i_dlpfc = strcmp(sub(sn).expt_name, 'dlpfc');
            i_vertex = strcmp(sub(sn).expt_name, 'vertex');
            
            
            idx = i_vertex;
            ind = (sub(sn).gameLength == GL(h)) & idx;
            time(h, t, 1, sn) = nanmean(sub(sn).RT(ind,t+4));
            
            idx = i_dlpfc;
            ind = (sub(sn).gameLength == GL(h)) & idx;
            time(h, t, 2, sn) = nanmean(sub(sn).RT(ind,t+4));
        end
    end
end

figure(1); clf;
set(gcf, 'position', [588   600   600   500]);
ax = easy_gridOfEqualFigures([0.2 0.18], [0.16 0.18]);

axes(ax(1)); hold on;
xlim([0.5 6.5])
M = nanmean(time,4);
S = nanstd(time,[],4)/sqrt(length(sub));
e = errorbar(M(:,:,1)', S(:,:,1)');
e(2,:) = errorbar(M(:,:,2)', S(:,:,2)');
xlabel('trial number')
ylabel('RT')
title('"Reaction Times"')
leg = legend(e(:), {'vertex, horizon 1' 'rDLPFC, horizon 1' 'vertex, horizon 6' 'rDLPFC, horizon 6'})
set(leg, 'position', [0.2658    0.7610    0.2425    0.1150])

set(e([1],1), 'color', arminn_pink)
set(e([2],1), 'color', arminn_pink*0.5+0.5*[1 1 1])
set(e([1],2), 'color', arminn)
set(e([2],2), 'color', arminn*0.5+0.5*[1 1 1], 'linestyle', '--')
set(e([1],:), 'marker', '+')
set(e([2],:), 'marker', 'x')
set(e, 'markersize', 10);%, 'marker', '.')
set(ax(1), 'ylim', [0.1 2], 'ytick', [0:0.2:1.6], 'xlim', [0.5 6.5], 'xtick', [1:6])

set(ax, 'tickdir', 'out')
