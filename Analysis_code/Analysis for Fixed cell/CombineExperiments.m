clear; clc; close all;
root='D:\Projects\4. p21p27\';
datadir=[root 'CombinedData\'];
resultdir=[root 'Result\'];
if ~exist(resultdir,'dir')
    mkdir(resultdir);
end
%% plot Data
datafiles=getFilenames(datadir,'pRb');
 
%datafiles=datafiles(boolRegExp(datafiles,cellLine));
%datafiles=datafiles(~boolRegExp(datafiles,'._results'));
 
for i=1:length(datafiles)
    load([datadir datafiles{i}]);
    %Idx_1 = find(contains(label,'DMSO'));
    %Idx_2 = find(contains(label,'CDK4i'));
    %Idx_3 = find(contains(label,'CDK2i'));
    %idx=[Idx_1 Idx_2 Idx_3];
    
    pRb_G1_wt(i,1:2)=data{1}(1:2);
    pRb_S_wt(i,1:2)=data{2}(1:2);
    pRb_G2_wt(i,1:2)=data{3}(1:2);
    
%     pRb_G1_p21(i,1:5)=data{1}(6:10);
%     pRb_S_p21(i,1:5)=data{2}(6:10);
%     pRb_G2_p21(i,1:5)=data{3}(6:10);
end
%% Percentage of CDK2 low cells
option=2;
 
if option==1
plotV_1=100-pRb_G1_wt*100; 
plotV_2=100-pRb_G1_p21*100; ymax=100;
elseif option==2
plotV_1=100-pRb_S_wt*100; 
plotV_2=100-pRb_S_p21*100; ymax=20;
elseif option==3
plotV_1=100-pRb_G2_wt*100;
plotV_2=100-pRb_G2_p21*100; ymax=20;
end
 
figure; hold on;
x=0:2:8;
y = nanmean(plotV_1);
err = nanstd(plotV_1);%./sqrt([length(control) length(NCS)]);
y2 = nanmean(plotV_2);
err2 = nanstd(plotV_2);%./sqrt([length(control) length(NCS)]);
 
errorbar(x,y,err); errorbar(x,y2,err2); 
xlim([-1 9]); ylim([0 ymax]);
cllegend({'wt','p21'});
 
%%
barwitherr(errY, y, 0.5);    % Plot with errorbars
H = notBoxPlot(plotV,[],0.3);
set([H.data] , 'markersize' , 6)
set(gcf,'color','w');
set([H.sdPtch],'FaceAlpha',0);
set([H.semPtch],'FaceAlpha',0);
set([H.mu],'visible','off');
set(gca , 'FontSize' , 10)
set(gca,'XTickLabel',name)
%title(['p value=' num2str(p) ', independent exp. n=' num2str(length(files))]);
ylabel('Hypo-Rb (%)'); ylim([0 ymax]);
%% Percentage of CDK2 low cells
option=3;
 
if option==1
plotV=pRb_G1_wt*100; ymax=100;
elseif option==2
plotV=pRb_S_wt*100; ymax=18;
elseif option==3
plotV=pRb_G2_wt*100; ymax=18;
end
 
figure;
y = nanmean(plotV);
errY = nanstd(plotV);%./sqrt([length(control) length(NCS)]);
barwitherr(errY, y, 0.5);    % Plot with errorbars
H = notBoxPlot(plotV,[],0.3);
set([H.data] , 'markersize' , 6)
set(gcf,'color','w');
set([H.sdPtch],'FaceAlpha',0);
set([H.semPtch],'FaceAlpha',0);
set([H.mu],'visible','off');
set(gca , 'FontSize' , 10)
set(gca,'XTickLabel',name)
%title(['p value=' num2str(p) ', independent exp. n=' num2str(length(files))]);
ylabel('Hypo-Rb (%)'); ylim([0 ymax]);
%% Normality test
for i=1:3
    [h(i),p(i)]=kstest(plotV(:,i)) %[h,p]=jbtest(plotV(:,i));
end
%% one-way ANOVA
[p,tbl,stats] = anova(plotV);
multcompare(stats)
%% repeated anova
[p,tbl,stats] = ranova(plotV);
multcompare(stats)
 
%% one-way ANOVA
[p,tbl,stats] = anova2(plotV);
multcompare(stats)
 
%% a paired t-test
clear stats
[h2(1),p2(1),ci(1,1:2),stats(1)]=ttest(plotV(:,1),plotV(:,2))
[h2(2),p2(2),ci(2,1:2),stats(2)]=ttest(plotV(:,1),plotV(:,3))
%%
clear stats
[p,tbl,stats]=anova1(plotV)
multcompare(stats)
