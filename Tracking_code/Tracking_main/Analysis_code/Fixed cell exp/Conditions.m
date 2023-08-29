%% Initializes and clears the workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
%% Setting up row and column numbers for each condition %%%%%%%%%%%%%%%%%%%

option=3; %1=MCF7; 2=T47D
siteV=1:32; 
Col_1=7:8; Col_2=11;
if option==1 %MCF7
    conditions={
        'MCF7 DMSO',2,Col_1,siteV; %
        'MCF7 1day',3,8,siteV; %
        'MCF7 2day',4,Col_1,siteV; %
        'MCF7 3day',5,Col_1,siteV; %
        'MCF7 4day',6,Col_1,siteV; %
        'MCF7 5day',7,Col_1,siteV; %
        'MCF7 6day',8,Col_1,siteV; %
        };
elseif option==2 %T47D
    conditions={
        'T47D DMSO',2,Col_2,siteV; %
        'T47D 1day',3,11:12,siteV; %
        'T47D 2day',4,Col_2,siteV; %
        'T47D 3day',5,Col_2,siteV; %
        'T47D 4day',6,Col_2,siteV; %
        'T47D 5day',7,Col_2,siteV; %
        'T47D 6day',8,Col_2,siteV; %     
        };
elseif option==3
    conditions={
        'MCF7 DMSO',2,Col_1,siteV; %
        'MCF7 1day',3,7,siteV; %
        'T47D DMSO',2,Col_2,siteV; %
        'T47D 1day',3,10,siteV; %
        };

end

%%
projectpath='D:\1. Projects\5. Breast cancer project\'; %% GENERAL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
experimentpath='2023-08-09 Day since palbo_MCF7,T47D_EdU,c-Myc\';
% combinepath='1. Combine\OptoFGFR (Duration, EdU)\';

resultdir=[projectpath,experimentpath,'Results\'];
% comresultdir=[projectpath,combinepath,'Data\'];
datadir=[projectpath,experimentpath,'Data\'];
if ~exist(resultdir,'dir')
    mkdir(resultdir)
end
% if ~exist(comresultdir,'dir')
%     mkdir(comresultdir)
% end
%%
close all 
saveoption=0;
Data=CellCycleGating(conditions,datadir,option);
label=conditions(:,1);
cellPhase={'G1','S','G2'};

figu(0.3,1.0);
if saveoption==1
    if option==1
        save([resultdir 'daysincepalbo_07312022'],'label','cellPhase','Data'); %
    elseif option==2
        save([resultdir 'daysincepalbo_08012022'],'label','cellPhase','Data'); %
    end
end
% end
% load('percentage.mat');

%% Hoechst vs c-Myc

close all 
saveoption=0;
Data=DScatter_Hoechst_cMyc(conditions,datadir,option);
label=conditions(:,1);
cellPhase={'G1','S','G2'};

figu(0.3,1.0);

%%
close all

label=conditions(:,1);
cellPhase={'G1','S','G2'};
saveoption=0;

cellS=[4]; %1: G1; 2: S; 3: G2, 4: all;
DisplayOption='Panel';

option2=2; % 2: without normalization
protein=5; protein2=6; %5: Nuc; 6: Cyto (CDK4/CDK6/CycD)

counter=0;
for i=cellS
    counter=counter+1;
    data{counter}=Stain_HistComparison_HW_GOI(conditions,datadir,resultdir,protein,protein2,DisplayOption,i,option,option2);
end

if saveoption
    if option==1
        save([resultdir 'pRb_MCF10A_20210522.mat'],'label','cellPhase','data');
    elseif option==2
        save([resultdir 'pRb_RPE1,HS68_20210522.mat'],'label','cellPhase','data');
    elseif option==3
        if option2==2
            save([resultdir 'Nuc CDK4_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        elseif option2==1
            save([resultdir 'Nuc_cyto CDK4_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        end
    elseif option==4
        if option2==2
            save([resultdir 'Nuc CDK6_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        elseif option2==1
            save([resultdir 'Nuc_cyto CDK6_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        end
    elseif option==5
        if option2==2
            save([resultdir 'Nuc CycD_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        elseif option2==1
            save([resultdir 'Nuc_cyto CycD_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        end
    end    
end
figu(0.2,1.0)
%% Hoechst vs c-Myc
close all

label=conditions(:,1);
cellPhase={'G1','S','G2'};
saveoption=0;

cellS=[4]; %1: G1; 2: S; 3: G2, 4: all;
DisplayOption='Panel';

option2=2; % 2: without normalization
protein=5; protein2=6; %5: Nuc; 6: Cyto (CDK4/CDK6/CycD)

counter=0;
for i=cellS
    counter=counter+1;
    data{counter}=Stain_HistComparison_HW_GOI2(conditions,datadir,resultdir,protein,protein2,DisplayOption,i,option,option2);
end

if saveoption
    if option==1
        save([resultdir 'pRb_MCF10A_20210522.mat'],'label','cellPhase','data');
    elseif option==2
        save([resultdir 'pRb_RPE1,HS68_20210522.mat'],'label','cellPhase','data');
    elseif option==3
        if option2==2
            save([resultdir 'Nuc CDK4_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        elseif option2==1
            save([resultdir 'Nuc_cyto CDK4_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        end
    elseif option==4
        if option2==2
            save([resultdir 'Nuc CDK6_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        elseif option2==1
            save([resultdir 'Nuc_cyto CDK6_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        end
    elseif option==5
        if option2==2
            save([resultdir 'Nuc CycD_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        elseif option2==1
            save([resultdir 'Nuc_cyto CycD_MCF10A,RPE1,HS68_20210522.mat'],'label','cellPhase','data');
        end
    end    
end
figu(0.2,1.0)


%%
close all;
clear Data;
numCond=size(data{1},2);
count=0;
for i=1:numCond
    count=count+1;
    Data{i}=data{1}{count};
    cellNum(i)=length(Data{i});
end

fig_position = [200 200 600 400]; % coordinates for figures
% make figure
f10 = figure('Position', fig_position);
% violin(Data,'xlabel',{'-light','+light','-light','+light','-light','+light','-light','+light','-light','+light','-light','+light','-light','+light','-light','+light'},'facecolor',[0 1 0;0 1 0;1 0 0;1 0 0;0 1 0;0 1 0;1 0 0;1 0 0;0 1 0;0 1 0;1 0 0;1 0 0;0 1 0;0 1 0;1 0 0;1 0 0],'edgecolor','b',...
%     'bw',0.3,...
%     'mc','k',...
%     'medc','b--')

violin(Data,'facecolor',[0 1 0;1 1 0;1 0 0;1 0 1;1 1 1;0 1 0;1 1 0;1 0 0;1 0 1;1 1 1;0 1 0;1 1 0;1 0 0;1 0 1;1 1 1;0 1 0;1 1 0;1 0 0;1 0 1;1 1 1;0 1 0;1 1 0;1 0 0;1 0 1;1 1 1;0 1 0;1 1 0;1 0 0;1 0 1;1 1 1;],'edgecolor','b',...
'bw',0.3,...
'mc','k',...
'medc','b--')
if option==1
    set(gca,'Xtick',[1:numCond],'XTickLabel',{'No','1','2','3','4','5','6'},'fontsize',10,'tickdir','out','linewidth',2);
else
    set(gca,'Xtick',[1:numCond],'XTickLabel',{'No','1','2','3','4','5','6'},'fontsize',10,'tickdir','out','linewidth',2);
end
ylabel('','FontSize',12)
if option==1
    ymax=13;
    ylim([8 ymax])
elseif option==2
    ymax=14;
    ylim([8 ymax])
elseif option==3
    ymax=14;
    ylim([8 ymax])    
end
ax=gca;
ax.TickLength(1)=0.01;
%legend({'BRAFi','BRAFi+MEKi'});

%ylabel('');
figu(0.5,0.5);
%%




%%
clear Data;
numCond=size(data{1},2);
count=0;
for i=1:2
   for j=1:2
       count=count+1;
Data{i,j}=data{1}{count};
end
end
%%
% color scheme
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);

fig_position = [200 200 600 400]; % coordinates for figures
% make figure
f10 = figure('Position', fig_position);
h   = rm_raincloud(Data, cl);
%set(gca, 'YLim', [-0.3 1.6]);
%title(['Figure M10' newline 'Repeated measures raincloud plot - some aesthetic options'])

% define new colour
new_cl = [0.2 0.2 0.2];

% change one subset to new colour and alter dot size
h.p{2, 2}.FaceColor         = new_cl;
h.s{2, 2}.MarkerFaceColor   = new_cl;
h.m(2, 2).MarkerEdgeColor   = 'none';
h.m(2, 2).MarkerFaceColor   = new_cl;
h.s{2, 2}.SizeData          = 300;

%% 



%%
close all;
figure;
plot([0 0.25 0.5 1 2 4],data{1}./data{1}(1))
%%
close all;
figure;
plot([0 0.25 0.5 1 2 4],data{1}./data{1}(:,1))





%%
map=1; %1: 3D, 2: density
Binarize3Dmap(conditions,datadir,resultdir,map)
%%
%  protein2=8;
%8; 53BP puncta, 11; gH2AX puncta
option=1;
for i=1
    if i==1
        protein=5; protein2=6; name='CDK2';
    else
        protein=7; protein2=8; name='CDK4';
    end
    xlabel('Time since mitogen release (hr)');
    legend(conditions{i});
end

%%
figure; hold on;
for i=5:6
    plot(1:4,plotData{i})
end
ylabel('mRNA copy number'); xlabel('Time since mitogen release (hr)');
legend(conditions{5:6})





% export_fig([resultdir filesep 'pRb histogram.eps'],'-eps','-transparent','-nocrop');