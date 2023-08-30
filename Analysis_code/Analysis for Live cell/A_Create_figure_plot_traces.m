%% Reorganize Data (select data folder)
% clear all;
clear; clc;
root='D:\1. Projects\2. Melanoma project\2020-05-22 (Melanoma_DHFR-EGFR_TMP titration)\';
datadir=[root 'Data\Fixed'];
resultdir=[root 'Results\Fixed'];
if ~exist(resultdir)
    mkdir(resultdir);
end
option=2;
reorganize_Data(datadir,resultdir,option);
%% Load combined data
clear; close all;
clear ret
label={'0','0.001','0.005','0.01','0.05','0.1','0.5','1','2','5','10','50'};
for option=1:2;
    %root='D:\1. Projects\2. Melanoma project\2020-04-07 (Melanoma_ERK-CDK2-CDK4-H2B; non-resist vs resist; WM858_Dabrafenib)\';
    %resultdir=[root 'Results\'];
    root='D:\1. Projects\2. Melanoma project\2020-05-22 (Melanoma_DHFR-EGFR_TMP titration)\';
    resultdir=[root 'Results\Fixed\'];
    
    load([resultdir 'CombinedData_' num2str(option)]);
    
    %% Select condition of experiments 31 frames - 120 frames
    numcondi=numel(ret);
    for condi=1:numcondi;
        IF=0; nodivision=1; drugFrame=31;
        condition=condi;
        data = ret{condition}.store.storeData;
        divFrame = ret{condition}.store.allParentDivFrame;
        if IF
            IFdata = ret{condition}.store.IFdata;
        else
            numcells=size(data,1);
            IFdata=nan(numcells,1);
        end
        D = sort_or_align_by_mitosis(data,divFrame,IF,IFdata,nodivision)
        %% #1 Clear bad traces and short traces
        %D.clean_divnum_frame(16,1000);%20,29 / 22 31
        % D.clean_remove_short_before_mitosis(20);
        %D.clean_remove_short_after_mitosis(20);
        %D.clean_leastnum_frames(200); % this will remove frames shorter than x.
        D.clean_leastnum_frames_5(1,31,20,1);
        D.clean_leastnum_frames_5(32,50,10,1);
        D.clean_leastnum_frames_5(50,150,50,1);
        %D.clean_leastnum_frames_5(30,51,21);
        %D.clean_remove_suddendropCDK2
        %D.clean_remove_CDK2_not_dropped
        D.clean_remove_by_lower_upper_limit(1, [0, 3])
        D.clean_remove_by_lower_upper_limit(2, [0, 3])
        D.clean_remove_by_lower_upper_limit(3, [0, 3])
        % D.clean_remove_by_lower_upper_limit(3, [-50, 10000])
        
        %% #4 Clear bad traces (Manually)
        %D.clean_manually_remove_data(1,0); % For CDK2 axis([-2 0 0 3]);
        % D.clean_manually_remove_data(4,1); % For ERK
        %% smoothing dataSortedByMitosis
        D.smoothDataSortedByMitosis
        %% plot heatmap
        % D.sortXTickTime = [0:4:100]; % show tick at 8h,16h,24h...........
        % D.ylabels = {'CDK2 activity','ERK activity'};
        % D.xlabels = {'time (h)'};
        % D.colorRange{1} = [0.5 2.2];
        % D.colorRange{2} = [0.95 1.3];
        % D.clean_remove_mean_outside_n_std(2,3); % 1st arg means we use second channel, 2nd arg means discard outside mean+-2.5*std
        %D.plot_heatmap_sort(0)% input; smoothoption
        % figu(0.5,0.5);
        %colormap(parula_gray);
        % export_fig([resultdir,'heatmap_',conditions{condi},'.eps'],'-eps','-transparent','-nocrop');
        data=D.dataSorted(:,:,1);
        cellNum=size(data,1);
        
        tres=zeros(cellNum,1);
        for j=1:cellNum
            tres(j)=ttest2(data(j,1:drugFrame-1),data(j,drugFrame+2:drugFrame+20),.05,'left');
        end
        gating=logical(tres);
        %%
        sensor(condition).gating=gating;
        sensor(condition).ERK=D.dataSorted(:,:,1);
        sensor(condition).CDK2=D.dataSortedSmooth(:,:,2);
        sensor(condition).CDK4=D.dataSortedSmooth(:,:,3);
        sensor(condition).ERK_corrected=ERK_activity_correction_allcells(sensor(condition).ERK,sensor(condition).CDK2);
        sensor(condition).CDK4_corrected=CDK4_activity_correction_allcells(sensor(condition).CDK4,sensor(condition).CDK2);
        %% Save data
        
    end
    if option==1
        drug='Vem';
    elseif option==2
        drug='Dab';
    end
    save([resultdir filesep 'EGFR_plot' num2str(option) '.mat'],'sensor','label','drug');
end
%% Plotting data
option=2;

root='D:\1. Projects\2. Melanoma project\2020-05-22 (Melanoma_DHFR-EGFR_TMP titration)\';
resultdir=[root 'Results\Fixed'];
load([resultdir filesep 'EGFR_plot' num2str(option) '.mat'],'sensor','label','drug');
%% Heatmap + Averaged value
gate=1;
yaxis_1=[0 1.2]; yMax2=[1.6]; %fisrtFrame=5; drug_frame=105-fisrtFrame; lastFrame=105+100;
ci = 0.95; alpha = 1 - ci;
firstFrame=1;drugFrame=31;
maxframe=drugFrame+120;
figure; condition=0;
for i=1:length(sensor)
    %     condition=condition+1;
    %
    %     ERK_plotV_heatmap=ERK_plot{i};
    %     ERK_plotV_heatmap(ERK_plotV_heatmap<yaxis_1(1))=yaxis_1(1)+0.01;
    %
    %     subplot(4,12,condition); imagesc(ERK_plotV_heatmap,[yaxis_1]); colorbar; vline(drugFrame); colormap(ERK_colormap);
    %     set(gca,'Xtick',[firstFrame:25:350],'XTickLabel',[-20:5:20],'fontsize',10,'tickdir','out','linewidth',2)
    %     ylabel('Cells','FontSize',10);xlabel('Time since stimulation (hr)','FontSize',10); title(drug,'fontsize',20);
    %     xlim([firstFrame maxframe])
    
     gating=sensor(i).gating;
    if gate==1
        data=sensor(i).ERK_corrected(gating,:);
    else
        data=sensor(i).ERK_corrected(:,:);
    end
    
    ERK_mean_1=nanmean(data);
    S=nanstd(data);
    N=sum(~isnan(data));%sum(obj.cellState);
    T_multiplier = tinv(1-alpha/2, N-1);
    ci95 = T_multiplier.*S./sqrt(N);
    
    condition=condition+1;
    subplot(4,12,condition); hold on;
    plot(ERK_mean_1); plot(ERK_mean_1-ci95,'g'); plot(ERK_mean_1+ci95,'g'); vline(drugFrame);
    set(gca,'Xtick',[firstFrame:30:350],'XTickLabel',[-6:6:300],'fontsize',24,'tickdir','out','linewidth',2);
    set(gca,'Ytick',[0.2:0.2:2],'YTickLabel',[0.2:0.2:2],'fontsize',10,'tickdir','in','linewidth',2)
    ylabel('Cells','FontSize',10);xlabel('Time since inhibition (hr)','FontSize',10);title(['ERK activity n=' num2str(size(sensor(i).ERK_corrected,1))],'fontsize',20);
    xlim([firstFrame maxframe])
    if ismember(i,[1:24])
        axis([1 maxframe 0.3 1.4]);
    else
        axis([1 maxframe 0.2 0.8]);
    end
end
figu(0.5,1);

%% making color-map
a=makeColorMap([0.18 0.39 1],[1 1 1],300);
b=makeColorMap([1 1 1],[1 0 0],300);
CDK_colormap=[0.8 0.8 0.8; a; b];

a=parula(600);
ERK_colormap=[0.8 0.8 0.8; a];

%%
drugFrame=31;
firstFrame=1;
numframe=size(sensor(1).ERK,2);
numcondi=numel(sensor);

for i=1:numcondi
    gating=sensor(i).gating;
    if gate==1
        ERK=sensor(i).ERK_corrected(gating,:);
    else
        ERK=sensor(i).ERK_corrected(:,:);
    end
    gating=sum(ERK<0,2);
    ERK_gate=ERK(~gating,:);
    numcell=size(ERK_gate,1);
    idx=randperm(numcell);
    ERK_plot{i}=ERK_gate(idx,:);
end
%% Heatmap + Averaged value

yMax2=[1.6]; %fisrtFrame=5; drug_frame=105-fisrtFrame; lastFrame=105+100;
ci = 0.95; alpha = 1 - ci;

firstFrame=1;drugFrame=31;
maxframe=drugFrame+120;
figure; condition=0;
for i=1:length(sensor)
        condition=condition+1;
        
        if ismember(i,[1:24])
            yaxis_1=[0 1.2];
        else
            yaxis_1=[0 1];  
        end
        ERK_plotV_heatmap=ERK_plot{i};
        ERK_plotV_heatmap(ERK_plotV_heatmap<yaxis_1(1))=yaxis_1(1)+0.01;
        
        subplot(4,12,condition); imagesc(ERK_plotV_heatmap,[yaxis_1]); 
        colorbar; vline(drugFrame); colormap(ERK_colormap);
        set(gca,'Xtick',[firstFrame:25:350],'XTickLabel',[-20:5:20],'fontsize',10,'tickdir','out','linewidth',2)
        ylabel('Cells','FontSize',10);xlabel('Time since stimulation (hr)','FontSize',10); title(drug,'fontsize',20);
        xlim([firstFrame maxframe])
 
end
figu(0.5,1);

%%
condi=1;
figure; plot(ERK_plot{condi}')
%%
figure;
subplot(1,2,1); plot(nanmean(ERK_plot{1})); ylim([0.4 1.6]);
subplot(1,2,2); plot(nanmean(ERK_plot{2})); ylim([0.4 1.6]);
%%
figure;
subplot(1,2,1); plot(nanmean(sensor(1).CDK4_corrected)); ylim([0.3 0.9]);
subplot(1,2,2); plot(nanmean(sensor(2).CDK4_corrected)); ylim([0.3 0.9]);
%%
figure;
subplot(1,2,1); plot(nanmean(sensor(1).CDK2)); ylim([0.3 1.1]);
subplot(1,2,2); plot(nanmean(sensor(2).CDK2)); ylim([0.3 1.1]);

%% Calculating T1/2 time
yaxis_1=[1.5 1.5]; yMax2=[1.4 1.4];
cellLine={'WM858';};

figure; condition=0;
for i=1:length(sensor)
    
    ERK_plotV=ERK_plot{i}(:,drugFrame-100:drugFrame+120);
    
    data=nanmean(ERK_plotV);
    MAXpoint=find(data==max(data));
    decayERK=data(MAXpoint:MAXpoint+100);
    
    minNorm=min(decayERK);
    decayERK=decayERK-minNorm;
    decayERK=decayERK/max(decayERK);
    
    condition=condition+1;
    subplot(1,3,condition); hold on;
    xvals=[0:0.2:0.2*(length(decayERK)-1)];
    fitpara=fit(xvals',decayERK','exp2');
    halfLife(i)=log(2)/abs((fitpara.b)); %convert to time
    
    plot(fitpara,xvals',decayERK');hold on;
    
    %     set(gca,'Xtick',[0:100:lastFrame-fisrtFrame],'XTickLabel',[0:20:300],'fontsize',24,'tickdir','out','linewidth',2);
    %     set(gca,'Ytick',[0.3:0.3:1.5],'YTickLabel',[0.3:0.3:1.5],'fontsize',24,'tickdir','in','linewidth',5)
    ylim([0 1]); xlim([0 20]);
    ylabel('ERK activity','FontSize',30);xlabel('Time (hr)','FontSize',30);title(['T 1/2=' num2str(halfLife(i)) ' hr'],'fontsize',30);
    figu(0.5,1);
end
%% Averaged values
numFrame=size(sensor(1).ERK,2);
ci = 0.95; alpha = 1 - ci;

Data_1=nanmean(sensor(1).ERK_corrected);
S=nanstd(sensor(1).ERK_corrected);
N=sum(~isnan(sensor(1).ERK_corrected));%sum(obj.cellState);
T_multiplier = tinv(1-alpha/2, N-1);
ci95_1 = T_multiplier.*S./sqrt(N);

Data_2=nanmean(sensor(2).ERK_corrected);
S=nanstd(sensor(2).ERK_corrected);
N=sum(~isnan(sensor(2).ERK_corrected));%sum(obj.cellState);
T_multiplier = tinv(1-alpha/2, N-1);
ci95_2 = T_multiplier.*S./sqrt(N);

Data_3=nanmean(sensor(3).ERK_corrected);
S=nanstd(sensor(3).ERK_corrected);
N=sum(~isnan(sensor(3).ERK_corrected));%sum(obj.cellState);
T_multiplier = tinv(1-alpha/2, N-1);
ci95_3 = T_multiplier.*S./sqrt(N);

figure;
subplot(1,3,1);hold on;
plot(Data_1,'r'); plot(Data_1+ci95,'r'); plot(Data_1-ci95,'r');
plot(Data_2,'g'); plot(Data_2+ci95,'g'); plot(Data_2-ci95,'g');
plot(Data_3,'b'); plot(Data_3+ci95,'b'); plot(Data_3-ci95,'b');
vline(drugFrame);
set(gca,'Xtick',[6:25:200],'XTickLabel',[-20:5:200],'fontsize',10,'tickdir','out','linewidth',2)
ylabel('ERK activity','FontSize',10);xlabel('Time since mitogen stimulation (hr)','FontSize',10); title(label{i},'fontsize',20);
ylim([0.3 1.4]); xlim([6 numFrame]);

Data_1=nanmean(sensor(1).CDK4_corrected);
S=nanstd(sensor(1).CDK4_corrected);
N=sum(~isnan(sensor(1).CDK4_corrected));%sum(obj.cellState);
T_multiplier = tinv(1-alpha/2, N-1);
ci95_1 = T_multiplier.*S./sqrt(N);

Data_2=nanmean(sensor(2).CDK4_corrected);
S=nanstd(sensor(2).CDK4_corrected);
N=sum(~isnan(sensor(2).CDK4_corrected));%sum(obj.cellState);
T_multiplier = tinv(1-alpha/2, N-1);
ci95_2 = T_multiplier.*S./sqrt(N);

Data_3=nanmean(sensor(3).CDK4_corrected);
S=nanstd(sensor(3).CDK4_corrected);
N=sum(~isnan(sensor(3).CDK4_corrected));%sum(obj.cellState);
T_multiplier = tinv(1-alpha/2, N-1);
ci95_3 = T_multiplier.*S./sqrt(N);

subplot(1,3,2);hold on;
plot(Data_1,'r'); plot(Data_1+ci95,'r'); plot(Data_1-ci95,'r');
plot(Data_2,'g'); plot(Data_2+ci95,'g'); plot(Data_2-ci95,'g');
plot(Data_3,'b'); plot(Data_3+ci95,'b'); plot(Data_3-ci95,'b');
vline(drugFrame);
set(gca,'Xtick',[6:25:200],'XTickLabel',[-20:5:200],'fontsize',10,'tickdir','out','linewidth',2)
ylabel('CDK4/6 activity','FontSize',10);xlabel('Time since mitogen stimulation (hr)','FontSize',10); title(label{i},'fontsize',20);
ylim([0.2 0.75]); xlim([6 numFrame]);

Data_1=nanmean(sensor(1).CDK2);
S=nanstd(sensor(1).CDK2);
N=sum(~isnan(sensor(1).CDK2));%sum(obj.cellState);
T_multiplier = tinv(1-alpha/2, N-1);
ci95_1 = T_multiplier.*S./sqrt(N);

Data_2=nanmean(sensor(2).CDK2);
S=nanstd(sensor(2).CDK2);
N=sum(~isnan(sensor(2).CDK2));%sum(obj.cellState);
T_multiplier = tinv(1-alpha/2, N-1);
ci95_2 = T_multiplier.*S./sqrt(N);

Data_3=nanmean(sensor(3).CDK2);
S=nanstd(sensor(3).CDK2);
N=sum(~isnan(sensor(3).CDK2));%sum(obj.cellState);
T_multiplier = tinv(1-alpha/2, N-1);
ci95_3 = T_multiplier.*S./sqrt(N);

subplot(1,3,3);hold on;
plot(Data_1,'r'); plot(Data_1+ci95,'r'); plot(Data_1-ci95,'r');
plot(Data_2,'g'); plot(Data_2+ci95,'g'); plot(Data_2-ci95,'g');
plot(Data_3,'b'); plot(Data_3+ci95,'b'); plot(Data_3-ci95,'b');
vline(drugFrame);
set(gca,'Xtick',[6:25:200],'XTickLabel',[-20:5:200],'fontsize',10,'tickdir','out','linewidth',2)
ylabel('CDK2 activity','FontSize',10);xlabel('Time since mitogen stimulation (hr)','FontSize',10); title(label{i},'fontsize',20);
ylim([0.25 1.1]); xlim([6 numFrame]);
figu(0.4,0.9);
%%





















%%
condi=1;
figure;
subplot(1,2,1); imagesc(sensor(condi).CDK4,[0 1.2]); colormap(parula_gray);
subplot(1,2,2); imagesc(sensor(condi).CDK4_corrected,[0 1.2]); colormap(parula_gray);

figure;
subplot(1,2,1); imagesc(sensor(condi).ERK,[0 1.2]); colormap(parula_gray);
subplot(1,2,2); imagesc(sensor(condi).ERK_corrected,[0 1.2]); colormap(parula_gray);
%%
numFrame=size(sensor(1).ERK,2);

figure;
subplot(1,3,1);hold on;
plot(nanmean(sensor(1).ERK_corrected));
plot(nanmean(sensor(2).ERK_corrected));
plot(nanmean(sensor(3).ERK_corrected));
vline(drugFrame);
set(gca,'Xtick',[0:25:numFrame],'XTickLabel',[0:5:200],'fontsize',24,'tickdir','out','linewidth',2);
ylabel('ERK activity','FontSize',30);xlabel('Time (hr)','FontSize',30);title('Corrected ERK activity','fontsize',30);
ylim([0.2 1.4]);
legend(label)

subplot(1,3,2);hold on;
plot(nanmean(sensor(1).CDK4_corrected));
plot(nanmean(sensor(2).CDK4_corrected));
plot(nanmean(sensor(3).CDK4_corrected));
vline(drugFrame);
set(gca,'Xtick',[0:25:numFrame],'XTickLabel',[0:5:200],'fontsize',24,'tickdir','out','linewidth',2);
ylabel('CDK4 activity','FontSize',30);xlabel('Time (hr)','FontSize',30);title('Corrected CDK4 activity','fontsize',30);
ylim([0.1 1]);

subplot(1,3,3);hold on;
plot(nanmean(sensor(1).CDK2));
plot(nanmean(sensor(2).CDK2));
plot(nanmean(sensor(3).CDK2));
vline(drugFrame);
set(gca,'Xtick',[0:25:numFrame],'XTickLabel',[0:5:200],'fontsize',24,'tickdir','out','linewidth',2);
ylabel('CDK2 activity','FontSize',30);xlabel('Time (hr)','FontSize',30);title('CDK2 activity','fontsize',30);
ylim([0.1 1.2]);
figu(0.4,0.9);
%%
yaxis_1=[1.2]; yMax2=[0.85]; %fisrtFrame=5; drug_frame=105-fisrtFrame; lastFrame=105+100;
cellLine={'WM858';};
ci = 0.95; alpha = 1 - ci;
figure; condition=0;
for i=2:length(sensor)
    condition=condition+1;
    subplot(2,2,condition); imagesc(sensor(i).ERK_corrected,[0 yaxis_1]); colorbar; vline(drugFrame); colormap(parula_gray);
    set(gca,'Xtick',[9:25:200],'XTickLabel',[-20:5:200],'fontsize',10,'tickdir','out','linewidth',2)
    ylabel('Cells','FontSize',10);xlabel('Time since inhibition (hr)','FontSize',10); title(label{i},'fontsize',20);
    xlim([9 209])
    
    ERK_mean_1=nanmean(sensor(i).ERK_corrected);
    S=nanstd(sensor(i).ERK_corrected);
    N=sum(~isnan(sensor(i).ERK_corrected));%sum(obj.cellState);
    T_multiplier = tinv(1-alpha/2, N-1);
    ci95 = T_multiplier.*S./sqrt(N);
    
    condition=condition+1;
    subplot(2,2,condition); hold on;
    plot(nanmean(sensor(i).ERK_corrected)); plot(nanmean(sensor(i).ERK_corrected)-ci95,'g'); plot(nanmean(sensor(i).ERK_corrected)+ci95,'g'); axis([0 315 0.25 yMax2]); vline(drugFrame);
    set(gca,'Xtick',[9:25:200],'XTickLabel',[-20:5:200],'fontsize',24,'tickdir','out','linewidth',2);
    set(gca,'Ytick',[0.2:0.2:1.2],'YTickLabel',[0.2:0.2:1.2],'fontsize',10,'tickdir','in','linewidth',2)
    ylabel('Cells','FontSize',10);xlabel('Time since inhibition (hr)','FontSize',10);title(['ERK activity n=' num2str(size(sensor(i).ERK_corrected,1))],'fontsize',20);
    xlim([9 209])
end
figu(1,0.5);
%export_fig([savedir filesep 'ERKplot.eps'],'-eps','-transparent','-nocrop');

%%
yaxis_1=[1.2]; yMax2=[0.8]; %fisrtFrame=5; drug_frame=105-fisrtFrame; lastFrame=105+100;
cellLine={'WM858';};
ci = 0.95; alpha = 1 - ci;
figure; condition=0;
for i=1:length(sensor)
    condition=condition+1;
    subplot(2,2,condition); imagesc(sensor(i).CDK4_corrected,[0 yaxis_1]); colorbar; vline(drugFrame); colormap(parula_gray);
    set(gca,'Xtick',[9:25:200],'XTickLabel',[-20:5:200],'fontsize',10,'tickdir','out','linewidth',2)
    ylabel('Cells','FontSize',10);xlabel('Time since inhibition (hr)','FontSize',10); title(label{i},'fontsize',20);
    xlim([9 209])
    
    ERK_mean_1=nanmean(sensor(i).CDK4_corrected);
    S=nanstd(sensor(i).ERK_corrected);
    N=sum(~isnan(sensor(i).ERK_corrected));%sum(obj.cellState);
    T_multiplier = tinv(1-alpha/2, N-1);
    ci95 = T_multiplier.*S./sqrt(N);
    
    condition=condition+1;
    subplot(2,2,condition); hold on;
    plot(nanmean(sensor(i).CDK4_corrected)); plot(nanmean(sensor(i).CDK4_corrected)-ci95,'g'); plot(nanmean(sensor(i).CDK4_corrected)+ci95,'g'); axis([0 315 0.2 yMax2]); vline(drugFrame);
    set(gca,'Xtick',[9:25:200],'XTickLabel',[-20:5:200],'fontsize',24,'tickdir','out','linewidth',2);
    set(gca,'Ytick',[0.2:0.2:1.2],'YTickLabel',[0.2:0.2:1.2],'fontsize',10,'tickdir','in','linewidth',2)
    ylabel('Cells','FontSize',10);xlabel('Time since inhibition (hr)','FontSize',10);title(['ERK activity n=' num2str(size(sensor(i).ERK_corrected,1))],'fontsize',20);
end
figu(1,0.5);



%% Select condition of experiments
for condi=1:3;
    IF=0; nodivision=0;
    condition=condi;
    conditions = {'non-resis';'non-resis_BRAFi';'resis_BRAFi';};
    data = ret{condition}.store.storeData;
    divFrame = ret{condition}.store.allParentDivFrame;
    if IF
        IFdata = ret{condition}.store.IFdata;
    else
        IFdata=nan;
    end
    D = sort_or_align_by_mitosis(data,divFrame,IF,IFdata,nodivision)
    %% #1 Clear bad traces and short traces
    drug_frame=0;
    %D.clean_leastnum_frames_5(drug_frame-10,drug_frame+10,10);
    %D.clean_leastnum_frames_5(1,50,10);
    %  D.clean_divnum_frame(0,600);%20,29 / 22 31
    % D.clean_remove_short_before_mitosis(10);
    % D.clean_remove_short_after_mitosis(20);
    %D.clean_leastnum_frames(30); % this will remove frames shorter than x.
    % D.clean_leastnum_frames_5(1,100,40);
    % D.clean_leastnum_frames_5(30,51,21);
    % D.clean_leastnum_frames_5(drug_frame-100,drug_frame+100,120);
    %D.clean_remove_suddendropCDK2
    %D.clean_remove_CDK2_not_dropped
    %D.clean_remove_by_lower_upper_limit(1, [1, 3]) %RAS
    %D.clean_remove_by_lower_upper_limit(2, [0, 3]) %ERK
    % D.clean_remove_by_lower_upper_limit(3, [-50, 10000]);
    %% 3 Clear bad traces (Manually)
    %D.clean_manually_remove_data(1,0); % For CDK2 axis([-2 0 0 3]);
    %D.clean_manually_remove_data(4,1); % For ERK
    
    %% smoothing dataSortedByMitosis
    D.smoothDataSortedByMitosis
    %%
    i=condi
    RAS{i}=D.dataSorted(:,:,1);
    RAS_smooth{i}=D.dataSortedSmooth(:,:,1);
    
    
    RAS_act{i}=nanmean(D.dataSorted(:,1:50,1)');
    numcellCondi=length(RAS_act{i});
    
    if option==1
        numCell=20;
    elseif option==2
        numCell=20;
    end
    smaple=randsample(1:numcellCondi,numCell);
    RAS_act_selected{i}=RAS_act{i}(smaple);
    
    
    %% Bleaching correction
    bleachCorrectionLinear=fitline(1:drug_frame,nanmean(RAS_smooth{i}(:,1:drug_frame),1)); %m: slope, b: y-intercept
    %bleachCorrectionExp=fitline(1:drug_frame,log(nanmean(RAS_smooth{i}(:,1:drug_frame),1)));
    
    b=bleachCorrectionLinear.b;
    bleachFactor=bleachCorrectionLinear.m./bleachCorrectionLinear.b;
    
    bleachCor_RAS{i}=RAS_smooth{i}./repmat(1+b*bleachFactor*(1:size(RAS_smooth{i},2)),[size(RAS_smooth{i},1) 1]);
end
%%
for i=1:3
    num(i)=size(RAS_act_selected{i},2);
    std(i)=nanstd(RAS_act_selected{i});
    plotV(i)=nanmean(RAS_act_selected{i});
    
end
x=1:3;

figure;
bar(x,plotV); hold on;
er=errorbar(1:3,plotV,std); er.Color = [0 0 0]; er.LineStyle = 'none';
ylim([1.5 2.5]); title(['non=',num2str(num(1)),' non+BRAFi=',num2str(num(2)),'res+BRAFi=',num2str(num(3))]);
%%
for i=1:3
    num(i)=size(RAS_act{i},2);
    std(i)=nanstd(RAS_act{i});
    plotV(i)=nanmean(RAS_act{i});
end
x=1:3;

figure;
bar(x,plotV); hold on;
er=errorbar(1:3,plotV,std); er.Color = [0 0 0]; er.LineStyle = 'none';
ylim([1.5 2.5]); title(['non=',num2str(num(1)),' non+BRAFi=',num2str(num(2)),'res+BRAFi=',num2str(num(3))]);
%%
save([savedir filesep 'Data_20180712.mat'],'RAS_act_selected','ERK_act_selected')
%%
for i=1:3;
    numFrame=size(ERK{i},2);
    if option==1
        yMax_RAS=3; yMax_erk=1.5;
    else
        yMax_RAS=3; yMax_erk=1.2;
    end
    % drug_frame=[6 121];
    figure(i);
    % subplot(2,2,1); imagesc(ERK,[0 1.5]); colorbar; vline([drug_frame]);
    % % set(gca,'YAxisLocation','right','YColor',[0 0.63 0],'YTick',0.6:0.4:1.8,'fontsize',24,'tickdir','out','linewidth',2);
    % set(gca,'Xtick',[0:100:numFrame],'XTickLabel',[0:20:300],'fontsize',24,'tickdir','out','linewidth',2)
    % ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title('Uncorrected ERK activity','fontsize',30);
    
    
    subplot(1,2,1); imagesc(RAS_smooth{i},[0 yMax_RAS]); colorbar; vline([drug_frame]);
    set(gca,'Xtick',[0:60:numFrame],'XTickLabel',[0:12:300],'fontsize',24,'tickdir','out','linewidth',2)
    ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title(['RAS activity n=' num2str(size(ERK{i},1))],'fontsize',30);
    colormap('parula_gray')
    
    subplot(1,2,2); imagesc(ERK_smooth{i},[0 yMax_erk]); colorbar; vline([drug_frame]);
    set(gca,'Xtick',[0:60:numFrame],'XTickLabel',[0:12:300],'fontsize',24,'tickdir','out','linewidth',2)
    ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title('ERK activity','fontsize',30);
    
    figu(1,1);
    %export_fig([savedir filesep 'Heatmap_' conditions{i} '.eps'],'-eps','-transparent','-nocrop');
end
%%
save([savedir filesep 'Data.mat'],'CDK2','CDK4','ERK','ERK_corrected','CDK4_corrected')


%%
option=1
cellLine='858';

savedir=[root 'Results\' cellLine];
load([savedir filesep 'Data.mat'],'CDK2','CDK4','ERK','ERK_corrected','CDK4_corrected')


%%
figure; hold on;
for i=1:3
    plot(nanmean(RAS_smooth{i}));
    %plot(nanmean(ERK{i}));
    %plot(nanmean(bleachCor_RAS{i}));
end
legend(conditions)
%%


%%
yaxis_1=[1.5]; yMax2=[1.2]; %fisrtFrame=5; drug_frame=105-fisrtFrame; lastFrame=105+100;
cellLine={'WM858'; 'WM938B';};
ci = 0.95; alpha = 1 - ci;
figure; condition=0;
for i=1:length(ERK_corrected)
    condition=condition+1;
    subplot(4,2,condition); imagesc(ERK_corrected{i},[0 yaxis_1]); colorbar; vline(drug_frame); colormap(parula_gray);
    set(gca,'Xtick',[2:40:350],'XTickLabel',[-16:8:300],'fontsize',10,'tickdir','out','linewidth',2)
    ylabel('Cells','FontSize',10);xlabel('Time since inhibition (hr)','FontSize',10); title(conditions{i},'fontsize',20);
    
    ERK_mean_1=nanmean(ERK_corrected{i});
    S=nanstd(ERK_corrected{i});
    N=sum(~isnan(ERK_corrected{i}));%sum(obj.cellState);
    T_multiplier = tinv(1-alpha/2, N-1);
    ci95 = T_multiplier.*S./sqrt(N);
    
    condition=condition+1;
    subplot(4,2,condition); hold on;
    plot(nanmean(ERK_corrected{i})); plot(nanmean(ERK_corrected{i})-ci95,'g'); plot(nanmean(ERK_corrected{i})+ci95,'g'); axis([0 315 0.2 yMax2]); vline(drug_frame);
    set(gca,'Xtick',[2:40:350],'XTickLabel',[-16:8:300],'fontsize',24,'tickdir','out','linewidth',2);
    set(gca,'Ytick',[0.2:0.2:1.2],'YTickLabel',[0.2:0.2:1.2],'fontsize',10,'tickdir','in','linewidth',2)
    ylabel('Cells','FontSize',10);xlabel('Time since inhibition (hr)','FontSize',10);title(['ERK activity n=' num2str(size(ERK_corrected{i},1))],'fontsize',20);
end
figu(1,0.5);
export_fig([savedir filesep 'ERKplot.eps'],'-eps','-transparent','-nocrop');
%% #3 smoothing dataSortedByMitosis
D.smoothDataSortedByMitosis
%%
if ismember(condi,[1,5])
    name='resistance_Control';
elseif ismember(condi,[3,7])
    name='resistance_CDK4i';
elseif ismember(condi,[4,7])
    name='resistance_MW';
end


%%
figure;
% subplot(2,2,1); plot(ERK'); axis([0 numFrame 0 2.5]); vline([145 drug_frame]);
% set(gca,'Xtick',[0:100:numFrame],'XTickLabel',[0:20:300],'fontsize',24,'tickdir','out','linewidth',2);
% ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title('Uncorrected ERK activity','fontsize',30);

subplot(1,3,1); plot(ERK_corrected'); axis([0 numFrame 0 2]); vline([drug_frame]);
set(gca,'Xtick',[0:100:numFrame],'XTickLabel',[0:20:300],'fontsize',24,'tickdir','out','linewidth',2);
ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title('Corrected ERK activity','fontsize',30);

subplot(1,3,2); plot(CDK2'); axis([0 numFrame 0 3]); vline([drug_frame]);
set(gca,'Xtick',[0:100:numFrame],'XTickLabel',[0:20:300],'fontsize',24,'tickdir','out','linewidth',2);
ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title('CDK2 activity','fontsize',30);

subplot(1,3,3); plot(CDK4_corrected'); axis([0 numFrame 0 1.5]); vline([drug_frame]);
set(gca,'Xtick',[0:100:numFrame],'XTickLabel',[0:20:300],'fontsize',24,'tickdir','out','linewidth',2);
ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title('CDK4 activity','fontsize',30);
figu(1,1)
% export_fig([savedir filesep 'Traces_resistance_' name '.eps'],'-eps','-transparent','-nocrop');

%%
% num=2;
% ERK_plot{num}=ERK_corrected;
% save([resultdir filesep 'ERK_plot.mat'],'ERK_plot');
% load([resultdir filesep 'ERK_plot.mat'],'ERK_plot');
yaxis_1=[1.5 1.2]; yMax2=[3 1]; fisrtFrame=13; drug_frame=[145 413]-fisrtFrame;
cellLine={'WM858'; 'WM938B';};
ci = 0.95; alpha = 1 - ci;
figure; condition=0;
for i=1:length(ERK_corrected)
    condition=condition+1;
    subplot(2,2,condition); imagesc(ERK_corrected{i}(:,fisrtFrame:end),[0 yaxis_1(i)]); colorbar; vline(drug_frame); colormap(parula_gray);
    set(gca,'Xtick',[0:100:numFrame],'XTickLabel',[0:20:300],'fontsize',24,'tickdir','out','linewidth',2)
    ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title(['ERK activity n=' num2str(size(ERK_corrected{i},1))],'fontsize',30);
    
    ERK_mean_1=nanmean(ERK_corrected{i}(:,fisrtFrame:end));
    S=nanstd(ERK_corrected{i}(:,fisrtFrame:end));
    N=sum(~isnan(ERK_corrected{i}(:,fisrtFrame:end)));%sum(obj.cellState);
    T_multiplier = tinv(1-alpha/2, N-1);
    ci95 = T_multiplier.*S./sqrt(N);
    
    condition=condition+1;
    subplot(2,2,condition); hold on;
    plot(nanmean(ERK_corrected{i}(:,fisrtFrame:end))); plot(nanmean(ERK_corrected{i}(:,fisrtFrame:end))-ci95,'g'); plot(nanmean(ERK_corrected{i}(:,fisrtFrame:end))+ci95,'g'); axis([fisrtFrame numFrame-fisrtFrame 0.2 yMax2(i)]); vline(drug_frame);
    set(gca,'Xtick',[0:100:numFrame-fisrtFrame],'XTickLabel',[0:20:300],'fontsize',24,'tickdir','out','linewidth',2);
    ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title(['ERK activity cell=' cellLine{i}],'fontsize',30);
end
figu(1,1);
% export_fig([resultdir filesep 'ERKplot_Control+Vem.eps'],'-eps','-transparent','-nocrop');



%% Section-2; Aligned By Mitosis
% #1 Smoothing dataAligned
D.smoothDataAlignedByMitosis
ERK_corrected=ERK_activity_correction_allcells(D.dataAligned(:,:,1),D.dataAlignedSmooth(:,:,2));
CDK4_corrected=CDK4_activity_correction_allcells(D.dataAlignedSmooth(:,:,3),D.dataAlignedSmooth(:,:,2));
%% Detectin CDK2 rising point
tracestats=nan(size(D.dataAlignedSmooth,1),2);
tracestats(:,1)=5001;
tracestats(:,2)=sum(~isnan(D.dataAligned(:,5001:5100,1)),2)+5000;


[risetime_cdk4 , cdk4_low , badtraces_cdk4] = getCdk2rise_yilin(D.dataAligned(:,:,3) , tracestats);
[risetime_cdk2 , cdk2_low , badtraces_cdk2] = getCdk2rise_yilin(D.dataAligned(:,:,2) , tracestats);
%%
figure; hold on;
for numCell=1:10
    
    cellplot=find(~isnan(risetime_cdk4));
    cell_to_plot=cellplot(numCell);
    plot(D.dataAlignedSmooth(cell_to_plot,:,3)');
    
    scatter(risetime_cdk4(cell_to_plot),D.dataAlignedSmooth(cell_to_plot,risetime_cdk4(cell_to_plot),3),50,'k','filled');
end
%%
figure; hold on;
for numCell=1:10
    
    cellplot=find(~isnan(risetime_cdk2));
    cell_to_plot=cellplot(numCell);
    plot(D.dataAlignedSmooth(cell_to_plot,:,3)');
    
    scatter(risetime_cdk2(cell_to_plot),D.dataAlignedSmooth(cell_to_plot,risetime_cdk2(cell_to_plot),2),50,'k','filled');
end

%% Adding data
condition=1;
pre_size=size(sensor(condition).ERK,1);
curr_size=size(D.dataAligned,1);
SF=4900; EF=5100;
sensor(condition).ERK((pre_size+1):(pre_size+curr_size),:)=D.dataAligned(:,SF:EF,1);
sensor(condition).ERK_corrected((pre_size+1):(pre_size+curr_size),:)=ERK_corrected(:,SF:EF,1);
sensor(condition).CDK2((pre_size+1):(pre_size+curr_size),:)=D.dataAlignedSmooth(:,SF:EF,2);
sensor(condition).CDK4((pre_size+1):(pre_size+curr_size),:)=D.dataAlignedSmooth(:,SF:EF,3);
sensor(condition).CDK4_corrected((pre_size+1):(pre_size+curr_size),:)=CDK4_corrected(:,SF:EF);
% sensor(condition).CDK4high=[sensor(condition).CDK4high; gating_1];
% sensor(condition).CDK4delay=[sensor(condition).CDK4delay; gating_3];
% sensor(condition).CDK4low=[sensor(condition).CDK4low; gating_2];
% sensor(condition).CDK2high=[sensor(condition).CDK2high; gating_4'];
% sensor(condition).CDK2delay=[sensor(condition).CDK2delay; gating_6];
% sensor(condition).CDK2low=[sensor(condition).CDK2low; gating_5'];
% sensor(condition).risetime_cdk4=[sensor(condition).risetime_cdk4; risetime_cdk4-SF-1];
% sensor(condition).risetime_cdk2=[sensor(condition).risetime_cdk2; risetime_cdk2-SF-1];
% sensor(condition).cdk4_low=[sensor(condition).cdk4_low; cdk4_low];
% sensor(condition).cdk2_low=[sensor(condition).cdk2_low; cdk2_low];
% sensor(condition).badtraces_cdk4=[sensor(condition).badtraces_cdk4; badtraces_cdk4];
% sensor(condition).badtraces_cdk2=[sensor(condition).badtraces_cdk2; badtraces_cdk2];




%%
condition=1;
SF=4900; EF=5100;
sensor(condition).ERK=D.dataAligned(:,SF:EF,1);
sensor(condition).ERK_corrected=ERK_corrected(:,SF:EF,1);
sensor(condition).CDK2=D.dataAlignedSmooth(:,SF:EF,2);
sensor(condition).CDK4=D.dataAlignedSmooth(:,SF:EF,3);
sensor(condition).CDK4_corrected=CDK4_corrected(:,SF:EF);
% sensor(condition).CDK4high=gating_1;
% sensor(condition).CDK4delay=gating_3;
% sensor(condition).CDK4low=gating_2;
% sensor(condition).CDK2high=gating_4;
% sensor(condition).CDK2delay=gating_6;
% sensor(condition).CDK2low=gating_5;
% sensor(condition).risetime_cdk4=risetime_cdk4-SF-1;
% sensor(condition).risetime_cdk2=risetime_cdk2-SF-1;
% sensor(condition).cdk4_low=cdk4_low;
% sensor(condition).cdk2_low=cdk2_low;
% sensor(condition).badtraces_cdk4=badtraces_cdk4;
% sensor(condition).badtraces_cdk2=badtraces_cdk2;



%%
cellLine='983B';
% cellLine='858';
savedir=[resultdir cellLine];
if ~exist(savedir,'dir')
    mkdir(savedir);
end
%%


save([savedir filesep 'Data_mitogen_replace.mat'],'sensor','CDK4high','CDK4low');
%%
load([savedir filesep 'Data.mat'],'sensor','CDK4high','CDK4low');
name={'Control';'resistance';'Control_V'};


%% CDK2 classification
gating_5=logical(sensor(i).CDK2low)';
gating_6=(~isnan(sensor(i).risetime_cdk2)' & nanmean(sensor(i).CDK2(:,100:120)')<0.7)';
gating_4=~gating_2 & ~gating_3 & (nanmean(sensor(i).CDK2(:,130:200)')>0.75)';%nansum(D.dataAlignedSmooth(:,5000:5020,3)')<14 & nansum(D.dataAlignedSmooth(:,5060:5100,2)')>14 & nansum(D.dataAlignedSmooth(:,5100:5200,2)')>75 & nansum(D.dataAlignedSmooth(:,5100:5200,3)')>50;

%% CDK4 classification
% gating_1=nanmean(CDK4_corrected(:,5000:5020,1)')>0.5 & nanmean(CDK4_corrected(:,5030:5060,1)')>0.55 & nanmean(CDK4_corrected(:,5030:5060,1)')>0.6 & nanmean(CDK4_corrected(:,5030:5100,1)')>0.65;%& nanmean(D.dataAlignedSmooth(:,5060:5100,3)')>18;
% gating_2=nanmean(CDK4_corrected(:,5000:5020,1)')<0.6 & nanmean(CDK4_corrected(:,5030:5060,1)')<0.6 & nanmean(CDK4_corrected(:,5030:5100,1)')<0.7 & nanmean(CDK4_corrected(:,5030:5100,1)')<0.75;% | gate_2;%D.dataAlignedSmooth(:,5040,2)<0.7 & D.dataAlignedSmooth(:,5060,2)<0.7 & D.dataAlignedSmooth(:,5080,2)<0.7;
risingOpt=1; condnum=length(sensor);
for i=1:condnum; %(mean(nanmean(sensor(i).CDK4(:,50:100)'))-0.4)
    tracestats=nan(size(sensor(i).CDK4,1),2);
    
    for j=1:size(sensor(i).CDK4,1)
        tracestats(j,1) = find(sum(~isnan(sensor(i).CDK4_corrected(j,101:end)),1) > 0, 1 ,'first')+100;
        tracestats(j,2) = find(sum(~isnan(sensor(i).CDK4_corrected(j,101:end)),1) > 0, 1 , 'last')+100;
        if tracestats(j,2)>150
            tracestats(j,2)-30;
        else
            tracestats(j,2)-20;
        end
        %tracestats(:,1)=41-sum(~isnan(sensor(i).CDK4(:,101:140,1)),2)+100;
        %tracestats(:,2)=sum(~isnan(sensor(i).CDK4(:,101:end,1)),2)+99;
    end
    
    if risingOpt==1
        %tracestats(:,1)=101;
        %tracestats(:,2)=sum(~isnan(sensor(i).CDK4(:,101:end,1)),2)+100-30;
        [risetime_cdk4 , cdk4_low , badtraces_cdk4] = getCdk2rise_yilin(sensor(i).CDK4_corrected,tracestats);
        % [risetime_cdk2 , cdk2_low , badtraces_cdk2] = getCdk2rise_yilin(D.dataAligned(:,:,2) , tracestats);
    else
        
        minlength=40; Drugoption=0; quiescentanalysis=0;
        [risetime_cdk4,badtraces_cdk4]=getCdk2features_drug(sensor(i).CDK4_corrected,tracestats,minlength,Drugoption,quiescentanalysis);
        cdk4_low=isnan(risetime_cdk4) & ~badtraces_cdk4;
        %[sensor(i).POI_cdk2,sensor(i).badtraces_signal2]=getCdk2features_steve_drug(sensor(i).cdk2_smooth,sensor(i).tracestats,minlength,Drugoption,quiescentanalysis);
        %[sensor(i).POI_geminin,sensor(i).badtraces_signal1]=getGemininfeatures_1(sensor(i).geminin,sensor(i).tracestats,50);
    end
    sensor(i).CDK4low=cdk4_low;
    sensor(i).risetime_cdk4=risetime_cdk4;
    
    gating_2=logical(sensor(i).CDK4low) & (nanmean(sensor(i).CDK4(:,101:120)')<0.65)' & (nanmean(sensor(i).CDK4(:,130:200)')<0.65)';
    % gating_3=(~isnan(sensor(i).risetime_cdk4)' & nanmean(sensor(i).CDK4(:,101:120)')<(mean(nanmean(sensor(i).CDK4(:,50:100)'))-0.4))' & (nanmean(sensor(i).CDK4(:,130:200)')>(mean(nanmean(sensor(i).CDK4(:,50:100)'))-0.3))';
    %gating_1=~gating_2 & ~gating_3 & (nanmean(sensor(i).CDK4(:,101:120)')>(mean(nanmean(sensor(i).CDK4(:,50:100)'))-0.4))';% & (nanmean(sensor(i).CDK4(:,130:200)')>0.8)';%nansum(D.dataAlignedSmooth(:,5000:5020,3)')<14 & nansum(D.dataAlignedSmooth(:,5060:5100,2)')>14 & nansum(D.dataAlignedSmooth(:,5100:5200,2)')>75 & nansum(D.dataAlignedSmooth(:,5100:5200,3)')>50;
    gating_1=~gating_2  & (nanmean(sensor(i).CDK4(:,101:120)')>0.55)' & (nanmean(sensor(i).CDK4(:,120:140)')>0.65)' & (nanmean(sensor(i).CDK4(:,140:200)')>0.7)'; %& ~sensor(i).badtraces_cdk4
    gating_1_2=(nanmean(sensor(i).CDK4(:,180:120)')>0.7)' | isnan(nanmean(sensor(i).CDK4_corrected(:,180:200)'));
    gating_3=gating_1 & gating_1_2;
    %gating_1=~gating_2  & (nanmean(sensor(i).CDK4(:,101:120)')>0.65)'; %& ~sensor(i).badtraces_cdk4
    %gating_3=gating_1;
    condition=i;
    sensor(condition).CDK4high=gating_1;
    %sensor(condition).CDK4delay=gating_3;
    sensor(condition).CDK4low=gating_2;
    %sensor(condition).CDK2high=gating_4;
    %sensor(condition).CDK2delay=gating_6;
    %sensor(condition).CDK2low=gating_5;
end
%%
figure; hold on;
i=1;
for numCell=9
    
    cellplot=find(~isnan(sensor(i).risetime_cdk4));
    cell_to_plot=cellplot(numCell);
    plot(sensor(i).CDK4_corrected(cell_to_plot,:)');
    
    scatter(sensor(i).risetime_cdk4(cell_to_plot),sensor(i).CDK4_corrected(cell_to_plot,sensor(i).risetime_cdk4(cell_to_plot)),50,'k','filled');
end
%%

for i=1:condnum
    CDK4high(i)=sum(sensor(i).CDK4high)/(sum(sensor(i).CDK4high)+sum(sensor(i).CDK4low));
    %     CDK4delay(i)=sum(sensor(i).CDK4delay)/(sum(sensor(i).CDK4high)+sum(sensor(i).CDK4delay)+sum(sensor(i).CDK4low))
    CDK4low(i)=sum(sensor(i).CDK4low)/(sum(sensor(i).CDK4high)+sum(sensor(i).CDK4low));
end
CDK4high
CDK4low
%%
condition=0;
figure; hold on;
for i=1:condnum
    condition=condition+1;
    subplot(1,condnum,condition),hold on
    plot(nanmean(sensor(i).CDK4_corrected(sensor(i).CDK4high | sensor(i).CDK4delay,100:end)),nanmean(sensor(i).CDK2(sensor(i).CDK4high | sensor(i).CDK4delay,100:end)),'g')
    axis([0.2 0.85 0.2 1.5]); xlabel('CDK4 activity'); ylabel('CDK2 activity');
    %  plot(nanmean(sensor(i).CDK4_corrected(sensor(i).CDK4delay,101:end)),nanmean(sensor(i).CDK2(sensor(i).CDK4delay,101:end)),'b')
    %   axis([0.2 0.85 0.2 1.5]); xlabel('CDK4 activity'); ylabel('CDK2 activity');
    plot(nanmean(sensor(i).CDK4_corrected(sensor(i).CDK4low,100:end)),nanmean(sensor(i).CDK2(sensor(i).CDK4low,100:end)),'r')
    axis([0.2 0.85 0.2 1.6]); xlabel('CDK4 activity'); ylabel('CDK2 activity');
    title(['control CDK4high=' num2str(CDK4high(i)) ' CDK4low=' num2str(CDK4low(i))]);
    legend({'CDK4high','CDK4low'});
end
figu(1,1);
% export_fig([resultdir filesep 'PhasePlot.eps'],'-eps','-transparent','-nocrop');
%%
condition=0;
figure;
name={'Control','resistance','Control+V'};
for i=[1:2]
    subplot(1,2,i),hold on
    %plot(nanmean(sensor(i).CDK4_corrected(sensor(i).CDK4high,100:end)),nanmean(sensor(i).CDK2(sensor(i).CDK4high,100:end)),'k')
    plot_colorcode(nanmean(sensor(i).CDK4_corrected(sensor(i).CDK4high,100:end)),nanmean(sensor(i).CDK2(sensor(i).CDK4high,100:end)),1:102,3);
    %axis([0.2 0.85 0.2 1.5]); xlabel('CDK4 activity'); ylabel('CDK2 activity');
    %axis([0.2 0.85 0.2 1.5]); xlabel('CDK4 activity'); ylabel('CDK2 activity');
    %title(['control CDK4high=' num2str(CDK4high(i)) ' CDK4low=' num2str(CDK4low(i))]);
    
    %subplot(1,2,2),hold on
    plot_colorcode(nanmean(sensor(i).CDK4_corrected(sensor(i).CDK4low,100:end)),nanmean(sensor(i).CDK2(sensor(i).CDK4low,100:end)),1:102,3)
    %plot(nanmean(sensor(i).CDK4_corrected(sensor(i).CDK4low,100:end)),nanmean(sensor(i).CDK2(sensor(i).CDK4low,100:end)),'Color',colors(i,:))
    axis([0.2 0.85 0.2 1.6]); xlabel('CDK4 activity'); ylabel('CDK2 activity');
    title(['control CDK4high=' num2str(CDK4high(i)) ' CDK4low=' num2str(CDK4low(i))]);
    %legend({'Control','resistance','Control+V'});
    title(name{i})
end
figu(0.25,0.4);
% export_fig([savedir filesep 'PhasePlot.eps'],'-eps','-transparent','-nocrop');
%%
ci = 0.95; alpha = 1 - ci;
figure;
for i=1:condnum;
    ERK_mean_1=nanmean(sensor(i).ERK_corrected(sensor(i).CDK4high,:));
    Shigh=nanstd(sensor(i).ERK_corrected(sensor(i).CDK4high,:));
    Nhigh=sum(~isnan(sensor(i).ERK_corrected(sensor(i).CDK4high,:)));%sum(obj.cellState);
    T_multiplier = tinv(1-alpha/2, Nhigh-1);
    ci95high = T_multiplier.*Shigh./sqrt(Nhigh);
    
    %ERK_median_2=nanmean(sensor(i).ERK_corrected(sensor(i).CDK4delay,:));
    ERK_median_3=nanmean(sensor(i).ERK_corrected(sensor(i).CDK4low,:));
    Slow=nanstd(sensor(i).ERK_corrected(sensor(i).CDK4low,:));
    Nlow=sum(~isnan(sensor(i).ERK_corrected(sensor(i).CDK4low,:)));%sum(obj.cellState);
    T_multiplier = tinv(1-alpha/2, Nlow-1);
    ci95low = T_multiplier.*Slow./sqrt(Nlow);
    
    subplot(1,condnum,i),hold on
    plot(1:201,ERK_mean_1,'g'); plot(1:201,ERK_mean_1-ci95high,'g'); plot(1:201,ERK_mean_1+ci95high,'g');
    plot(1:201,ERK_median_3,'r'); plot(1:201,ERK_median_3-ci95low,'r'); plot(1:201,ERK_median_3+ci95low,'r');
    %shadedErrorBar(1:201,ERK_median_1,ci95high,'g',1);shadedErrorBar(1:201,ERK_median_3,ci95low,'r',1);
    ylim([0.3 1.2]); xlim([50 200]); vline([95 105]); %plot(ERK_median_2,'Color',[0.6 0.6 0]);
    %set(hline2,'color',[0 0.63 0],'linewidth',10);
    set(gca,'YAxisLocation','left','YColor',[0 0 0],'YTick',0.3:0.3:1.2,'fontsize',24,'tickdir','out','linewidth',5);
    set(gca,'Xtick',[50:25:200],'XTickLabel',[-10:5:20],'fontsize',24,'tickdir','out','linewidth',5)
    ylabel('ERK activity','FontSize',30); xlabel('Time since mitosis (hr)','FontSize',30);
    title(['control CDK4high=' num2str(CDK4high(i)) ' CDK4low=' num2str(CDK4low(i))]);
    %     legend({'CDK4high','CDK4low'});
end
figu(1,1)
% export_fig([savedir filesep 'ERKPlot (control, resistance).eps'],'-eps','-transparent','-nocrop');
%%
numToPlot=10;
condition=0; figure;
for i=1:condnum
    numCDK4high=round(CDK4high(i)*numToPlot);
    numCDK4low=round(CDK4low(i)*numToPlot);
    
    CDK4highIdx=find((sensor(i).CDK4high)==1);
    CDK4lowIdx=find(sensor(i).CDK4low==1);
    
    cell_plot_CDK4high=randsample(CDK4highIdx,numCDK4high);
    cell_plot_CDK4low=randsample(CDK4lowIdx,numCDK4low);
    
    condition=condition+1;
    subplot(condnum,3,condition); hold on;
    plot(sensor(i).ERK_corrected(cell_plot_CDK4high,:)','g');
    plot(sensor(i).ERK_corrected(cell_plot_CDK4low,:)','r'); axis([50 200 0 2]);
    set(gca,'YAxisLocation','left','YColor',[0 0 0],'YTick',0:0.5:2.1,'fontsize',24,'tickdir','out','linewidth',5);
    set(gca,'Xtick',[50:25:200],'XTickLabel',[-10:5:20],'fontsize',24,'tickdir','out','linewidth',5)
    ylabel('ERK activity','FontSize',30); xlabel('Time since mitosis (hr)','FontSize',30);
    %     legend({'CDK4high','CDK4low'});
    
    condition=condition+1;
    subplot(condnum,3,condition); hold on;
    plot(sensor(i).CDK2(cell_plot_CDK4high,:)','g');
    plot(sensor(i).CDK2(cell_plot_CDK4low,:)','r'); axis([50 200 0 3]);
    set(gca,'YAxisLocation','left','YColor',[0 0 0],'YTick',0:0.5:3,'fontsize',24,'tickdir','out','linewidth',5);
    set(gca,'Xtick',[50:25:200],'XTickLabel',[-10:5:20],'fontsize',24,'tickdir','out','linewidth',5)
    ylabel('CDK2 activity','FontSize',30); xlabel('Time since mitosis (hr)','FontSize',30);
    %     legend({'CDK4high','CDK4low'});
    
    condition=condition+1;
    subplot(condnum,3,condition); hold on;
    plot(sensor(i).CDK4_corrected(cell_plot_CDK4high,:)','g');
    plot(sensor(i).CDK4_corrected(cell_plot_CDK4low,:)','r'); axis([50 200 0 1.5]);
    set(gca,'YAxisLocation','left','YColor',[0 0 0],'YTick',0:0.5:1.5,'fontsize',24,'tickdir','out','linewidth',5);
    set(gca,'Xtick',[50:25:200],'XTickLabel',[-10:5:20],'fontsize',24,'tickdir','out','linewidth',5)
    ylabel('CDK4 activity','FontSize',30); xlabel('Time since mitosis (hr)','FontSize',30);
    title(['control CDK4high=' num2str(CDK4high(i)) ' CDK4low=' num2str(CDK4low(i))]);
    %     legend({'CDK4high','CDK4low'});
end
figu(1,1);
% export_fig([savedir filesep 'Traces (Mitogen replace).eps'],'-eps','-transparent','-nocrop');








%%
cellLine='983B';
% cellLine='858';
savedir=[resultdir filesep cellLine];
name={'Control+Vem', 'Resistance+Vem'};
%%
num=2;
ERK_corrected{num}=ERK_corrected(:,drug_frame-99:drug_frame+100);
condnum=length(ERK_corrected);
%%
% 1.2, 1.5

figure;
for i=1:2
    subplot(1,condnum,i);
    imagesc(ERK_corrected{i},[0. 1.2]); colorbar; vline([100]);
    set(gca,'Xtick',[0:25:300],'XTickLabel',[0:5:300],'fontsize',24,'tickdir','out','linewidth',2)
    ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title(['ERK activity; ' name{i}],'fontsize',30);
end
colormap(parula_gray);
figu(1,1);
export_fig([savedir filesep 'ERK_Heatmap (control, resistance).eps'],'-eps','-transparent','-nocrop');
%%
ci = 0.95; alpha = 1 - ci;
figure;
for i=1:2
    ERK_mean_1=nanmean(ERK_corrected{i});
    S=nanstd(ERK_corrected{i});
    N=sum(~isnan(ERK_corrected{i}));%sum(obj.cellState);
    T_multiplier = tinv(1-alpha/2, N-1);
    ci95 = T_multiplier.*S./sqrt(N);
    
    subplot(1,condnum,i),hold on
    plot(ERK_mean_1,'g'); plot(ERK_mean_1-ci95,'g'); plot(ERK_mean_1+ci95,'g');
    
    ylim([0.35 2.5]); xlim([0 200]); vline([100]);%vline([95 105]); %plot(ERK_median_2,'Color',[0.6 0.6 0]);
    set(gca,'Xtick',[0:25:200],'XTickLabel',[0:5:300],'fontsize',24,'tickdir','out','linewidth',2)
    ylabel('Cells','FontSize',30);xlabel('Time (hr)','FontSize',30);title(['ERK activity; ' name{i}],'fontsize',30);
end
figu(1,1)
export_fig([savedir filesep 'ERK_Traces (control, resistance).eps'],'-eps','-transparent','-nocrop');

%%















figure;
gating_1=CDK2(:,end)>1;
subplot(2,3,1); plot(CDK2(gating_1,:)'); subplot(2,3,2); plot(CDK4_corrected(gating_1,:)'); subplot(2,3,3); plot(ERK_corrected(gating_1,:)'); ylim([-0.4 1]);
subplot(2,3,4); plot(CDK2(~gating_1,:)'); subplot(2,3,5); plot(CDK4_corrected(~gating_1,:)'); subplot(2,3,6); plot(ERK_corrected(~gating_1,:)'); ylim([-0.4 1]);



%%
% clear nucVScyto; clear Nuc; clear Nuc_2; clear Cyto_1; clear Cyto_2;
%%
D.Cdk2Classification %parameter for classification

%% #3 Define rising point of signal and classify CDK2ind ,CDKdelay, and CDKlow
D.judgecdk2_detect_rising_pts(2,0)
D.judgeCDK2IncByTimeAfterMitosis(5,1,2); %time,smoothing Option
D.plot_compare_risingpoint
% export_fig([resultdir filesep 'CDK classification and ERK.eps'],'-eps','-transparent','-nocrop');
%% #2 Classify highERK and lowERK cells
% clf
% D=sort_or_align_by_mitosis([resultdir 'CONTROL.mat']);
D.judgeERKInc(5,10,2); %G2, percentile, cellcycle (1; before mitosis, 2; after mitosis, 3; both)
%% plot data
cdk2Inc=sum(D.cellState)/(sum(D.cellState)+sum(D.cellDelay)+sum(D.cellLow))
cdk2Delay=sum(D.cellDelay)/(sum(D.cellState)+sum(D.cellDelay)+sum(D.cellLow))
cdk2Low=sum(D.cellLow)/(sum(D.cellState)+sum(D.cellDelay)+sum(D.cellLow))
sample=100/(sum(D.cellState)+sum(D.cellDelay)+sum(D.cellLow))
%% smoothing dataSortedByMitosis
D.smoothDataSortedByMitosis
%% clf
gate='G2';
align=1; individual=1; smooth=1; sample=nan;
D.xlabels = {'time relative to mitosis (h)'}
D.timeInterval=1/30;
D.timeToPlotBeforeAfterMitosis = 25;
D.plot_timecourse_align(1,10,smooth,60,individual,1,gate,align,sample);
title([num2str(cdk2Inc) '% ' num2str(cdk2Delay) '% ' num2str(cdk2Low) '% ']);
%#1 plotType(1,;timelaps, 2;hist), #2 timepoint for hist (hr after mitosis), #3 smoothOptoion (1), #4 wellnumber (for SEM), #5 individual trace (1), Classification by CDK(1) or ERK(2)
% figu(0.5,1) %MitoWithdraw Control
% export_fig([resultdir filesep conditions{condi} '(G2).eps'],'-eps','-transparent','-nocrop');



%% Save figure
export_fig([resultdir filesep 'mitoWithdraw(G2 5hrs)-ave.eps'],'-eps','-transparent','-nocrop');
%% plot heatmap (aligned by mitosis)
clf
D.plot_heatmap_alignByMitosis
figu(0.5,1)
% clearvars -except D condi conditions resultdir
% exportfig(resultdir,['heatmap_align ',conditions{condi}]);



%%
sensor(condi).cdk2=D.dataSorted(:,:,1);
sensor(condi).cdk4=D.dataSorted(:,:,2);
sensor(condi).geminin=D.dataSorted(:,:,3);
sensor(condi).POI=D.divFrameSorted;

%%
condi=1;
sensor(condi).cdk2=D.dataSorted(:,1:100,1);
sensor(condi).cdk4=D.dataSorted(:,1:100,2);
sensor(condi).geminin=D.dataSorted(:,1:100,3);
sensor(condi).POI=D.divFrameSorted;



%%
condnum=length(sensor);
minlength=100; Drugoption=0; quiescentanalysis=1;
for i=1%:condnum;
    
    %sensor(i).geminin_area=sensor(i).geminin.*sensor(i).nucleus_area;
    %sensor(i).geminin_norm=normalizeMyTracesGeminin_alt2(sensor(i).geminin_area,0.01);
    sensor(i).geminin_norm=normalizeMyTracesGeminin_alt2(sensor(i).geminin,0.01);
    sensor(i).APC=E3_activity_from_trace(sensor(i).geminin_norm,0.01);
    
    sensor(i).cdk4_smooth=[];sensor(i).cdk2_smooth=[];
    for j=1:size(sensor(i).cdk2,1)
        sensor(i).cdk2_smooth(j,:)=nansmooth(sensor(i).cdk2(j,:),5);
        sensor(i).cdk4_smooth(j,:)=nansmooth(sensor(i).cdk4(j,:),5);
        
        sensor(i).tracestats(j,1) = find(sum(~isnan(sensor(i).cdk2(j,:)),1) > 0, 1 ,'first');
        sensor(i).tracestats(j,2) = find(sum(~isnan(sensor(i).cdk2(j,:)),1) > 0, 1 , 'last');
    end
    %sensor(i).cdk4_corrected=cdk4_activity_correction_allcells(sensor(i).cdk4_smooth,sensor(i).cdk2_smooth);
    
    sensor(i).cdk4_median=nanmedian(sensor(i).cdk4_smooth);
    %sensor(i).cdk4_corrected_median=nanmedian(sensor(i).cdk4_corrected);
    sensor(i).cdk2_median=nanmedian(sensor(i).cdk2_smooth);
    [sensor(i).POI_cdk4,sensor(i).badtraces_signal4]=getCdk2features_steve_drug(sensor(i).cdk4_smooth,sensor(i).tracestats,minlength,Drugoption,quiescentanalysis);
    [sensor(i).POI_cdk2,sensor(i).badtraces_signal2]=getCdk2features_steve_drug(sensor(i).cdk2_smooth,sensor(i).tracestats,minlength,Drugoption,quiescentanalysis);
    [sensor(i).POI_geminin,sensor(i).badtraces_signal1]=getGemininfeatures_1(sensor(i).geminin,sensor(i).tracestats,50);
end


%% Cell selection
plotCond=1;
xaxis=1:100;
for cell_to_plot=200:220;
    figure(cell_to_plot); hold on;
    plot(xaxis,sensor(plotCond).cdk2_smooth(cell_to_plot,:),'k','linewidth',2);
    [haxes,hline1,hline2] = plotyy(xaxis,sensor(plotCond).cdk4_smooth(cell_to_plot,:),xaxis,sensor(plotCond).geminin(cell_to_plot,:));
    axes(haxes(1));%axis([0 numframes 0.2 2]);
    set(gca,'YAxisLocation','left','YColor',[0.83 0 0],'YTick',0:0.5:2.5,'fontsize',16,'tickdir','out','linewidth',1.5);
    set(hline1,'color','r','linewidth',2);
    %ylim([0 1]); xlim([0 numframes]);
    %h=vline(frame_drug_added,'k:');
    %title(['Cell ',num2str(randcells(kk))]);
    %     axis([minX,maxX,minY,maxY]);
    xlim([0 100]); ylim([0 2.5])
    set(gca,'Xtick',[0:25:201],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
    ylabel('CDK2 and CDK4 activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('a single CDK2 trace','fontsize',18);
    
    if ~isnan(sensor(plotCond).POI_cdk4(cell_to_plot)) && ~isnan(sensor(plotCond).POI_cdk2(cell_to_plot))
        scatter(sensor(plotCond).POI_cdk4(cell_to_plot),sensor(plotCond).cdk4_smooth(cell_to_plot,sensor(plotCond).POI_cdk4(cell_to_plot)),100,'r','filled');
        text([sensor(plotCond).POI_cdk4(cell_to_plot)+5],[sensor(plotCond).cdk4_smooth(cell_to_plot,sensor(plotCond).POI_cdk4(cell_to_plot))],{'\leftarrow CDK4 activation'},'color',[1 0 0],'fontsize',18);
        scatter(sensor(plotCond).POI_cdk2(cell_to_plot),sensor(plotCond).cdk2_smooth(cell_to_plot,sensor(plotCond).POI_cdk2(cell_to_plot)),100,'k','filled');
        text([sensor(plotCond).POI_cdk2(cell_to_plot)+5],[sensor(plotCond).cdk2_smooth(cell_to_plot,sensor(plotCond).POI_cdk2(cell_to_plot))],{'\leftarrow CDK2 activation'},'color',[0 0 0],'fontsize',18);
        
    end
    
    axes(haxes(2));%axis([0 numframes 0 1]);hold on
    xlim([0 100]); ylim([0 1000]);  hold on;
    set(hline2,'color',[0 0.63 0],'linewidth',2);
    set(gca,'YAxisLocation','right','YColor',[0 0.63 0],'fontsize',16,'tickdir','out','linewidth',1.5);
    ylabel('APC activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('a single CDK2 trace','fontsize',18);
    %set(gca,'YAxisLocation','right','YColor',[0 0.63 0],'YTick',0:0.2:1);
    if ~isnan(sensor(plotCond).POI_geminin(cell_to_plot))
        scatter(sensor(plotCond).POI_geminin(cell_to_plot),sensor(plotCond).geminin(cell_to_plot,sensor(plotCond).POI_geminin(cell_to_plot)),100,'filled','MarkerEdgeColor',[0 0.63 0],'MarkerFaceColor',[0 0.63 0]);
        text([sensor(plotCond).POI_geminin(cell_to_plot)+5],[sensor(plotCond).geminin(cell_to_plot,sensor(plotCond).POI_geminin(cell_to_plot))],{'\leftarrow APC inactivation'},'color',[0 0.63 0],'fontsize',18);
    end
end
%%
for plotCond=1
    gating_1=~isnan(sensor(plotCond).POI_cdk4);
    % subplot(2,6,plotCond); hold on;
    subplot(2,1,plotCond); hold on;
    plot(sensor(plotCond).cdk4_smooth(gating_1,:)','k');
    plot(sensor(plotCond).cdk2_smooth(gating_1,:)','r');
    ylim([0 2]); xlim([0 100]); title(conditions{plotCond})
    
    subplot(2,1,plotCond+1); hold on;
    plot(sensor(plotCond).geminin(gating_1,:)','k');
    
    scatter(sensor(plotCond).POI_geminin(gating_1),sensor(plotCond).geminin(gating_1,sensor(plotCond).POI_geminin(gating_1))',100,'k','filled');
    % plot(nanmedian(sensor(plotCond).cdk4_smooth(gating_1,:)),'k');
    % plot(nanmedian(sensor(plotCond).cdk4_smooth(~gating_1,:)),'r');
    % plot(sensor(plotCond).cdk2_median,'r');
    ylim([-50 2000]); xlim([0 100]); title(conditions{plotCond})
    
end
%% Plot Cdk2 traces
figure(99)
CDK2_low=CDK2(:,NumOfFrame)<1 & CDK2(:,140)<0.8 & CDK2(:,100)<0.8;

subplot(1,2,1); hold on
plot(xaxis,CDK2(CDK2_low,:),'r');
plot(xaxis,CDK2(~CDK2_low,:),'g');
subplot(1,2,2); hold on
plot(xaxis,ERK(CDK2_low,:),'r');
plot(xaxis,ERK(~CDK2_low,:),'g');

plot(xaxis,nanmean(ERK(CDK2_low,:)),'r');
plot(xaxis,nanmean(ERK(~CDK2_low,:)),'k');

figure(100)
early_CDK2on=POI_rising<80;
subplot(1,2,1); hold on
plot(xaxis,CDK2(early_CDK2on,:),'r');
plot(xaxis,CDK2(CDK2_low,:),'g');
subplot(1,2,2); hold on
plot(xaxis,ERK(early_CDK2on,:),'r');
plot(xaxis,ERK(CDK2_low,:),'g');

plot(xaxis,nanmean(ERK(early_CDK2on,:)),'r');
plot(xaxis,nanmean(ERK(CDK2_low,:)),'k');




%%
save([resultdir filesep 'Data.mat'],'sensor','conditions');
%%
root='E:\2. Data\2. Development of CDK4 sensor\2017-12-19 (MCF10A; Geminin-CDK2-CDK4-H2B; Serum starvation)\';
resultdir=[root 'Results\'];
load([resultdir filesep 'Data.mat'],'sensor','conditions');
condnum=length(sensor);
%%
for plotCond=1
    gating_1=~isnan(sensor(plotCond).POI_cdk4);
    % subplot(2,6,plotCond); hold on;
    subplot(2,1,plotCond); hold on;
    plot(sensor(plotCond).cdk4_smooth(gating_1,:)','k');
    plot(sensor(plotCond).cdk2_smooth(gating_1,:)','r');
    ylim([0 2]); xlim([0 100]); title(conditions{plotCond})
    
    subplot(2,1,plotCond+1); hold on;
    plot(sensor(plotCond).geminin(gating_1,:)','k');
    % plot(nanmedian(sensor(plotCond).cdk4_smooth(gating_1,:)),'k');
    % plot(nanmedian(sensor(plotCond).cdk4_smooth(~gating_1,:)),'r');
    % plot(sensor(plotCond).cdk2_median,'r');
    ylim([-50 2000]); xlim([0 100]); title(conditions{plotCond})
end



%%
figure; hold on;
for i=1:condnum;
    subplot(2,6,i); hold on;
    plot(sensor(i).cdk4_median,'k');
    plot(sensor(i).cdk2_median,'r');
    ylim([0 2]); xlim([0 100]); title(conditions{i})
end

figure; hold on;
for i=1:condnum;
    subplot(2,6,i); hold on;
    plot(sensor(i).cdk4_median(1:100),sensor(i).cdk2_median(1:100));
    axis([0 2 0 2]); title(conditions{i})
end
%%
cond=[1 5:8]
figure; hold on;
for i=cond
    cdfplot_2(sensor(i).POI_cdk4)
    ylim([0 1])
end
legend(conditions{cond})
%%
cond=[1 8:9]
for i=cond
    figure(i); hold on;
    cdfplot_2(sensor(i).POI_cdk4)
    cdfplot_2(sensor(i).POI_cdk2)
    ylim([0 1])
end
legend({'Cdk4'; 'Cdk2'})
%%
figure;
plotCond=12;
subplot(1,3,1); hold on;
plot(sensor(plotCond).cdk2_smooth'); xlim([1 100]);
prcCellcycle(plotCond)=sum(sensor(plotCond).cdk2_smooth(:,100)>0.8)/size(sensor(plotCond).cdk2_smooth,1);
title(num2str(prcCellcycle(plotCond)))
subplot(1,3,2); hold on;
plot(sensor(plotCond).cdk4_smooth');xlim([1 100]);
prcCellcycle(plotCond)=sum(sensor(plotCond).cdk4_smooth(:,100)>0.8)/size(sensor(plotCond).cdk4_smooth,1);
title(num2str(prcCellcycle(plotCond)))
subplot(1,3,3); hold on;
plot(sensor(plotCond).geminin');xlim([1 100]);
prcCellcycle(plotCond)=sum(sensor(plotCond).geminin(:,100)>0.8)/size(sensor(plotCond).geminin,1);
title(num2str(prcCellcycle(plotCond)))
%%
plotCond=10;
gating_1=sensor(plotCond).cdk4_smooth(:,50)>1 & nanmean(sensor(plotCond).geminin(:,1:30)')'<100;
gating_2=nanmean(sensor(plotCond).geminin(:,100-20:100)')>100;
gating_3=gating_1 & gating_2';
gating_4=gating_1 & ~gating_2';

figure; hold on;
subplot(1,3,1); hold on;
plot(sensor(plotCond).cdk2_smooth(gating_3,:)','r'); xlim([1 100]);
plot(sensor(plotCond).cdk2_smooth(gating_4,:)','k'); xlim([1 100]);

subplot(1,3,2); hold on;
plot(sensor(plotCond).cdk4_smooth(gating_3,:)','r'); xlim([1 100]);
plot(sensor(plotCond).cdk4_smooth(gating_4,:)','k'); xlim([1 100]);

subplot(1,3,3); hold on;
plot(sensor(plotCond).geminin(gating_3,:)','r'); xlim([1 100]);
plot(sensor(plotCond).geminin(gating_4,:)','k'); xlim([1 100]);

%%
prcCellcycle(plotCond)=sum(sensor(plotCond).cdk2_smooth(:,100)>0.8)/size(sensor(plotCond).cdk2_smooth,1);
title(num2str(prcCellcycle(plotCond)))
subplot(1,3,2); hold on;
plot(sensor(plotCond).cdk4_smooth');xlim([1 100]);
prcCellcycle(plotCond)=sum(sensor(plotCond).cdk4_smooth(:,100)>0.8)/size(sensor(plotCond).cdk4_smooth,1);
title(num2str(prcCellcycle(plotCond)))
subplot(1,3,3); hold on;
plot(sensor(plotCond).geminin');xlim([1 100]);
prcCellcycle(plotCond)=sum(sensor(plotCond).geminin(:,100)>0.8)/size(sensor(plotCond).geminin,1);
title(num2str(prcCellcycle(plotCond)))


%%
cond=3;
numframes=size(sensor(1).cdk2,2);
frame_drug_added=0;
startcell=1;
endcell=size(sensor(cond).geminin_norm,1);%
endcell=startcell+40;
xaxis=1:numframes;
numcells=100;
randcells=randsample(size(sensor(i).APC,1),numcells)';
for kk=1:numcells%startcell:endcell;%size(geminin_goodcells,1);
    figure(ceil((kk-startcell+1)/20));hold on
    subaxis(4,5,mod(kk-1,20)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);hold on
    [haxes,hline1,hline2] = plotyy(xaxis,sensor(cond).geminin_norm(randcells(kk),:),xaxis,sensor(cond).cdk4(randcells(kk),:));
    hold on;
    axes(haxes(1));%axis([0 numframes 0.2 2]);
    set(gca,'YAxisLocation','left','YColor',[0.83 0 0],'YTick',0:0.2:1);
    set(hline1,'color','r','linewidth',1.1);
    ylim([0 1]); xlim([0 numframes]);
    h=vline(frame_drug_added,'k:');
    title(['Cell ',num2str(randcells(kk))]);
    
    axes(haxes(2));%axis([0 numframes 0 1]);hold on
    ylim([0 3]);%max(sensor(cond).cdk4_corrected(randcells(kk),:))]);
    set(gca,'YAxisLocation','right','YColor',[0 0.63 0]);
    %set(gca,'YAxisLocation','right','YColor',[0 0.63 0],'YTick',0:0.2:1);
    xlim([0 numframes]);    %ylim([0 1]);
    set(hline2,'color','g','linewidth',1.1);
end









%%
%Input, matrix for Geminin and CDK2 time course

APC=Geminin(:,:);
TimeInterval=0.2;  %time interval in hours
APCthresh=120;  %Degron amplitude threshold used to find cells that inactivate APC for first time
APCactivityThr= -.2; %log10 APC activity threshold to determine time point for APC inactivation

minReporterLevel=0.02; %minimal Geminin Degron reporter level to prevent negative expression values
minAPCact=.004; %minimal APC activity to prevent negative activity values

close all
a0=length(CDK2(:,1)); % number of time points
c0=length(CDK2(1,:)); % number of cells tracked
a(1:a0)=0;

for i=1:a0
    APC(i,2:c0-1)=0.33333*(APC(i,1:c0-2)+APC(i,2:c0-1)+APC(i,3:c0));%Averaging three time points
    CDK2(i,2:c0-1)=0.33333*(CDK2(i,1:c0-2)+CDK2(i,2:c0-1)+CDK2(i,3:c0));
    k=APC(i,40:322)>APCthresh; %only checks for a degron increase 8 hours into the stimulation protocol
    for j=1:260
        kk(j)=sum(k(j:j+19)); %finds the first time point where 20 APC values are greater than maxval
    end
    aa=find(kk==20,1,'first')+39;
    if aa<285 & APC(i,30)>0 & APC(i,aa+20)>0
        a(i)=aa;
    end
end
indAPC=find(a>0);

b(1:a0)=0; %Determines maximal slope
for i=indAPC
    x=(1:9)*TimeInterval; xav=mean(x);
    y=APC(i,a(i)+5:a(i)+13); yav=mean(y);
    b(i)=sum((x-xav).*(y-yav))/sum((x-xav).^2); %linear regression to measure maximal degron slope
end

sel=find(b<10);   %eliminates cells that have shallow slopes, likely reflecting cells with low expression of degron reporter
indAPC=setdiff(indAPC,sel);

for i=indAPC
    APC(i,:)=APC(i,:)/b(i);    %normalizes degron levels to the maximal slope after APC inactivation
    APC(i,APC(i,:)<minReporterLevel)=minReporterLevel; %ensures that reporter levels do not go negative
end
%figure,plot(APC(indAPC,:)') %control plot to see how corrected APC data looks like

clear Slope
Slope(1:a0,1:c0)=0;  %calculates slopes of the normalized reporter for all time points
for i=indAPC
    for j=3:c0-2
        x=(1:5)*TimeInterval; xav=mean(x);
        y=APC(i,j-2:j+2); yav=mean(y);
        Slope(i,j)=sum((x-xav).*(y-yav))/sum((x-xav).^2); %linear regression to measure slope as function of time
    end
end

%---------------------------------------------------------------------------
%Calculates APC activity parameter
APCact(1:a0,1:c0)=0;
alog(1:a0)=0;
figure,hold on
sel=[];
for i=indAPC
    ind=find(~isnan(Slope(i,1:a(i)+20))); %only analyzes time points that have real Slope values
    
    APCact(i,ind)=(1-Slope(i,ind))./APC(i,ind); %Main equation to derive APC activity
    
    APCact(i,APCact(i,:)<minAPCact)=minAPCact; %ensures that APC activity is not negative
    plot(ind*TimeInterval,log10(APCact(i,ind)))
    
    %Threshold for APCact is used to call the time for APC inactivation
    aa=find((log10(APCact(i,a(i)-20:a(i)+8))+log10(APCact(i,a(i)-19:a(i)+9))+log10(APCact(i,a(i)-18:a(i)+10)))<3*APCactivityThr,1,'first')+a(i)-20;
    if ~isnan(aa) %selects cells that pass below threshold, marked by line in plot
        sel=[sel i];
        alog(i)=aa;
    end
end
line([0 60],[-.2 -.2],'Color','k')
title('Degradation rate of degron reporter (APC activity). Black line marks time of inactivation')
axis([0 300/5 -1.8 1.5])
xlabel('Time after mitogen stimulation (hours)')
ylabel('APC acivity, Log10 (1/hour)')

indAPC=sel;  %limits analyzed cells to those that had APC inactivated below threshold

%--------------------------------------------------------------------------
%Some plots to evaluate APC inactivation and the corresponding CDK2 activation

%Next plot shows the log10 APC activity aligned by the APC inactivation time
figure,hold on
for i=indAPC
    if alog(i)>0
        plot(([1:30]-21)*TimeInterval,log10(APCact(i,alog(i)-20:alog(i)+9)))
    end
    
end
line([-21/5 6-21/5],[-.2 -.2],'Color','k')
axis([-21/5 6-21/5 -1.8 1.5])
title('APC activity changes aligned by time of APC inactivation')
xlabel('Time before APC inactivation (hour)')
ylabel('APC activity, Log10 (1/hour)')

%Next plot shows the derived half-life of the degron reporter
figure,hold on
for i=indAPC
    if alog(i)>0
        plot(([1:21]-19)*TimeInterval,1./(APCact(i,alog(i)-18:alog(i)+2))*log(2))
    end
end
title('Half-life of degron reporter from ~10 min in G1 to >30 hours in S/G2')
xlabel('Time before APC inactivation (hour)')
ylabel('Derived APC degron reporter half-life (hour)')
axis([-18/5 5/5 0 1.5])

%Next plot shows the degron reporter expression change aligned by APC
%inactivation time (normalized by slope but otherwise raw data)
n=0;
figure,hold on
z=[];
for i=indAPC
    if alog(i)>0
        if mean(APC(i,alog(i)-18:alog(i)-6))<0.3
            plot(([1:21]-19)*TimeInterval,APC(i,alog(i)-18:alog(i)+2), 'b')
            n=n+1;
            z(n,1:21)=APC(i,alog(i)-18:alog(i)+2);
        end
    end
end
title('Rel. increase in degron reporter intensity aligned by APC inactivation')
xlabel('Time before APC inactivation (hour)')
ylabel('APC degron reporter intensity values, mean in red')
axis([-18/5 2/5 0 .8])
hold on
plot(([1:21]-19)*TimeInterval,mean(z,1),'Linewidth', 3,'Color','r')

%next plot shows the CDK2 activity reporter data co-aligned by the time of
%APC inactivation, red shows the relative APC reporter increase shifted by 4
%time points to correct for folding time
figure,hold on
for i=indAPC
    if alog(i)-30>0
        plot(([1:19]-19)*TimeInterval,(CDK2(i,alog(i)-22:alog(i)-4)))
    end
end
axis([-19/5 0 0.2 1])
ylabel('CDK2 activity')
xlabel('Time before APC inactivation (hour)')
title('CDK2 activity aligned by APC inactivation & corrected for folding, red marks APC reporter increase')
plot(([1:19]-19)*TimeInterval,1.4*mean(z(:,1:19),1),'Linewidth', 3,'Color','r')

%next is a phase plot analysis comparing the increase in CDK2 and APC
%reporter, shifted by folding time
clear x y
x(1:a0,1:23)=0;
y(1:a0,1:23)=0;
isel=[];
for i=indAPC
    if alog(i)>40 & sum(isnan(CDK2(i,alog(i)-20:alog(i)-4)))==0 & sum(isnan(log10(APCact(i,alog(i)-16:alog(i)))))==0
        isel=[isel i];
        x(i,1:23)=CDK2(i,alog(i)-20:alog(i)+2);
        y(i,1:23)=log10(APCact(i,alog(i)-16:alog(i)+6));
    end
end
indAPC=isel;

figure,hold on
for i=indAPC
    plot(x(i,6:23),y(i,6:23))
end
axis([0 1.5 -2.3 1])

xlabel('CDK2 activity')
ylabel('APC activity, Log10 (1/hour)')
title('Initial APC inactivation by CDK2 before steep inactivation switch is triggered')
hold on,plot(mean(x(indAPC,6:22),1),mean(y(indAPC,6:22),1),'Linewidth', 3,'Color','r')
%%
APC=Geminin(:,:);
TimeInterval=0.2;  %time interval in hours
maxval=120;  %initial minimal Degron amplitude parameter to find cells that inactivated APC

minReporterLevel=0.1; %minimal Geminin Degron reporter level to prevent negative values
minAPCact=.004; %minimal APC activity to prevent negative values
RemAPC=.004; %estimate of residual APC activity in S-phase, assumes a 50 hour life-time of Degron reporter

close all
a0=length(CDK4(:,1));
c0=length(CDK4(1,:));
a(1:a0)=0;

for i=1:a0
    APC(i,2:c0-1)=0.33333*(APC(i,1:c0-2)+APC(i,2:c0-1)+APC(i,3:c0));%Averaging three time points
    CDK4(i,2:c0-1)=0.3333333*(CDK4(i,1:c0-2)+CDK4(i,2:c0-1)+CDK4(i,3:c0));
    CDK2(i,3:c0-2)=0.2*(CDK2(i,1:c0-4)+CDK2(i,2:c0-3)+CDK2(i,3:c0-2)+CDK2(i,4:c0-1)+CDK2(i,5:c0));
    k=APC(i,40:322)>maxval;
    for j=1:260
        kk(j)=sum(k(j:j+19));
    end
    aa=find(kk==20,1,'first')+39;
    if aa<285 & APC(i,30)>0 & APC(i,aa+20)>0
        a(i)=aa;
    end
end
indAPC=find(a>0);

b(1:a0)=0;
sel=[];
for i=indAPC
    x=1:9; xav=mean(x);
    y=APC(i,a(i)+2:a(i)+10); yav=mean(y);
    b(i)=sum((x-xav).*(y-yav))/sum((x-xav).^2); %linear regression to measure maximal degron slope
end
%figure,hist(b,0:1:50)

sel=find(b<10);   %eliminates cells that express low levels of degron reporter
indAPC=setdiff(indAPC,sel);

for i=indAPC
    APC(i,:)=APC(i,:)/b(i);    %normalizes degron level to maximal slope after APC inactivation
    APC(i,APC(i,:)<minReporterLevel)=minReporterLevel; %makes sure that reporter level is not negative
end
figure,plot(APC(indAPC,:)')

clear Slope
Slope(1:a0,1:322)=0;
for i=indAPC
    for j=3:322-2
        x=1:5; xav=mean(x);
        y=APC(i,j-2:j+2); yav=mean(y);
        Slope(i,j)=sum((x-xav).*(y-yav))/sum((x-xav).^2); %linear regression to measure slope as function of time
    end
end

APCact(1:a0,1:322)=0;
alog(1:a0)=0;
figure,hold on
sel=[];
for i=indAPC
    ind=find(~isnan(Slope(i,1:a(i)+20)));
    APCact(i,ind)=(1-Slope(i,ind)+RemAPC)./APC(i,ind)/TimeInterval; %Main equation for APC activity
    APCact(i,APCact(i,:)<minAPCact)=minAPCact; %ensures that APC activity is not negative
    plot(ind*TimeInterval,log10(APCact(i,ind)))
    aa=find((log10(APCact(i,a(i)-20:a(i)+8))+log10(APCact(i,a(i)-19:a(i)+9))+log10(APCact(i,a(i)-18:a(i)+10)))<-.6,1,'first')+a(i)-20;
    if ~isnan(aa) & ~isnan(APCact(i,a(i)-20))
        sel=[sel i];
        alog(i)=aa;
    end
end
indAPC=sel;
axis([0 300/5 -1.8 1.5])
xlabel('Time (hours)')
ylabel('Rate of degradation by APC, Log10 (1/hours)')

figure,hold on
for i=indAPC
    if alog(i)>0
        plot([1:30]*TimeInterval,log10(APCact(i,alog(i)-20:alog(i)+9)))
    end
    
end
axis([0 6 -1.8 1.5])

figure,hold on
for i=indAPC
    if alog(i)>0
        plot([1:27]*TimeInterval,(CDK2(i,alog(i)-30:alog(i)-4)))
    end
end
axis([0 27/5 0 2])

figure,hold on
for i=indAPC
    if alog(i)>0
        plot([1:27]*TimeInterval,(CDK4(i,alog(i)-30:alog(i)-4)))
    end
end
axis([0 27/5 0 2])

figure,hold on
for i=indAPC
    if alog(i)>0
        plot([1:27]*TimeInterval,log10(APCact(i,alog(i)-26:alog(i))))
    end
end
axis([0 27/5 0 2])

clear x y z
x(1:a0,1:37)=0;
y(1:a0,1:37)=0;
z(1:a0,1:37)=0;
figure,hold on
isel=[];
for i=indAPC
    if alog(i)>40 & sum(isnan(CDK2(i,alog(i)-40:alog(i)-4)))==0 & sum(isnan(log10(APCact(i,alog(i)-36:alog(i)))))==0
        isel=[isel i];
        x(i,1:37)=CDK2(i,alog(i)-40:alog(i)-4);
        z(i,1:37)=CDK4(i,alog(i)-40:alog(i)-4)-0.2*CDK2(i,alog(i)-40:alog(i)-4);
        y(i,1:37)=log10(APCact(i,alog(i)-36:alog(i)));
        plot(CDK2(i,alog(i)-30:alog(i)-4),log10(APCact(i,alog(i)-26:alog(i))))
    end
end
indAPC=isel;
axis([0 1.5 -2 2])
figure,plot(mean(x(indAPC,1:37),1),mean(y(indAPC,1:37),1),'b-')
xlabel('CDK2 activity')
ylabel('log10 initial APC inactivation by CDK2 before switch is triggered')

indAPC1=indAPC(alog(indAPC)<80);
indAPC2=indAPC(alog(indAPC)>110);


figure, plot([-36:0]*TimeInterval,mean(0.5*(y(indAPC1,:)+1.1),1),'k-')
hold on, plot([-36:0]*TimeInterval,mean(0.5*(y(indAPC2,:)+1.1),1),'k--')

hold on,plot([-36:0]*TimeInterval,mean(x(indAPC1,:),1),'r-')
hold on,plot([-36:0]*TimeInterval,mean(x(indAPC2,:),1),'r--')

hold on,plot([-36:0]*TimeInterval,mean(z(indAPC1,:),1),'b-')
hold on,plot([-36:0]*TimeInterval,mean(z(indAPC2,:),1),'b--')

axis([-37/5 0 0.2 1.2])
xlabel('hour before APC inactivation, corrected for YFP folding 4x12 minutes')
ylabel('rel amplitude APC, CDK4, CDK2')
title('CDK4 blue, CDK2 red, log10 & shifted APC black, solid early APC inact <16 h, dashed late >22 h')























checking=0;
n=0; n1=0;
if condi>5
    j=condi-5;
else
    j=condi;
end
for i=1:size(D.dataAlignedSmooth,1)
    if sum(isnan(D.dataAlignedSmooth(i,4920:4997,1)))>10
        continue;
    end
    
    %av=.2*(D.dataAlignedSmooth(i,4941:4994,1) + D.dataAlignedSmooth(i,4942:4995,1) + D.dataAlignedSmooth(i,4943:4996,1) + D.dataAlignedSmooth(i,4944:4997,1) + D.dataAlignedSmooth(i,4945:4998,1));
    val=(D.dataAlignedSmooth(i,4920:4994,1));
    [p1 vin]=max(val);
    ind=find(val>0);
    p2=val(ind(1));
    normVal=(val-p2)/(p1-p2);
    [c index] = min(abs(normVal-0.5));
    
    if length(normVal)<index+10
        %nucVScyto(n,1:length(normVal(index-10:end)))=normVal(index-10:end);
    else
        if nanmean(normVal(index+5:index+7))>0.65 && nanmean(normVal(index-5:index-3))<0.4 && (index-10)>0
            n=n+1;
            nucVScyto{j}(n,1:21)=normVal(index-10:index+10);
            slop=diff(nucVScyto{j}(n,1:21),1,2);
            tempRising=slop>0.02;
            rising = find(tempRising==1,1,'first');
            risingTime{j}(n)=10-rising;
        end
    end
    
    clear normVal;
    val=(D.dataAlignedSmooth(i,4920:5060,2));
    [p1 vin]=max(val);
    ind=find(val>0);
    p2=val(ind(1));
    normVal=(val-p2)/(p1-p2);
    [c, index] = min(abs(normVal-0.5));
    if length(normVal)<index+20
        %Nuc(n,1:length(normVal(index-10:end)))=normVal(index-10:end);
    else
        if index-10>0
            if nanmean(normVal(index-10:index-5))<0.4 & nanmean(normVal(index+5:index+7))>0.65 & nanmean(normVal(index-8:index-5))<0.6
                n1=n1+1;
                if index+60>size(normVal,2)
                    Nuc{j}(n1,:)=nan(1,71);
                    Nuc{j}(n1,1:size(normVal(index-10:end),2))=normVal(index-10:end);
                else
                    Nuc{j}(n1,1:71)=normVal(index-10:index+60);
                end
                
                slop=diff(Nuc{j}(n1,:),1,2);
                tempRising=slop>0.02;
                risingPoint{j}(n1) = find(tempRising==1,1,'first');
                timeMax{j}(n1)=find(Nuc{j}(n1,:)==1);
                
                temp=[];
                if risingPoint{j}(n1)-10<=0
                    addV=nan(1,abs(risingPoint{j}(n1)-10));
                    temp=[addV Nuc{j}(n1,1:(risingPoint{j}(n1)+60+risingPoint{j}(n1)-10))];
                else
                    temp=Nuc{j}(n1,risingPoint{j}(n1)-10:risingPoint{j}(n1)+60);
                end
                Nuc_2{j}(n1,1:size(temp,2))=temp; %Nuc2; aligned by rising time
                
                %             if n1==25
                %                 keyboard;
                %             end
                %smoothSlop=smooth(slop);
                tempDecreasing=slop<-0.02;
                decreasingPoint{j}(n1) = find(tempDecreasing==1,1,'first')+1;
                
                temp=[];
                if decreasingPoint{j}(n1)-20<=0
                    addV=nan(1,abs(decreasingPoint{j}(n1)-20));
                    temp=[addV Nuc{j}(n1,decreasingPoint{j}(n1)-19+length(addV):decreasingPoint{j}(n1)+30)];
                elseif decreasingPoint{j}(n1)+30>length(Nuc{j}(n1,:))
                    addV=nan(1,abs(decreasingPoint{j}(n1)+30-length(Nuc{j}(n1,:))));
                    temp=[Nuc{j}(n1,decreasingPoint{j}(n1)-20:decreasingPoint{j}(n1)+30-length(addV)) addV];
                else
                    temp=Nuc{j}(n1,decreasingPoint{j}(n1)-20:decreasingPoint{j}(n1)+30);
                end
                Nuc_3{j}(n1,1:size(temp,2))=temp; %Nuc2; aligned by decreasing time
                
                inx=1;
                while slop(decreasingPoint{j}(n1)+3)>-0.02 | slop(decreasingPoint{j}(n1)+4)>-0.02 %| slop(decreasingPoint(n1)+7)>-0.02
                    inx=inx+1;
                    if inx >= 66
                        break
                    end
                    decreasingPoint{j}(n1) = max(find(tempDecreasing==1,inx))+1;
                end
                
                Cyto_1{j}(n1,:)=(D.dataAlignedSmooth(i,4700:5060,2));
                risingPointForCyto{j}(n1)=risingPoint{j}(n1)+index-10+(4920-4700);
                
                if checking
                    figure(n1);
                    subplot(2,1,1); hold on;
                    plot(Nuc(n1,:),'k');
                    scatter(risingPoint{j}(n1),Nuc(n1,risingPoint{j}(n1)),'r','filled');
                    scatter(decreasingPoint{j}(n1),Nuc(n1,decreasingPoint{j}(n1)),'g','filled');
                    
                    subplot(2,1,2); hold on;
                    plot(Cyto_1{j}(n1,:),'k');
                    scatter(risingPointForCyto{j}(n1),Cyto_1{j}(n1,risingPointForCyto{j}(n1)),'r');
                end
                
                temp=nan(1,350);
                alignPoint=350-risingPointForCyto{j}(n1);
                temp(1+alignPoint:risingPointForCyto{j}(n1)+alignPoint)=Cyto_1{j}(n1,1:risingPointForCyto{j}(n1));
                Cyto_2{j}(n1,:)=temp;
                
            end
        end
    end
    
    %     plot(1:41,Nuc(n1,:),'-p','MarkerIndices',[risingPoint(n1) decreasingPoint(n1)],...
    %     'MarkerFaceColor','red','MarkerSize',15)
    
    %     kk=find(av>0.5);
    %     avnorm=av(kk(1)-11:kk(1)+8);
    %     av0=mean(avnorm(1:3));
    %     avnorm=(avnorm-av0)/(1-av0);
    %     vvv(n,1:23,1)=avnorm;
    % figure,plot(1:23,avnorm,'b')
    %
    %
    % av=.2*(MEKi_Cyto(i,4949:5008) + MEKi_Cyto(i,4950:5009) + MEKi_Cyto(i,4951:5010) + MEKi_Cyto(i,4952:5011) + MEKi_Cyto(i,4953:5012));
    %     [p1 vin]=max(av);
    %     kk=find(av>0);
    %     p2=av(kk(1));
    %     av=(MEKi_Cyto(i,4951:5010)-p2)/(p1-p2);
    %     kk=find(av>0.5);
    %     avnorm=av(kk(1)-11:kk(1)+11);
    %     av0=mean(avnorm(1:3));
    %     avnorm=(avnorm-av0)/(1-av0);
    %     vvv(n,1:23,2)=avnorm;
    % hold on,plot(1:23,avnorm,'r')
    % axis([0 25 -.1 1.1])
end
%%
figure;
for i=1:4;%1:4
    %     if condi>5
    %         j=i-5;
    %     else
    j=i;
    %     end
    subplot(1,4,j); hold on;
    Cyto_cyclinB{j}=nanmedian(Cyto_2{j});
    plot(Cyto_2{j}','k'); plot(Cyto_cyclinB{j},'r','Linewidth',5);
    axis([100 350 0 1500]); xlabel('time relative to rising point (hr)'); ylabel('Cyto cyclin B1 (FI)');
    set(gca,'XTick',[110:30:350]); set(gca,'XTickLabel',[-8:1:0]);
    title(conditions{j});
end
figu(0.5,1);
%%
figure; hold on;
colorMap={'k','r','g','b'};
for i=1:4
    plot(Cyto_cyclinB{i},colorMap{i},'Linewidth',5);
    Cyto_slope{i}=diff(smooth(Cyto_cyclinB{i}(260:350))',1,2);
    median_Cyto_slop(i)=nanmedian(Cyto_slope{i})
end
axis([100 350 0 700]); xlabel('time relative to rising point (hr)'); ylabel('Cyto cyclin B1 (FI)');
set(gca,'XTick',[110:30:350]); set(gca,'XTickLabel',[-8:1:0]);
figu(1,1); legend(conditions);
%%
figure; hold on;
for condi=1:4;
    % clear time_start_to_anaphase;
    timeInterval=2;
    time_start_to_anaphase{condi}=decreasingPoint{condi}-risingPoint{condi};
    time_start_to_max{condi}=timeMax{condi}-risingPoint{condi};
    % histogram(time_start_to_anaphase{condi})
    histogram(time_start_to_anaphase{condi})
    
    % figure,plot(1:23,mean(vvv(1:n,:,1),1),'b');
    % hold on , plot(1:23,mean(vvv(1:n,:,2),1),'r');
    % axis([0 25 -.1 1.1])
end
legend(conditions);

condition=0; figure;
for i=1:4;
    condition=condition+1;
    Ratio_cyclinB{i}=nanmedian(nucVScyto{i});
    subplot(4,3,condition); hold on;
    plot(nucVScyto{i}','k'); plot(Ratio_cyclinB{i},'r','linewidth',4)
    axis([1 21 0 1]); xlabel('time relative to half maximum (min)'); ylabel('Norm Nuc/Cyto');
    set(gca,'XTick',[1:2:21]); set(gca,'XTickLabel',[-20:4:20]);
    
    condition=condition+1;
    subplot(4,3,condition); hold on;
    Nuc_cyclinB{i}=nanmedian(Nuc_2{i});
    
    [val ind]=max(Nuc_cyclinB{i});
    Nuc_slope_1{i}=diff(Nuc_cyclinB{i}(10:ind),1,2);
    median_Nuc_slop_1(i)=nanmedian(Nuc_slope_1{i})
    
    plot(Nuc_2{i}','k'); plot(Nuc_cyclinB{i},'r','linewidth',4) %Nuc_2{condi}; time relative to rising point (min), Nuc; time relative to half maximum (min)
    scatter(10,Nuc_cyclinB{i}(10),'b','filled');
    scatter(ind,Nuc_cyclinB{i}(ind),'g','filled');
    axis([5 61 -0.1 1]); xlabel('time relative to rising point (min)'); ylabel('Norm Nuc');
    set(gca,'XTick',[5:5:61]); set(gca,'XTickLabel',[-10:10:100]);
    title(conditions{i});
    
    Nuc_cyclinB_2{i}=nanmedian(Nuc_3{i});
    idx_1=0;%find(diff(Nuc_cyclinB{i}(20:25),1,2)<0);
    idx_2=find(diff(Nuc_cyclinB_2{i}(20:30),1,2)>-0.03);
    Nuc_slope_2{i}=diff(Nuc_cyclinB_2{i}(idx_1(1)+20:idx_2(1)+25),1,2);
    median_Nuc_slop_2(i)=nanmedian(Nuc_slope_2{i})
    
    condition=condition+1;
    subplot(4,3,condition); hold on;
    plot(Nuc_3{i}','k'); plot(Nuc_cyclinB_2{i},'r','linewidth',4) %Nuc_2{condi}; time relative to rising point (min), Nuc; time relative to half maximum (min)
    scatter(idx_1(1)+20,Nuc_cyclinB_2{i}(idx_1(1)+20),'b','filled');
    scatter(idx_2(1)+20,Nuc_cyclinB_2{i}(idx_2(1)+20),'g','filled');
    axis([1 31 -0.1 1]); xlabel('time relative to rising point (min)'); ylabel('Norm Nuc');
    set(gca,'XTick',[0:5:31]); set(gca,'XTickLabel',[-40:10:20]);
end
figu(1,1);
% ratio=nucVScyto; nuc_intensity=Nuc_2{condi};
% save([resultdir filesep conditions{condi} '_20170316.mat'],'timeInterval','time_start_to_max','time_start_to_anaphase','nuc_intensity','ratio');
% export_fig([resultdir filesep conditions{condi} '_CyclinB translocation.eps'],'-eps','-transparent','-nocrop');
% export_fig([resultdir filesep conditions{condi} '_CyclinB translocation.png'],'-nocrop');
%%
colorMap={'k','r','g','b'};

figure; hold on;
for i=1:4
    plot(Ratio_cyclinB{i},colorMap{i},'Linewidth',5);
    axis([1 21 0 1]); xlabel('time relative to half maximum (min)'); ylabel('Norm Nuc/Cyto');
    set(gca,'XTick',[1:2:21]); set(gca,'XTickLabel',[-20:4:20]);
end
legend(conditions);

figure; hold on;
legend(conditions);
for i=1:4
    plot(Nuc_cyclinB{i},colorMap{i},'Linewidth',5);
    axis([5 61 -0.1 1]); xlabel('time relative to rising point (min)'); ylabel('Norm Nuc');
    set(gca,'XTick',[5:5:61]); set(gca,'XTickLabel',[-10:10:100]);
end
legend(conditions);

figure; hold on;
for i=1:4
    plot(Nuc_cyclinB_2{i},colorMap{i},'Linewidth',5);
    axis([1 31 -0.1 1]); xlabel('time relative to rising point (min)'); ylabel('Norm Nuc');
    set(gca,'XTick',[0:5:31]); set(gca,'XTickLabel',[-40:10:20]);
end
legend(conditions);
%%
mkdir(resultdir);
save([resultdir 'Data.mat'],'nucVScyto','Ratio_cyclinB','Cyto_2','Cyto_cyclinB','conditions','Cyto_slope','median_Cyto_slop','time_start_to_anaphase',.....
    'risingPoint','decreasingPoint','Nuc_2','Nuc_cyclinB','median_Nuc_slop_1','Nuc_3','Nuc_cyclinB_2','median_Nuc_slop_2');
%%





















%%
clear time_M;
time_M=decreasingPoint-risingPoint;
%figure; histogram(time_M*2)
% figure,plot(1:23,mean(vvv(1:n,:,1),1),'b');
% hold on , plot(1:23,mean(vvv(1:n,:,2),1),'r');
% axis([0 25 -.1 1.1])
figure; subplot(1,2,1); hold on;
plot(nucVScyto{condi}','k'); plot(nanmedian(nucVScyto{condi}),'r','linewidth',4)
axis([1 21 0 1]); xlabel('time relative to half maximum (min)'); ylabel('Norm Nuc/Cyto');
set(gca,'XTick',[1:2:21]); set(gca,'XTickLabel',[-10:2:10]);

subplot(1,2,2); hold on;
plot(Nuc_2{condi}','k'); plot(nanmedian(Nuc_2{condi}),'r','linewidth',4) %Nuc_2{condi}; time relative to rising point (min), Nuc; time relative to half maximum (min)
axis([6 61 -0.1 1]); xlabel('time relative to rising point (min)'); ylabel('Norm Nuc');
set(gca,'XTick',[6:2:61]); set(gca,'XTickLabel',[-8:4:100]);
figu(0.5,1);

% save([resultdir filesep conditions{condi} '_20170330.mat'],'risingTime','time_M','nucVScyto','Nuc_2');
% export_fig([resultdir filesep conditions{condi} '_CyclinB translocation.eps'],'-eps','-transparent','-nocrop');
% export_fig([resultdir filesep conditions{condi} '_CyclinB translocation.png'],'-nocrop');
%% Save data
clearvars -except D condi conditions resultdir gate
exportData(resultdir,conditions{condi});

%%
figure; hold on;
histogram(Control,'Normalization','probability');
histogram(Wee1,'Normalization','probability');
histogram(Cdc25,'Normalization','probability');
histogram(Aura,'Normalization','probability');




%% Section-1; SortByMitosis

%% plot heatmap
D.sortXTickTime = [5 10 15 20]; % show tick at 8h,16h,24h...........
D.ylabels = {'CDK2 activity','ERK activity'};
D.xlabels = {'time (h)'};
D.colorRange{1} = [0.5 2.2];
% D.colorRange{2} = [0.95 1.3];
% D.clean_remove_mean_outside_n_std(2,3); % 1st arg means we use second channel, 2nd arg means discard outside mean+-2.5*std
D.plot_heatmap_sortByMitosis(1)% input; smoothoption
% figu(0.5,0.5);
colormap(parula);
export_fig([resultdir,'heatmap_',conditions{condi},'.eps'],'-eps','-transparent','-nocrop');


%% CDK and ERK activity; drug added before or after mitosis
D.sortXTickTime = [0,8,16,24,32,40]; % show tick at 8h,16h,24h...........
D.ylabels = {'CDK2 activity','ERK activity'};
D.xlabels = {'time (h)'};
D.plot_timecourse_align_Drug(40,2,0); %#1 Time of drug addition, #2 1;Mitosis before drug, 2;Mitosis after drug, #3 1; individual, 2; average
ylim([1.2 1.8]);
%% Save figure
export_fig([resultdir filesep 'CDK and ERK activity (NoMitogen)_drug added before mitosis.eps'],'-eps','-transparent','-nocrop'); %1 after 2 before






%% Section-3 Box plot of ERK activity based on cell phase
clf
D.boxplot_median_erk(5, 10, 5); % duration of G2, S, G1

% does G2 duration include mitotic pulse as well? or, do you want to
% include it? Right now G2 is -8 to -2h relative mitosis. If you want to
% change this, we have to change a line 17 of the function.

%% plot cycling and quiescent cells separately... but there are not many quiescent cells in this condition
D.boxplot_median_erk_cyc_qui_separated(5, 10, 5); % duration of G2, S, G1
%% plot cycling and quiescent cells separately... but there are not many quiescent cells in this condition
D.bargraph_mean_erk_cyc_qui_separated(5, 10, 5); % duration of G2, S, G1
%% plot CDKfast, CDKslow, and CDKlow
D.boxplot_median_erk_separated(5,5); %duration of G2 and G1
%%
D.boxplot_mean_erk(5,5); %duration of G2 and G1






%% #2 Clear bad traces
D.clean_remove_by_lower_upper_limit(1, [0, 3.5]) % for channel 1, remove something higher than 3.5.
% D.clean_remove_by_lower_upper_limit(2, [0.4, 3]) % for channel 2, remove something lower than 0.5 and higher than 3.
%% #3 Clear bad traces by indicating channel, frame and range to remove (added 9/12/14)
D.clean_remove_dataaligned_indicate_timeandrange
switch gate
    case 'G1'
        D.clean_remove_dataaligned_indicate_timeandrange(1,4990:4993,[0 1.2])
        D.clean_remove_dataaligned_indicate_timeandrange(1,5010:5050,[1.6 4])
    case 'G2'
        D.clean_remove_dataaligned_indicate_timeandrange(1,4970:4980,[0 0.85])
        D.clean_remove_dataaligned_indicate_timeandrange(1,4980:4985,[0 1])
        D.clean_remove_dataaligned_indicate_timeandrange(1,5010:5025,[2 4])
end
%%
D.clean_remove_dataaligned_indicate_timeandrange(1,4990:4993,[0.5 1])
% D.clean_remove_dataaligned_indicate_timeandrange(1,5050:5055,[0.5 1])

%% Save data
clearvars -except D condi conditions resultdir
exportfig(resultdir,['_',conditions{condi}]);
%% Check individual traces (CDK2)
for i=[1051 1331]%900:1400%[558]%34]% 79]%81:100%length(D.dataSorted(:,:,1))
    if max(smooth(D.dataSorted(i,:,1)))<=D.colorRange{1}(2) && min(smooth(D.dataSorted(i,:,1)))>=D.colorRange{1}(1)
        figure(i);
        plot(smooth(D.dataSorted(i,:,1))); ylim(D.colorRange{1}); xlim([0 241]);
        %export_fig([resultdir,'MitogenWithdraw_CDK2_' num2str(i) '_plot'],'-eps','-transparent','-nocrop');
        imagesc(smooth(D.dataSorted(i,:,1))',D.colorRange{1});
        %export_fig([resultdir,'MitogenWithdraw_CDK2_' num2str(i) '_heatmap'],'-eps','-transparent','-nocrop');
    end
end
%% Check individual traces (ERK)
for i=907 %Mitogen withdraw %[60 92] %DMSO   %[802 869 1081 1097] %800:1200%[16 46]%1:50%[663 710 802 1182]%1:length(D.dataSorted(:,:,1));%[340 398]%310:510%100%length(D.dataSorted(:,:,1))
    if max(D.dataSorted(i,:,2))<=D.colorRange{2}(2) && min(D.dataSorted(i,:,2))>=D.colorRange{2}(1)
        figure(i);
        %plot(smooth(D.dataSorted(i,:,2))); ylim(D.colorRange{2});
        plot(D.dataSorted(i,:,2)); ylim(D.colorRange{2}); xlim([0 241]);
        %export_fig([resultdir,'MitogenWithdraw_ERK_' num2str(i) '_plot'],'-eps','-transparent','-nocrop');
        %imagesc(smooth(D.dataSorted(i,:,2))',D.colorRange{2});
        imagesc(D.dataSorted(i,:,2),D.colorRange{2});
        %export_fig([resultdir,'MitogenWithdraw_ERK_' num2str(i) '_heatmap'],'-eps','-transparent','-nocrop');
    end
end


%% Section-4; Calculate probability of quiescence
%% Load combined data
root='J:\2. Data\IXmicro (New Axon)\1. Project\1. ERK and Cell cycle\2014-05-16 (MCF10A; ERK-EV-nls, hDHB-mCherry; 12min, 48hr)\';
resultdir=[root 'Results\'];
load([resultdir 'CombinedData']);
%% Select condition of experiments
condi = 9;
condition=condi;
conditions = {'EGF20Serum5';'EGF10Serum2.5';'EGF5Serum1.25';'MitogenWithdraw';'DMSO';'MEKi';'CDK1_2i';'CDK2i';'CDK4i';'MEKiCDK1_2i';'MEKi CDK2i';'MEKi CDK4i'};
data = ret{condition}.store.storeData;
divFrame = ret{condition}.store.allParentDivFrame;
D = calculate_next_cell_cycle_entry_prob(data,divFrame) %Cell probability
D.clean_leastnum_frames(50);%(241); % this will remove frames shorter than 241.
D.clean_remove_suddendropCDK2
D.clean_remove_CDK2_not_dropped
D.clean_remove_by_lower_upper_limit(1, [-Inf, 3.5]) % for channel 1, remove something higher than 3.5.
%% Classify CDK2inc and CDK2low cells
D.clean_remove_CDK2_above_thres_at_mitosis(1.4)
D.judgeCDK2Inc(0.01) %parameter for classification
%% Classify CDK2inc and CDK2low cells (CDK2 activity at 10hr after mitosis < 1)
D.clean_remove_CDK2_above_thres_at_mitosis(1.4)
D.judgeCDK2IncByTimeAfterMitosis(10,1); %time,1 or 2
%% Calculate probability of quiescence
D.timeWindow = 2;
D.inhibitionTime = 8;
D.calc_prob;
D.plot_line_probability_cell_divide('MEKi',resultdir);