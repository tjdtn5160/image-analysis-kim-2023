%% Initializes and clears the workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
%% GENERAL SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='D:\1. Projects\5. Breast cancer project\';
experimentpath='2023-04-27 Extended EdU_Palbo,Fulv_2day\';
combinepath='1. Combine\EdU_MCF7,T47D,CAMA1_palbo,fulv 2day\';

datadir=[projectpath,experimentpath,'Data\'];
resultdir=[projectpath,experimentpath,'Results\'];
combinedir=[projectpath,combinepath,'Data\'];
if ~exist(resultdir,'dir')
    mkdir(resultdir)
end

%%
rows=1:8;
cols=1:6;%[3 4 7 8 11 12];%; %[1 2 5 6 9 10];

numrows=length(rows);
numcols=length(cols);
shots=numrows*numcols;

for idx=1:shots
    colidx=mod(ceil(idx),numcols);
    if colidx==0
        colidx=numcols;
    end
    col=cols(colidx);
    rowidx=ceil(idx/(numcols));
    row=rows(rowidx);

%     if row==1 & ismember(col,3:6)
%         sites=11:32;
%     else
        sites=1:32;
%     end

    numsites=length(sites);

    Alldata=[];
    for site=sites
        %Zeiss
        %shot=wellnum2str(row,col,site);
        shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
        
        %Nikon
        %[rowstring,colstring,sitestring]=wellnum2strRCS_3(row,col,site);
        %shot=[rowstring,'_',colstring,'_',sitestring];
        if exist([datadir,'IF_',shot,'.mat'])
            load([datadir,'IF_',shot,'.mat'],'IFdata');
            Alldata=[Alldata;IFdata];
        else
            continue
        end
    end
    if isempty(Alldata)
        continue
    end
    Hoechstval=Alldata(:,3).*Alldata(:,4);
%     if ismember(col,[5:6 11:12])
        posval_hoechst=Hoechstval>1000000 & Hoechstval<20000000;
%     elseif ismember(col,[])
%         posval_hoechst=Hoechstval>500000 & Hoechstval<10000000;
%     else
%         posval_hoechst=Hoechstval>500000 & Hoechstval<10000000;
%     end
    %figure;histogram(Hoechstval);vline(1000000);
    %keyboard;

%     Rb=Alldata(:,5);%./Alldata(:,protein);
%     posval_protein_1=Rb>=0; Rb(Rb<1)=1; Rb=log2(Rb);

    EdU_channel=5;
    
    EdUval=Alldata(:,EdU_channel);
    EdUval(EdUval<1)=1; EdUval=log2(EdUval);
    %figure;histogram(EdUval);vline(4);vline(14);
    posval_EdU=EdUval>1 & EdUval<16;    
    gating=posval_hoechst & posval_EdU; %& posval_protein_1;    
    
    numCell=sum(gating);
    numCellData(row,col)=numCell;    
    
%     Rb=Rb(gating);
%     RbRaw{row,col}=Rb; RbData(row,col)=nanmean(Rb);
        
%         %figure;histogram(pRb);vline(4);vline(14);
%         threshpRb=10;
%         numpRb=sum(pRb>threshpRb);
%         perpRb=numpRb/numCell*100;
        
%         pRbData(row,col)=perpRb;
    
    
%     if ismember(col,[2:7])
%         CDK2=CDK2(gating);
%         
%         CDK2Raw{row,col}=CDK2;
%         CDK2Data(row,col)=nanmean(CDK2);
%     end
    
    Hoechstval=Hoechstval(gating);
    HoechstRaw{row,col}=Hoechstval;
    
    EdUval=EdUval(gating);
    EdURaw{row,col}=EdUval;

%     if ismember(col,1:6) & row==1      
        threshEdU=8.5;
%     elseif ismember(col,[6 12]) & ismember(row,7:8)
%         threshEdU=5;
%     elseif col==12 & ismember(row,3:4)
%         threshEdU=5;
%     else
%         threshEdU=8;
%     end
    numS=sum(EdUval>threshEdU);
    perS=numS/numCell*100;
    
    SData(row,col)=perS;
    
    %CDK4=Alldata(:,8)./Alldata(:,7);
    %posval_CDK4=CDK4>=0;
    
    %CDK4_stain=Alldata(:,11);%./Alldata(:,protein);
    %posval_CDK4_stainig=CDK4_stain>=0; CDK4_stain(CDK4_stain<1)=1; CDK4_stain=log2(CDK4_stain);
    
    %CDK4=CDK4(gating);
    %CDK4_stain=CDK4_stain(gating);
    
    %keyboard;
    %%%%%% visualize Hoechst vs EdU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    figure;dscatter(Hoechstval,EdUval);
    figure;dscatter(CDK2,EdUval);
    figure;dscatter(Rb,EdUval);
    figure;dscatter(CDK4_stain,EdUval);
    figure;histogram(EdUval,100,'Normalization','probability'); vline(threshEdU); xlim([3 14.5]);
    %}
    
    
    
    
    %CDK4Raw{row,col}=CDK4_activity_correction_allcells_MCF7(CDK4,CDK2);
    %CDK4_stainRaw{row,col}=CDK4_stain;
    
    
    %CDK4_stainData(row,col)=nanmean(CDK4_stain);
    %CDK4Data(row,col)=nanmean(CDK4Raw{row,col});
    
end
%%
% conditions={'DMSO'	'Bafilomycin A1 (10然)'	'Leupeptin (10然)' 'Chloroquine (10然)' 'MG-132 (5 然)'	'Borezomib (100 nM)' 'Carfilzomib (100 nM)'	'Ixazomib (100 nM)'	'Rapamycin (100nM)'	'Everolimus (1 ?M)'	'AZD-8055 (100nM)'	'LY293002 (50 然)' 'Tram (100 nM)' 'MK2206' 'Calpe' 'PD150' 'QX77' 'Lys05'};
%%
save([resultdir 'Data_20230427'],'HoechstRaw','EdURaw','numCellData','SData');
save([combinedir 'Data_20230427'],'HoechstRaw','EdURaw','numCellData','SData');
%%
load([resultdir 'Data_20230412'],'HoechstRaw','EdURaw','numCellData','SData');


%%
yMax=20000;
figure;imagesc(numCellData(1:end,1:end),[0 yMax]); colorbar; title('Number of Cells');
set(gca,'Xtick',[1:12],'XTickLabel',[1:12],'fontsize',12,'tickdir','out','linewidth',2);
colormap('cyan_yellow');
%colormap('blue_yellow')
%%
figure;imagesc(SData(1:end,1:end),[0 80]); colorbar; title('%S-phase intensity');
%set(gca,'Xtick',[1:10],'XTickLabel',[1:6],'fontsize',12,'tickdir','out','linewidth',2);
colormap('cyan_yellow');
%colormap('blue_yellow')
%%
% figure;imagesc(CDK2Data(2:end,1:end),[0.2 0.6]); colorbar; title('CDK2 activity');
% set(gca,'Xtick',[1:12],'XTickLabel',[1:12],'fontsize',12,'tickdir','out','linewidth',2);
% colormap('cyan_yellow');
%colormap('blue_yellow')
%%
% yMax=14;
% 
% figure;imagesc(RbData(1:end,2:end),[8 yMax]); colorbar; title('Rb intensity'); colormap('cyan_yellow');
% %colormap('blue_yellow')
%set(gca,'Xtick',[1:12],'XTickLabel',[0 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 1 2],'fontsize',12,'tickdir','out','linewidth',2);
%set(gca,'Ytick',[1:8],'YTickLabel',{'Bortezomib','Ixazomib','Celastrol','Carfilzomib','Epoxomicin','ONX0914','Delanzomib','Oprozomib'},'fontsize',12,'tickdir','out','linewidth',2);
%%
% yMax=105;
% 
% figure;imagesc(pRbData(2:end,1:end),[8 yMax]); colorbar; title('pRb intensity'); colormap('cyan_yellow');
%colormap('blue_yellow')
%set(gca,'Xtick',[1:12],'XTickLabel',[0 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 1 2],'fontsize',12,'tickdir','out','linewidth',2);
%set(gca,'Ytick',[1:8],'YTickLabel',{'Bortezomib','Ixazomib','Celastrol','Carfilzomib','Epoxomicin','ONX0914','Delanzomib','Oprozomib'},'fontsize',12,'tickdir','out','linewidth',2);
%%
% yMax=13;
% 
% figure;imagesc(CDK4Data(6,1:end),[8 yMax]); colorbar; title('CDK4 intensity'); colormap('cyan_yellow');
%colormap('blue_yellow')
%set(gca,'Xtick',[1:12],'XTickLabel',[0 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 1 2],'fontsize',12,'tickdir','out','linewidth',2);
%set(gca,'Ytick',[1:8],'YTickLabel',{'Bortezomib','Ixazomib','Celastrol','Carfilzomib','Epoxomicin','ONX0914','Delanzomib','Oprozomib'},'fontsize',12,'tickdir','out','linewidth',2);
%%
% yMin=8;
% yMax=14;
% 
% figure;imagesc(CDK6Data(7,1:end),[yMin yMax]); colorbar; title('CDK6 intensity'); colormap('cyan_yellow');
%%
% yMin=8;
% yMax=14;
% 
% figure;imagesc(CycDData(8,1:end),[yMin yMax]); colorbar; title('CycD intensity'); colormap('cyan_yellow');
%%
% yMax=10;
% 
% figure;imagesc(RbData(2:end,1:end),[6 yMax]); colorbar; title('Rb intensity');
% colormap('cyan_yellow');
%%
% figure; hold on;

% histogram(RbRaw{1,8});
% histogram(RbRaw{3,8});
% histogram(RbRaw{6,5});
% histogram(RbRaw{3,2});
% histogram(RbRaw{5,2});
% histogram(RbRaw{2,8});
% histogram(RbRaw{3,8});
% histogram(RbRaw{5,8});

%%
% figure;imagesc(CDK4Data(2:end,2:end),[0.05 0.5]); colorbar; title('Number of Cells');
% set(gca,'Xtick',[1:12],'XTickLabel',[1:12],'fontsize',12,'tickdir','out','linewidth',2);
% colormap('cyan_yellow');
% %colormap('blue_yellow')



%%
% figure;imagesc(pRbData(2:end,2:end),[0 100]); colorbar; title('%pRB intensity');
% set(gca,'Xtick',[1:10],'XTickLabel',[1:6],'fontsize',12,'tickdir','out','linewidth',2);
% colormap('cyan_yellow');
% %colormap('blue_yellow')


%%
% figure;imagesc(numCellData(2:end,2:end),[0 10000]); colorbar; title('Number of Cells');
% set(gca,'Xtick',[1:12],'XTickLabel',[1:12],'fontsize',12,'tickdir','out','linewidth',2);
% colormap('cyan_yellow');
%colormap('blue_yellow')

%%
% figure;imagesc(CDK4_stainData(2:end,2:end)); colorbar; %title('Rb intensity');
% colormap('cyan_yellow');
%colormap('blue_yellow')
%set(gca,'Xtick',[1:12],'XTickLabel',[0 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 1 2],'fontsize',12,'tickdir','out','linewidth',2);
%set(gca,'Ytick',[1:8],'YTickLabel',{'Bortezomib','Ixazomib','Celastrol','Carfilzomib','Epoxomicin','ONX0914','Delanzomib','Oprozomib'},'fontsize',12,'tickdir','out','linewidth',2);
%%
% figure;imagesc(numCellData(2:end,2:end),[0 10000]); colorbar; title('Number of Cells');
% set(gca,'Xtick',[1:12],'XTickLabel',[1:12],'fontsize',12,'tickdir','out','linewidth',2);
% colormap('cyan_yellow');



%%
clear hillCoeff ec50
count=0;
%dose=[0.001 0.01 0.1 0.5 1 3 5 8 10 50 100 500];
dose=[10 30 100 300 1000 3000 6000 10000]; % nM
for i=[1:4] %
    count=count+1;
    doseresponse=SData(1:8,i);
    response{count}=reshape(doseresponse,[1 8]);
    figure; hold on;
    [hillCoeff ec50(count)]=doseResponse(dose,response{count});
    
%     com_data{count}(plateNum,1:9)=response;
%     if i==2
%         ylim([0 65])
%     end
end
xlabel('Dose (nM)'); ylabel('S-phase cells (%)');
%% Normalization
clear hillCoeff ec50
count=0;
%dose=[0.001 0.01 0.1 0.5 1 3 5 8 10 50 100 500];
dose=[10 30 100 300 1000 3000 6000 10000]; % nM
for i=[1:2:8] %
    count=count+1;
    doseresponse=mean(SData(1:8,i:i+1)')/mean(SData(1,i:i+1)');
    
    
    response{count}=reshape(doseresponse,[1 8]);
    figure; hold on; 
    [hillCoeff ec50(count)]=doseResponse(dose,response{count});
    ylim([0 1.2])
end
xlabel('Dose (nM)'); ylabel('S-phase cells (%)');

%% Normalization
clear hillCoeff EC50
count=0;
%dose=[0.001 0.01 0.1 0.5 1 3 5 8 10 50 100 500];
conc=[10 30 100 300 1000 3000 6000 10000]; % nM
for i=[1:2:11] %
    count=count+1;
    if i==11;
        responses=SData(:,i)';
    else
        responses=mean(SData(1:8,i:i+1)')/mean(SData(1,i:i+1)');
    end
%     responses{count}=reshape(doseresponse,[1 8]);
    figure; hold on; 
    [hillCoeff EC50(count)]=ec50(conc',responses');
    ylim([0 1.2])
end
xlabel('Dose (nM)'); ylabel('S-phase cells (%)');
%%
% dose=[0.01 0.03 0.1 0.3 1 3 6 10];
% for i=1:3
%     response=mean(com_data{i});
%     figure; hold on;
%     [hillCoeff ec50(i)]=doseResponse(dose,response);
%     if i==1
%         ylim([0 65])
%     end
% end
% xlabel('Dose (nM)'); ylabel('S-phase cells (%)');

%% IC50

figure; count=0;
for i=[1:2:12]   
    
        doseresponse=mean(SData(1:8,i:i+1)')/mean(SData(1,i:i+1)'); %normalization by control condition
        if ismember(i,1:8)
            conc=[10 30 100 300 1000 3000 6000 10000];
        else
            conc=[10 30 100 300 1000 3000 10000 30000];
        end
    
    %MAXpoint=find(doseresponse==max(doseresponse));
    %MINpoint=find(doseresponse==min(doseresponse));
    %decreaseS=doseresponse(MAXpoint:MINpoint);        
    %minNorm=min(decreaseS);
    %decreaseS=decreaseS-minNorm;
    %decreaseS=decreaseS/max(decreaseS);
    
    count=count+1;
    subplot(1,6,count); hold on;
    xvals=(conc);
    fitpara=fit(xvals',doseresponse','exp1');
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time
    
    plot(fitpara,xvals',doseresponse');hold on;
    
    %     set(gca,'Xtick',[0:100:lastFrame-fisrtFrame],'XTickLabel',[0:20:300],'fontsize',24,'tickdir','out','linewidth',2);
    %     set(gca,'Ytick',[0.3:0.3:1.5],'YTickLabel',[0.3:0.3:1.5],'fontsize',24,'tickdir','in','linewidth',5)
%     ylim([0 1.2]); xlim([0 end]);
    ylabel('Response','FontSize',20);xlabel('Palbociclib','FontSize',20);title(['IC50=' num2str(halfLife(i))],'fontsize',20);
    figu(0.5,1);
end

%%

conc=[10 30 100 300 1000 3000 6000 10000];
count=0;
% Response=[];
for i=[1:2:8]    
    if ismember(i,1:10)
        doseresponse=mean(SData(1:8,i:i+1)')/mean(SData(1,i:i+1)'); %normalization by control condition
    else
        doseresponse=SData(1:8,i)'/SData(1,i)';    
    end  
    Response{i}=doseresponse;
    count=count+1;    
    xvals=(conc);
    fitpara=fit(xvals',doseresponse','exp1');
    halfLife(i)=(log(2)/abs((fitpara.b))); %convert to time    
end
%%
save([resultdir 'Data_20210623'],'HoechstRaw','EdURaw','numCellData','SData','Response');


%%
figure;imagesc(RbData(4,2:end),[8 10]); colorbar; title('Rb intensity');
colormap('cyan_yellow');
%colormap('blue_yellow')
set(gca,'Xtick',[1:10],'XTickLabel',[0 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 1 2],'fontsize',12,'tickdir','out','linewidth',2);
%set(gca,'Ytick',[1:8],'YTickLabel',{'Bortezomib','Ixazomib','Celastrol','Carfilzomib','Epoxomicin','ONX0914','Delanzomib','Oprozomib'},'fontsize',12,'tickdir','out','linewidth',2);



%%
count=1;
for i=[2:7]
    count=count+1;
    xvals=[0.00001 0.001 5 10 50];
    val=plotV(i,2:6);
    % fitpara=fit(xvals',val','exp2');
    %     halfLife(i)=log(2)/abs((fitpara.b)); %convert to time
    figure; hold on;
    [hillCoeff ec50(count)]=doseResponse(xvals,val)
    
    
    %plot(fitpara,xvals',val');hold on;
end
%%
figure;imagesc(RbData(2:8,2:11),[7 11]); colorbar; title('Rb intensity');
colormap('cyan_yellow');
%colormap('blue_yellow')
%set(gca,'Xtick',[1:12],'XTickLabel',[0 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 1 2],'fontsize',12,'tickdir','out','linewidth',2);
%set(gca,'Ytick',[1:8],'YTickLabel',{'Bortezomib','Ixazomib','Celastrol','Carfilzomib','Epoxomicin','ONX0914','Delanzomib','Oprozomib'},'fontsize',12,'tickdir','out','linewidth',2);
%%
figure;
plotV=(SData(2:7,1:12)./SData(2:7,1));
xaxis=[0.001	0.01	0.05	0.08	0.1	0.3	0.5	0.8	1	3	5	10];
plot(xaxis,plotV')
%%
figure;
plotV=(SData([1:4,7,8],2:11)./SData([1:4,7,8],2));
xaxis=[0 0.001 5 10 50 100 500 1000 5000 10000];
plot([0 1 5 10],plotV(:,2:5)')
%%
figure;
plotV=(SData([4:5 6:7],1:12));
xaxis=[0.001	0.01	0.05	0.08	0.1	0.3	0.5	0.8	1	3	5	10];
plot(xaxis,plotV')
%%
figure;
plotV=(SData(:,2:11));
xaxis=[0 0.001 5 10 50 100 500 1000 5000 10000];
plot(log([0.1 1 5 10]),plotV(:,2:5)')
%%
figure; hold on;
data_1=RbRaw{2,2};
data_2=RbRaw{2,4};
data_3=RbRaw{2,6};
data_3=RbRaw{2,8};
cellNum_1=size(data_1,1);
cellNum_2=size(data_2,1);
cellNum_3=size(data_3,1);
histogram(data_1,100,'Normalization','probability');
histogram(data_2,100,'Normalization','probability');
histogram(data_3,100,'Normalization','probability');
%histogram(data_4,100,'Normalization','probability');
xlim([9 13]);
xlabel('Rb expression (log2)'); ylabel('Fraction')
legend({'resist','drug-naive','Rb KO'})
%%
figure; hold on;
histogram(RbRaw{8,2},100,'Normalization','probability');
histogram(RbRaw{8,6},100,'Normalization','probability');
histogram(RbRaw{5,2},100,'Normalization','probability');
xlim([6 13]);
%%
figure;dscatter(RbRaw{2,2},EdURaw{2,2})
%%
row=2;
col=2;

color=flipud(pink(1000));

RbData=RbRaw{row,col};
EdUData=EdURaw{row,col};
numCell=size(RbData,1);
gating=RbData>0 & EdUData>0;
RbData=RbData(gating);
EdUData=EdUData(gating);
figure;dscatter_gray(RbData,EdUData);%colormap('parula_white');
vline(10);hline(10);
pop_1=RbData>6 & RbData<10 & EdUData<10;
pop_2=RbData>6 & RbData<10 & EdUData>10;
pop_3=RbData>10 & EdUData<10;
pop_4=RbData>10 & EdUData>10;

numPop_1=sum(pop_1);
numPop_2=sum(pop_2);
numPop_3=sum(pop_3);
numPop_4=sum(pop_4);

total=numPop_1+numPop_2+numPop_3+numPop_4;
per_1=numPop_1/total %EdU low / Rb low
per_2=numPop_2/total %EdU high / Rb low
per_3=numPop_3/total %EdU low / Rb high
per_4=numPop_4/total %EdU high / Rb high
%ylim([0 2]);
% axis([7 13 4 15]);
xlabel('Rb expression (log2)'); ylabel('EdU (log2)'); colorbar;
%%
i=7;
%i=1;
color=flipud(pink(1000));

RbData=RbRaw{i,2};
HoechstData=HoechstRaw{i,2};
numCell=size(RbData,1);
gating=RbData>6 & HoechstData<6000000;
RbData=RbData(gating);
HoechstData=HoechstData(gating);
figure;dscatter_gray(RbData,HoechstData);%colormap('parula_white');
vline(10);hline(10);
pop_1=RbData>6 & RbData<10 & HoechstData<3000000;
pop_2=RbData>6 & RbData<10 & HoechstData>3000000;
pop_3=RbData>10 & HoechstData<3000000;
pop_4=RbData>10 & HoechstData>3000000;

numPop_1=sum(pop_1);
numPop_2=sum(pop_2);
numPop_3=sum(pop_3);
numPop_4=sum(pop_4);

total=numPop_1+numPop_2+numPop_3+numPop_4;
per_1=numPop_1/total %EdU low / Rb low
per_2=numPop_2/total %EdU high / Rb low
per_3=numPop_3/total %EdU low / Rb high
per_4=numPop_4/total %EdU high / Rb high
%axis([7 13 4 15]);
xlabel('Rb expression (log2)'); ylabel('Hoechst (log2)'); colorbar;
%%
figure;dscatter(HoechstRaw{2,2},RbRaw{2,2})
%%
figure;dscatter(RbRaw{8,2},HoechstRaw{8,2})

%%
figure; hold on;
histogram(EdURaw{6,2},100,'Normalization','probability');
histogram(EdURaw{6,4},100,'Normalization','probability');
xlim([3 16]); vline(10);
legend({'resist','resist+Bort'})
%%
figure;imagesc(SData,[0 35]); colorbar; title('%S-phase intensity');
set(gca,'Xtick',[1:12],'XTickLabel',[0 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 1 2],'fontsize',12,'tickdir','out','linewidth',2);
colormap('cyan_yellow');
%colormap('blue_yellow')

%%
figure;imagesc(CleavedCas3Data(2:end,2:end),[0 100]); colorbar; title('CleavedCaspase3 intensity');
set(gca,'Xtick',[1:12],'XTickLabel',[0 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 1 2],'fontsize',12,'tickdir','out','linewidth',2);
colormap('cyan_yellow');
%colormap('blue_yellow')
%%
figure; hold on;
histogram(CleavedCas3Raw{2,2},100,'Normalization','probability');
histogram(CleavedCas3Raw{2,10},100,'Normalization','probability');
% xlim([6 18]);
legend({'resist','resist+Bort'})



%%
figure;
RbDataPlot=[RbData(1:2,1:10) RbData(3:5,1:10)]';
imagesc(RbDataPlot,[12 13]); colorbar; title('Rb intensity (log2)');
colormap('cyan_yellow');
%colormap('blue_yellow')
set(gca,'Ytick',[1:size(RbDataPlot,1)],'YTickLabel',{'Omarigliptin','PDSF','Bortezomib','Maribavir','Ixazomib citrate','Ixazomib','Nafamostat mesylate','Celastrol','N-Ethylmaleimide','Carfilzomib','Camostat mesilate','Leupeptin Hemisulfate','Epoxomicin','MG101','Loxistatin Acid','ONX0914','PI1840','Aloxistatin','VR23','Paritaprevir','Delanzomib','Oprozomib','Z-FA-FMK','MG132'},'fontsize',12,'tickdir','out','linewidth',2);
ylabel('drugs')
set(gca,'Xtick',[1:4],'XTickLabel',[0 2 5 10],'fontsize',12,'tickdir','out','linewidth',2);
xlabel('drug concentration (然)');
figu(0.8,0.4)

%%
figure;
SDataPlot=[SData(1:4,:) SData(5:8,:)]';
imagesc(SDataPlot,[0 45]); colorbar; title('S phase cells (%)');
colormap('cyan_yellow');
%colormap('blue_yellow')
set(gca,'Ytick',[1:size(SDataPlot,1)],'YTickLabel',{'Omarigliptin','PDSF','Bortezomib','Maribavir','Ixazomib citrate','Ixazomib','Nafamostat mesylate','Celastrol','N-Ethylmaleimide','Carfilzomib','Camostat mesilate','Leupeptin Hemisulfate','Epoxomicin','MG101','Loxistatin Acid','ONX0914','PI1840','Aloxistatin','VR23','Paritaprevir','Delanzomib','Oprozomib','Z-FA-FMK','MG132'},'fontsize',12,'tickdir','out','linewidth',2);
ylabel('drugs')
set(gca,'Xtick',[1:4],'XTickLabel',[0 2 5 10],'fontsize',12,'tickdir','out','linewidth',2);
xlabel('drug concentration (然)');
figu(0.8,0.4)

%%
figure;
Cas9DataPlot=[CleavedCas3Data(1:4,:) CleavedCas3Data(5:8,:)]';
imagesc(Cas9DataPlot,[0 95]); colorbar; title('Apoptotic cells (%)');
colormap('cyan_yellow');
%colormap('blue_yellow')
set(gca,'Ytick',[1:size(Cas9DataPlot,1)],'YTickLabel',{'Omarigliptin','PDSF','Bortezomib','Maribavir','Ixazomib citrate','Ixazomib','Nafamostat mesylate','Celastrol','N-Ethylmaleimide','Carfilzomib','Camostat mesilate','Leupeptin Hemisulfate','Epoxomicin','MG101','Loxistatin Acid','ONX0914','PI1840','Aloxistatin','VR23','Paritaprevir','Delanzomib','Oprozomib','Z-FA-FMK','MG132'},'fontsize',12,'tickdir','out','linewidth',2);
ylabel('drugs')
set(gca,'Xtick',[1:4],'XTickLabel',[0 2 5 10],'fontsize',12,'tickdir','out','linewidth',2);
xlabel('drug concentration (然)');
figu(0.8,0.4)



%%
figure(1);
RbDataPlot=[RbData(1,[1 4:6 8])];

imagesc(RbDataPlot,[12 13]); colorbar; title('Rb intensity (log2)');
colormap('cyan_yellow');
%colormap('blue_yellow')
set(gca,'Ytick',[1:size(RbDataPlot,1)],'YTickLabel',{'Bortezomib'},'fontsize',12,'tickdir','out','linewidth',2);
ylabel('drugs')
set(gca,'Xtick',[1:5],'XTickLabel',[0 0.05 0.1 0.5 5],'fontsize',12,'tickdir','out','linewidth',2);
xlabel('drug concentration (然)');
figu(0.2,0.4)


%%
figure(1);
nuCellDataPlot=[numCellData(1,[1 4:6 8])];

imagesc(nuCellDataPlot); colorbar; title('Rb intensity (log2)');
colormap('cyan_yellow');
%colormap('blue_yellow')
set(gca,'Ytick',[1:size(nuCellDataPlot,1)],'YTickLabel',{'Bortezomib'},'fontsize',12,'tickdir','out','linewidth',2);
ylabel('drugs')
set(gca,'Xtick',[1:5],'XTickLabel',[0 0.05 0.1 0.5 5],'fontsize',12,'tickdir','out','linewidth',2);
xlabel('drug concentration (然)');
figu(0.2,0.4)




%% only 3 FDA-approved proteasome inhibitors
figure;
RbDataPlot=[RbData(1:4,:) RbData(5:8,:)]';
RbDataPlot_2=RbDataPlot([3 6 10],:);
imagesc(RbDataPlot_2,[10 12]); colorbar; title('Rb intensity (log2)');
colormap('cyan_yellow');
%colormap('blue_yellow')
set(gca,'Ytick',[1:size(RbDataPlot_2,1)],'YTickLabel',{'Bortezomib','Ixazomib','Carfilzomib'},'fontsize',12,'tickdir','out','linewidth',2);
ylabel('drugs')
set(gca,'Xtick',[1:4],'XTickLabel',[0 2 5 10],'fontsize',12,'tickdir','out','linewidth',2);
xlabel('drug concentration (然)');
figu(0.8,0.4)

%%
figure;
SDataPlot=[SData(1:4,:) SData(5:8,:)]';
SDataPlot_2=SDataPlot([3 6 10],:);
imagesc(SDataPlot_2,[0 45]); colorbar; title('S phase cells (%)');
colormap('cyan_yellow');
%colormap('blue_yellow')
set(gca,'Ytick',[1:size(SDataPlot_2,1)],'YTickLabel',{'Bortezomib','Ixazomib','Carfilzomib'},'fontsize',12,'tickdir','out','linewidth',2);
ylabel('drugs')
set(gca,'Xtick',[1:4],'XTickLabel',[0 2 5 10],'fontsize',12,'tickdir','out','linewidth',2);
xlabel('drug concentration (然)');
figu(0.8,0.4)
%%
numCellPlot=[numCellData(1:4,:) numCellData(5:8,:)]';
numCellPlot_2=numCellPlot([3 6 10],:);

%%
figure(2);imagesc(sigData); colorbar; title('Fluorescent intensity');
%%
figure(3);imagesc(densityData); colorbar; title(['Cell Density ExpNum ' num2str(plate)]);
colormap('cyan_yellow')

writematrix(densityData,[resultdir 'Data.xls'],'Sheet','Data','Range','A1');
%%
GC376=flip(densityData(:,9:12)',2);
compund18=flip(densityData(:,5:8)',2);
DMSO=flip(densityData(1:7,1:4)',2);

subplot(3,1,1);
imagesc(GC376,[0 35]); colormap('cyan_yellow'); title('GC376'); colorbar
subplot(3,1,2);
imagesc(compund18,[0 35]); colormap('cyan_yellow'); title('Compound18'); colorbar
subplot(3,1,3);
imagesc(DMSO,[0 35]); colormap('cyan_yellow'); title('DMSO'); colorbar
%%
conc=log10([0.0001 0.001 0.01 0.1 1 10 50 100]);
% DMSO_V=nanmean(

for i=1:12
    if ismember(i,[1:4:20])
        figure; hold on;
    end
    if ismember(i,[1:4])
        plot(conc(2:end),densityData(7:-1:1,i)); ylim([10 40]);
    else
        plot(conc,densityData(8:-1:1,i)); ylim([10 40]);
    end
    if i==4
        title('Control')
    elseif i==8
        title('compount18')
    elseif i==12
        title('GC376');
    end
    xlabel('drug concentration (log10)');
    ylabel('density (%)');
end
%%
clear plotV; clear STD;
conc=log10([0.0001 0.001 0.01 0.1 1 10 50 100]);
% DMSO_V=nanmean(
label={'control','GC376'};

plotV{1}=nanmean(densityData(7:-1:1,1:4)');
STD{1}=nanstd(densityData(7:-1:1,1:4)');

plotV{2}=nanmean(densityData(8:-1:1,9:12)');
STD{2}=nanstd(densityData(8:-1:1,9:12)');

errorbar(conc,plotV{2},STD{2})

%plot(conc,1);
ylim([10 40]); xlim([-4.5 2.5]);
xlabel('drug concentration (log10)');
ylabel('density (%)');
%set(gca,'XTickLabel',conc);
%%
conc=log10([0.0001 0.001 0.01 0.1 1 10 50 100]);
% DMSO_V=nanmean(

for i=1:12
    if ismember(i,[1:4:20])
        figure; hold on;
    end
    if ismember(i,[1:4])
        plot(conc(2:end),areaData(7:-1:1,i)); ylim([4000000 14000000]);
    else
        plot(conc,areaData(8:-1:1,i)); ylim([4000000 14000000]);
    end
    if i==4
        title('Control')
    elseif i==8
        title('compount18')
    elseif i==12
        title('GC376');
    end
    xlabel('drug concentration (log10)');
    ylabel('cell area');
end