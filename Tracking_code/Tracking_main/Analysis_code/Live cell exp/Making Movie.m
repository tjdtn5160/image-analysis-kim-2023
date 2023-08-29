%% Make Alignment movie
root='E:\2. Data\2. Development of CDK4 sensor\2017-09-29 (MCF10A, BJ5ta; Geminin-CDK2-CDK4-H2B; Serum starvation)\';
moviedir=[root 'Results\'];
% load([moviedir 'Control']);

videoname = [moviedir,'Movie_v1']; % Ignore error
videoname2 = [moviedir,'Movie_v2'];
%% Gathering data
badCDK2=sensor(1).badtraces_signal2(:,1)>0.7;
badCDK4=sensor(1).badtraces_signal4(:,1)>0.7;
badtraces=sensor(1).badtraces_signal2 | sensor(1).badtraces_signal4 | badCDK2 | badCDK4;
CDK2=sensor(1).cdk2_smooth(~badtraces,:);
CDK4=sensor(1).cdk4_smooth(~badtraces,:);
APC=APCact(~badtraces,:);
% APC=sensor(1).APC(~badtraces,:);

POI_cdk4rising=sensor(1).POI_cdk4(~badtraces);
POI_cdk2rising=sensor(1).POI_cdk2(~badtraces);
POI_APCinactivation=alog(~badtraces);
%%
NumOfFrame=size(CDK2,2);
numFramesPerMovie=75;
numcells_to_plot=200;%size(CDK2,1);
minX=0;
maxX=NumOfFrame;
minY=0;
maxY=2.5;
xaxis=1:NumOfFrame;
%% Plot Cdk2 traces
figure(99),hold on
for i=1:size(CDK2,1);
    plot(xaxis,CDK2(i,1:NumOfFrame),'r','LineWidth',1);
    axis([minX,maxX,minY,maxY]);
    hold on
end
 plot(xaxis,nanmean(CDK2),'k','LineWidth',5);
%% Plot Cdk2 traces aligned by Cdk2 on 
figure(100),hold on
for i=1:size(CDK2,1);
    plot(xaxis-POI_cdk4rising(i),CDK4(i,1:NumOfFrame),'g','LineWidth',1);
    plot(xaxis-POI_cdk4rising(i),CDK2(i,1:NumOfFrame),'r','LineWidth',1);
    axis([-20,100,minY,maxY]);
    hold on
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
for i=1:size(CDK2,1);
    plot(xaxis,CDK2(i,1:NumOfFrame),'r','LineWidth',1);
    axis([minX,maxX,minY,maxY]);
    hold on
end
 plot(xaxis,nanmean(CDK2),'k','LineWidth',5);
%% Cell selection
% xtime=1:321;
xtime=1:322;
for cell_to_plot=117;
    figure(cell_to_plot); hold on;
    plot(xaxis,CDK2(cell_to_plot,:),'k','linewidth',2);
    [haxes,hline1,hline2] = plotyy(xtime,CDK4(cell_to_plot,:),xtime,log10(APC(cell_to_plot,:)));
% [haxes,hline1,hline2] = plotyy(xtime,CDK4(cell_to_plot,1:(end-1)),xtime,APC(cell_to_plot,:));
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
    if isnan(POI_cdk4rising(cell_to_plot)) | isnan(POI_cdk2rising(cell_to_plot))
        continue
    end
    scatter(POI_cdk4rising(cell_to_plot),CDK4(cell_to_plot,POI_cdk4rising(cell_to_plot)),100,'r','filled');
    text([POI_cdk4rising(cell_to_plot)+5],[CDK4(cell_to_plot,POI_cdk4rising(cell_to_plot))],{'\leftarrow CDK4 activation'},'color',[1 0 0],'fontsize',18);
    scatter(POI_cdk2rising(cell_to_plot),CDK2(cell_to_plot,POI_cdk2rising(cell_to_plot)),100,'k','filled');
    text([POI_cdk2rising(cell_to_plot)+5],[CDK2(cell_to_plot,POI_cdk2rising(cell_to_plot))],{'\leftarrow CDK2 activation'},'color',[0 0 0],'fontsize',18);
    if isnan(POI_APCinactivation(cell_to_plot)) | POI_APCinactivation(cell_to_plot)==0
        continue
    end
    
    axes(haxes(2));%axis([0 numframes 0 1]);hold on
    xlim([0 100]); ylim([-3 2]); hold on;
    set(hline2,'color',[0 0.63 0],'linewidth',2);
    set(gca,'YAxisLocation','right','YColor',[0 0.63 0],'fontsize',16,'tickdir','out','linewidth',1.5);
    ylabel('APC activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('a single CDK2 trace','fontsize',18);
    %set(gca,'YAxisLocation','right','YColor',[0 0.63 0],'YTick',0:0.2:1);
    scatter(POI_APCinactivation(cell_to_plot),log10(APC(cell_to_plot,POI_APCinactivation(cell_to_plot))),100,'filled','MarkerEdgeColor',[0 0.63 0],'MarkerFaceColor',[0 0.63 0]);
    text([POI_APCinactivation(cell_to_plot)+5],[log10(APC(cell_to_plot,POI_APCinactivation(cell_to_plot)))],{'\leftarrow APC inactivation'},'color',[0 0.63 0],'fontsize',18);
end

%% Detecting cell cycle entering point in single cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cell_to_plot=117;

aviobj=VideoWriter([videoname,'.avi']);
aviobj.FrameRate=5;
open(aviobj);
hh=figure(2);
set(gcf,'Position',[200 450 600 400],'Color','w');
figu(0.9,0.7)
NumOfFrame=100;
maxX=NumOfFrame;
xaxis=1:NumOfFrame;
for f=1:NumOfFrame;
    p1=plot(xaxis(1:f),CDK2(cell_to_plot,1:f),'k','linewidth',10); hold on
    scatter(f,CDK2(cell_to_plot,f),300,'ko','filled');
    [haxes,hline1,hline2] = plotyy(xaxis(1:f),CDK4(cell_to_plot,1:f),xaxis(1:f),log10(APC(cell_to_plot,1:f)));
    scatter(f,CDK4(cell_to_plot,f),300,'r','filled');
    axes(haxes(1));%axis([0 numframes 0.2 2]);
    set(gca,'YAxisLocation','left','YColor',[0.83 0 0],'YTick',0:0.5:2.5,'fontsize',24,'tickdir','out','linewidth',5);
    set(hline1,'color','r','linewidth',10);
    xlim([0 100]); ylim([0 2.5])
    set(gca,'Xtick',[0:25:201],'XTickLabel',[0:5:40],'fontsize',24,'tickdir','out','linewidth',5)
    ylabel('CDK2 and CDK4 activity','FontSize',30);xlabel('Time (hr)','FontSize',30);
    
    
    axes(haxes(2));%axis([0 numframes 0 1]);hold on
    xlim([0 100]); ylim([-3 2]); hold on;
    set(hline2,'color',[0 0.63 0],'linewidth',10);
    set(gca,'YAxisLocation','right','YColor',[0 0.63 0],'YTick',-3:1:2,'fontsize',24,'tickdir','out','linewidth',5);
    ylabel('APC activity','FontSize',30);title('a single cell trace','fontsize',30);
    scatter(f,log10(APC(cell_to_plot,f)),300,'filled','MarkerEdgeColor',[0 0.63 0],'MarkerFaceColor',[0 0.63 0]);
    hold off
    
    
    if f>=POI_cdk4rising(cell_to_plot)
        hold on
        axes(haxes(1));
        scatter(POI_cdk4rising(cell_to_plot),CDK4(cell_to_plot,POI_cdk4rising(cell_to_plot)),500,'r','filled');
        text([POI_cdk4rising(cell_to_plot)+3],[CDK4(cell_to_plot,POI_cdk4rising(cell_to_plot))],{'\leftarrow CDK4 activation'},'color',[1 0 0],'fontsize',30);
        axes(haxes(2));
        hold off
    end
    
    if f>=POI_cdk2rising(cell_to_plot)
        hold on
        axes(haxes(1));
        scatter(POI_cdk2rising(cell_to_plot),CDK2(cell_to_plot,POI_cdk2rising(cell_to_plot)),500,'k','filled');
        text([POI_cdk2rising(cell_to_plot)+3],[CDK2(cell_to_plot,POI_cdk2rising(cell_to_plot))],{'\leftarrow CDK2 activation'},'color',[0 0 0],'fontsize',30);
        axes(haxes(2));
        hold off
    end
    
    if f>=POI_APCinactivation(cell_to_plot)
        hold on
        axes(haxes(2));
        scatter(POI_APCinactivation(cell_to_plot),log10(APC(cell_to_plot,POI_APCinactivation(cell_to_plot))),500,'filled','MarkerEdgeColor',[0 0.63 0],'MarkerFaceColor',[0 0.63 0]);
        text([POI_APCinactivation(cell_to_plot)+3],[log10(APC(cell_to_plot,POI_APCinactivation(cell_to_plot)))],{'\leftarrow APC inactivation'},'color',[0 0.63 0],'fontsize',30);
        hold off
    end
    lgd=legend([p1 hline1 hline2],'CDK2','CDK4','APC');
    lgd.FontSize = 24;
    lgd.TextColor = 'black';
    box(haxes(1),'off'); %box(hline1,'off'); box(hline2,'off');
    
    frame=getframe(gcf);
    writeVideo(aviobj,frame);
end
close(aviobj);





%%
for i=1:40
    writeVideo(aviobj,frame);
end
%% Add a couple of traces on the top of a single trace
cells_to_plot2=[1:5:size(CDK2,1)];
CDK2_plot2=CDK2(cells_to_plot2,:);
CDK2_low_plot2=CDK2_low(cells_to_plot2); 

colors_to_plot=makeColorMap([0 0 0],[0.5 0.5 0.5],length(cells_to_plot2));
for f=1:NumOfFrame;
    plot(xaxis(1:NumOfFrame),CDK2(cell_to_plot,1:NumOfFrame),'k','linewidth',1);hold on
    for numcell=1:length(cells_to_plot2)
        plot(xaxis(1:f),CDK2(cells_to_plot2(numcell),1:f),'color',colors_to_plot(numcell,:),'linewidth',1);
        hold on
        scatter(f,CDK2(cells_to_plot2(numcell),f),100,'o','filled','markerfacecolor',colors_to_plot(numcell,:));
        
        axis([minX,maxX,minY,maxY]);
        set(gca,'Xtick',[0:25:201],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
        ylabel('Proliferation signal (CDK2)','FontSize',18);xlabel('Time since mitogen stimulation (hr)','FontSize',18);title('Heterogeneous population','fontsize',18);
    end
    box off
    hold off
    frame=getframe(gcf);
    writeVideo(aviobj,frame);
end

for i=1:20
    writeVideo(aviobj,frame);
end


plot(xaxis,CDK2_plot2(~CDK2_low_plot2,:),'Color',[0.2 1 0.2]);hold on;
plot(xaxis,CDK2_plot2(CDK2_low_plot2,:),'r');
text([2],[2],{'Proliferation'},'color',[0.2 0.8 0.2],'fontsize',18);
text([2],[0.2],{'Quiescence'},'color',[1 0 0],'fontsize',18);
axis([minX,maxX,minY,maxY]); box off
set(gca,'Xtick',[0:25:201],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
ylabel('Proliferation signal (CDK2)','FontSize',18);xlabel('Time since mitogen stimulation (hr)','FontSize',18);title('Heterogeneous population','fontsize',18);
frame=getframe(gcf);
writeVideo(aviobj,frame);
for i=1:50
    writeVideo(aviobj,frame);
end

plot(xaxis,nanmean(CDK2),'k','LineWidth',5);
text([100],[1.6],{'Averaged'},'color',[0 0 0],'fontsize',18,'FontWeight','bold');
frame=getframe(gcf);
writeVideo(aviobj,frame);
for i=1:30
    writeVideo(aviobj,frame);
end

plot(xaxis,nanmedian(CDK2_plot2(~CDK2_low_plot2,:)),'Color',[0 0.5 0],'LineWidth',5);
plot(xaxis,nanmedian(CDK2_plot2(CDK2_low_plot2,:)),'Color',[0.5 0 0],'LineWidth',5);
frame=getframe(gcf);
writeVideo(aviobj,frame);
for i=1:50
    writeVideo(aviobj,frame);
end
close(aviobj);
















%% Select specific windows
clf
plot(xaxis(1:NumOfFrame),CDK2(cells_to_plot2,1:NumOfFrame),'k','linewidth',1);hold on
for numcell=1:length(cells_to_plot2)
    scatter(POI_mitosis(cells_to_plot2(numcell),1),CDK2(cells_to_plot2(numcell),POI_mitosis(cells_to_plot2(numcell),1)),100,'r','filled');
    if ~isnan(POI_rising(cells_to_plot2(numcell),1))
        scatter(POI_rising(cells_to_plot2(numcell),1),CDK2(cells_to_plot2(numcell),POI_rising(cells_to_plot2(numcell),1)),100,'g','filled');
    end
end
axis([minX,maxX,minY,maxY]);
set(gca,'Xtick',[0:25:201],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
ylabel('CDK2 activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('\it in silico \rm \bf selection','fontsize',18);
box off
frame=getframe(gcf);
for i=1:20
    writeVideo(aviobj,frame);
end

vline([100 125],'r:'); 
frame=getframe(gcf);
for i=1:40
    writeVideo(aviobj,frame);
end

clf
selection=POI_mitosis(cells_to_plot2,1)>100 & POI_mitosis(cells_to_plot2,1)<125;
cell_to_plot3=cells_to_plot2(selection);
plot(xaxis(1:NumOfFrame),CDK2(cell_to_plot3,1:NumOfFrame),'k','linewidth',1);hold on
for numcell=1:length(cell_to_plot3)
    scatter(POI_mitosis(cell_to_plot3(numcell),1),CDK2(cell_to_plot3(numcell),POI_mitosis(cell_to_plot3(numcell),1)),100,'r','filled');
    if ~isnan(POI_rising(cell_to_plot3(numcell),1))
        scatter(POI_rising(cell_to_plot3(numcell),1),CDK2(cell_to_plot3(numcell),POI_rising(cell_to_plot3(numcell),1)),100,'g','filled');
    end
end
vline([100 125],'r:'); 
axis([minX,maxX,minY,maxY]);
set(gca,'Xtick',[0:25:201],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
ylabel('CDK2 activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('\it in silico \rm \bf selection','fontsize',18);
box off
frame=getframe(gcf);
for i=1:70
    writeVideo(aviobj,frame);
end

clf
plot(xaxis(1:NumOfFrame),CDK2(cells_to_plot2,1:NumOfFrame),'k','linewidth',1);hold on
for numcell=1:length(cells_to_plot2)
    scatter(POI_mitosis(cells_to_plot2(numcell),1),CDK2(cells_to_plot2(numcell),POI_mitosis(cells_to_plot2(numcell),1)),100,'r','filled');
    if ~isnan(POI_rising(cells_to_plot2(numcell),1))
        scatter(POI_rising(cells_to_plot2(numcell),1),CDK2(cells_to_plot2(numcell),POI_rising(cells_to_plot2(numcell),1)),100,'g','filled');
    end
end
axis([minX,maxX,minY,maxY]);
set(gca,'Xtick',[0:25:201],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
ylabel('CDK2 activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('\it in silico \rm \bf selection','fontsize',18);
box off
vline([100 125],'g:'); 
frame=getframe(gcf);
for i=1:40
    writeVideo(aviobj,frame);
end

clf
selection2=POI_rising(cells_to_plot2,1)>100 & POI_rising(cells_to_plot2,1)<125;
cell_to_plot3=cells_to_plot2(selection2);
plot(xaxis(1:NumOfFrame),CDK2(cell_to_plot3,1:NumOfFrame),'k','linewidth',1);hold on
for numcell=1:length(cell_to_plot3)
    scatter(POI_mitosis(cell_to_plot3(numcell),1),CDK2(cell_to_plot3(numcell),POI_mitosis(cell_to_plot3(numcell),1)),100,'r','filled');
    if ~isnan(POI_rising(cell_to_plot3(numcell),1))
        scatter(POI_rising(cell_to_plot3(numcell),1),CDK2(cell_to_plot3(numcell),POI_rising(cell_to_plot3(numcell),1)),100,'g','filled');
    end
end
vline([100 125],'g:'); 
axis([minX,maxX,minY,maxY]);
set(gca,'Xtick',[0:25:201],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
ylabel('CDK2 activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('\it in silico \rm \bf selection','fontsize',18);
box off
frame=getframe(gcf);
for i=1:70
    writeVideo(aviobj,frame);
end
%% alignment movie
clf
CDK2_new=CDK2(cells_to_plot2,1:NumOfFrame);
POI_new=POI_mitosis(cells_to_plot2,:);
alignVec=POI_new(:,1);
alignStepVec = (alignVec-80)./numFramesPerMovie;
colors_to_plot_all_cells=makeColorMap([0 0 0],[0.5 0.5 0.5],length(cells_to_plot2));

for numcell=1:length(cells_to_plot2)
    hold on
    plot(xaxis,CDK2_new(numcell,:),'color',colors_to_plot_all_cells(numcell,:),'LineWidth',1);
    scatter(POI_new(numcell,1),CDK2_new(numcell,POI_new(numcell,1)),100,'r','filled');
    axis([minX,maxX,minY,maxY]);
    set(gca,'Xtick',[0:25:200],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
    ylabel('CDK2 activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
    box off
end
frame=getframe(gcf);
for i=1:40
    writeVideo(aviobj,frame);
end

for (numFrame=1:numFramesPerMovie)
    hold off;
    thisFrameBias = numFrame.*alignStepVec;
    for (biasIndex=1:length(cells_to_plot2));%length(alignVec))
        thisXVec = xaxis-thisFrameBias(biasIndex);
        plot(thisXVec,CDK2_new(biasIndex,:),'color',colors_to_plot_all_cells(biasIndex,:),'LineWidth',1)
        hold on
        scatter(POI_new(biasIndex,1)-thisFrameBias(biasIndex),CDK2_new(biasIndex,POI_new(biasIndex,1)),100,'r','filled');
    end
    box off
    axis([minX,maxX,minY,maxY]);
    if numFrame==numFramesPerMovie
        set(gca,'Xtick',[0:20:200],'Xticklabel',[-16:4:24],'fontsize',16,'tickdir','out','linewidth',1.5);
        xlabel('Time (hr)','FontSize',18);
    else
        set(gca,'XTick',[]);
        %set(gca,'Xtick',[0:25:200],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
    end
    ylabel('CDK2 activity','FontSize',18);xlabel('Time since mitosis (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
  
    frame=getframe(gcf);
    writeVideo(aviobj,frame);
end
frame=getframe(gcf);
for i=1:50;
    writeVideo(aviobj,frame);
end

clf
%indCDK2inc=cells_to_plot2(D.cellState(cells_to_plot2));
indCDK2inc=cells_to_plot2(cells_to_plot2<=Cdk2inc_end);
for numcell=1:length(indCDK2inc)
    hold on
    plot(xaxis-POI_mitosis(indCDK2inc(numcell)),CDK2(indCDK2inc(numcell),1:NumOfFrame),'color',[0 0 1],'LineWidth',1);
end

indCDK2delay=cells_to_plot2(Cdk2inc_end<cells_to_plot2 & cells_to_plot2<=Cdk2delay_end);
for numcell=1:length(indCDK2delay)
    plot(xaxis-POI_mitosis(indCDK2delay(numcell)),CDK2(indCDK2delay(numcell),1:NumOfFrame),'color',[1 0.5 0],'LineWidth',1);
end

indCDK2Low=cells_to_plot2(Cdk2delay_end<cells_to_plot2 & cells_to_plot2<=Cdk2low_end);
for numcell=1:length(indCDK2Low)
    plot(xaxis-POI_mitosis(indCDK2Low(numcell)),CDK2(indCDK2Low(numcell),1:NumOfFrame),'color',[1 0 0],'LineWidth',1);
end

vline(0,'k:'); box off
axis([-80,120,minY,maxY]);
set(gca,'Xtick',[-80:20:120],'Xticklabel',[-16:4:24] ,'fontsize',16,'tickdir','out','linewidth',1.5)
ylabel('CDK2 activity','FontSize',18);xlabel('Time since mitosis (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
text(3,2.7,'Cdk2^{high}','color',[0 0 1],'fontsize',18);
text(3,2.4,'Cdk2^{delay}','color',[1 0.5 0],'fontsize',18);
text(3,2.1,'Cdk2^{low}','color',[1 0 0],'fontsize',18);
frame=getframe(gcf);
for i=1:70
    writeVideo(aviobj,frame);
end
%% IF align movie
for numcell=1:length(indCDK2inc)
    tempAxis=xaxis-POI_mitosis(indCDK2inc(numcell));
    XvalCDK2inc(numcell)=tempAxis(end);
    scatter(XvalCDK2inc(numcell),CDK2(indCDK2inc(numcell),NumOfFrame),100,'filled','MarkerFaceColor',[0.4 0.4 1]);
end

for numcell=1:length(indCDK2delay)
    tempAxis=xaxis-POI_mitosis(indCDK2delay(numcell));
    XvalCDK2delay(numcell)=tempAxis(end);
    scatter(XvalCDK2delay(numcell),CDK2(indCDK2delay(numcell),NumOfFrame),100,'filled','MarkerFaceColor',[1 0.5 0]);
end

for numcell=1:length(indCDK2Low)
    tempAxis=xaxis-POI_mitosis(indCDK2Low(numcell));
    XvalCDK2low(numcell)=tempAxis(end);
    scatter(XvalCDK2low(numcell),CDK2(indCDK2Low(numcell),NumOfFrame),100,'filled','MarkerFaceColor',[1 0.4 0.4]);
end

frame=getframe(gcf);
for i=1:60
    writeVideo(aviobj,frame);
end

IFCDK2inc_1=CDK2(indCDK2inc,NumOfFrame); IFCDK2inc_2=CDK2(indCDK2inc,NumOfFrame)*0.1;
IFCDK2delay_1=CDK2(indCDK2delay,NumOfFrame); IFCDK2delay_2=CDK2(indCDK2delay,NumOfFrame)*3;
IFCDK2low_1=CDK2(indCDK2Low,NumOfFrame); IFCDK2low_2=CDK2(indCDK2Low,NumOfFrame)*2;

numFramesPerMovie=40;
alignStepVec_inc = abs(IFCDK2inc_2-IFCDK2inc_1)./numFramesPerMovie;
alignStepVec_delay = abs(IFCDK2delay_2-IFCDK2delay_1)./numFramesPerMovie;
alignStepVec_low = abs(IFCDK2low_2-IFCDK2low_1)./numFramesPerMovie;

for (numFrame=1:numFramesPerMovie)
    clf; hold on;
    thisFrameBias_inc = numFrame.*alignStepVec_inc;
    thisYVec_inc = IFCDK2inc_1+thisFrameBias_inc;
    scatter(XvalCDK2inc,thisYVec_inc,100,'filled','MarkerFaceColor',[0.4 0.4 1]);
    
    thisFrameBias_delay = numFrame.*alignStepVec_delay;
    thisYVec_delay = IFCDK2delay_1+thisFrameBias_delay;
    scatter(XvalCDK2delay,thisYVec_delay,100,'filled','MarkerFaceColor',[1 0.5 0]);
    
    thisFrameBias_low = numFrame.*alignStepVec_low;
    thisYVec_low = IFCDK2low_1+thisFrameBias_low;
    scatter(XvalCDK2low,thisYVec_low,100,'filled','MarkerFaceColor',[1 0.4 0.4]);
    
    axis([-80,120,0,8]);
    vline(0,'k:'); box off
    set(gca,'Xtick',[-80:20:120],'Xticklabel',[-16:4:24] ,'fontsize',16,'tickdir','out','linewidth',1.5)
    ylabel('IF (log2)','FontSize',18);xlabel('Time since mitosis (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
    text(3,6,'Cdk2^{high}','color',[0.4 0.4 1],'fontsize',18);
    text(3,5,'Cdk2^{delay}','color',[1 0.5 0],'fontsize',18);
    text(3,4,'Cdk2^{low}','color',[1 0.4 0.4],'fontsize',18);

    frame=getframe(gcf);
    writeVideo(aviobj,frame);
end
frame=getframe(gcf);
for i=1:50;
    writeVideo(aviobj,frame);
end
close(aviobj);




%% Make alignment movie-1
aviobj=VideoWriter([videoname2,'.avi']);
aviobj.FrameRate=15;
open(aviobj);
hh=figure(3);
set(gcf,'Position',[200 100 600 800],'Color','w');

random_cells=randperm(numcells_to_plot);
% random_cells=randi([1,size(CDK2,1)],1,numcells_to_plot);
CDK2_new=CDK2(random_cells,1:NumOfFrame);
ERK_new=ERK(random_cells,1:NumOfFrame);
POI_new=POI_mitosis(random_cells,:);
alignVec=POI_new(:,1);
alignStepVec = (alignVec-80)./numFramesPerMovie;
colors_to_plot_all_cells=makeColorMap([0 0 0],[0.5 0.5 0.5],size(CDK2,1));
subplot(2,1,1);
for numcell=1:numcells_to_plot
    hold on
    plot(xaxis,CDK2_new(numcell,:),'color',colors_to_plot_all_cells(numcell,:),'LineWidth',1);
    scatter(POI_new(numcell,1),CDK2_new(numcell,POI_new(numcell,1)),100,'r','filled');
    axis([minX,maxX,minY,maxY]);
    set(gca,'Xtick',[0:25:200],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
    ylabel('CDK2 activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
    box off
end

subplot(2,1,2);
for numcell=1:numcells_to_plot
    hold on
    plot(xaxis,ERK_new(numcell,:),'color',colors_to_plot_all_cells(numcell,:),'LineWidth',1);
    scatter(POI_mitosis(cell_to_plot,1),CDK2(cell_to_plot,POI_mitosis(cell_to_plot,1)),100,'r','filled');
    axis([minX,maxX,1.2,2.2]);
    set(gca,'Xtick',[0:25:200],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
    ylabel('ERK activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
    box off
end

frame=getframe(gcf);
for i=1:20
    writeVideo(aviobj,frame);
end
%%
for (numFrame=1:numFramesPerMovie)
    subplot(2,1,1);
    hold off;
    thisFrameBias = numFrame.*alignStepVec;
    for (biasIndex=1:numcells_to_plot);%length(alignVec))
        thisXVec = xaxis-thisFrameBias(biasIndex);
        plot(thisXVec,CDK2_new(biasIndex,:),'color',colors_to_plot_all_cells(biasIndex,:),'LineWidth',1)
        hold on
        scatter(POI_new(biasIndex,1)-thisFrameBias(biasIndex),CDK2_new(biasIndex,POI_new(biasIndex,1)),100,'r','filled');
    end
    box off
    axis([minX,maxX,minY,maxY]);
    if numFrame==numFramesPerMovie
        set(gca,'Xtick',[0:20:200],'Xticklabel',[-16:4:24],'fontsize',16,'tickdir','out','linewidth',1.5);
        xlabel('Time (hr)','FontSize',18);
    else
        set(gca,'XTick',[]);
        %set(gca,'Xtick',[0:25:200],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
    end
    ylabel('CDK2 activity','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
    
    subplot(2,1,2);
    hold off;
    for (biasIndex=1:numcells_to_plot);%length(alignVec))
        thisXVec = xaxis-thisFrameBias(biasIndex);
        plot(thisXVec,ERK_new(biasIndex,:),'color',colors_to_plot_all_cells(biasIndex,:),'LineWidth',1)
        hold on
    end
    box off
    axis([minX,maxX,1.2,2.2]);
    if numFrame==numFramesPerMovie
        set(gca,'Xtick',[0:20:200],'Xticklabel',[-16:4:24],'fontsize',16,'tickdir','out','linewidth',1.5);
        xlabel('Time (hr)','FontSize',18);
    else
        set(gca,'XTick',[]);
        %set(gca,'Xtick',[0:25:200],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
    end
    ylabel('ERK activity','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
    
    frame=getframe(gcf);
    writeVideo(aviobj,frame);
end
frame=getframe(gcf);
for i=1:40;
    writeVideo(aviobj,frame);
end
% %%
% hold off
% for numcell=1:size(APC_off)
%    hold on
%    plot(xaxis-POI_off(numcell)+10,APC_off(numcell,:),'r','LineWidth',1);
%    set(gca,'XTickLabel',[-2,0,2,4,6,8,10],'fontsize',16,'tickdir','out','linewidth',1.5);xlabel('Time Relative to Mitosis (hr)','FontSize',18);
%    box off
% end
% for numcell=1:size(APC_on)
%    hold on
%    plot(xaxis-POI_on(numcell)+10,APC_on(numcell,:),'color',[0 0.4 1],'LineWidth',1);
%    set(gca,'XTickLabel',[-2,0,2,4,6,8,10],'fontsize',16,'tickdir','out','linewidth',1.5);xlabel('Time Relative to Mitosis (hr)','FontSize',18);
%    box off
% end
% text(38,0.95,'Quiescent Cells','color',[0 0.4 1],'fontsize',18);
% text(36,0.3,'Cycling Cells','color',[0.9 0 0],'fontsize',18);
% frame=getframe(gcf);
% for i=1:30;
%     writeVideo(aviobj,frame);
% end
% %%
% plot(xaxis,APCallignedMean_off,'color',[0.3 0 0],'LineWidth',6);set(gca,'XTickLabel',[-2,0,2,4,6,8,10],'fontsize',16,'tickdir','out','linewidth',1.5);xlabel('Time Relative to Mitosis (hr)','FontSize',18);
% plot(xaxis,APCallignedMean_on,'color',[0 0 0.5],'LineWidth',6);set(gca,'XTickLabel',[,-2,0,2,4,6,8,10],'fontsize',16,'tickdir','out','linewidth',1.5);xlabel('Time Relative to Mitosis (hr)','FontSize',18);
% box off
% frame=getframe(gcf);
% writeVideo(aviobj,frame);
%%
clf
subplot(2,1,1);
for numcell=1:size(Cdk2inc)
    hold on
    plot(xaxis-POIinc_mitosis(numcell),Cdk2inc(numcell,1:NumOfFrame),'color',[0 0 1],'LineWidth',1);
end
for numcell=1:size(Cdk2delay)
    plot(xaxis-POIdelay_mitosis(numcell),Cdk2delay(numcell,1:NumOfFrame),'color',[1 0.5 0],'LineWidth',1);
end
for numcell=1:size(Cdk2low)
    plot(xaxis-POIlow_mitosis(numcell),Cdk2low(numcell,1:NumOfFrame),'color',[1 0 0],'LineWidth',1);
end
vline(0,'k:'); box off
axis([-80,120,minY,maxY]);
set(gca,'Xtick',[-80:20:120],'Xticklabel',[-16:4:24] ,'fontsize',16,'tickdir','out','linewidth',1.5)
ylabel('CDK2 activity','FontSize',18);xlabel('Time since mitosis (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
text(0,2.9,'Cdk2^{high}','color',[0 0 1],'fontsize',18);
text(0,2.6,'Cdk2^{delay}','color',[1 0.5 0],'fontsize',18);
text(0,2.3,'Cdk2^{low}','color',[1 0 0],'fontsize',18);

temp=subplot(2,1,2);
for numcell=1:size(ERKinc)
    hold on
    plot(xaxis-POIinc_mitosis(numcell),ERKinc(numcell,1:NumOfFrame),'color',[0 0 1],'LineWidth',1);
end
for numcell=1:size(ERKdelay)
    plot(xaxis-POIdelay_mitosis(numcell),ERKdelay(numcell,1:NumOfFrame),'color',[1 0.5 0],'LineWidth',1);
end
for numcell=1:size(ERKlow)
    plot(xaxis-POIlow_mitosis(numcell),ERKlow(numcell,1:NumOfFrame),'color',[1 0 0],'LineWidth',1);
end
vline(0,'k:'); box off
axis([-80,120,1.2,2.2]);
set(gca,'Xtick',[-80:20:120],'Xticklabel',[-16:4:24] ,'fontsize',16,'tickdir','out','linewidth',1.5)
ylabel('ERK activity','FontSize',18);xlabel('Time since mitosis (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);

frame=getframe(gcf);
for i=1:40
    writeVideo(aviobj,frame);
end

numFramesToPlot = 35./0.2;
frameToPlot = 5000-numFramesToPlot+1:5000+numFramesToPlot;
timeToPlot = linspace(-25,25,numFramesToPlot*2);
delete(temp)
subplot(2,1,2);
hold on

% Plot mean +/- 95% CI
ci = 0.95;
alpha = 1 - ci;
Sinc=nanstd(ERKinc_align(:,frameToPlot));
Sdelay=nanstd(ERKdelay_align(:,frameToPlot));
Slow=nanstd(ERKlow_align(:,frameToPlot));

Ninc=sum(~isnan(ERKinc_align(:,frameToPlot)));
T_multiplier = tinv(1-alpha/2, Ninc-1);
ci95inc = T_multiplier.*Sinc./sqrt(Ninc);

Ndelay=sum(~isnan(ERKdelay_align(:,frameToPlot)));
T_multiplier = tinv(1-alpha/2, Ndelay-1);
ci95delay = T_multiplier.*Sdelay./sqrt(Ndelay);

Nlow=sum(~isnan(ERKlow_align(:,frameToPlot)));
T_multiplier = tinv(1-alpha/2, Nlow-1);
ci95low = T_multiplier.*Slow./sqrt(Nlow);

shadedErrorBar(timeToPlot,nanmean(ERKinc_align(:,frameToPlot)),ci95inc,{'Color',[0 0 1],'linewidth',2},1);
shadedErrorBar(timeToPlot,nanmean(ERKdelay_align(:,frameToPlot)),ci95delay,{'Color',[1 0.5 0],'linewidth',2},1);
shadedErrorBar(timeToPlot,nanmean(ERKlow_align(:,frameToPlot)),ci95low,{'Color',[1 0 0],'linewidth',2},1);
vline(0,'k:');
axis([-16,24,1.35,2]);
set(gca,'Xtick',[-16:4:24],'Xticklabel',[-16:4:24] ,'fontsize',16,'tickdir','out','linewidth',1.5)
ylabel('ERK activity','FontSize',18);xlabel('Time relative to mitosis (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);

frame=getframe(gcf);
for i=1:40
    writeVideo(aviobj,frame);
end


% rectangle('Position',[9,1.65,24,0.13],'edgecolor',[0 0.2 1],'FaceColor','w');
% rectangle('Position',[59,0.55,29,0.13],'edgecolor',[0.9 0 0],'FaceColor','w');
% text(10,1.7,'Cycling Cells','color',[0 0.2 1],'fontsize',18);
% text(60,0.6,'Quiescent Cells','color',[0.9 0 0],'fontsize',18);
%% write movie
close(aviobj);



%% Make alignment movie-2
aviobj=VideoWriter([videoname2,'.avi']);
aviobj.FrameRate=15;
open(aviobj);
hh=figure(3);
set(gcf,'Position',[200 100 600 500],'Color','w');

cells_to_plot3=[1:5:size(CDK2,1)];
% numcells_to_plot=200
% cells_to_plot3=randperm(numcells_to_plot);
% random_cells=randi([1,size(CDK2,1)],1,numcells_to_plot);
ERK_new=ERK(cells_to_plot3,1:NumOfFrame);
POI_new=POI_mitosis(cells_to_plot3,:);
alignVec=POI_new(:,1);
alignStepVec = (alignVec-80)./numFramesPerMovie;
colors_to_plot_all_cells=makeColorMap([0 0 0],[0.5 0.5 0.5],size(CDK2,1));
minX=30; maxX=155;

for numcell=1:length(cells_to_plot3)
    hold on
    plot(xaxis,ERK_new(numcell,:),'color',colors_to_plot_all_cells(numcell,:),'LineWidth',1);
    scatter(POI_new(numcell,1),ERK_new(numcell,POI_new(numcell,1)),100,'r','filled');
    box off
end
axis([minX,maxX,1.2,2.2]);
set(gca,'Xtick',[minX:25:200],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5);
set(gca,'Ytick',[1.2:0.2:2.2],'fontsize',16,'tickdir','out','linewidth',1.5);
ylabel('ERK activity','FontSize',18);xlabel('Time (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
    
frame=getframe(gcf);
for i=1:20
    writeVideo(aviobj,frame);
end
%%
for (numFrame=1:numFramesPerMovie)
    hold off;
    thisFrameBias = numFrame.*alignStepVec;

    for (biasIndex=1:length(cells_to_plot3));%length(alignVec))
        thisXVec = xaxis-thisFrameBias(biasIndex);
        plot(thisXVec,ERK_new(biasIndex,:),'color',colors_to_plot_all_cells(biasIndex,:),'LineWidth',1)
        hold on
        scatter(POI_new(biasIndex,1)-thisFrameBias(biasIndex),ERK_new(biasIndex,POI_new(biasIndex,1)),100,'r','filled');
    end
    box off
    axis([minX,maxX,1.2,2.2]);
    if numFrame==numFramesPerMovie
        set(gca,'Xtick',[minX:25:200],'Xticklabel',[-10:5:15],'fontsize',16,'tickdir','out','linewidth',1.5);
        xlabel('Time since mitosis (hr)','FontSize',18);
    else
        set(gca,'XTick',[]);
        %set(gca,'Xtick',[0:25:200],'XTickLabel',[0:5:40],'fontsize',16,'tickdir','out','linewidth',1.5)
        xlabel('Time (hr)','FontSize',18);
    end
    set(gca,'Ytick',[1.2:0.2:2.2],'fontsize',16,'tickdir','out','linewidth',1.5);
    ylabel('ERK activity','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
    
    frame=getframe(gcf);
    writeVideo(aviobj,frame);
end
frame=getframe(gcf);
for i=1:40;
    writeVideo(aviobj,frame);
end
%%
clf
indCDK2inc=cells_to_plot3(cells_to_plot3<=Cdk2inc_end);
indCDK2delay=cells_to_plot3(Cdk2inc_end<cells_to_plot3 & cells_to_plot3<=Cdk2delay_end);
indCDK2Low=cells_to_plot3(Cdk2delay_end<cells_to_plot3 & cells_to_plot3<=Cdk2low_end);
for numcell=1:length(indCDK2inc)
    hold on
    plot(xaxis-POI_mitosis(indCDK2inc(numcell)),ERK(indCDK2inc(numcell),1:NumOfFrame),'color',[0 0 1],'LineWidth',1);
end
for numcell=1:length(indCDK2delay)
    plot(xaxis-POI_mitosis(indCDK2delay(numcell)),ERK(indCDK2delay(numcell),1:NumOfFrame),'color',[1 0.5 0],'LineWidth',1);
end
for numcell=1:length(indCDK2Low)
    plot(xaxis-POI_mitosis(indCDK2Low(numcell)),ERK(indCDK2Low(numcell),1:NumOfFrame),'color',[1 0 0],'LineWidth',1);
end
vline(0,'k:'); box off
axis([-50,75,1.2,2.2]);
set(gca,'Xtick',[-50:25:75],'Xticklabel',[-10:5:15] ,'fontsize',16,'tickdir','out','linewidth',1.5)
set(gca,'Ytick',[1.2:0.2:2.2],'fontsize',16,'tickdir','out','linewidth',1.5);
ylabel('ERK activity','FontSize',18);xlabel('Time since mitosis (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
text(12,2.14,'Cdk2^{high}','color',[0 0 1],'fontsize',18);
text(12,2.07,'Cdk2^{delay}','color',[1 0.5 0],'fontsize',18);
text(12,2,'Cdk2^{low}','color',[1 0 0],'fontsize',18);

frame=getframe(gcf);
for i=1:40
    writeVideo(aviobj,frame);
end

numFramesToPlot = 35./0.2;
frameToPlot = 5000-numFramesToPlot+1:5000+numFramesToPlot;
timeToPlot = linspace(-25,25,numFramesToPlot*2);

clf
hold on
% Plot mean +/- 95% CI
ci = 0.95;
alpha = 1 - ci;
Sinc=nanstd(ERKinc_align(:,frameToPlot));
Sdelay=nanstd(ERKdelay_align(:,frameToPlot));
Slow=nanstd(ERKlow_align(:,frameToPlot));

Ninc=sum(~isnan(ERKinc_align(:,frameToPlot)));
T_multiplier = tinv(1-alpha/2, Ninc-1);
ci95inc = T_multiplier.*Sinc./sqrt(Ninc);

Ndelay=sum(~isnan(ERKdelay_align(:,frameToPlot)));
T_multiplier = tinv(1-alpha/2, Ndelay-1);
ci95delay = T_multiplier.*Sdelay./sqrt(Ndelay);

Nlow=sum(~isnan(ERKlow_align(:,frameToPlot)));
T_multiplier = tinv(1-alpha/2, Nlow-1);
ci95low = T_multiplier.*Slow./sqrt(Nlow);

shadedErrorBar(timeToPlot,nanmean(ERKinc_align(:,frameToPlot)),ci95inc,{'Color',[0 0 1],'linewidth',2},1);
shadedErrorBar(timeToPlot,nanmean(ERKdelay_align(:,frameToPlot)),ci95delay,{'Color',[1 0.5 0],'linewidth',2},1);
shadedErrorBar(timeToPlot,nanmean(ERKlow_align(:,frameToPlot)),ci95low,{'Color',[1 0 0],'linewidth',2},1);
vline(0,'k:');
axis([-10,15,1.35,2]);
set(gca,'Xtick',[-10:5:15],'Xticklabel',[-10:5:15] ,'fontsize',16,'tickdir','out','linewidth',1.5)
set(gca,'Ytick',[1.4:0.2:2.0],'fontsize',16,'tickdir','out','linewidth',1.5);
ylabel('ERK activity','FontSize',18);xlabel('Time relative to mitosis (hr)','FontSize',18);title('\it in silico \rm \bf alignment','fontsize',18);
text(2,1.95,'Cdk2^{high}','color',[0 0 1],'fontsize',18);
text(2,1.9,'Cdk2^{delay}','color',[1 0.5 0],'fontsize',18);
text(2,1.85,'Cdk2^{low}','color',[1 0 0],'fontsize',18);


frame=getframe(gcf);
for i=1:40
    writeVideo(aviobj,frame);
end


% rectangle('Position',[9,1.65,24,0.13],'edgecolor',[0 0.2 1],'FaceColor','w');
% rectangle('Position',[59,0.55,29,0.13],'edgecolor',[0.9 0 0],'FaceColor','w');
% text(10,1.7,'Cycling Cells','color',[0 0.2 1],'fontsize',18);
% text(60,0.6,'Quiescent Cells','color',[0.9 0 0],'fontsize',18);
%% write movie
close(aviobj);