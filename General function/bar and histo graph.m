%% Make bar graph

y=[mean(datMat{1}(:,1)),mean(datMat{3}(:,1))];
ey=[std(datMat{1}(:,1))/sqrt(length(datMat{1}(:,1))),std(datMat{3}(:,1))/sqrt(length(datMat{3}(:,1)))];
% subplot(1,3,1);
bar(1:2,y,'facecolor','k')
hold on;errorbar(1:2,y,ey,'.k')
axis([0,3,1,1.15])
nanmean(datMat{1}(:,1)')


%% Make hist graph
% data set %%%%%%%%%%%%%%%%%%%%%%%%%% 
none=datMat{1}(:,1);
mock=datMat{2}(:,1);
LY29=datMat{3}(:,1);


colors={'Cyan','blue','red'};
figure(13), hold on
for i=2:3

% subplot(1,3,1), hold on  
%     xi=linspace(0.4,2,50);
%     f=ksdensity(none,xi);
%     plot(xi,f);
%    subplot(1,3,i), hold on
    bins=0.4:0.03:2;
    x1=hist(datMat{i}(:,1),bins);

    percent_total_cells=100*x1/sum(x1);
    plot(bins,smooth(percent_total_cells,5),'Color',colors{i});xlim([0.4,2]); % bin, smoothing graph
    h=vline(nanmedian(datMat{i}(:,1)),'k:');
end
    
    subplot(1,3,2), hold on
    bins=-0.6:0.02:0;
    x1=hist(LY29,bins);