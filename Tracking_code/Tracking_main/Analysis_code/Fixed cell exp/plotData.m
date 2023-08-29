function [plotCDK4V plotCDK2V plotRb]=plotData(conditions,plotMatCDK2,plotMatCDK4,pRb_pos,inhibitor,conc)
count=1;
norm=1.15;
figure; hold on;
color=jet(length(inhibitor));
for i=1:length(inhibitor)
    %     figure; hold on;
    for j=1:length(conc)
        val=plotMatCDK4(count,:)*norm;
        plotCDK4V{i}(j)=mean(val);
        stdardErrorCDK4{i}(j)=std(val);
        count=count+1;
    end
    x=log10(conc*1000);
    x=[0 x(2:end)];
    errorbar(x,plotCDK4V{i},stdardErrorCDK4{i},'.-','MarkerSize',15);
    ylim([0.5 1.05]); xlim([-0.5 5]);
    xlabel('Inhibitor (log10, nM)');
    ylabel('Mean Cyt/Nuc RBCterm');
end
legend(inhibitor);

count=1;
figure; hold on;
color=jet(length(inhibitor));
for i=1:length(inhibitor)
    %     figure; hold on;
    for j=1:length(conc)
        val=plotMatCDK2(count,:)*norm;
        plotCDK2V{i}(j)=mean(val);
        stdardErrorCDK2{i}(j)=std(val);
        count=count+1;
    end
    x=log10(conc*1000);
    x=[0 x(2:end)];
    errorbar(x,plotCDK2V{i},stdardErrorCDK2{i},'.-','MarkerSize',15);
    ylim([0.5 1.05]); xlim([-0.5 5]);
    xlabel('Inhibitor (log10, nM)');
    ylabel('Mean Cyt/Nuc DHB');
end
legend(inhibitor);

count=1;
figure; hold on;
color=jet(length(inhibitor));
for i=1:length(inhibitor)
    %     figure; hold on;
    for j=1:length(conc)
        val=pRb_pos(count,:)*norm;
        plotRb{i}(j)=mean(val);
        stdardErrorRB{i}(j)=std(val);
        count=count+1;
    end
    x=log10(conc*1000);
    x=[0 x(2:end)];
    errorbar(x,plotRb{i},stdardErrorRB{i},'.-','MarkerSize',15);
    ylim([0 100]); xlim([-0.5 5]);
    xlabel('Inhibitor (log10, nM)');
    ylabel('% p-Rb');
end
legend(inhibitor);