function CDK2vspRb(conditions,datadir)
%condnum=size(conditions,1);
allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);
sigchoice='p21';
statarray=cell(uniquecondnum,1);

G1minH=200000; G1maxH=350000; G1minE=3; G1maxE=6.5;
SminH=200000; SmaxH=600000; SminE=8; SmaxE=12;
G2minH=500000; G2maxH=700000; G2minE=3; G2maxE=6.5;
xylim=[100000 800000 1 13];

for i=1:uniquecondnum
    Alldata=[];
    condrow=find(ismember(conditions(:,1),uniquenames{i}));
    wellstats=[];
    for c=condrow'
        rowmat=cell2mat(conditions(c,2));
        colmat=cell2mat(conditions(c,3));
        sitemat=cell2mat(conditions(c,4));
        for row=rowmat
            for col=colmat
                for site=sitemat
                    %shot=wellnum2str(row,col,site);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    load([datadir,shot,'.mat'],'IFdata');
                    Alldata=[Alldata;IFdata];
                end
                Hoechstval=Alldata(:,3).*Alldata(:,4);
                EdUval=Alldata(:,7);
                EdUval(EdUval<1)=1;
                EdUval=log2(EdUval);
                G1cells=Hoechstval>G1minH & Hoechstval<G1maxH & EdUval>G1minE & EdUval<G1maxE;
                %%% Gate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                p21=Alldata(:,5); %global:5 tophat:8
                gatep21=p21>=1; p21(p21<1)=1; p21=log2(p21);
                cycD=Alldata(:,6); %global:6 tophat:9
                gatecycD=cycD>=1; cycD(cycD<1)=1; cycD=log2(cycD);
                Cellgating=G1cells;
                Gating=Cellgating & gatep21 & gatecycD;
                p21=p21(Gating);
                cycD=cycD(Gating);
                switch sigchoice
                    case 'cycD'
                        IFcalc=cycD;
                    case 'p21'
                        IFcalc=p21;
                    case 'cycD-p21'
                        IFcalc=cycD-p21;
                end
                wellstats=[wellstats;nanmean(IFcalc)];
            end
        end
    end
    statarray{i}=wellstats;
end

statmean=ones(uniquecondnum,1)*NaN; statstd=statmean;
for i=1:uniquecondnum
    statmean(i)=mean(statarray{i});
    statstd(i)=std(statarray{i});
    p=ranksum(statarray{1},statarray{i});
    fprintf('pval = %0.4f\n',p);
end

%%% Bar graph & Error bars %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
bar(1:uniquecondnum,statmean,'w');
%set(gca,'XTickLabel',1:uniquecondnum);
%rotateXLabels(gca,45);
hold on;
errorbar(statmean,statstd,'.','linewidth',2,'color','k');
ylim([0 10]);
set(gca,'TickDir','out','box','off');
%set(gca,'TickDir','out','box','off','fontsize',18);
set(gcf,'color','w','PaperPosition',[0 0 5 2]); %[5 5]
saveas(gcf,'D:\Downloads\Bar.jpg');
saveas(gcf,'D:\Downloads\Bar.eps');

end