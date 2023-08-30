function [Data CDK2_Rbneg CDK2_Rbpos CDK4_Rbneg CDK4_Rbpos]=Stain_HistComparison_HW_GOI2(conditions,datadir,resultdir,protein,protein2,DisplayOption,cellS,option,option2)
EdUMeasured=0;
showylabel=0;
%condnum=size(conditions,1);
allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);
bmin=4; %cycD=4, p21=6, pRb=8
bmax=30; %cycD=13, p21=13, pRb=13
bstep=(bmax-bmin)/40; %pRb/tRb:100
bin=bmin:bstep:bmax;
[~,idx1]=min(abs(bin-1));
numbins=numel(bin);
binfill=[bin fliplr(bin)];
namecount=cell(uniquecondnum,1);
colorcode='krgcmb';
colors=colorcode;
fontsizevar=8; %2col:8 4col:6
if showylabel
    mlvar=0.25; %default 0.25
    %pphvar=2; %4col per slide
    %pphvar=3;
    pphvar=4; %2col per slide
    %phvar=6; %1col per slide
else
    mlvar=0.1;
    %pphvar=1.75; %4col per slide
    pphvar=2.5;
    %pphvar=3.5; %2col per slide
end

%%% Varying cell cycle gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if option==1
    G1minH=[000000,200000,200000,200000,200000,200000,200000,200000,200000];
        G1maxH=[2100000,1700000,1700000,1200000,1200000,1200000,1200000,1200000];
        G1minE=[2,2,2,2,2,2,2,2]; G1maxE=[8,8,8,8,8,8,8,8];
        SminH=[200000,200000,200000,200000,200000,200000,200000,200000,200000];
        SmaxH=[4000000,3200000,3200000,2400000,2400000,2400000,2400000,2400000];
        SminE=[8,8,8,8,8,8,8,8];
        SmaxE=[16,16,16,16,16,16,16,16];
        G2minH=[2100000,1700000,1700000,1200000,1200000,1200000,1200000,1200000];
        G2maxH=[4000000,3200000,3200000,2400000,2400000,2400000,2400000,2400000];
        G2minE=[2,2,2,2,2,2,2,2]; G2maxE=[8,8,8,8,8,8,8,8];
        xylim=[200000 4000000 2 16];
elseif option==2
    G1minH=[000000,000000,000000,000000,000000,000000,000000,000000,000000];
    G1maxH=[1500000,1500000,1100000,1100000,1100000,1100000,1100000,1500000];
    G1minE=[2,2,2,2,2,2,2,2]; G1maxE=[8,8,8,8,8,8,8,8];
    SminH=[000000,000000,000000,000000,000000,000000,000000,000000,000000];
    SmaxH=[3000000,3000000,2500000,2500000,2500000,2500000,2500000,3000000];
    SminE=[8,8,8,8,8,8,8,8];
    SmaxE=[16,16,16,16,16,16,16,16];
    G2minH=[1500000,1500000,1100000,1100000,1100000,1100000,1100000,1500000,1500000];
    G2maxH=[3000000,3000000,2500000,2500000,2500000,2500000,2500000,3000000,3000000];
    G2minE=[2,2,2,2,2,2,2,2]; G2maxE=[8,8,8,8,8,8,8,8];
    xylim=[00000 3000000 2 16];
    elseif option==3
        G1minH=[000000,200000,000000,000000,000000,000000,000000,000000,000000];
        G1maxH=[2100000,1700000,3000000,3000000,1100000,1100000,1100000,1500000];
        G1minE=[2,2,2,2,2,2,2,2]; G1maxE=[8,8,8,8,8,8,8,8];
        SminH=[000000,200000,000000,000000,000000,000000,000000,000000,000000];
        SmaxH=[4000000,3200000,3000000,3000000,2500000,2500000,2500000,3000000];
        SminE=[8,8,8,8,8,8,8,8];
        SmaxE=[16,16,16,16,16,16,16,16];
        G2minH=[2100000,1700000,1500000,1500000,1100000,1100000,1100000,1500000,1500000];
        G2maxH=[4000000,3200000,3000000,3000000,2500000,2500000,2500000,3000000,3000000];
        G2minE=[2,2,2,2,2,2,2,2]; G2maxE=[8,8,8,8,8,8,8,8];
        xylim=[00000 4000000 2 16];

end

switch DisplayOption
    case {'Panel'}
        figure('Position',[100 100 250 900]);hold on;
    case 'Scatter'
        
    case 'boxplot'
        figure('Position',[100 100 800 400]);hold on;
    case 'plot'
        figure; hold on;
    case '3D_hist'
        figure; hold on;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pmat=cell(uniquecondnum,1);
for i=1:uniquecondnum
    Alldata=[];
    condrow=find(ismember(conditions(:,1),uniquenames{i}));
    for c=condrow'
        rowmat=cell2mat(conditions(c,2));
        colmat=cell2mat(conditions(c,3));
        sitemat=cell2mat(conditions(c,4));
        for row=rowmat
            for col=colmat
                for site=sitemat
                    %shot=wellnum2str(row,col,site);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    
                    %[rowstring,colstring,sitestring]=wellnum2strRCS_3(row,col,site);
                    %shot=[rowstring,colstring,'_',sitestring];
                    if exist([datadir,'IF_',shot,'.mat'])
                        %load([datadir,shot,'.mat'],'IFdata');
                        load([datadir,'IF_',shot,'.mat'],'IFdata');
                        Alldata=[Alldata;IFdata];
                    else
                        continue
                    end
                end
            end
        end
    end
    Hoechstval=Alldata(:,3).*Alldata(:,4);
    %histogram(Hoechstval,100); vline([G1minH(i) G2maxH(i)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if EdUMeasured==1
        %EdUval=log2(Alldata(:,3).*Alldata(:,6));
        %%%%%%%%%%%%%%change%%%%%%%%%%%%%%%%%%%%%%%%
        EdUval=Alldata(:,6);
        EdUval(EdUval<1)=1;
        EdUval=log2(EdUval);
        %%%%%% visualize Hoechst vs EdU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %figure;dscatter(Hoechstval,EdUval);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %drawHoechstEdUgates(Hoechstval,EdUval,xylim,G1minH,G1maxH,G1minE,G1maxE);
        %drawHoechstEdUgates(Hoechstval,EdUval,xylim,SminH,SmaxH,SminE,SmaxE);
        %drawHoechstEdUgates(Hoechstval,EdUval,xylim,G2minH,G2maxH,G2minE,G2maxE);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Allcells=Hoechstval>xylim(1) & Hoechstval<xylim(2) & EdUval>xylim(3) & EdUval<xylim(4);
        G1cells=Hoechstval>G1minH(i) & Hoechstval<G1maxH(i) & EdUval>G1minE(i) & EdUval<G1maxE(i);
        Scells=Hoechstval>SminH(i) & Hoechstval<SmaxH(i) & EdUval>SminE(i) & EdUval<SmaxE(i);
        G2cells=Hoechstval>G2minH(i) & Hoechstval<G2maxH(i) & EdUval>G2minE(i) & EdUval<G2maxE(i);
    else
        %figure;histogram(Hoechstval); vline(G1minH(i),'r'); vline(G2maxH(i),'k');
        Allcells=Hoechstval>G1minH(1) & Hoechstval<G2maxH(1);
        G1cells=Hoechstval>G1minH(1) & Hoechstval<G1maxH(1);
        G2cells=Hoechstval>G1maxH(1) & Hoechstval<G2maxH(1);
    end
    %%double values
    %Abvalsx=Alldata(:,5);
    %Abvals=Alldata(:,9)./Alldata(:,8); %pRb/tRb
    %posval=Abvals>=0;
    %posvalx=Abvalsx>=0; %!!
    
    switch DisplayOption %(median) 5; p53, 6; p21, 7; EdU
        case {'Panel','plot','boxplot','3D_hist'}
            %%single value
            if option2==1
                Abvals=Alldata(:,protein)./Alldata(:,protein2);
                posval=Abvals>=0;
            elseif option2==2
                Abvals=Alldata(:,protein);%./Alldata(:,protein);
                posval=Abvals>=0; Abvals(Abvals<1)=1; Abvals=log2(Abvals);
            elseif option2==3
                Abvals=Alldata(:,protein2);
                posval=Abvals>=0;
                elseif option2==4
                Abvals=Alldata(:,protein2);
                posval=Abvals>=0;
                
                p21=Alldata(:,protein);%./Alldata(:,protein);
                posval=p21>=0; p21(p21<1)=1; p21=log2(p21);
                
                low_p21=p21<7.5;
                middle_p21=p21>8 & p21<10;
                high_p21=p21>11;
            end
%             %correct?
%                         CDK2=Alldata(:,5)./Alldata(:,4);
%                         CDK2gate=Alldata(:,4)>0;
%             
%                         CDK4=Alldata(:,7)./Alldata(:,6);
%                         CDK4gate=Alldata(:,6)>0;
%             
%                         CDK4_corr=CDK4_activity_correction_allcells(CDK4,CDK2);
        case 'Scatter'
            if protein==21
                Abvalsx=Alldata(:,3).*Alldata(:,4); %Hoechst
                posvalx=Abvalsx>=0;
                
                Abvals=Alldata(:,protein2);%./Alldata(:,protein);
                posval=Abvals>=0; Abvals(Abvals<1)=1; Abvals=log2(Abvals);
            elseif option2==1 & protein==11
                Abvals=Alldata(:,protein2)./Alldata(:,protein);
                posval=Abvals>=0;
            elseif option2==2 & protein==11
                Abvalsx=Alldata(:,6)./Alldata(:,5);
                posvalx=Abvalsx>=0 & Alldata(:,5)>250;
            elseif option2==3 & protein==11
                Abvals=Alldata(:,protein2)./Alldata(:,protein);
                posval=Abvals>=0;
                
                Abvalsx=Alldata(:,6)./Alldata(:,5);
                posvalx=Abvalsx>=0 & Alldata(:,5)>250;
            end
    end
    %%% compare subpopulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cellS==1
        Cellstage=G1cells;
        CellPhase='G1';
    elseif cellS==2
        Cellstage=Scells;
        CellPhase='S';
    elseif cellS==3
        Cellstage=G2cells;
        CellPhase='G2';
    elseif cellS==4
        Cellstage=Allcells;
        CellPhase='all';
    elseif cellS==5
        Cellstage=G1cells | G2cells;
        CellPhase='G1_G2';
    end
    
    switch DisplayOption
        case {'Panel','plot','boxplot'}
            if option2==4
                gating_1=Cellstage & posval & low_p21;% & CDK2gate & CDK4gate;
                gating_2=Cellstage & posval & middle_p21;
                gating_3=Cellstage & posval & high_p21;
                
                cycD_1=Abvals(gating_1);
                cycD_2=Abvals(gating_2);
                cycD_3=Abvals(gating_3);
            else
                gating=Cellstage & posval;% & CDK2gate & CDK4gate;
                Abvals=Abvals(gating); %panel
            end
            %CDK2=CDK2(gating);
            %CDK4=CDK4(gating);
            %CDK4_corr=CDK4_corr(gating);
            %pRb=pRb(Cellstage & posval);
        case 'Scatter'
            Abvals=Abvals(Cellstage & posval & posvalx); %scatter Abvals=Abvals(Cellstage & posval & posvalx);
            Abvalsx=Abvalsx(Cellstage & posval & posvalx); %scatter Abvalsx=Abvalsx(Cellstage & posval & posvalx);
            [R P]=corrcoef(Abvals',Abvalsx');
            Rval{i}=R(1,2);
            Pval{i}=P(1,2);
    end
    pdfvals=histc(Abvals,bin);
    pdfvals=100*pdfvals/sum(pdfvals);
    cdfvals=cumsum(pdfvals);
    cdfvals=100*cdfvals/max(cdfvals);
    cdfvals(cdfvals==max(cdfvals))=99;
    plotMat{i}=Abvals;
    plotV(i)=nanmean(plotMat{i},1);
    plotV2(i)=nanmedian(plotMat{i},1);
    errV(i)=nanstd(plotMat{i});%Standard deviation    /sqrt(size(plotMat{1},1));
    %SEMV(i)=errV(i)/sqrt(72);
    cellNum(i)=size(plotMat{i},1);
    pmat{i}=Abvals;
    %Data=Abvals;
    
    switch DisplayOption
        case 'Panel'
            %figure(99)
            %subaxis(1,uniquecondnum,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02); %MB 0.15
            %subaxis(uniquecondnum,1,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02); %vertical
            
            %%%%%%%%%%%%STOP BELOW AND TEST THIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %histogram(Abvals,100);
            
            if option2==1
                if ismember(option,[1 3])
                    bmin=0; bmax=8; ymax=0.1;
                    thresh=4;
                elseif ismember(option,[5 7])
                    bmin=0; bmax=6; ymax=0.1;
                    thresh=4;
                elseif ismember(option,[2])
                    bmin=1; bmax=15; ymax=0.1;
                    thresh=4;
                    
                   elseif ismember(option,[4 8])
                    bmin=2; bmax=50; ymax=0.1;
                    thresh=4; 
                   elseif ismember(option,[6])
                    bmin=1; bmax=25; ymax=0.1;
                    thresh=4; 
                end

            elseif option2==2
                if option==1
                    bmin=7; bmax=13;  ymax=0.07; thresh=10;
                elseif option==2
                    bmin=7; bmax=15;  ymax=0.1; thresh=11;
                end
            end
            if cellS==2
                binNum=40;
            else
                binNum=60;
            end
            subplot(1,uniquecondnum,i);
            [h]=histogram(Abvals,binNum,'Normalization','probability','BinLimits',[bmin,bmax]);
            if ismember(option,[1 2]);
                cMyclow(i)=sum(h.Values(h.BinEdges>=0 & h.BinEdges<=thresh));
                cMychigh(i)=sum(h.Values(h.BinEdges>thresh & h.BinEdges<h.BinLimits(2)));

                vline(thresh); title(['c-Myc low=',num2str(round(cMyclow(i),2)),' c-Myc high=',num2str(round(cMychigh(i),2))],'FontSize',9); %' cellNum=',num2str(cellNum(i))],
                xlim([bmin bmax]);
                ylim([0 ymax]);
            end
            xlim([bmin,bmax])
            %%%%%%%%%%%%%%%%%%%%%%%%%change below%%%%%%%%%%%%%%%%%%%%%%%%%

                if ismember(option,[11]);
                    Data(i)=cMychigh(i);
                else
                    %Data(i)=nanmean(Abvals);
                    Data{i}=Abvals;
                end

            %             ylim([0 ymax]); xlim([bmin bmax]);
            
            %gateRb=Abvals>thresh;
            %CDK4_pRbpos(i)=nanmean(CDK4_corr(gateRb));
            %CDK4_pRbneg(i)=nanmean(CDK4_corr(~gateRb));
            
            %CDK2_pRbpos(i)=nanmean(CDK2(gateRb));
            %CDK2_pRbneg(i)=nanmean(CDK2(~gateRb));
            
            %mean_Rbpos(i)=mean(Abvals(gateRb));
            %mean_pRbneg(i)=mean(Abvals(~gateRb));
            %disp('hi');
            
            %             figure(100); hold on;
            %             subaxis(uniquecondnum,1,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02); %vertical
            %             histogram(CDK4_pRbpos,40,'Normalization','probability','BinLimits',[0,1]);
            %             histogram(CDK4_pRbneg,40,'Normalization','probability','BinLimits',[0,1]);
            %             xlim([0 1])
            
            %             figure(101); hold on;
            %             subaxis(uniquecondnum,1,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02); %vertical
            %             histogram(CDK2_pRbpos,40,'Normalization','probability','BinLimits',[0.2,1.3]);
            %             histogram(CDK2_pRbneg,40,'Normalization','probability','BinLimits',[0.2,1.3]);
            %             xlim([0.2 1.3])
            
            %%% draw 2N DNA content boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %hold on;
            %line([G1maxH(i) G1maxH(i)],[0 10],'color','k','linewidth',2,'linestyle','--');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if showylabel
                y=ylabel(char(uniquenames(i)),'rot',0,'fontsize',fontsizevar);
                set(y,'Units','Normalized','Position',[-0.2,0.4,0]);
            end
            
            if option==1000
                Data=0;
                CDK2=Alldata(:,6)./Alldata(:,5);
                gatex=CDK2>=0 & Alldata(:,5)>250;
                CDK4=Alldata(:,8)./Alldata(:,7);
                gatey=CDK4>=0 & Alldata(:,7)>100;
                
                pRb=Alldata(:,10)./Alldata(:,11);
                gateRb=pRb>=0; %pRb(pRb<1)=1; pRb=log2(pRb);
                %thresh=9.8;
                thresh=0.55;
                
                gating=Cellstage & gatex & gatey & gateRb;
                
                CDK2=CDK2(gating); CDK4=CDK4(gating); pRb=pRb(gating);
                CDK4_corr=CDK4_activity_correction_allcells(CDK4,CDK2);
                
                Rb_neg=pRb<thresh; Rb_pos=pRb>thresh;
                
                CDK2_Rbneg{i}=CDK2(Rb_neg); CDK2_Rbpos{i}=CDK2(Rb_pos);
                CDK4_Rbneg{i}=CDK4_corr(Rb_neg); CDK4_Rbpos{i}=CDK4_corr(Rb_pos);
                
                figure; subplot(1,2,1);
                if length(CDK2_Rbneg{i})>5
                    %dscatter(CDK2_Rbneg{i},CDK4_Rbneg{i},'MSIZE',8); %colormap(parula_gradwhite) % 60
                    scatter(CDK2_Rbneg{i},CDK4_Rbneg{i});
                    axis([0 1.5 0 1]);
                end
                subplot(1,2,2);
                %dscatter(CDK2_Rbpos{i},CDK4_Rbpos{i},'MSIZE',8); colormap(parula_gradwhite) % 60
                scatter(CDK2_Rbpos{i},CDK4_Rbpos{i});
                %figure; dscatter(CDK2,CDK4_corr);
                axis([0 1.8 0 1]);
                %xlabel(xname); ylabel('p-Rb (S807/811)');
            end
        case 'Scatter'
            %figure('Position',[100 100 800 400]);hold on;
            Data(i)=nanmean(Abvals);
            
            subaxis(1,uniquecondnum,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02); %MB 0.15
            %subaxis(1,uniquecondnum,i,'ML',0.1,'MR',0.03,'MT',0.02,'MB',0.1,'SH',0.02); %MB 0.15
            if protein==21
                xname='DNS content';
            elseif protein==7
                xname='CDK4 activity';
                xname='CDK2 activity';
            end
            
            
            binNum=30; xrange=[0.2 1.6]; xthresh=0.75;
            
            %yrange=[0 1.2]; ythresh=0.5;
            %subplot(3,2,1);
            %[h]=histogram(Abvalsx,binNum,'Normalization','probability','BinLimits',xrange); xlim(xrange); ymax=ceil(max(h.Values)*100)/100; ylim([0 ymax]);
            %high(i)=sum(h.Values(h.BinEdges>=0 & h.BinEdges<=xthresh));
            %low(i)=sum(h.Values(h.BinEdges>xthresh & h.BinEdges<h.BinLimits(2)));
            %vline(xthresh); title(['high=',num2str(round(high(i),2)),' low=',num2str(round(low(i),2)),' cellNum=',num2str(cellNum(i))],'FontSize',9);
            
            %subplot(3,2,2);
            %[h]=histogram(Abvals,binNum,'Normalization','probability','BinLimits',yrange); xlim(yrange); ymax=ceil(max(h.Values)*100)/100; ylim([0 ymax]);
            %high(i)=sum(h.Values(h.BinEdges>=0 & h.BinEdges<=ythresh));
            %low(i)=sum(h.Values(h.BinEdges>ythresh & h.BinEdges<h.BinLimits(2)));
            %vline(ythresh); title(['high=',num2str(round(high(i),2)),' low=',num2str(round(low(i),2)),' cellNum=',num2str(cellNum(i))],'FontSize',9);
            
            %subplot(3,2,[3:6]); hold on;
            dscatter(Abvalsx,Abvals,'MSIZE',8); %colormap(parula_gradwhite) % 60
            xlabel(xname); ylabel('cyclin D1');
            ylim([4 15])
            %%% annotate correlation between two data%%%%%%%%%%%%%%%%%%%%%%
            %             textpos=[0.1+((i-1)*0.3) 0.9 0.1 0.1];
            %             annotation('textbox',textpos,'String',['R=',num2str(Rval{i}),'  p=',num2str(Pval{i})]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            title(char(uniquenames(i)));
            
            
        case 'boxplot'
            subaxis(1,uniquecondnum,i,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02); %MB 0.15
            boxplot(Abvals,'notch','on','symbol','o');
            %hyperRb=pRb>11;
            %hypoRb=pRb>2 & pRb<8;
            %subaxis(1,2,1,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02);
            %boxplot(Abvals(hypoRb),'notch','on','symbol','o'); ylim([0 100]);
            %subaxis(1,2,2,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02);
            %boxplot(Abvals(hyperRb),'notch','on','symbol','o'); ylim([0 100]);
            %             set(gca,'XTickLabel',{char(uniquenames(i))});
            
            %             binNum=50;
            %             subaxis(1,2,1,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02);
            %             histogram(Abvals(hypoRb),binNum,'Normalization','probability');
            %             subaxis(1,2,2,'ML',0.1,'MR',0.03,'MT',0.1,'MB',0.15,'SH',0.02);
            %             histogram(Abvals(hyperRb),binNum,'Normalization','probability');
    end
    set(gca,'fontsize',fontsizevar);
    namecount{i}=char(uniquenames(i));
end
% p=ranksum(pmat{1},pmat{7}); fprintf('p value=%8.10f\n',p);
% p2=ranksum(pmat{1},pmat{10}); fprintf('p value=%8.10f\n',p2);
%p(1)=ranksum(pmat{1},pmat{2}); p(2)=ranksum(pmat{1},pmat{3}); p(3)=ranksum(pmat{1},pmat{4}); p(4)=ranksum(pmat{1},pmat{5});
%fprintf('mean(set1)=%8.2f  mean(set2)=%8.2f\n',mean(pmat{1}),mean(pmat{2}));

% condi={'pRb (807/811)' 'CyclinD1' 'rH2AX' 'p53' 'p21'};

switch DisplayOption
    case '3D_hist'
        if ismember(protein,[5 7])
            bmin=0; bmax=1.5;
        elseif protein==10
            bmin=0; bmax=100;
            %         elseif protein==4
            %             bmin=0; bmax=6000;
        else
            bmin=3; bmax=15;
        end
        num_bins=50; offset_increment=0.1; color1=[1 0 0]; color2=[0 0 1]; xlims=[]; mthd='histogram';
        histogram_3D_v3(plotMat,bmin,bmax,num_bins,offset_increment,color1,color2,xlims,mthd)
        
    case 'Panel'
        %        if option==3
        %save([resultdir CellPhase '_pRb_tRb.mat'],'plotMat','uniquenames')
        %             xlabel('p21 log2(medianRFU)');
        %             hxl=get(gca,'xlabel');
        %             set(hxl,'fontsize',fontsizevar);
        %             set(gcf,'color','w','PaperPosition',[0 0 pphvar 4]);
        %             export_fig([resultdir 'CyclinD1_' CellPhase '_.eps'],'-eps','-transparent','-nocrop');
        %             export_fig([resultdir 'CyclinD1_' CellPhase '_.png']);
        %        end
    case 'Scatter'
        set(gcf,'color','w','PaperPosition',[0 0 15 3]);
        %         save([resultdir 'G2_Corr#6-CyclinD.mat'],'Rval','Pval');
        %         if staining==1
        %             if protein==18
        %                 xlabel('Number of nuclear 53BP1 folci'); ylabel('Number of nuclear p-rH2AX folci');
        %                 export_fig([resultdir '53BP1foci-rH2AXfoci_' CellPhase '_Correlation.eps'],'-eps','-transparent','-nocrop');
        %             elseif protein==6
        %                 xlabel('Intensity of nuclear 53BP1'); ylabel('Number of nuclear p-rH2AX folci');
        %                 export_fig([resultdir '53BP1intensity-rH2AXfoci_' CellPhase '_Correlation.eps'],'-eps','-transparent','-nocrop');
        %             end
        %         elseif staining==3
        %             if protein==6
        %                 xlabel('Intensity of nuclear 53BP1'); ylabel('p21');
        %                 export_fig([resultdir '53BP1intensity-p21_' CellPhase '_Correlation.eps'],'-eps','-transparent','-nocrop');
        %             elseif protein==12
        %                 xlabel('Number of nuclear 53BP1 folci'); ylabel('p21');
        %                 export_fig([resultdir '53BP1foci-p21_' CellPhase '_Correlation.eps'],'-eps','-transparent','-nocrop');
        %             end
        %         elseif staining==5
        %             if protein==6
        %                 xlabel('Intensity of nuclear 53BP1'); ylabel('p53');
        %
        %             elseif protein==12
        %                 xlabel('Number of nuclear 53BP1 folci'); ylabel('p53');
        %                 export_fig([resultdir '53BP1foci-p53_' CellPhase '_Correlation.eps'],'-eps','-transparent','-nocrop');
        %             end
        %         end
        
        % xlabel('Hoechst'); ylabel('p21');
        % export_fig([resultdir 'HoechstVSp21.eps'],'-eps','-transparent','-nocrop');
    case 'plot'
        colors={[0 0 0],[0.8 0 0],[0 0.7 0],[0 0 1],[0.5 0.5 0],[0 0.5 0.5]};
        timeVals=[1:6];
        %         errorbar(timeVals,plotV,SEMV,'Color',colors{1});%errV
        plot(timeVals,plotV);
        ylabel('p21 mRNA'); %ylim([0 1.1]);
        
        xlabel('inhibtion');
        hxl=get(gca,'ylabel');
        set(hxl,'fontsize',fontsizevar);
        set(gcf,'color','w','PaperPosition',[0 0 pphvar 4]);
        
        if option==1
            if protein==8
                save([resultdir CellPhase '_53BP1_area.mat'],'plotMat','uniquenames')
            elseif protein==11
                save([resultdir CellPhase '_gH2AX_area.mat'],'plotMat','uniquenames')
            end
        elseif option==2
            if protein==5
                save([resultdir CellPhase '_p53.mat'],'plotMat','uniquenames')
            elseif protein==6
                save([resultdir CellPhase '_p21.mat'],'plotMat','uniquenames')
            end
        elseif option==4
            if protein==5
                save([resultdir CellPhase '_pChk1.mat'],'plotMat','uniquenames')
            end
        elseif option==5
            if protein==5
                save([resultdir CellPhase '_pChk2.mat'],'plotMat','uniquenames')
            end
        elseif option==6
            if protein==5
                save([resultdir CellPhase '_cycD1.mat'],'plotMat','uniquenames')
            end
        elseif option==7
            if protein==5
                save([resultdir CellPhase '_p27.mat'],'plotMat','uniquenames')
            end
        elseif option==8
            if protein==5
                save([resultdir CellPhase '_p57.mat'],'plotMat','uniquenames')
            end
        end
    case 'boxplot'
        
end
%textpos=[0.15 0.55 0.3 0.15];
%annotation('textbox',textpos,'String',['p-value = ',num2str(p)]);
%saveas(gcf,'h:\Downloads\Fig.jpg');
%%%
end