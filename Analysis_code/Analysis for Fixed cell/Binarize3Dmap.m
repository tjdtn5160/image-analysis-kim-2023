function Binarize3Dmap(conditions,datadir,resultdir,option,map)
panel=0; %0:single image
if panel==0
    titlefont=8;
else
    titlefont=16;
end
%condnum=size(conditions,1);
allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);
%colorcode=distributecolormap(jet,condnum);

        G1minH=[900000,1000000,900000,900000,900000,900000,900000,900000];
        G1maxH=[1750000,2000000,2000000,2000000,2000000,2000000,2000000,2500000];
        G1minE=[4.5,4.5,4.5,4.5,4.5,4.5,4.5,4.5]; G1maxE=[9.5,9.5,9.5,9.5,9.5,9.5,9.5,9.5];
        SminH=[1000000,1000000,1000000,1000000,1000000,1000000,1000000,1000000];
        SmaxH=[3000000,3000000,3000000,3000000,3000000,3000000,3000000,3000000];
        SminE=[10.5,10.5,10.5,10.5,10.5,10.5,10.5,10.5,10.5];
        SmaxE=[16,16,16,16,16,16,16,16,16];
        G2minH=[2300000,2500000,2500000,2500000,2200000,2500000,2500000,2800000];
        G2maxH=[3300000,4000000,4000000,4000000,4000000,4000000,4000000,4500000];
        G2minE=[6,6,6,6,6,6,6,6,6,6]; G2maxE=[10,10,10,10,10,10,10,10];
        xylim=[200000 5000000 3 17];
    
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
                    load([datadir,shot,'.mat'],'IFdata');
                    Alldata=[Alldata;IFdata];
                end
            end
        end
    end
    figure;
    Hoechstval=Alldata(:,3).*Alldata(:,4);
    EdUval=Alldata(:,7); EdUval(EdUval<1)=1; EdUval=log2(EdUval);
    
    Allcells=Hoechstval>xylim(1) & Hoechstval<xylim(2) & EdUval>xylim(3) & EdUval<xylim(4);
        G1cells=Hoechstval>G1minH(i) & Hoechstval<G1maxH(i) & EdUval>G1minE(i) & EdUval<G1maxE(i);
        Scells=Hoechstval>SminH(i) & Hoechstval<SmaxH(i) & EdUval>SminE(i) & EdUval<SmaxE(i);
        G2cells=Hoechstval>G2minH(i) & Hoechstval<G2maxH(i) & EdUval>G2minE(i) & EdUval<G2maxE(i);
    
        %CDK2=Alldata(:,6)./Alldata(:,5);
        %        gatepCDK2=CDK2>=0 & Alldata(:,5)>250;
        %        
        %        CDK4=Alldata(:,8)./Alldata(:,7);
        %        gateCDK4=CDK4>=0 & Alldata(:,7)>100;
    
     pRb=Alldata(:,5)./Alldata(:,6);
                gateRb=pRb>=0; %pRb(pRb<1)=1; pRb=log2(pRb);
                %thresh=9.8;
                pRbthresh=0.55;
                
                
    %gate_p21=log2(CDK2);
    %[kval,xval]=ksdensity(gate_p21);
    %MAX_peak=xval(find(kval==max(kval),1,'first'));
    %MAX_peak=2000;
    %histogram(p21,120); xlim([0 10000]); ylim([0 2000]); xlabel('p21 intensity'); ylabel('Number of cells');
    
    %minGate=10000;%(MAX_peak-0.5)
    %maxGate=10000000;%(MAX_peak+0.5)
    %gateMEAN=minGate<CDK2 & maxGate>CDK2;
    %gateMEAN=p21+cycD>500 & p21+cycD<10000;
    
    
    %APC=Alldata(:,4);
    %gateAPC=APC>=0; APC(APC<1)=1; APC=log2(APC);
    
    Cellgating=G2cells;
    Gating=Cellgating & gatepCDK2 & gateCDK4 & gateRb & gateAPC;
    
    %pRbthresh=10; %histogram(pRb); vline(pRbthresh);
    %APCthresh=7.5;
    %CDK2range=[0 1.4];
    %CDK4range=[0 1.2];
    
    %CDK2=CDK2(Gating);
    %CDK4=CDK4(Gating);
    pRb=pRb(Gating);
    APC=APC(Gating);
    %CDK4=CDK4_activity_correction_allcells(CDK4,CDK2);
    
    if option==1
        protein=pRb;
        thresh=pRbthresh;
    elseif option==2
        protein=APC;
        thresh=APCthresh;
    elseif option==3
        protein=mRNA;
    end
    
    if panel==1
        subaxis(1,5,i,'ML',0.02,'MR',0.02,'MT',0.05,'MB',0.05,'SV',0.1,'SH',0.05);
    end
    
    %     cycDrange=[6 14];
    plotoption='prob';
    
    %     figure;histogram(p21,35,'Normalization','probability','BinLimits',[8 15]);
    %     figure;histogram(cycD,35,'Normalization','probability','BinLimits',[3 15]);%ylim([0 0.15])
    %     figure;histogram(pRb);vline(pRbthresh);
    
    switch plotoption
        case '3d'
            pRbmin=5; pRbanchorlow=5; pRbanchorhigh=10; pRbmax=10;
            badvals=pRb<pRbmin | pRb>pRbmax;
            CDK2=CDK2(~badvals);
            CDK4=CDK4(~badvals);
            pRb=pRb(~badvals);
            make3dmap(CDK2,CDK4,pRb,pRbmin,pRbanchorlow,pRbthresh,pRbanchorhigh,pRbmax,CDK2range,CDK4range);
        case 'prob'
            if option==3
                makemeanbmap(CDK2,CDK2range,CDK4,CDK4range,protein)
            else
                makeprobmap(CDK2,CDK2range,CDK4,CDK4range,protein,thresh,map);
            end
        case 'histo'
            pRBgate=pRb>pRbthresh;
            cycDpos=CDK4(pRBgate);
            cycDneg=CDK4(~pRBgate);
            
            p21pos=CDK2(pRBgate);
            p21neg=CDK2(~pRBgate);
            %{
            figure;
             histogram(p21pos,'Normalization','probability')
hold on
histogram(p21neg,'Normalization','probability')
% xlim([5 12])
            xlim([8 14])
xlabel('blue Rb positive Red Rb negative')
ylabel('p21 in G1')
            
            %}
        case 'stoi'
            stoich=CDK4./CDK2;
            binarytransform(stoich,pRb,pRbthresh);
    end
    %     title(uniquenames{i});
end

if panel==0
    paperwidth=4; paperheight=3;
else
    paperwidth=15; paperheight=8;
end
% set(gcf,'color','w','PaperPosition',[0 0 paperwidth paperheight]); %6 5 (small panel: 4 3.3)
% saveas(gcf,'D:\Downloads\Fig.jpg');
% print('-depsc','-tiff','-r300','D:\Downloads\Fig.eps');
end

function make3dmap(valuesx,valuesy,valuesz,minz,anchorlowz,midz,anchorhighz,maxz,rangex,rangey)
rangez=maxz-minz;
normvaluesz=(valuesz-minz)/rangez;
cmapz=normvaluesz*64;
cmapz=round(cmapz); cmapz(cmapz<1)=1; cmapz(cmapz>64)=64;
lowcode=[0 0 1]; highcode=[1 1 0]; midcode=mean([lowcode;highcode]); %yellow:[1 1 0]
anchorlow=round(64*(anchorlowz-minz)/rangez)+1;
anchormid=round(64*(midz-minz)/rangez);
anchorhigh=round(64*(anchorhighz-minz)/rangez);
cmap=makecmap(lowcode,midcode,highcode,anchorlow,anchormid,anchorhigh);
colorsz=cmap(cmapz,:);
scatter(valuesx,valuesy,20,colorsz,'fill'); %default demo:10 panel:20
axis([rangex rangey]);
set(gca,'box','off');
%%% optional: add line %%%
hold on;
plot([5.5 12.5],[5 10.5],'k--');
end

function makeprobmap(valuesx,xrange,valuesy,yrange,valuesz,midz,map)
stepsize=0.05; %default 0.2
halfstep=0.5*stepsize;
minz=min(valuesz);
lowerleft=[min(xrange) min(yrange) minz];
upperright=[max(xrange) max(yrange) minz];
valuesx=[valuesx;lowerleft(1);upperright(1)];
valuesy=[valuesy;lowerleft(2);upperright(2)];
valuesz=[valuesz;minz;minz];
binx=min(valuesx)+halfstep:stepsize:max(valuesx)+halfstep;
biny=min(valuesy)+halfstep:stepsize:max(valuesy)+halfstep;
numbinx=length(binx); numbiny=length(biny);
countmap=ones(numbiny,numbinx)*-1;
densitymap=ones(numbiny,numbinx)*-1;
vbinx=zeros(numel(countmap),1); vbiny=vbinx; vbinz=vbinx;
cc=0;
for i=1:numbinx
    samplex=find(valuesx>=binx(i)-halfstep & valuesx<binx(i)+halfstep);
    for j=1:numbiny
        cc=cc+1;
        vbinx(cc)=binx(i);
        vbiny(cc)=biny(j);
        if isempty(samplex)
            continue;
        end
        samplebinidx=valuesy(samplex)>=biny(j)-halfstep & valuesy(samplex)<biny(j)+halfstep;
        samplebin=samplex(samplebinidx);
        samplecount=numel(samplebin);
        samplepos=sum(valuesz(samplebin)>midz);
        if samplecount<2 %default <1
            continue;
        end
        countmap(j,i)=100*samplepos/samplecount;
        densitymap(j,i)=samplecount;
    end
end
if map==1
    countmap(countmap<3 & countmap>=0)=3;
    countmap(countmap<0)=1; %anywhere there are no datapoints set to background
    zrange=[0 100];
    imagesc(binx,biny,countmap,zrange);
    set(gca,'YDir','normal');
    cmap=summer;
    cmap(1,:)=[0.75 0.75 0.75]; %gray background
    colormap(cmap);
else
    maxCellNum=max(densitymap(:));
    densitymap(densitymap<0)=0; %anywhere there are no datapoints set to background
    densitymap=densitymap/maxCellNum;
    
    imagesc(binx,biny,densitymap);
    set(gca,'YDir','normal');
    dmap=parula;
    dmap(1,:)=[0.75 0.75 0.75];
    colormap(dmap)
end
% cmap(1,:)=[1 1 1];
colorbar;
% xlim([5.5 13]);
axis([xrange yrange]);
set(gca,'tickdir','out');
xlabel('cycE/A-CDK activity'); ylabel('CDK4/6 activity');
% export_fig([resultdir filesep 'FISH_',shot,'_processed.eps'],'-eps','-transparent','-nocrop');
end

function makemeanbmap(valuesx,xrange,valuesy,yrange,valuesz)
stepsize=0.015; %default 0.2
halfstep=0.5*stepsize;
minz=min(valuesz);
lowerleft=[min(xrange) min(yrange) minz];
upperright=[max(xrange) max(yrange) minz];
valuesx=[valuesx;lowerleft(1);upperright(1)];
valuesy=[valuesy;lowerleft(2);upperright(2)];
valuesz=[valuesz;minz;minz];
binx=min(valuesx)+halfstep:stepsize:max(valuesx)+halfstep;
biny=min(valuesy)+halfstep:stepsize:max(valuesy)+halfstep;
numbinx=length(binx); numbiny=length(biny);
countmap=ones(numbiny,numbinx)*-1;
vbinx=zeros(numel(countmap),1); vbiny=vbinx; vbinz=vbinx;
cc=0;
for i=1:numbinx
    samplex=find(valuesx>=binx(i)-halfstep & valuesx<binx(i)+halfstep);
    for j=1:numbiny
        cc=cc+1;
        vbinx(cc)=binx(i);
        vbiny(cc)=biny(j);
        if isempty(samplex)
            continue;
        end
        samplebinidx=valuesy(samplex)>=biny(j)-halfstep & valuesy(samplex)<biny(j)+halfstep;
        samplebin=samplex(samplebinidx);
        samplecount=numel(samplebin);
        samplemean=mean(valuesz(samplebin));
        if samplecount<10 %default <1
            continue;
        end
        countmap(j,i)=samplemean;
    end
end
countmap(countmap<3 & countmap>=0)=3;
countmap(countmap<0)=1; %anywhere there are no datapoints set to background
zrange=[10 55];
imagesc(binx,biny,countmap,zrange);
set(gca,'YDir','normal');
cmap=summer;
cmap(1,:)=[0.75 0.75 0.75]; %gray background
% cmap(1,:)=[1 1 1];
colormap(cmap);
colorbar;
% xlim([5.5 13]);
axis([xrange yrange]);
set(gca,'tickdir','out');
xlabel('cycE/A-CDK activity'); ylabel('CDK4/6 activity');
% export_fig([resultdir filesep 'FISH_',shot,'_processed.eps'],'-eps','-transparent','-nocrop');
end

function binarytransform(input,output,midpoint)
numbins=100;
bmin=0.2;%min(input);
bmax=10;%max(input);
bstep=(bmax-bmin)/numbins;
bin=bmin:bstep:bmax;
bcalc=ones(numbins,1)*NaN;
for i=1:numbins
    inputidx=input>bin(i) & input<bin(i+1);
    ovals=output(inputidx);
    bcalc(i)=sum(ovals>midpoint)/length(ovals);
end
truebin=bin+bstep/2;
scatter(truebin(1:numbins),bcalc,100,'Marker','o','markeredgecolor','k','markerfacecolor',[.49 1 .63],'SizeData',21);
hold on;
fo = fitoptions('Method','NonlinearLeastSquares',...
    'Lower',[0,0],...
    'StartPoint',[0.6,2]);
ft = fittype('x^n/(a^n+x^n)','options',fo);
nanID=~isnan(bcalc);
truebin=truebin(nanID);
bcalc=bcalc(nanID);
[fitData,gof] = fit(truebin',bcalc,ft);
plot(fitData,'m')

% [estimated_params]=sigm_fit(truebin(1:numbins),bcalc(1:numbins)') %[min, max, x50, slope]
% title(['min=' num2str(estimated_params(1)) ' max=' num2str(estimated_params(2)) ' x50=' num2str(estimated_params(3)) ' slope' num2str(estimated_params(4))]);
title([' n=' num2str(fitData.n) ' a=' num2str(fitData.a)]);
axis([bmin bmax -0 1]);
xlabel('cycD1/p21 stoichiometry');
ylabel('% pRb-positive')
set(gcf,'color','w','PaperPosition',[0 0 4 4]);
% saveas(gcf,'h:\Downloads\Fig.jpg');
end