function [plotMatCDK2 plotMatCDK4 pRb_pos]=CDK2CDK4analysis(conditions,datadir,cellS)
showylabel=1;
%condnum=size(conditions,1);
allnames=conditions(:,1);
[~,uidx]=unique(allnames,'first');
uniquenames=allnames(sort(uidx));
uniquecondnum=numel(uniquenames);
namecount=cell(uniquecondnum,1);
colorcode='krgcmb';
colors=colorcode;
fontsizevar=18; %boxplots:12 histogram:12
if showylabel
    mlvar=0.25; %default 0.25
    %pphvar=2; %4col per slide
    %pphvar=2.5;
    pphvar=4; %2col per slide
    %phvar=6; %1col per slide
else
    mlvar=0.1;
    %pphvar=1.75; %4col per slide
    %pphvar=2;
    pphvar=3.5; %2col per slide
end
pmat=cell(uniquecondnum,1);
boxplotdata=[];
offset=[];
for i=1:uniquecondnum
    colcount=1;
    Alldata=[];
    condrow=find(ismember(conditions(:,1),uniquenames{i}));
    for c=condrow'
        rowmat=cell2mat(conditions(c,2));
        colmat=cell2mat(conditions(c,3));
        sitemat=cell2mat(conditions(c,4));
        for col=colmat
            for row=rowmat
                for site=sitemat
                    %shot=wellnum2str(row,col,site);
                    shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
                    load([datadir,shot,'.mat'],'IFdata');
                    %load([datadir,'IF_',shot,'.mat'],'IFdata');
                    Alldata=[Alldata;IFdata];
                end
            end
            Hoechstval=Alldata(:,3).*Alldata(:,4);
            %histogram(Hoechstval,100);
            %Hoechstval=Alldata(:,3);
            %EdUval=Alldata(:,7);
            %EdUval(EdUval<1)=1;
            %EdUval=log2(EdUval);
            %%% Varying cell cycle gating %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            c=1;
            G1minH=100000; G1maxH=500000; G1minE=3; G1maxE=6.5;
            SminH=200000; SmaxH=600000; SminE=8; SmaxE=12;
            G2minH=500000; G2maxH=1000000; G2minE=3; G2maxE=6.5;
            xylim=[100000 1000000 1 13];
            
            %%% Report Cell-cycle percentages
            G1cells=Hoechstval>G1minH(c) & Hoechstval<G1maxH(c);
            G2cells=Hoechstval>G2minH(c) & Hoechstval<G2maxH(c);
            Allcells=Hoechstval>G1minH(c) & Hoechstval<G2maxH(c);
            
            %%% CDK2 and CDK4
            CDK2vals=Alldata(:,6)./Alldata(:,5); %CDK2
            posvalCDK2=Alldata(:,5)>70 & CDK2vals>0 & CDK2vals<=2;
            
            CDK4vals=Alldata(:,8)./Alldata(:,7); %CDK4
            posvalCDK4=Alldata(:,7)>120 & CDK4vals>0 & CDK4vals<=2;
            
            %figure; subplot(1,2,1);histogram(Alldata(:,5),100); subplot(1,2,2);histogram(Alldata(:,7),100);
            
            IFvals=Alldata(:,9);
            posval=IFvals>=0; IFvals(IFvals<1)=1; IFvals=log2(IFvals);
            %%% compare subpopulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if cellS==1
                Cellstage=G1cells;
                Cell='G1';
            elseif cellS==2
                Cellstage=Scells;
                Cell='S';
            elseif cellS==3
                Cellstage=G2cells;
                Cell='G2';
            elseif cellS==4
                Cellstage=Allcells;
                Cell='all';
            end
            CDK2vals=CDK2vals(Cellstage & posvalCDK2);
            CDK4vals=CDK4vals(Cellstage & posvalCDK4);
            IFvals=IFvals(Cellstage & posval);
            
            Rbthresh=6;
            plotMatCDK2(i,colcount)=nanmean(CDK2vals);
            plotMatCDK4(i,colcount)=nanmean(CDK4vals);
            pRb_pos(i,colcount)=sum(IFvals>Rbthresh)/length(IFvals)*100;
            colcount=colcount+1;
        end
    end
end