%%% INSTRUCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Set file paths, cmosoffsetfile, separatedirectories option, file root
%    names, ringcalc option, and StartFrame & EndFrame.
% 2. Set synchronization under tracking parameters.
% 3. Choose firstsegmethod & secondsegmethod. Default: concavity
% 4. Visually confirm acceptable segmentation.
% 5. Visually confirm bias correction & background subtraction:
%       a. In command line, type imtool(imreal1,[]), e.g. for 1st channel.
%       b. Click the contrast icon (black and white circle at top).
%       c. Drag down upper saturation limit to see how well the bias has
%          been corrected by assessing the background uniformity.
%       d. Visually confirm that the background is approximately zero.
%       e. Repeat for all channels and frames of interest.
% 6. Visually confirm tracking.
% 7. Set the desired parameterization of signals (e.g. mean vs median)
%    under the 'compile data' section.
% 8. If 'ringcalc' was chosen, set the ringthresh. Default=0.
% 9. Create a Data folder and run a few frames to make sure the final
%    dataset gets properly saved.
%%% TROUBLESHOOTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Segmentation errors:
%       a. Visually confirm segmentation after each line in the code that
%          is commented "Troubleshoot segmentation". For the first frame,
%          if errors arise after:
%           i. first selected algorithm: try another algorithm.
%           ii. bwareaopen: increase debrisarea.
%           iii. bridgedeflections: adjust nuclear radius nucr.
%           iv. excludelargeandwarped: adjust boulderarea &
%               eccentricitythresh (0.95-strict, 0.85-loose).
% 2. Inaccurate illumination bias correction (non-uniform background):
%       a. If you see an overcorrection at the edge, this could be
%          minimized by increasing tilenum in the biascalc script.
% 3. False-positive badframe calls (determined retrospectively by watching
%    movie:
%       a. Increase blurthreshhigh, decrease blurthreshlow, increase
%          numthresh.
function Process_2_Timelapse(row,col,site)
%row=6;col=7;site=4;

FileOption='tif';
ImageOption='multisite';

binning=2;
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='D:\';
imagepath='E:\';
experimentpath='Folder\';

%shot=wellnum2str(row,col,site);
saveshot=[num2str(row),'_',num2str(col),'_',num2str(site)];
biasdir=[imagepath,experimentpath,'Bias\'];
datadir=[projectpath,experimentpath,'Data\'];
if ~exist(datadir,'dir')
    mkdir(datadir);
end
separatedirectories=1;
if separatedirectories==1
    rawdir=[imagepath,experimentpath,'Raw\',saveshot,'\'];
    maskdir=[imagepath,experimentpath,'Mask\',saveshot];
else
    rawdir=[imagepath,experimentpath,'Raw\',shot,'_'];
    maskdir=[imagepath,experimentpath,'Mask'];
end
maskwrite=0;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end
if separatedirectories==0
    maskdir=[maskdir,'\',shot,'_'];
end

% if exist([datadir,'tracedata_',shot,'.mat'])
%     return
% end
%%% general settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SF=1;
pre_frame=100;
EF=pre_frame+105; %1:39
ringcalc=1; %set to zero if ring is unneeded
name1='H2B_'; %nuclear channel
name2='ERK_';
name3='CDK2_';
name4='CDK4_';
namenucedge='nucedge_';
%%% segmentation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=13; %WM983B
eccentricitythresh=0.90; 
debrisarea=200;
boulderarea=1500;
blobthreshold=-0.03;
%%% tracking parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxjump=nucr*3; %4
masschangethreshold=0.20;
areachangethreshold=0.60;
daughtervariance=0.10;
trackparams={nucr,maxjump,debrisarea,masschangethreshold,areachangethreshold,daughtervariance};
%%% measurement parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
punctathresh=250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frames=SF:EF;
totalframes=numel(frames);
badframes=ones(EF,1)*NaN;
if SF>1
    badframes(1:SF-1)=0;
end
jitters=zeros(EF,2);
blocksize=10000;
maxcellnum=blocksize;
parameternum=11;
jitter_size=1;
tracedata=ones(maxcellnum,EF,parameternum)*NaN;
tracking=ones(maxcellnum,5)*NaN;
timetotal=tic;

%rawdir = [imagepath,experimentpath,'Raw\1\'];
%[firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_Nikon(row,col,site,rawdir,signame{1},frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
[firstgoodindex,blurthreshhigh,blurthreshlow,numthresh,badframes,height,width]=timelapsesetup_4(rawdir,name1,frames,nucr,blobthreshold,debrisarea,badframes,maskwrite);
regheight=1:jitter_size*height; regwidth=1:jitter_size*width;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load(cmosoffsetfile,'cmosoffset'); %load CMOS offset
load([biasdir,'H2B.mat']); bias1=bias;
load([biasdir,'ERK.mat']); bias2=bias;
load([biasdir,'CDK2.mat']); bias3=bias;
%load([biasdir,name4,num2str(site),'.mat']); bias4=bias;
for i=firstgoodindex:totalframes
    f=frames(i); fprintf('frame %0.0f\n',f);
    timeframe=tic;
    %%% read images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist([rawdir,name1,num2str(f),'.tif'],'file')
        rawImage{1}=double(imread([rawdir,name1,num2str(f),'.tif'])); rawImage{1}=rawImage{1}./bias1;
        rawImage{2}=double(imread([rawdir,name2,num2str(f),'.tif'])); rawImage{2}=rawImage{2}./bias2;
        rawImage{3}=double(imread([rawdir,name3,num2str(f),'.tif'])); rawImage{3}=rawImage{3}./bias3;
        rawImage{4}=double(imread([rawdir,name4,num2str(f),'.tif'])); %raw4=(raw4-cmosoffset)./bias4;
    else
        fprintf('badframe: frame %0.0f\n',f);
        badframes(f)=1;
        continue;
    end
    
    %%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load([datadir,'bleedthroughrate_TxRedtoCy5.mat'],'bleedthroughrate');
    %bleedthroughrate(2)=0.03;
    %rawImage{1}bleedthrough=rawImage{4}*bleedthroughrate(2)+bleedthroughrate(1);
    %rawImage{1}=rawImage{1}-rawImage{1}bleedthrough;
    %%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==firstgoodindex
        firstsegmethod='log'; %single
        switch firstsegmethod
            case 'concavity'
                blurradius=1;
                blur=imfilter(rawImage{1},fspecial('gaussian',blurradius),'symmetric');
                nuc_mask=ThreshImage(blur);
                nuc_mask=markershed(nuc_mask,round(nucr*2/3)); %nucr*2/3
            case 'log'
                nuc_mask=blobdetector_4(log(rawImage{1}),nucr,blobthreshold,debrisarea);
            case 'single'
                blurradius=3;
                nuc_mask=threshmask(rawImage{1},blurradius);
                nuc_mask=markershed(nuc_mask,round(nucr*2/3));
            case 'double'
                blurradius=3;
                nuc_mask=threshmask(rawImage{1},blurradius);
                nuc_mask=markershed(nuc_mask,round(nucr*2/3));
                nuc_mask=secondthresh(rawImage{1},blurradius,nuc_mask,boulderarea*2);
        end
        foreground=nuc_mask;
        nuc_mask=bwareaopen(nuc_mask,debrisarea);
        %%% Deflection-Bridging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nuc_mask=segmentdeflections_bwboundaries(nuc_mask,nucr,debrisarea);
        nuc_mask=excludelargeandwarped_3(nuc_mask,boulderarea,eccentricitythresh);
    else
        nuc_mask=threshmask(rawImage{1},1);
        %%% calculate jitter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lastgoodframe=find(badframes==0,1,'last');
        [reljitx,reljity]=registerimages(imfill(extractmask(regheight,regwidth),'holes'),nuc_mask(regheight,regwidth));
        jitters(f,:)=jitters(lastgoodframe,:)+[reljitx,reljity];
        secondsegmethod='log';
        switch secondsegmethod
            case 'log'
                nuc_mask=blobdetector_4(log(rawImage{1}),nucr,blobthreshold,debrisarea);
            case 'double'
                nuc_mask=threshmask(rawImage{1},blurradius);
                nuc_mask=markershed(nuc_mask,round(nucr*2/3));
                nuc_mask=secondthresh(rawImage{1},blurradius,nuc_mask,boulderarea*2);
            case 'apriori'
                %nuc_mask=blobdetector_4(log(rawImage{1}),nucr,blobthreshold,debrisarea);
                [nuc_mask,marker_mask]=apriori_markermask(nuc_mask,nuc_center,jitters(f,:));
        end
        foreground=nuc_mask;
        nuc_mask=bwareaopen(nuc_mask,debrisarea);
    end
    %%% remove border objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_mask=nuc_mask;
    antiborder_mask=imclearborder(nuc_mask);
    border_mask=nuc_mask-antiborder_mask;
    nuc_mask=logical(all_mask-border_mask);
    %%% background subtract: masked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    compression=4;
    nanmask=imdilate(foreground,strel('disk',nucr));
    nanmaskcyto=imdilate(foreground,strel('disk',nucr*2));
    blur1=imfilter(rawImage{1},fspecial('disk',3),'symmetric');
    blur2=imfilter(rawImage{2},fspecial('disk',3),'symmetric');
    blur3=imfilter(rawImage{3},fspecial('disk',3),'symmetric');
    blur4=imfilter(rawImage{4},fspecial('disk',3),'symmetric');
    
    [real1 bg badframe]=bgsubmasked_global_2(blur1,nanmask,1,compression,50);
    if badframe
        fprintf('badframe: frame %0.0\n',f);
        badframes(f)=1;
        continue;
    end
    real2=bgsubmasked_global_2(blur2,nanmaskcyto,1,compression,50);
    real3=bgsubmasked_global_2(blur3,nanmaskcyto,1,compression,50);
    real4=bgsubmasked_global_2(blur4,nanmaskcyto,1,compression,50);
    
    %true1=bgsubmasked_global_2(rawImage{1},nanmask,1,compression,50);
    %true2=bgsubmasked_global_2(rawImage{2},nanmask,1,compression,50);
    %real2top=imtophat(rawImage{2},strel('disk',2,0));
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nuc_label,numcells]=bwlabel(nuc_mask);
    nuc_info=struct2cell(regionprops(nuc_mask,real1,'Area','Centroid','MeanIntensity')');
    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
    %%%%%% detect bad frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mednuc=median(nuc_area);
    if i>firstgoodindex+1 && (numcells<numthresh || mednuc>blurthreshhigh || mednuc<blurthreshlow)
        fprintf('badframe: frame %0.0f/n',f);
        badframes(f)=1;
        badextractmask=bwmorph(nuc_mask,'remove');
        if maskwrite
            imwrite(uint16(extractmask),[maskdir,'\',namenucedge,num2str(f),'.tif']);
            imwrite(uint16(real1),[maskdir,'\',name,num2str(f),'.tif']);
            imwrite(uint16(real2),[maskdir,'\',name,num2str(f),'.tif']);
        end
        continue;
    end
    blurthreshhigh=1.2*mednuc;
    blurthreshlow=0.7*mednuc;
    numthresh=0.5*numcells;
    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
    nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
    %%% calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nuc_mass=nuc_density.*nuc_area;
    curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i>firstgoodindex
        nuc_center(:,1)=nuc_center(:,1)+jitters(f,1);
        nuc_center(:,2)=nuc_center(:,2)+jitters(f,2);
        %%% temporarily store values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        curdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass];
        debugpackage={extractmask,jitters(lastgoodframe,:),[reljitx,reljity]};
        %%% track & correct merges (update centers, masses and labels) %%%%
        [tracedata,curdata,tracking,nuc_label]=adaptivetrack_9(f,lastgoodframe,f,tracedata,curdata,tracking,real1,nuc_label,jitters(f,:),trackparams,debugpackage);
        badframes(f)=0;
    end
    %%% visualize extractions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    extractmask=bwmorph(nuc_label,'remove');
    if maskwrite
        imwrite(uint16(extractmask),[maskdir,'\',namenucedge,num2str(f),'.tif']);
        imwrite(uint16(true1),[maskdir,'\',name,num2str(f),'.tif']);
        imwrite(uint16(true2),[maskdir,'\',name,num2str(f),'.tif']);
    end
    %%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellid=find(~isnan(curdata(:,1)));
    numlivecells=numel(cellid);
    curdata=curdata(cellid,:);
    nuc_center=curdata(:,[1 2]);
    nuc_area=curdata(:,3);
    nuc_mass=curdata(:,4);
    nuc_info=regionprops(nuc_label,'PixelIdxList');
    %puncta_info=getpuncta(nuc_label,real2top,punctathresh);
    %puncta_mask=real2>punctathresh;
    %%% account for all objects, including border cells %%%%%%%%%%%%%%%%%%%
    maxtrackedID=max(nuc_label(:));
    border_label=bwlabel(border_mask);
    border_label=border_label+maxtrackedID;
    border_label(border_label==maxtrackedID)=0;
    all_label=nuc_label+border_label;
    %%% compile data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nanvec=ones(numlivecells,1)*NaN; sig1=nanvec; sig2=nanvec; sig3=nanvec; sig4=nanvec;%sig2area=nanvec; sig2puncta=nanvec;
    for n=1:numlivecells
        cc=cellid(n);
        sig1(n)=median(real1(nuc_info(cc).PixelIdxList));
        sig2(n)=median(real2(nuc_info(cc).PixelIdxList));
        sig3(n)=median(real3(nuc_info(cc).PixelIdxList));
        sig4(n)=median(real4(nuc_info(cc).PixelIdxList));
        %if (cc<=numel(puncta_info)), sig2area(n)=puncta_info(cc).Area; end
    end
    
    if ringcalc==1
        innerrad=1; outerrad=5; %10xB1|20xB2: 1/5
        %showImagesWithLinkedAxes({real2,real3,real4})
        ringthresh=20;
        sig2ring=sigringcalc(cellid,all_label,innerrad,outerrad,numlivecells,real2,ringthresh);
        sig3ring=sigringcalc(cellid,all_label,innerrad,outerrad,numlivecells,real3,ringthresh);
        sig4ring=sigringcalc(cellid,all_label,innerrad,outerrad,numlivecells,real4,ringthresh);
    end
    %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tracedata(cellid,f,:)=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring,sig3,sig3ring,sig4,sig4ring];
    if maxcellnum-max(cellid)<blocksize
        tempdata=ones(blocksize,EF,parameternum)*NaN;
        temptrack=ones(blocksize,5)*NaN;
        tracedata=[tracedata;tempdata];
        tracking=[tracking;temptrack];
        maxcellnum=maxcellnum+blocksize;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    toc(timeframe);
end
[tracedata,genealogy,jitters]=postprocessing_nolinking(tracedata,cellid,jitters,badframes,tracking,maxcellnum,nucr);
%%% save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([datadir,'tracedata_',saveshot,'.mat'],'tracedata','genealogy','jitters');
toc(timetotal);
clear all; clear mex;
end

function puncta_info=getpuncta(nuclabel,raw,thresh)
puncta_mask=raw>thresh;
puncta_label=nuclabel.*puncta_mask;
puncta_info=regionprops(puncta_label,'Area');
%imtool(raw)
%showImagesWithLinkedAxes({raw,puncta_mask})
end

%{
%%% debugging: view images %%%%%%%%%%
RGB=zeros(height,width,3);
RGB(:,:,1)=imadjust(mat2gray(rawImage{1}));
RGB(:,:,2)=bwperim(nuc_mask); %all_mask
RGB(RGB>1)=1; RGB(RGB<0)=0;
figure,imshow(RGB);

areas=cell2mat(struct2cell(regionprops(nuc_mask,'area'))');
histogram(areas(:),200);
big_mask=bwareaopen(nuc_mask,4000);
small_mask=nuc_mask-big_mask;
%}
