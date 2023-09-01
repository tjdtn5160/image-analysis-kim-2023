%%% INSTRUCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Set file paths, cmosoffsetfile, file root names, separatedirectories
%    option, and ringcalc option.
% 2. Choose segmethod. Default: concavity
% 3. For smear detection, pause before removesmears() and:
%       a. In command line, type imtool(rawImage{1},[]).
%       b. Click the contrast icon (black and white circle at top).
%       c. Reduce upper saturation limit to reveal all foreground signal.
%       d. Approximate an absolute threshold to constitute a smear, and set
%          foregroundthresh to this value.
%       e. Repeat imtool visualization of 'smear-corrected' rawImage{1} to
%          to confirm proper execution.
% 3. Visually confirm acceptable segmentation.
% 4. Visually confirm bias correction & background subtraction (same
%    instructions as Timelapse.m).
% 5. Set the desired parameterization of signals (e.g. mean vs median)
%    under the 'compile data' section.
% 6. If 'ringcalc' was chosen, set the ringthresh. Default=0.
% 7. If multiple rounds of imaging are being combined, uncomment the
%    relevant sections of code.
%%% TROUBLESHOOTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Segmentation and bias solutions are same as Timelapse.m.
% 2. Inaccurate background subtraction (background far from zero):
%       a. Save bias-corrected image (e.g. blurA1 to realA1).
%       b. Save estimated global background: uncomment save command at end.
function Process_3_Immunostain(row,col,site)
% row=6;col=11;site=13;

timetotal=tic;
secondImaging=0;
binning=2; %20x objective
puncta=0;
FileOption='nd2'; %nd2 tif
ImageOption='multisite';

ringcalc=1;
%%% set paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
projectpath='D:\';
imagepath='E:\';
experimentpath='Breast Cancer\';
saveshot=[num2str(row),'_',num2str(col),'_',num2str(site)];

datadir=[projectpath,experimentpath,'Data\'];
if ~exist(datadir,'dir')
    mkdir(datadir);
end
biasdir=[imagepath,experimentpath,'Bias\'];
separatedirectories=0;
if separatedirectories
    maskdir = [imagepath,experimentpath,'Mask\',saveshot];
else
    rawdir = [imagepath,experimentpath,'Raw\'];
    maskdir = [imagepath,experimentpath,'Mask\'];
end
maskwrite=0;
if ~exist(maskdir,'dir') && maskwrite
    mkdir(maskdir);
end

%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucedgename='nucedge';
%%% segmentation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nucr=12; %
eccentricitythresh=0.90;
debrisarea=250;
boulderarea=1200;

nameA1='Hoechst';
nameA2='EGFP';
nameA3='CY5';

%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load(cmosoffsetfile,'cmosoffset'); %load CMOS offset bgcmos=shadingcorrection;
% [height,width]=size(cmosoffset);
% load([biasdir,nameA1,'_',num2str(site),'.mat']); biasA1=bias;
% load([biasdir,nameA2,'_',num2str(site),'.mat']); biasA2=bias;
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAXintensity=2^16-1; %16 bit image
sigThresh=10^4;

switch FileOption
    case 'nd2'
        [rowstring,colstring,sitestring]=wellnum2strRCS_3(row,col,site);
        rawdir = [imagepath,experimentpath,'Raw\'];
        shot=[rowstring,colstring];
        FileName=getFilenames(rawdir,['^Well',shot]);
        Image=BioformatsImage([rawdir,FileName{1}]);
        numChannel=numel(Image.channelNames);
        for idxChannel=1:numChannel
            switch ImageOption
                case 'multisite'
                    rawA{idxChannel}=double(Image.getXYplane(idxChannel, site, 1)); %channel, XYloc, frame
                case 'tile'
                    rawA{idxChannel}=getTile(Image, [1 idxChannel 1], [5 5], site);%getTile(bfObj, [Z, C, T], [numRows, numColumns], tileIndex)
            end
            %showImagesWithLinkedAxes({rawA{1},rawA{2},rawA{3}});%,rawImage{4}})
        end
    case 'tif'
        shot=saveshot;
        rawdir = [imagepath,experimentpath,'Raw\Seq_1\',shot,'\'];
        rawA{1}=double(imread([rawdir,nameA1,'.tif'])); %rawA1=(rawA1-bgcmos)./biasA1;
        rawA{2}=double(imread([rawdir,nameA2,'.tif'])); %rawA2=(rawA2-bgcmos)./biasA2;
        %rawImage_A{3}=double(imread([rawdir,nameA3,'.tif'])); %rawA3=(rawA3-bgcmos)./biasA3;
        %rawImage{4}=double(imread([rawdir,nameA4,'.tif'])); %rawA4=(rawA4-bgcmos)./biasA4;
        %rawImage{1}(rawImage{1} == MAXintensity) = NaN;
        %shot=[rowstring,colstring,'_',sitestring];
        if secondImaging
            rawdir = [imagepath,experimentpath,'Raw\Seq_2\',shot,'\'];
            rawImage_B{1}=double(imread([rawdir,nameB1,'.tif'])); %rawA1=(rawA1-bgcmos)./biasA1;
            rawImage_B{2}=double(imread([rawdir,nameB2,'.tif'])); %rawA2=(rawA2-bgcmos)./biasA2;
            rawImage_B{3}=double(imread([rawdir,nameB3,'.tif'])); %rawA3=(rawA3-bgcmos)./biasA3;
            if ismember(row,2:5)
                rawImage_B{4}=double(imread([rawdir,nameB4,'.tif'])); %rawA4=(rawA4-bgcmos)./biasA4;
            end
        end
        
end

if separatedirectories
    NEfile=[maskdir,'\',nucedgename,'_stain.tif'];
else
    NEfile=[maskdir,'\',shot,'_',nucedgename,'_stain.tif'];
end
%%% detect and remove smears %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%foregroundthresh=2000;
%areathresh=20000;
%bgprctile=20;
%[rawImage{1},stainflag]=removesmears_2(rawImage{1},foregroundthresh,areathresh,bgprctile);
%rawImage{1}(rawImage{1}<1)=1;

%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
segmethod='log';
switch segmethod
    case 'concavity'
        blurradius=5;
        blur=imfilter(rawA{1},fspecial('gaussian',blurradius),'symmetric');
        nuc_maskA=ThreshImage(blur);
    case 'log'
        blobthreshold=-0.03;
        nuc_maskA=blobdetector_5(log(rawA{1}),nucr,blobthreshold);
    case 'single'
        blurradius=1;
        nuc_maskA=threshmask(rawA{1},blurradius);
    case 'double'
        blurradius=3;
        nuc_maskA=threshmask(rawA{1},blurradius);
        nuc_maskA=markershed(nuc_maskA,round(nucr*2/3));
        nuc_maskA=secondthresh(rawA{1},blurradius,nuc_maskA,boulderarea*2);
end
nuc_maskA=logical(imfill(nuc_maskA,'holes')); %Troubleshoot segmentation
foreground=nuc_maskA;
nuc_maskA=bwareaopen(nuc_maskA,debrisarea); %Troubleshoot segmentation
nuc_maskA=bridgedeflections(nuc_maskA,nucr); %Troubleshoot segmentation
%[nuc_maskA,deflections,bridges]=bridgedeflections_illustrate(nuc_maskA,nucr); %show features
nuc_maskA=bwareaopen(nuc_maskA,debrisarea); %Troubleshoot segmentation
nuc_maskA=logical(nuc_maskA);
%%% segment second nuclear image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if secondImaging
%     nuc_maskB=threshmask(rawImage_B{1},3);
%     [reljitxAB,reljityAB]=registerimages(nuc_maskA,nuc_maskB);
%     jitmatx=[reljitxAB];jitmaty=[reljityAB];
%%% align images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [height,width]=size(rawA{1});
% cropcoors=getcropcoors([height width],jitmatx,jitmaty);
% foreground=foreground(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
% nuc_maskA=nuc_maskA(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
%
% rawA{1}=rawA{1}(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
% rawA{2}=rawA{2}(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
% rawA{3}=rawA{3}(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
% if ismember(row,2:4)
% rawA{4}=rawA{4}(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
% end
%     if ismember(row,2:4)
%         rawImage_B{1}=rawImage_B{1}(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
%         rawImage_B{2}=rawImage_B{2}(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
%         rawImage_B{3}=rawImage_B{3}(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
%         rawImage_B{4}=rawImage_B{4}(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
%     else
%         rawImage_B{1}=rawImage_B{1}(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
%         rawImage_B{2}=rawImage_B{2}(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
%         rawImage_B{3}=rawImage_B{3}(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
%     end
% end
%%% Note border\large\warped objects for later removal %%%%%%%%%%%%%%%%%%%%
antiborder_mask=imclearborder(nuc_maskA);
border_mask=nuc_maskA-antiborder_mask;
antilargewarped_mask=excludelargeandwarped_3(nuc_maskA,boulderarea,eccentricitythresh);
largewarped_mask=nuc_maskA-antilargewarped_mask;
badmask=border_mask | largewarped_mask;
goodmask=logical(nuc_maskA-badmask);
%%% Quality Control: Visually confirm acceptable segmentation %%%%%%%%%%%%%
%{
[height width]=size(goodmask);
extractmask=bwmorph(goodmask,'remove');
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(rawA{1}));
tempframe(:,:,2)=extractmask;
figure,imshow(tempframe);

% Alignment
nuc_maskB=nuc_maskB(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
% nuc_maskC=nuc_maskC(cropcoors(3,1):cropcoors(3,2),cropcoors(3,3):cropcoors(3,4));
% nuc_maskD=nuc_maskD(cropcoors(4,1):cropcoors(4,2),cropcoors(4,3):cropcoors(4,4));
[height width]=size(nuc_maskA);
extractmask_A=bwmorph(nuc_maskA,'remove');
extractmask_B=bwmorph(nuc_maskB,'remove');
% extractmask_C=bwmorph(nuc_maskC,'remove');
% extractmask_D=bwmorph(nuc_maskD,'remove');

%Align A-B
tempframe=zeros(height,width,3);
tempframe(:,:,1)=extractmask_B;
tempframe(:,:,2)=extractmask_A;
tempframe(:,:,3)=imadjust(mat2gray(rawImage_A{1}));
figure,imshow(tempframe);

%Align A-C
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(rawImage{1}));
tempframe(:,:,2)=extractmask_A;
tempframe(:,:,3)=extractmask_C;
figure,imshow(tempframe);

%Align A-D
tempframe=zeros(height,width,3);
tempframe(:,:,1)=imadjust(mat2gray(rawImage{1}));
tempframe(:,:,2)=extractmask_A;
tempframe(:,:,3)=extractmask_D;
figure,imshow(tempframe);
%}
%%% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nanmask=imdilate(foreground,strel('disk',nucr/2));
nanmaskcyto=imdilate(foreground,strel('disk',nucr*2));
blurA1=imfilter(rawA{1},fspecial('disk',3),'symmetric');
blurA2=imfilter(rawA{2},fspecial('disk',3),'symmetric');
blurA3=imfilter(rawA{3},fspecial('disk',3),'symmetric');

%%% general parameter for Back-ground Substraction %%%%%%%%%%%%%%%%%%%%%%%%
compression=1;
numblocks_1=4;
numblocks_2=1;
sampleprctile=50;

[realA1,bgA1]=bgsubmasked_global_2(blurA1,nanmask,numblocks_1,compression,sampleprctile);
[realA2,bgA2]=bgsubmasked_global_2(blurA2,nanmaskcyto,numblocks_1,compression,sampleprctile);
[realA3,bgA3]=bgsubmasked_global_2(blurA3,nanmaskcyto,numblocks_1,compression,sampleprctile);

%figure; imagesc(realA3,[0 2000])
%figure; imagesc(rawImage{3},[0 10000])
%keyboard

% [realB1,bgB1]=bgsubmasked_global_2(blurB1,nanmask,numblocks_2,compression,sampleprctile);
% [realB2,bgB2]=bgsubmasked_global_2(blurB2,nanmask,numblocks_2,compression,sampleprctile);
% [realB3,bgB3]=bgsubmasked_global_2(blurB3,nanmask,numblocks_2,compression,sampleprctile);
% if ismember(row,2:4)
%     [realB4,bgB4]=bgsubmasked_global_2(blurB4,nanmask,numblocks_2,compression,sampleprctile);
% end
%figure; showImagesWithLinkedAxes({rawImage{1},realB1,realB2,realB3,realB4})
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_568to647.mat'],'bleedthroughrate');
% real3bleedthrough=real3*bleedthroughrate(2)+bleedthroughrate(1);
% real4=real4-real3bleedthrough;
%%% label both good and bad nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nuc_label_good,numgood]=bwlabel(goodmask);
nuc_label_bad=bwlabel(badmask);
nuc_label_bad=nuc_label_bad+numgood;
nuc_label_bad(nuc_label_bad==numgood)=0;
nuc_label_all=nuc_label_good+nuc_label_bad;
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_info=regionprops(nuc_label_all,'Centroid','Area','PixelIdxList');
if puncta
    %thickenradius=2*nucr;
    %outer_label=labelthicken(nuc_label_all,thickenradius);    
    %imtool(realA3)
    punctathresh_low=5000; punctathresh_high=10000;
    [puncta_info_low_A2,punctacount_low_A2]=getpuncta(nuc_label_all,realA2,punctathresh_low);
    [puncta_info_high_A2,punctacount_high_A2]=getpuncta(nuc_label_all,realA2,punctathresh_high);    
end
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nanvec=ones(numgood,1)*NaN;
xcoor=nanvec; ycoor=nanvec; nuc_area=nanvec;
sigA1=nanvec; sigA2=nanvec; sigA3=nanvec; 

if puncta
    zerovec=zeros(numgood,1); %needed for possible zero counts for puncta
    sigA2arealow=zerovec; sigA2areahigh=zerovec;
    sigA2sumlow=zerovec; sigA2sumhigh=zerovec;
    sigA2countlow=zerovec; sigA2counthigh=zerovec;
    
    sigA3arealow=zerovec; sigA3areahigh=zerovec;
    sigA3sumlow=zerovec; sigA3sumhigh=zerovec;
    sigA3countlow=zerovec; sigA3counthigh=zerovec;
end

for cc=1:numgood
    xcoor(cc)=nuc_info(cc).Centroid(1);
    ycoor(cc)=nuc_info(cc).Centroid(2);
    nuc_area(cc)=nuc_info(cc).Area;
    sigA1(cc)=median(realA1(nuc_info(cc).PixelIdxList));
    sigA2(cc)=median(realA2(nuc_info(cc).PixelIdxList));
    sigA3(cc)=median(realA3(nuc_info(cc).PixelIdxList));      
end

ringcalc==1;
innerrad=1; outerrad=5; %10xB1|20xB2: 1\5
ringthresh=20;
MAXthresh=100000; %100000
%showImagesWithLinkedAxes({realA1,realA2,realA3})
% sigA2ring=sigringcalc([],nuc_label_all,innerrad,outerrad,numgood,realA2,ringthresh);
% sigA3ring=sigringcalc([],nuc_label_all,innerrad,outerrad,numgood,realA3,ringthresh);

%%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IFdata=[xcoor,ycoor,nuc_area,sigA1,sigA2,sigA3];

save([datadir,'IF_',saveshot,'.mat'],'IFdata');
toc(timetotal);
end

function [puncta_info,puncta_count]=getpuncta(nuclabel,raw,thresh)
%imtool(raw)
puncta_mask=raw>thresh;
puncta_label=nuclabel.*puncta_mask;
puncta_info=regionprops(puncta_label,'PixelIdxList','Area');
numcells=max(nuclabel(:));
puncta_count=zeros(numcells,1);
for i=1:numcells
    tempmask=puncta_label==i;
    [~,tempcount]=bwlabel(tempmask,4);
    puncta_count(i)=tempcount;
end
% showImagesWithLinkedAxes({puncta_mask,raw})
end

%{
%%% debugging: view images %%%%%%%%%%
extractmask=bwmorph(nuc_maskA,'remove');
rawImage{1}=rawImage{1}(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
tempframe=imadjust(mat2gray(raw{1}));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);

extractmask=bwmorph(nuc_maskA,'remove');
% rawB1=rawB1(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
tempframe=imadjust(mat2gray(rawB1));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);

nuc_info=struct2cell(regionprops(smear_mask,'Area')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
hist(nuc_area,100);

anti_mask=bwareaopen(nuc_mask,debrisarea);
temp_mask=nuc_mask-anti_mask;
extractmask=bwmorph(temp_mask,'remove');

anti_mask=bwareaopen(smear_mask,10000);
extractmask=bwmorph(anti_mask,'remove');
%}