row=3;col=6;site=15;
ROIxmin=14; ROIxmax=ROIxmin+400; ROIwidth=ROIxmax-ROIxmin+1;
ROIymin=300; ROIymax=ROIymin+400; ROIheight=ROIymax-ROIymin+1;

projectpath='\Volumes\HW_Data\2. Data\2. Development of CDK4 sensor\';
imagepath='\Volumes\CDK4 study\Data for CDK4 sensor\'
%shadingpath='D:\Images\ShadingImages\20140425 DCYTC 10x\';
cmosoffsetfile='\Volumes\HW_Data\3. Processing data\IXmicro\CameraNoise\cmosoffset_bin2.mat';
experimentpath='2018-09-17 (MCF10A; Geminin-CDK2-CDK4-H2B; multi p-Rb 12hr)\';
shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
rawdir = [imagepath,experimentpath,'Raw\',shot,'\'];
maskdir = [imagepath,experimentpath,'Mask\',shot];
%%% setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nameA1='H2B'; %nuc
nameA2='APC';
nameA3='CDK2';
nameA4='CDK4';

nameB1='Hoechst-1';
nameB2='A488-1';

nameC1='Hoechst-2';
nameC2='A488-2';
nameC3='A647-2';

nameD1='Hoechst-3';
nameD2='A488-3';
%%% segmentation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nucr=12;
debrisarea=100; %150
boulderarea=2000; %1000
thickenradius=2*nucr;
timetotal=tic;
%%% shading correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(cmosoffsetfile,'cmosoffset'); %load CMOS offset bgcmos=shadingcorrection;
[height,width]=size(cmosoffset);
%%% load images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawA1=double(imread([rawdir,nameA1,'_stain.tif'])); %rawA2=(rawA2-cmosoffset).\bias2;
rawA2=double(imread([rawdir,nameA2,'_stain.tif'])); %rawA2=(rawA2-cmosoffset).\bias2;
rawA3=double(imread([rawdir,nameA3,'_stain.tif'])); %rawA3=(rawA3-cmosoffset).\bias3;
rawA4=double(imread([rawdir,nameA4,'_stain.tif'])); %rawA2=(rawA2-cmosoffset).\bias2;

rawB1=double(imread([rawdir,nameB1,'_stain.tif'])); %raw1=(raw1-bgcmos).\bias1; 
rawB2=double(imread([rawdir,nameB2,'_stain.tif'])); %raw2=(raw2-bgcmos).\bias2; %S807\811

rawC1=double(imread([rawdir,nameC1,'_stain.tif'])); %raw1=(raw1-bgcmos).\bias1; 
rawC2=double(imread([rawdir,nameC2,'_stain.tif'])); %raw1=(raw1-bgcmos).\bias1; 
rawC3=double(imread([rawdir,nameC3,'_stain.tif'])); %raw1=(raw1-bgcmos).\bias1; 

rawD1=double(imread([rawdir,nameD1,'_stain.tif'])); %raw1=(raw1-bgcmos).\bias1; 
rawD2=double(imread([rawdir,nameD2,'_stain.tif'])); %raw1=(raw1-bgcmos).\bias1; 
% %%% remove smears %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
foregroundthresh=1000;
areathresh=20000;
bgprctile=20;
% [rawA1,stainflag]=removesmears_2(rawA1,foregroundthresh,areathresh,bgprctile);
rawB1(rawB1<1)=1;
%%% segment nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
segmethod='single';
switch segmethod
    case 'concavity'
        blurradius=5;
        blur=imfilter(rawB1,fspecial('gaussian',blurradius),'symmetric');
        nuc_maskB=ThreshImage(blur);
        nuc_maskB=markershed(nuc_maskB,round(nucr*2\3));
    case 'log'
        blobthreshold=-0.03;
        nuc_maskB=blobdetector_3(log(rawB1),nucr,blobthreshold,debrisarea);
    case 'single'
        blurradius=3;
        nuc_maskB=threshmask(rawB1,blurradius);
        nuc_maskB=markershed(nuc_maskB,round(nucr*2\3));
    case 'double'
        blurradius=3;
        nuc_maskB=threshmask(rawB1,blurradius);
        nuc_maskB=markershed(nuc_maskB,round(nucr*2\3));
        nuc_maskB=secondthresh(rawB1,blurradius,nuc_maskB,boulderarea*2);
end
nuc_maskB=imfill(nuc_maskB,'holes');
nuc_maskB=bwareaopen(nuc_maskB,debrisarea);
nuc_maskB=segmentdeflections_bwboundaries(nuc_maskB,nucr,debrisarea);
nuc_maskB=logical(nuc_maskB);
%%% segment second nuclear image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_maskA=threshmask(rawA1,3);
nuc_maskA=markershed(nuc_maskA,nucr*2\3);
%%% align images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[reljitxAB,reljityAB]=registerimages(nuc_maskA,nuc_maskB);

nuc_maskC=threshmask(rawC1,3);
[reljitxAC,reljityAC]=registerimages(nuc_maskA,nuc_maskC);

nuc_maskD=threshmask(rawD1,3);
[reljitxAD,reljityAD]=registerimages(nuc_maskA,nuc_maskD);
jitmatx=[reljitxAB;reljitxAC;reljitxAD]; jitmaty=[reljityAB;reljityAC;reljityAD];

cropcoors=getcropcoors([height width],jitmatx,jitmaty);
rawA1=rawA1(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
rawA2=rawA2(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
rawA3=rawA3(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));
rawA4=rawA4(cropcoors(1,1):cropcoors(1,2),cropcoors(1,3):cropcoors(1,4));

rawB1=rawB1(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
rawB2=rawB2(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));
nuc_maskB=nuc_maskB(cropcoors(2,1):cropcoors(2,2),cropcoors(2,3):cropcoors(2,4));

rawC1=rawC1(cropcoors(3,1):cropcoors(3,2),cropcoors(3,3):cropcoors(3,4));
rawC2=rawC2(cropcoors(3,1):cropcoors(3,2),cropcoors(3,3):cropcoors(3,4));
rawC3=rawC3(cropcoors(3,1):cropcoors(3,2),cropcoors(3,3):cropcoors(3,4));

rawD1=rawD1(cropcoors(4,1):cropcoors(4,2),cropcoors(4,3):cropcoors(4,4));
rawD2=rawD2(cropcoors(4,1):cropcoors(4,2),cropcoors(4,3):cropcoors(4,4));

% showImagesWithLinkedAxes({rawB1,rawB2,rawC2,rawC3,rawD2})
%%% calculate background: semi-local %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nanmaskcyto=imdilate(nuc_maskB,strel('disk',nucr*2));
compression=1;

blurA4=medfilt2(rawA4,[3 3]);
% showImagesWithLinkedAxes({rawA4,blurA4})
realA2=bgsubmasked_global_2(rawA2,nanmaskcyto,1,compression,50);
realA3=bgsubmasked_global_2(rawA3,nanmaskcyto,1,compression,50);
realA4=bgsubmasked_global_2(rawA4,nanmaskcyto,1,compression,50);
%realA4=imresize(realA4,0.5,'nearest');

realB1=bgsubmasked_global_2(rawB1,nanmaskcyto,1,compression,50);
% realB2top=imtophat(rawB2,strel('disk',2,0));
realB2=bgsubmasked_global_2(rawB2,nanmaskcyto,1,compression,50);

realC2=bgsubmasked_global_2(rawC2,nanmaskcyto,1,compression,50);
realC3=bgsubmasked_global_2(rawC3,nanmaskcyto,1,compression,50);

realD2=bgsubmasked_global_2(rawD2,nanmaskcyto,1,compression,50);

%showImagesWithLinkedAxes({realB1,realB2,realC2,realC3,realD2,realA4})

%%% Note border\large\warped objects for later removal %%%%%%%%%%%%%%%%%%%%
antiborder_mask=imclearborder(nuc_maskB);
border_mask=nuc_maskB-antiborder_mask;
antilargewarped_mask=excludelargeandwarped_3(nuc_maskB,boulderarea,0.95);
largewarped_mask=nuc_maskB-antilargewarped_mask;
badmask=border_mask | largewarped_mask;
goodmask=logical(nuc_maskB-badmask);
%%% correct IF for bleedthrough %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load([datadir,'bleedthroughrate_568to647.mat'],'bleedthroughrate');
% real4bleedthrough=real3*bleedthroughrate(2)+bleedthroughrate(1);
% real4=real4-real4bleedthrough;
%%% subtract background and calculate masses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% label both good and bad nuclei %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nuc_label_good,numgood]=bwlabel(goodmask);
nuc_label_bad=bwlabel(badmask);
nuc_label_bad=nuc_label_bad+numgood;
nuc_label_bad(nuc_label_bad==numgood)=0;
nuc_label_all=nuc_label_good+nuc_label_bad;
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nuc_info=regionprops(nuc_label_all,'Centroid','Area','PixelIdxList');
%%% extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nanvec=ones(numgood,1)*NaN; zerovec=zeros(numgood,1);
% Hoechst=nanvec; pRb=nanvec;
% sigB2area=zerovec; sigB2mean=zerovec; sigB2count=zerovec;
% for cc=1:numgood
%     Hoechst(cc)=mean(realB1(nuc_info(cc).PixelIdxList));
%     pRb(cc)=median(realB3(nuc_info(cc).PixelIdxList));
% end

% showImagesWithLinkedAxes({realB1,realB3,puncta_mask,realA4,realA3,realA2})

%%%%%% pRb
pRb(pRb<1)=1;
pRb=log2(pRb);
%histogram(pRb)
pRbthresh=10.5;
%%% Generate Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imscale=2;
pixlength=31; %20x bin2; 0.65um
pixwidth=3;
scalebarmargin=10;
pixdims=[pixwidth pixlength];
startpos=[ROIheight-pixwidth-scalebarmargin ROIwidth-pixlength-scalebarmargin];

%%%%%% Hoechst
crop_Hoechst=realB1(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
RGB(:,:,1)=imadjust(mat2gray(crop_Hoechst));
RGB(:,:,2)=imadjust(mat2gray(crop_Hoechst));
RGB(:,:,3)=imadjust(mat2gray(crop_Hoechst));
Hoechst_image=RGB;

crop_pRb=realB2(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
RGB(:,:,1)=imadjust(mat2gray(crop_pRb),[0.03 0.4]);
RGB(:,:,2)=imadjust(mat2gray(crop_pRb),[0.03 0.4]);
RGB(:,:,3)=imadjust(mat2gray(crop_pRb),[0.03 0.4]);
pRb807_image=RGB;

crop_pRb=realC2(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
RGB(:,:,1)=imadjust(mat2gray(crop_pRb),[0.1 0.5]);
RGB(:,:,2)=imadjust(mat2gray(crop_pRb),[0.1 0.5]);
RGB(:,:,3)=imadjust(mat2gray(crop_pRb),[0.1 0.5]);
pRbS373_image=RGB;

crop_pRb=realC3(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
RGB(:,:,1)=imadjust(mat2gray(crop_pRb),[0.08 0.35]);
RGB(:,:,2)=imadjust(mat2gray(crop_pRb),[0.08 0.35]);
RGB(:,:,3)=imadjust(mat2gray(crop_pRb),[0.08 0.35]);
pRbT780_image=RGB;

crop_pRb=realD2(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
RGB(:,:,1)=imadjust(mat2gray(crop_pRb),[0.1 0.3]);
RGB(:,:,2)=imadjust(mat2gray(crop_pRb),[0.1 0.3]);
RGB(:,:,3)=imadjust(mat2gray(crop_pRb),[0.1 0.3]);
pRbT608_image=RGB;

crop_CDK4=realA4(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
RGB(:,:,1)=imadjust(mat2gray(crop_CDK4));
RGB(:,:,2)=imadjust(mat2gray(crop_CDK4));
RGB(:,:,3)=imadjust(mat2gray(crop_CDK4));
CDK4_image=RGB;

crop_CDK2=realA3(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
RGB(:,:,1)=imadjust(mat2gray(crop_CDK2),[0.05 0.5]);
RGB(:,:,2)=imadjust(mat2gray(crop_CDK2),[0.05 0.5]);
RGB(:,:,3)=imadjust(mat2gray(crop_CDK2),[0.05 0.5]);
CDK2_image=RGB;

crop_APC=realA2(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
RGB(:,:,1)=imadjust(mat2gray(crop_APC),[0.15 0.5]);
RGB(:,:,2)=imadjust(mat2gray(crop_APC),[0.15 0.5]);
RGB(:,:,3)=imadjust(mat2gray(crop_APC),[0.15 0.5]);
RGB=addscalebar(RGB,startpos,pixdims);
APC_image=RGB;

showImagesWithLinkedAxes({Hoechst_image,pRb807_image,pRbS373_image,pRbT780_image,pRbT608_image,CDK4_image,CDK2_image,APC_image})

% %%% remove smears %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%seg_nuc_label=bwlabel(seg_nuc);
% pRbposID=[1,2,4,5,6,8,10,16,18];
% pRbnegID=[3,7,9,11,12,13,14,15,17];
pRbposID=find(pRb>pRbthresh);
pRbnegID=find(pRb<=pRbthresh);
%seg_pRbPos=ismember(seg_nuc_label,pRbposID);
%seg_pRbNeg=ismember(seg_nuc_label,pRbnegID);
seg_pRbPos=ismember(nuc_label_all,pRbposID);
seg_pRbNeg=ismember(nuc_label_all,pRbnegID);
crop_pRbPos=seg_pRbPos(ROIymin:ROIymax,ROIxmin:ROIxmax);
crop_pRbNeg=seg_pRbNeg(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
greenchannel=crop_pRbPos | crop_pRbNeg;
RGB(:,:,1)=crop_pRbPos;
RGB(:,:,2)=greenchannel;
RGB(:,:,3)=crop_pRbNeg;
RGB=addscalebar(RGB,startpos,pixdims);
pRbClass_image=RGB;

%%%%% FISH: pRb-pos:1 pRb-neg:2
threshFISH=500;
FISH_mask=realB2top>threshFISH;
pRb_label=seg_pRbPos+2*seg_pRbNeg;
outer_label=labelthicken(pRb_label,thickenradius);
outer_ring=outer_label-imerode(outer_label,strel('disk',1,0));
FISH_label=FISH_mask.*outer_label;
seg_FISHPos=outer_ring==1 | FISH_label==1;
seg_FISHNeg=outer_ring==2 | FISH_label==2;
crop_FISHPos=seg_FISHPos(ROIymin:ROIymax,ROIxmin:ROIxmax);
crop_FISHNeg=seg_FISHNeg(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
greenchannel=crop_FISHPos | crop_FISHNeg;
RGB(:,:,1)=crop_FISHPos;
RGB(:,:,2)=greenchannel;
RGB(:,:,3)=crop_FISHNeg;
E2F_image=RGB;


FISH_mask_Im=FISH_mask(ROIymin:ROIymax,ROIxmin:ROIxmax);
showImagesWithLinkedAxes({Hoechst_image,pRb_image,CDK4_image,CDK2_image,APC_image,FISH_mask_Im,E2F_image,pRbClass_image})

%% Processed FISH image
outer_label=labelthicken(nuc_label_all,thickenradius);
badIDs=unique(nuc_label_all.*nuc_label_bad);
badIDs(1)=[];
numobs=max(outer_label(:));
badIDs(badIDs>numobs)=[];
outer_info=regionprops(outer_label,'PixelIdxList');
puncta_label=outer_label.*FISH_mask;
puncta_label(1,1)=numobs+1; %ensure that pucta_info has at least nummem slots
puncta_info=regionprops(puncta_label,'PixelIdxList','Area');
punctacount=ones(numobs,1)*NaN;
image_number=zeros(size(outer_label));
image_area=zeros(size(outer_label));
for i=1:numobs
    if ismember(i,badIDs)
        image_number(outer_label==i)=-100;
        image_area(outer_label==i)=-100;
    else
        tempmask=puncta_label==i;
        [~,tempcount]=bwlabel(tempmask,4);
        punctacount(i)=tempcount;
        image_number(outer_label==i)=tempcount;
        image_area(outer_label==i)=puncta_info(i).Area;
    end
end
image_area=image_area(ROIymin:ROIymax,ROIxmin:ROIxmax);

figure; imagesc(image_area,[-1 50]); colormap(hot_gray(51));
axis image;  set(gca,'XTick',[]);set(gca,'YTick',[]); colorbar; title('FISH image','FontSize',12);

% image_G=mat2gray(realB2top,[0 5000]);
crop_FISH=realB2top(ROIymin:ROIymax,ROIxmin:ROIxmax);
RGB=zeros(ROIheight,ROIwidth,3);
RGB(:,:,1)=imadjust(mat2gray(crop_FISH));
RGB(:,:,2)=imadjust(mat2gray(crop_FISH));
RGB(:,:,3)=imadjust(mat2gray(crop_FISH));
FISH_image=RGB;
figure; imagesc(FISH_image);
axis image;  set(gca,'XTick',[]);set(gca,'YTick',[]); title('Processed FISH_area','FontSize',12);

% crop_H2B=rawB1(ROIymin:ROIymax,ROIxmin:ROIxmax);
% RGB=zeros(ROIheight,ROIwidth,3);
% RGB(:,:,1)=imadjust(mat2gray(crop_H2B));
% RGB(:,:,2)=imadjust(mat2gray(crop_H2B));
% RGB(:,:,3)=imadjust(mat2gray(crop_H2B));
% H2B_image=RGB;





RGB=addscalebar(RGB,startpos,pixdims);
imwrite(imresize(RGB,imscale),[presentationdir,'Seg_',shot,'_FISH.tif']);




imwrite(imresize(RGB,imscale),[presentationdir,'Seg_',shot,'_pRb.tif']);


%%% Generate Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;imshow(RGB);
%%%%%% Hoechst
RGB=imresize(zeros(ROIheight,ROIwidth,3),2);
RGB(:,:,3)=imresize(imadjust(mat2gray(real1(ROIymin:ROIymax,ROIxmin:ROIxmax))),2);
imwrite(double(RGB),[presentationdir,'Raw_6_5_5_Hoechst.tif']);
%%%%%% DHB
RGB=imresize(zeros(ROIheight,ROIwidth,3),2);
RGB(:,:,2)=imresize(imadjust(mat2gray(real2(ROIymin:ROIymax,ROIxmin:ROIxmax))),2);
imwrite(double(RGB),[presentationdir,'Raw_6_5_5_DHB.tif']);
%%%%%% pRb
RGB=imresize(zeros(ROIheight,ROIwidth,3),2);
RGB(:,:,1)=imresize(imadjust(mat2gray(real4(ROIymin:ROIymax,ROIxmin:ROIxmax))),2);
imwrite(double(RGB),[presentationdir,'Raw_6_5_5_pRb.tif']);
%%%%%% FISH
RGB=imresize(zeros(ROIheight,ROIwidth,3),2);
RGB(:,:,1)=imresize(imadjust(mat2gray(real3(ROIymin:ROIymax,ROIxmin:ROIxmax))),2);
RGB(:,:,2)=imresize(imadjust(mat2gray(real3(ROIymin:ROIymax,ROIxmin:ROIxmax))),2);
imwrite(double(RGB),[presentationdir,'Raw_6_5_5_FISH.tif']);







