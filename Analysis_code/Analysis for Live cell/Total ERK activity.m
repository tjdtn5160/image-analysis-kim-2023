%% Load combined data
%clear; close all;
root='E:\2. Data\4. Drug resistance\2018-03-26 (Melanoma; ERK-CDK2-CDK4-H2B)\';
option=3;
resultdir=[root 'Results\'];
load([resultdir 'CombinedData_' num2str(option)]);
%% Total ERK activity (first 100 frames)
condi=5; firstFram=6; lastFrame=51;
IF=0; inhibitionFrame=0; nodivision=1;
if IF
    IFdata = ret{condition}.store.IFdata;
else
    IFdata=nan;
end
condition=condi;
data = ret{condition}.store.storeData;
divFrame = ret{condition}.store.allParentDivFrame;
D = sort_or_align_by_mitosis(data,divFrame,IF,IFdata,nodivision)
D.clean_leastnum_frames_5(firstFram,lastFrame,15);
D.clean_remove_by_lower_upper_limit(1, [0, 15])
D.clean_remove_by_lower_upper_limit(2, [0, 3])
D.clean_remove_by_lower_upper_limit(3, [0, 3])
D.clean_remove_dataaligned_indicate_timeandrange(2,4980:4985,[0 1])
D.clean_remove_dataaligned_indicate_timeandrange(2,5010:5025,[1 4])
D.clean_remove_dataaligned_indicate_timeandrange(2,4990:4993,[0.5 1])
D.smoothDataSortedByMitosis
ERK_corrected=ERK_activity_correction_allcells(D.dataSorted(:,:,1),D.dataSortedSmooth(:,:,2));
%%
% D.clean_manually_remove_data(1,0);
%%
label={'WM858_Control';'WM983_Control';'WM858_Vem';'WM983_Vem';'WM858_resis_Vem';'WM983_resis_Vem'};
num=6;

ERK{num}=ERK_corrected(:,firstFram:lastFrame);
%%
figure; hold on;
for i=1:length(ERK)
    plot(nanmean(ERK{i}));
end
legend(label)
%%
for i=1:length(ERK)
    dataERK{i}=ERK{i};
end

figure; hold on;
for i=1:length(dataERK)
    plot(nanmean(dataERK{i}));
end
legend(label)
%% save Data
save([resultdir 'total_ERKdata_20180326-2'],'dataERK','label');
