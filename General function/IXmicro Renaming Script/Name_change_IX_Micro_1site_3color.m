%% Init
clear
close all
clc

time1=tic;  %Starts Timer
%% Load Directory with images
dir='F:\Hee Won\20130803-MCF10A-after';     %'/Users/scappell/Documents/Meyer_Lab/Data/2012-07-04/32';
alldir=getSubdirectories(dir); %Creates list of all the subfolders

%% Designate Channels
channels={'CFP','FRET','RFP'};  %Change this to match the correct channels used in the experiment  %e.g. {'Channel1','Channel2','Channel3'}  

for i=1:length(alldir);
    disp(i)                                                                               %Displays what timepoint you are on
    time2=tic;                                                                            %Starts timer to record time of each timepoint
    filenames=getFilenames([dir,'\TimePoint_',num2str(i)]);                           %Returns list of all file names inside the folder
    filenames=filenames(boolRegExp(filenames,'\.tif$') & ~boolRegExp(filenames,'thumb')); %removes accessory files and thumbnails from the file list
    
    %% extract wellname, site, and channel from original file name
    for j=1:length(filenames);
        wellname=char(getTokens(filenames(j),'^[^_]+_([A-Z][0-9][0-9])_'));  %Finds well name in format wellColumn (eg B04)
        row_column=wellNameToRowColumn(wellname);                            %Converts well name to format row_column (eg 2_4)
        %site=char(getTokens(filenames(j),'^[^_]+_[A-Z][0-9][0-9]_s([0-9])_')); %Finds site image was taken (eg 1)
        channel=char(getTokens(filenames(j),'^[^_]+_[A-Z][0-9][0-9]_w([1-9])'));  %Finds channel number (eg 1)      
        
        if i==1;
            if ~exist([dir,'\',wellname],'dir');
                mkdir([dir,'\',wellname]);
            end
%             if ~exist([dir,'/', wellname,'/site' num2str(site)],'dir')
%                 mkdir([dir,'/', wellname,'/site_' num2str(site)]);
%             end
        end
        
        %% move file
         newFileName=[row_column,'_1_',channels{str2double(channel)},'_',num2str(i+105),'.tif'];  %designates new name (eg 2_4_1_EYFP_1) (i+120)
         oldFileName=['TimePoint_',num2str(i),'\',filenames{j}];        %Designates old file name to be used in old file path
%          movefile([dir,'\',oldFileName],[dir,'\',newFileName],'f'); %Moves file from old folder to new folder and changes the name
         movefile([dir,'\',oldFileName],[dir,'\',wellname,'\',newFileName],'f'); %Moves file from old folder to new folder and changes the name. Use this line if you only have 1 site.

    end
    toc(time2) %Displays time it takes to rename all the files in each Time point
end
toc(time1)  %Displays total time to rename all files
