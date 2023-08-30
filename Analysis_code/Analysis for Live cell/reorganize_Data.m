function [ret] = reorganize_Data(outputPath,resultdir,option)
IFoption=0;
%addpath([pwd,'\functions'])
if ~exist(resultdir,'dir')
    mkdir(resultdir);
end
siteV=1:4;
if option==1
    conditionCell = {
        {1,1,siteV}; %
        {1,2,siteV}; %
        {1,3,siteV}; %
        {1,4,siteV}; %
        {1,5,siteV}; %
        {1,6,siteV}; %
        {1,7,siteV}; %
        {1,8,siteV}; %
        {1,9,siteV}; %
        {1,10,siteV}; %
        {1,11,siteV}; %
        {1,12,siteV}; %
        
        {3,1,siteV}; %
        {3,2,siteV}; %
        {3,3,siteV}; %
        {3,4,siteV}; %
        {3,5,siteV}; %
        {3,6,siteV}; %
        {3,7,siteV}; %
        {3,8,siteV}; %
        {3,9,siteV}; %
        {3,10,siteV}; %
        {3,11,siteV}; %
        {3,12,siteV}; %
        
        {2,1,siteV}; %
        {2,2,siteV}; %
        {2,3,siteV}; %
        {2,4,siteV}; %
        {2,5,siteV}; %
        {2,6,siteV}; %
        {2,7,siteV}; %
        {2,8,siteV}; %
        {2,9,siteV}; %
        {2,10,siteV}; %
        {2,11,siteV}; %
        {2,12,siteV}; %
        
        {4,1,siteV}; %
        {4,2,siteV}; %
        {4,3,siteV}; %
        {4,4,siteV}; %
        {4,5,siteV}; %
        {4,6,siteV}; %
        {4,7,siteV}; %
        {4,8,siteV}; %
        {4,9,siteV}; %
        {4,10,siteV}; %
        {4,11,siteV}; %
        {4,12,siteV}; %
        };
elseif option==2
    conditionCell = {
        {5,1,siteV}; %
        {5,2,siteV}; %
        {5,3,siteV}; %
        {5,4,siteV}; %
        {5,5,siteV}; %
        {5,6,siteV}; %
        {5,7,siteV}; %
        {5,8,siteV}; %
        {5,9,siteV}; %
        {5,10,siteV}; %
        {5,11,siteV}; %
        {5,12,siteV}; %
        
        {7,1,siteV}; %
        {7,2,siteV}; %
        {7,3,siteV}; %
        {7,4,siteV}; %
        {7,5,siteV}; %
        {7,6,siteV}; %
        {7,7,siteV}; %
        {7,8,siteV}; %
        {7,9,siteV}; %
        {7,10,siteV}; %
        {7,11,siteV}; %
        {7,12,siteV}; %
        
        {6,1,siteV}; %
        {6,2,siteV}; %
        {6,3,siteV}; %
        {6,4,siteV}; %
        {6,5,siteV}; %
        {6,6,siteV}; %
        {6,7,siteV}; %
        {6,8,siteV}; %
        {6,9,siteV}; %
        {6,10,siteV}; %
        {6,11,siteV}; %
        {6,12,siteV}; %
        
        {8,1,siteV}; %
        {8,2,siteV}; %
        {8,3,siteV}; %
        {8,4,siteV}; %
        {8,5,siteV}; %
        {8,6,siteV}; %
        {8,7,siteV}; %
        {8,8,siteV}; %
        {8,9,siteV}; %
        {8,10,siteV}; %
        {8,11,siteV}; %
        {8,12,siteV}; %
        };
    
end

operation = 'data(:,:,1) = tracedata(:,:,7)./tracedata(:,:,6); data(:,:,2) = tracedata(:,:,9)./tracedata(:,:,8); data(:,:,3) = tracedata(:,:,11)./tracedata(:,:,10);'

conditions = {'wt'};
%outputPath = uigetdir;
ret = reconstruct_tracedata(conditionCell,outputPath,operation,IFoption);
save([resultdir filesep 'CombinedData_' num2str(option) '.mat'],'ret','conditions','-v7.3');
end