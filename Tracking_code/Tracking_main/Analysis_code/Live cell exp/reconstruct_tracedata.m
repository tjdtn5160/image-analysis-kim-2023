function [ret] = reconstruct_tracedata(conditionCell,outputpath,operation,leastFrameNum)
if nargin==3
    leastFrameNum=0;%default
end
for condi = 1:length(conditionCell)
    store.data = [];
    store.storeData = [];
    store.firstParentDivFrame = [];
    store.parentStat = [];
    store.allParentDivFrame = [];
    for coli = conditionCell{condi}{2}
        for rowi = conditionCell{condi}{1}
            for sitei = conditionCell{condi}{3}
                if exist([outputpath,filesep,'tracedata_',num2str(rowi),'_',num2str(coli),'_',num2str(sitei),'.mat'])
                    load([outputpath,filesep,'tracedata_',num2str(rowi),'_',num2str(coli),'_',num2str(sitei),'.mat']);
                    
                    %                 data(:,:,1) = tracedata(:,:,9)./tracedata(:,:,8);%CDK2
                    %                 data(:,:,2) = tracedata(:,:,7)./tracedata(:,:,6);%FRET
                    eval(operation);
                    
                    allnanidx = find(sum(~isnan(data),2)==0);
                    genealogy(ismember(genealogy, allnanidx)) = NaN;
                    
                    %data = interpolate_value_ifInf(data);
                    parentStat = parentstat(data,genealogy);
                    %leastFrameNum = 150; % cells that have more than 'leastFrameNum' frames
                    [storeData,storeFirstParentDivFrame,allParentDivFrame] = combine_parentanddaughter(data,genealogy,parentStat,leastFrameNum);
                    
                    [data,storeData] = find_outlierframe(data,storeData);
                    
                    store.data = [store.data;data];
                    store.storeData = [store.storeData;storeData];
                    store.parentStat = [store.parentStat;parentStat + size(store.parentStat,1)];
                    
                    store.firstParentDivFrame = [store.firstParentDivFrame;storeFirstParentDivFrame];
                    store.allParentDivFrame = [store.allParentDivFrame;allParentDivFrame];
                    
                    clearvars -except coli rowi cond store ret rowCond colCond siteCond condi conditionCell outputpath operation leastFrameNum
                else
                    clearvars -except coli rowi cond store ret rowCond colCond siteCond condi conditionCell outputpath operation leastFrameNum
                end
            end
        end
        
    end
    ret{condi}.store = store;
end
end

function [parentStat] = parentstat(data,genealogy)
% parentStat(:,1) shows the parentId. If parentStat(122,1) is 45, that
% means 122th cell is a daughter of 45th cell.
%parentStat(:,2) indicates when the parent arises, and parentStat(:,3)
%indicates when the parent disappers. parentStat(:,4) shows how long (how
%many frames) parent exists.

%get parentStat-----------------------------------------------
parentStat = ones(length(genealogy),4).*NaN;parentStat(:,1) = genealogy;
for i = size(genealogy,1):-1:1
    if ~isnan(genealogy(i))
        startParent = find(~isnan(data(genealogy(i),:,1)),1,'first');
        lastParent = find(~isnan(data(genealogy(i),:,1)),1,'last');
        parentStat(i,2) = startParent;
        parentStat(i,3) = lastParent;
        parentStat(i,4) = lastParent - startParent + 1;
    end
end
end

function [data] = interpolate_value_ifInf(data)
%interpolate value if Inf--------------------------------
data(data==Inf) = 0;
[infRow,infCol] = find(sum(data==0,3)>0);
for i = 1:length(infRow)
    if infCol(i) == size(data,2)
        data(infRow(i),infCol(i),:) = data(infRow(i),infCol(i)-1,:);
    elseif infCol(i)>1
        data(infRow(i),infCol(i),:) = mean([data(infRow(i),infCol(i)-1,:),data(infRow(i),infCol(i)+1,:)]);
    elseif infCol(i)==1
        data(infRow(i),infCol(i),:) = data(infRow(i),infCol(i)+1,:);
    end
end
end

function [data,storeData] = find_outlierframe(data,storeData)

%find outliers in patterns and interpolate it. first do in cdk2 and then
%fret. actually this is artificial smoothing for some drastic changes. not
%really good for some type of stimulation?
for chan = 1:size(storeData,3)
    x = storeData(:,:,chan) - repmat(nanmin(storeData(:,:,chan),[],2),1,size(storeData(:,:,chan),2));
    x = x./repmat(max(x,[],2),1,size(x,2));
    artifactFrame = mean(diff(x,[],2));
    artifactFrame = [mean(artifactFrame),artifactFrame];%just to adjust the number
    outlierBoundaryLow = mean(artifactFrame) - 2*std(artifactFrame);%2sigma
    outlierBoundaryHigh = mean(artifactFrame) + 2*std(artifactFrame);%2sigma
    artifactFrame = outlierBoundaryLow>artifactFrame|outlierBoundaryHigh<artifactFrame;
    artifactFrame = find(artifactFrame==1);
    for i = 1:length(artifactFrame)
        if artifactFrame(i)==1
            storeData(:,artifactFrame(i),:) = [storeData(:,artifactFrame(i)+1,:)];
            data(:,artifactFrame(i),:) = [data(:,artifactFrame(i)+1,:)];
        elseif artifactFrame(i)==size(storeData,2)
            storeData(:,artifactFrame(i),:) = [storeData(:,artifactFrame(i)-1,:)];
            data(:,artifactFrame(i),:) = [data(:,artifactFrame(i)-1,:)];
        else
            storeData(:,artifactFrame(i),:) = mean([storeData(:,artifactFrame(i)-1,:),storeData(:,artifactFrame(i)+1,:)],2);
            data(:,artifactFrame(i),:) = mean([data(:,artifactFrame(i)-1,:),data(:,artifactFrame(i)+1,:)],2);
        end
    end
end
end

function [storeData, storeFirstParentDivFrame,storeAllDivFrame] = combine_parentanddaughter(data,genealogy,parentStat,leastFrameNum)

%combine parent data to daughter's and extract that has profile throughout
%all the timecourse
%parentStat(:,1) = genealogy; parentStat(i,2) = startParent; parentStat(i,3) = lastParent; parentStat(i,4) = lastParent - startParent + 1;
storeData = ones(size(data)).*NaN;
storeFirstParentDivFrame = ones(size(data,1),1).*NaN;
storeAllDivFrame = ones(size(data,1),20).*NaN;

% remove duplicated mother traces
u=unique(genealogy(~isnan(genealogy)));
n=histc(genealogy,u);
d = u(n > 1);
for i=1:length(d)
    out=find(ismember(genealogy,d(i)));
    genealogy(out(2:end))=nan;
end

% combine mother and daughter traces
for i = size(genealogy,1):-1:1
    
    temp = data(i,:,:);num = i;firstParentDivFrame = NaN;
    allDivFrame = [];% added
    while ~isnan(genealogy(num))
        num_temp=num;
        temp(:,parentStat(num,2):parentStat(num,3),:) = data(parentStat(num,1),parentStat(num,2):parentStat(num,3),:);
        firstParentDivFrame = parentStat(num,3);
        allDivFrame = [allDivFrame,parentStat(num,3)];% added
        num = genealogy(num);
        if num_temp==num
            break
        end
    end
    storeAllDivFrame(i,1:size(allDivFrame,2)) = fliplr(allDivFrame);% added
    storeFirstParentDivFrame(i) = firstParentDivFrame;
    storeData(i,:,:) = temp;
    
end

storeFirstParentDivFrame = storeFirstParentDivFrame(sum(~isnan(storeData(:,:,1)),2)>=leastFrameNum,:);
storeAllDivFrame = storeAllDivFrame(sum(~isnan(storeData(:,:,1)),2)>=leastFrameNum,:);% added

storeData = storeData(sum(~isnan(storeData(:,:,1)),2)>=leastFrameNum,:,:);
end