classdef sort_or_align_by_mitosis < handle
    
    properties
        rawData
        rawDivFrame
        IFdata
        rawSorted
        IFSorted
        dataSortedSmooth
        dataAlignedSmooth
        rawDivFrameSorted
        sortByDivision = 1;
        goodTrace
        risingPoint
        
        timeInterval = 12/60; % default interval in hour
        timeToPlotBeforeAfterMitosis = 20; % default: plot 6 hours before and after mitosis
        cellState % for example, Cdk2 Inc and Cdk2 Low
        cellLow
        cellDelay
        cellERKinc
        cellERKlow
        sortXTickTime = [1,2,3,4];
        colorRange
        titles
        xlabels = {'time relative to division (h)'};
        ylabels
        mRNAarea
        mRNAintensity
        mRNAcount
        Protein
        nodivision
        hyperpRb
    end
    
    properties (Access = private)
        ifor % needed to pass info to add_plot_info
    end
    
    properties (Dependent = true)
        divFrameSorted
        dataSorted
        IFdataSorted
        dataAligned
        IFdataAligned
    end
    properties (Dependent = true, Access=private)
        divFrameAfterAlignment
    end
    
    
    methods
        function obj = sort_or_align_by_mitosis(rawData,rawDivFrame,IF,IFdata,nodivision)
            obj.rawData = rawData;obj.rawDivFrame = rawDivFrame;obj.nodivision=nodivision; if IF; obj.IFdata = IFdata; end
            
            if nodivision
                
            else
                rawData = rawData(~isnan(rawDivFrame(:,obj.sortByDivision)),:,:);
                if IF
                    IFdata = IFdata(~isnan(rawDivFrame(:,obj.sortByDivision)),:);
                end
                rawDivFrame = rawDivFrame(~isnan(rawDivFrame(:,obj.sortByDivision)),:);
            end
            
            [~,order] = sort(rawDivFrame(:,obj.sortByDivision));
            rawDivFrameSorted = ones(size(rawDivFrame)).*NaN;
            rawSorted = ones(size(rawData)).*NaN;
            if IF
                IFSorted = ones(size(IFdata)).*NaN;
            end
            for i = 1:size(rawData,1)
                rawSorted(i,:,:) = rawData(order(i),:,:);
                if IF
                    IFSorted(i,:) = IFdata(order(i),:);
                end
                rawDivFrameSorted(i,:) = rawDivFrame(order(i),:);
            end
            obj.rawSorted = rawSorted;
            if IF
                obj.IFSorted = IFSorted;
            end
            obj.rawDivFrameSorted = rawDivFrameSorted;
            
            % set default color range
            for i = 1:size(rawData,3)
                temp = rawData(:,:,i);
                colorRange{i,:} = [prctile(temp(:),5),prctile(temp(:),95)];
            end
            obj.colorRange = colorRange;
            obj.goodTrace = ones(size(rawSorted,1),1);
        end
        %% get and set functions
        function value = get.dataSorted(obj)
            value = obj.rawSorted(logical(obj.goodTrace),:,:);
        end
        
        function value = get.divFrameSorted(obj)
            value = obj.rawDivFrameSorted(logical(obj.goodTrace),:);
        end
        
        function value = get.IFdataSorted(obj)
            if IF
                value = obj.IFSorted(logical(obj.goodTrace),:,:);
            end
        end
        
        function value = get.dataAligned(obj)
            divFrameSorted = obj.divFrameSorted;dataSorted = obj.dataSorted;
            dataAligned = ones(size(dataSorted,1),10000,size(dataSorted,3))*NaN;
            for j = 1:size(dataSorted,1)
                temp = dataSorted(j,:,:);
                if obj.nodivision
                    dataAligned(j,1:size(temp,2),:) = temp;
                else
                    dataAligned(j,5000-divFrameSorted(j,obj.sortByDivision):5000-divFrameSorted(j,obj.sortByDivision)+size(temp,2)-1,:) = temp;
                    %                    divFrameAfterAlignment(j,:)  = divFrameAfterAlignment(j,:)-tempdiv(j); % division frame relative to alignment position
                end
            end
            value = dataAligned;
        end
        
        function value = get.divFrameAfterAlignment(obj)
            divFrameSorted = obj.divFrameSorted;
            dataSorted = obj.dataSorted;
            divFrameAfterAlignment = divFrameSorted;
            for j = 1:size(dataSorted,1)
                temp = dataSorted(j,:,:);
                divFrameAfterAlignment(j,:)  = divFrameAfterAlignment(j,:)-divFrameSorted(j,:); % division frame relative to alignment position
            end
            value = divFrameAfterAlignment;
        end
        
        function obj = set.sortByDivision(obj,varargin)
            obj.sortByDivision = varargin{1};
            rawData = obj.rawData;rawDivFrame = obj.rawDivFrame;
            rawData = rawData(~isnan(rawDivFrame(:,obj.sortByDivision)),:,:);
            rawDivFrame = rawDivFrame(~isnan(rawDivFrame(:,obj.sortByDivision)),:);
            [~,order] = sort(rawDivFrame(:,obj.sortByDivision));
            rawDivFrameSorted = ones(size(rawDivFrame)).*NaN;
            rawSorted = ones(size(rawData)).*NaN;
            for i = 1:size(rawData,1)
                rawSorted(i,:,:) = rawData(order(i),:,:);
                rawDivFrameSorted(i,:) = rawDivFrame(order(i),:);
            end
            obj.rawSorted = rawSorted;
            obj.rawDivFrameSorted = rawDivFrameSorted;
            
            % set default color range
            for i = 1:size(rawData,3)
                temp = rawData(:,:,i);
                colorRange{i,:} = [prctile(temp(:),5),prctile(temp(:),95)];
            end
            obj.colorRange = colorRange;
            obj.goodTrace = ones(size(rawSorted,1),1);
        end
        
        function obj = set.titles(obj, varargin)
            obj.titles = varargin{:};
        end
        function obj = set.xlabels(obj, varargin)
            obj.xlabels = varargin{:};
        end
        
        function obj = set.ylabels(obj, varargin)
            obj.ylabels = varargin{:};
        end
        
        function obj = set.timeInterval(obj, varargin)
            obj.timeInterval = varargin{:};
        end
        
        function obj = set.timeToPlotBeforeAfterMitosis(obj, varargin)
            obj.timeToPlotBeforeAfterMitosis = varargin{:};
        end
        %% plotting function main
        function obj = smoothDataSortedByMitosis(obj,varargin)
            dataSortedSmooth=NaN(size(obj.dataSorted));
            for ifor = 1:size(obj.dataSorted,3)
                for i = 1:size(obj.dataSorted,1)
                    dataSortedSmooth(i,:,ifor)=nansmooth(obj.dataSorted(i,:,ifor));
                end
            end
            obj.dataSortedSmooth = dataSortedSmooth;
        end
        
        function obj = smoothDataAlignedByMitosis(obj,varargin)
            dataAlignedByMitosisSmooth=NaN(size(obj.dataAligned));
            for ifor = 2:3%:size(obj.dataAligned,3)
                %dataAlignedByMitosisSmooth(:,frameToPlot,ifor)=smooth(obj.dataAlignedByMitosis(:,frameToPlot,ifor));
                %dataAlignedByMitosisSmooth(:,:,ifor)=nansmooth(obj.dataAligned(:,:,ifor));
                
                for i = 1:size(obj.dataAligned,1)
                   %ind=~isnan(obj.dataAligned(i,:,ifor));
                    %dataAlignedByMitosisSmooth(i,ind,ifor)=smooth(obj.dataAligned(i,ind,ifor));
                    dataAlignedByMitosisSmooth(i,:,ifor)=nansmooth(obj.dataAligned(i,:,ifor),5);
                end
            end
            obj.dataAlignedSmooth = dataAlignedByMitosisSmooth;
            %obj.dataAlignedSmooth(:,5010:end,:) = dataAlignedByMitosisSmooth_2(:,5010:end,:);
        end
        
        function obj = plot_heatmap_sort(obj,varargin)
            smoothoption=varargin{1};
            if smoothoption
                for ifor = 1:size(obj.dataSorted,3)
                    for i = 1:size(obj.dataSorted,1)
                        dataSortedSmooth(i,:,ifor)=smooth(obj.dataSorted(i,:,ifor));
                    end
                end
            end
            for ifor = 1:size(obj.dataSorted,3)
                subplot(1,size(obj.dataSorted,3),ifor)
                if smoothoption
                    imagesc(dataSortedSmooth(:,:,ifor),obj.colorRange{ifor})
                else
                    imagesc(obj.dataSorted(:,:,ifor),obj.colorRange{ifor})
                end
                time = [0:size(obj.dataSorted,2)-1]*obj.timeInterval;
                for i = 1:length(obj.sortXTickTime)
                    tickFrame(i,:) = find(time==obj.sortXTickTime(i));
                    set(gca,'XTick',tickFrame);
                    set(gca,'XTickLabel',obj.sortXTickTime);
                end
                colorbar
                obj.ifor = ifor;
                add_plot_info(obj)
            end
        end
        
        function obj = plot_timecourse_align(obj,varargin)
            
            plotType=varargin{1};
            timepoint=varargin{2};
            smoothoption=varargin{3};
            wellNum=varargin{4};
            individual=varargin{5};
            CDK_ERK=varargin{6};
            gating=varargin{7};
            align=varargin{8};
            sample=varargin{9};
            inhibitionFrame=varargin{10};
            
            samplePercent=sample; %control; 34, Nutlin; 19, Tenovin; 30
            samplenum=round(sum(obj.cellState)*(samplePercent));
            orgidx=find(obj.cellState==1);
            cellState=randsample(orgidx,samplenum);
            samplenum=round(sum(obj.cellDelay)*(samplePercent));
            orgidx=find(obj.cellDelay==1);
            cellDelay=randsample(orgidx,samplenum);
            samplenum=round(sum(obj.cellLow)*(samplePercent));
            orgidx=find(obj.cellLow==1);
            cellLow=randsample(orgidx,samplenum);            
            
            figure; hold on;
            if plotType==1
                numFramesToPlot = obj.timeToPlotBeforeAfterMitosis./obj.timeInterval;
                timeToPlot = linspace(-obj.timeToPlotBeforeAfterMitosis,obj.timeToPlotBeforeAfterMitosis,numFramesToPlot*2);
                frameToPlot = 5000-numFramesToPlot+1:5000+numFramesToPlot;
                if CDK_ERK==1
                    if any(obj.cellState) | any(obj.cellLow)
                        if smoothoption && ~individual
                            for ifor = 1:size(obj.dataAlignedSmooth,3)
                                subplot(1,size(obj.dataAlignedSmooth,3),ifor)
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAlignedSmooth(cellState,frameToPlot,ifor)),nanstd(obj.dataAlignedSmooth(obj.cellState,frameToPlot,ifor))/sqrt(wellNum),'b',1);hold on
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAlignedSmooth(cellLow,frameToPlot,ifor)),nanstd(obj.dataAlignedSmooth(obj.cellLow,frameToPlot,ifor))/sqrt(wellNum),'r',1);
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAlignedSmooth(cellDelay,frameToPlot,ifor)),nanstd(obj.dataAlignedSmooth(obj.cellDelay,frameToPlot,ifor))/sqrt(wellNum),'g',1);
                                obj.ifor = ifor;add_plot_info(obj)
                            end
                        elseif smoothoption && individual
                            for ifor = 1:size(obj.dataAligned,3)
                                subplot(1,size(obj.dataAligned,3),ifor)
                                    if align
                                        if any(obj.cellState)
                                            plot(timeToPlot,obj.dataAlignedSmooth(cellState,frameToPlot,1),'Color','b');hold on
                                        end
                                        plot(timeToPlot,obj.dataAlignedSmooth(cellDelay,frameToPlot,1),'Color','g');hold on
                                        plot(timeToPlot,obj.dataAlignedSmooth(cellLow,frameToPlot,1),'Color','r');
                                    else
                                        if any(obj.cellState)
                                            plot(obj.dataSortedSmooth(cellState,:,1)','Color','b');hold on
                                        end
                                        plot(obj.dataSortedSmooth(cellDelay,:,1)','Color','g');hold on
                                        plot(obj.dataSortedSmooth(cellLow,:,1)','Color','r');
                                    end
                                    ylim([0 3]);
                                    xlim([inhibitionFrame-29 105]);
                                    vline([inhibitionFrame-27 inhibitionFrame-25],'g:'); 
                                    vline(inhibitionFrame);
                                    set(gca,'XTick',[inhibitionFrame-25:25:150]);
                                    set(gca,'XTickLabel',[-5:5:30]);
                                    xlabel('Time since drug addition (hr)');
                                    ylabel('Cdk2 activity');
                                    obj.ifor = ifor;
                                
                            end
                        elseif individual && ~smoothoption
                            for ifor = 1:size(obj.dataAligned,3)
                                subplot(1,size(obj.dataAligned,3),ifor)
                                plot(timeToPlot,obj.dataAligned(obj.cellState,frameToPlot,ifor),'Color','b');hold on
                                plot(timeToPlot,obj.dataAligned(obj.cellLow,frameToPlot,ifor),'Color','r');
                                plot(timeToPlot,obj.dataAligned(obj.cellDelay,frameToPlot,ifor),'Color','g');
                                obj.ifor = ifor;add_plot_info(obj)
                                if ifor==1
                                    ylim([0 3])
                                    xlim([-5 15]);
                                else
                                    ylim([1.1 2])
                                    xlim([-5 15]);
                                end
                            end
                        else
                            for ifor = 1:size(obj.dataAligned,3)
                                subplot(1,size(obj.dataAligned,3),ifor)
                                %shadedErrorBar(timeToPlot,nanmean(obj.dataAligned(obj.cellState,frameToPlot,ifor)),nanstd(obj.dataAligned(obj.cellState,frameToPlot,ifor))/sqrt(wellNum),'b',1);hold on
                                %shadedErrorBar(timeToPlot,nanmean(obj.dataAligned(obj.cellLow,frameToPlot,ifor)),nanstd(obj.dataAligned(~obj.cellState,frameToPlot,ifor))/sqrt(wellNum),'r',1);
                                %shadedErrorBar(timeToPlot,nanmean(obj.dataAligned(obj.cellDelay,frameToPlot,ifor)),nanstd(obj.dataAligned(~obj.cellState,frameToPlot,ifor))/sqrt(wellNum),'g',1);
                                plot(timeToPlot,nanmean(obj.dataAligned(cellState,frameToPlot,ifor)),'b');hold on
                                plot(timeToPlot,nanmean(obj.dataAligned(cellState,frameToPlot,ifor))-nanstd(obj.dataAligned(obj.cellState,frameToPlot,ifor))/sqrt(wellNum),'b:');
                                plot(timeToPlot,nanmean(obj.dataAligned(cellState,frameToPlot,ifor))+nanstd(obj.dataAligned(obj.cellState,frameToPlot,ifor))/sqrt(wellNum),'b:');
                                plot(timeToPlot,nanmean(obj.dataAligned(cellLow,frameToPlot,ifor)),'r');
                                plot(timeToPlot,nanmean(obj.dataAligned(cellLow,frameToPlot,ifor))-nanstd(obj.dataAligned(obj.cellLow,frameToPlot,ifor))/sqrt(wellNum),'r:');
                                plot(timeToPlot,nanmean(obj.dataAligned(cellLow,frameToPlot,ifor))+nanstd(obj.dataAligned(obj.cellLow,frameToPlot,ifor))/sqrt(wellNum),'r:');
                                plot(timeToPlot,nanmean(obj.dataAligned(cellDelay,frameToPlot,ifor)),'g');
                                plot(timeToPlot,nanmean(obj.dataAligned(cellDelay,frameToPlot,ifor))-nanstd(obj.dataAligned(obj.cellDelay,frameToPlot,ifor))/sqrt(wellNum),'g:');
                                plot(timeToPlot,nanmean(obj.dataAligned(cellDelay,frameToPlot,ifor))+nanstd(obj.dataAligned(obj.cellDelay,frameToPlot,ifor))/sqrt(wellNum),'g:');
                                obj.ifor = ifor;add_plot_info(obj)
                                switch gating
                                        case 'G1'
                                            xlim([-3 20]);
                                        case 'G2'
                                            xlim([-10 20]);
                                    end
                                if ifor==2
                                    ylim([1 1.6])
                                end
                            end
                        end
                    else
                        if smoothoption && ~individual
                            for ifor = 1:size(obj.dataAlignedSmooth,3)
                                subplot(1,size(obj.dataAlignedSmooth,3),ifor)
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAlignedSmooth(:,frameToPlot,ifor)),nanstd(obj.dataAlignedSmooth(:,frameToPlot,ifor))/sqrt(wellNum),'b',1);
                                obj.ifor = ifor;add_plot_info(obj)
                            end
                        elseif smoothoption && individual
                            for ifor = 1:size(obj.dataAlignedSmooth,3)
                                subplot(1,size(obj.dataAlignedSmooth,3),ifor)
                                %subplot(1,size(obj.dataAligned,3),ifor)
                                if ifor==1
                                    plot(timeToPlot,obj.dataAlignedSmooth(:,frameToPlot,1),'Color','b');
                                else
                                    plot(timeToPlot,obj.dataAligned(:,frameToPlot,2),'Color','b');
                                end
                                obj.ifor = ifor;add_plot_info(obj)
                            end
                        else
                            for ifor = 1:size(obj.dataAligned,3)
                                subplot(1,size(obj.dataAligned,3),ifor)
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAligned(:,frameToPlot,ifor)),nanstd(obj.dataAligned(:,frameToPlot,ifor))/sqrt(wellNum),'b',1);
                                obj.ifor = ifor;add_plot_info(obj)
                            end
                        end
                    end
                elseif CDK_ERK==2
                    counter=0;
                    if any(obj.cellERKinc)
                        if smoothoption && ~individual
                            for ifor = 1:size(obj.dataAlignedSmooth,3)
                                subplot(1,size(obj.dataAlignedSmooth,3),ifor)
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAlignedSmooth(obj.cellERKinc,frameToPlot,ifor)),nanstd(obj.dataAlignedSmooth(obj.cellERKinc,frameToPlot,ifor))/sqrt(wellNum),'b',1);hold on
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAlignedSmooth(obj.cellERKlow,frameToPlot,ifor)),nanstd(obj.dataAlignedSmooth(obj.cellERKlow,frameToPlot,ifor))/sqrt(wellNum),'r',1);
                                obj.ifor = ifor;add_plot_info(obj)
                            end
                        elseif smoothoption && individual
                            for ifor = 1%:size(obj.dataAligned,3)
                                figure;
                                counter=counter+1;
                                subplot(1,2,counter);
                                plot(timeToPlot,obj.dataAlignedSmooth(obj.cellERKlow,frameToPlot,ifor),'Color','r');
                                counter=counter+1;
                                subplot(1,2,counter);
                                plot(timeToPlot,obj.dataAlignedSmooth(obj.cellERKinc,frameToPlot,ifor),'Color','b');
                                
                                obj.ifor = ifor;add_plot_info(obj)
                            end
                        elseif individual && ~smoothoption
                            for ifor = 1%1:size(obj.dataAligned,3)
                                %                         subplot(1,size(obj.dataAligned,3),ifor)
                                %                         plot(timeToPlot,obj.dataAligned(obj.cellERKinc,frameToPlot,ifor),'Color','b');hold on
                                %                         plot(timeToPlot,obj.dataAligned(obj.cellERKlow,frameToPlot,ifor),'Color','r'); %hold on
                                figure;
                                counter=counter+1;
                                subplot(1,2,counter);
                                plot(timeToPlot,obj.dataAligned(obj.cellERKlow,frameToPlot,ifor),'Color','r');
                                %                                 ylim([1.1 2.3]);
                                counter=counter+1;
                                subplot(1,2,counter);
                                plot(timeToPlot,obj.dataAligned(obj.cellERKinc,frameToPlot,ifor),'Color','b');
                                %                                 ylim([1.1 2.3]);
                                obj.ifor = ifor;add_plot_info(obj)
                            end
                        else
                            for ifor = 1:size(obj.dataAligned,3)
                                subplot(1,size(obj.dataAligned,3),ifor)
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAligned(obj.cellERKinc,frameToPlot,ifor)),nanstd(obj.dataAligned(obj.cellERKinc,frameToPlot,ifor))/sqrt(wellNum),'b',1);hold on
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAligned(obj.cellERKlow,frameToPlot,ifor)),nanstd(obj.dataAligned(obj.cellERKlow,frameToPlot,ifor))/sqrt(wellNum),'r',1);
                                obj.ifor = ifor;add_plot_info(obj)
                            end
                        end
                    else
                        if smoothoption
                            for ifor = 1:size(obj.dataAlignedSmooth,3)
                                subplot(1,size(obj.dataAlignedSmooth,3),ifor)
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAlignedSmooth(:,frameToPlot,ifor)),nanstd(obj.dataAlignedSmooth(:,frameToPlot,ifor))/sqrt(wellNum),'b',1);
                                obj.ifor = ifor;add_plot_info(obj)
                            end
                        else
                            for ifor = 1:size(obj.dataAligned,3)
                                subplot(1,size(obj.dataAligned,3),ifor)
                                shadedErrorBar(timeToPlot,nanmean(obj.dataAligned(:,frameToPlot,ifor)),nanstd(obj.dataAligned(:,frameToPlot,ifor))/sqrt(wellNum),'b',1);
                                obj.ifor = ifor;add_plot_info(obj)
                            end
                        end
                    end
                end
            else
                frameToPlot = 5000 + timepoint/obj.timeInterval;
                for ifor = 1%:size(obj.dataAlignedSmooth,3)
                    if ifor==1
                        bmin=0; bmax=2.5; bstep=(bmax-bmin)/28; bin=bmin:bstep:bmax; ymax=20; bi=12;
                    else
                        bmin=1; bmax=2; bstep=(bmax-bmin)/30; bin=bmin:bstep:bmax; ymax=20;
                    end
                    figure(ifor);
                    subplot(1,3,1);
                    Abvals=mean(obj.dataAlignedSmooth(:,frameToPlot-1:frameToPlot+1,ifor),2);
                    pdfvals=histc(Abvals,bin);
                    pdfvals=100*pdfvals/sum(pdfvals);
                    cdfvals=cumsum(pdfvals);
                    cdfvals=100*cdfvals/max(cdfvals);
                    cdfvals(cdfvals==max(cdfvals))=99;
                    line(bin,pdfvals,'linewidth',2,'color','r'); hold on;
                    line(bin(bi),5:10,'color','k');
                    ylim([0 ymax]);
                    thresh=round(cdfvals(bi));
                    title(['fraction (' num2str(bi) ') hypo=',num2str(thresh),' hyper=',num2str(100-thresh)],'FontSize',9);
                    
                    subplot(1,3,2);
                    Abvals=mean(obj.dataAlignedSmooth(obj.cellERKinc,frameToPlot-1:frameToPlot+1,ifor),2);
                    pdfvals=histc(Abvals,bin);
                    pdfvals=100*pdfvals/sum(pdfvals);
                    cdfvals=cumsum(pdfvals);
                    cdfvals=100*cdfvals/max(cdfvals);
                    cdfvals(cdfvals==max(cdfvals))=99;
                    line(bin,pdfvals,'linewidth',2,'color','r'); hold on;
                    line(bin(bi),5:10,'color','k');
                    ylim([0 ymax]);
                    thresh=round(cdfvals(bi));
                    title(['fraction (' num2str(bi) ') hypo=',num2str(thresh),' hyper=',num2str(100-thresh)],'FontSize',9);
                    
                    subplot(1,3,3);
                    Abvals=mean(obj.dataAlignedSmooth(obj.cellERKlow,frameToPlot-1:frameToPlot+1,ifor),2);
                    pdfvals=histc(Abvals,bin);
                    pdfvals=100*pdfvals/sum(pdfvals);
                    cdfvals=cumsum(pdfvals);
                    cdfvals=100*cdfvals/max(cdfvals);
                    cdfvals(cdfvals==max(cdfvals))=99;
                    line(bin,pdfvals,'linewidth',2,'color','r'); hold on;
                    line(bin(bi),5:10,'color','k');
                    ylim([0 ymax]);
                    thresh=round(cdfvals(bi));
                    title(['fraction (' num2str(bi) ') hypo=',num2str(thresh),' hyper=',num2str(100-thresh)],'FontSize',9);
                end
            end
        end
        
        function obj = plot_heatmap_align(obj)
            numFramesToPlot = obj.timeToPlotBeforeAfterMitosis./obj.timeInterval;
            %             timeToPlot = linspace(-obj.timeToPlotBeforeAfterMitosis,obj.timeToPlotBeforeAfterMitosis,numFramesToPlot*2);
            frameToPlot = 5000-numFramesToPlot+1:5000+numFramesToPlot;
            for ifor = 1:size(obj.dataAligned,3)
                subplot(1,size(obj.dataAligned,3),ifor)
                imagesc([obj.dataAligned(:,frameToPlot,ifor)],obj.colorRange{ifor});
                %                 colormap('jet_black');
                set(gca,'XTick',[1,ceil(length(frameToPlot)/2),length(frameToPlot)]);
                set(gca,'XTickLabel',[-obj.timeToPlotBeforeAfterMitosis,0,obj.timeToPlotBeforeAfterMitosis]);
                obj.ifor = ifor;add_plot_info(obj)
                %                 xlim([-20 10]);
                colorbar
            end
        end
        
        
        
        %% assorted
        function obj = judgeERK(obj,varargin)
            if nargin==4
                time = varargin{1};
                percentile = varargin{2};
                cellcycle = varargin{3};
            else
                sprintf('enter duration')
            end
            
            if cellcycle==1
                cycleFrame = [5000-round((time+2)/obj.timeInterval),5000-11];
            elseif cellcycle==2
                cycleFrame = [5001,5000+round(time/obj.timeInterval)];
            else cellcycle==3
                cycleFrame = [5000-round((time+2)/obj.timeInterval),5000-11,5000,5000+round(time/obj.timeInterval)];
            end
            percentilethresh=[100-percentile percentile];
            
            if cellcycle==3
                dataCycPhase(:,1) = nanmean(obj.dataAligned(:,[cycleFrame(1,1):cycleFrame(1,2),cycleFrame(1,3):cycleFrame(1,4)],1),2);
            else
                dataCycPhase(:,1) = nanmean(obj.dataAligned(:,cycleFrame(1,1):cycleFrame(1,2),1),2);
            end
            obj.cellERKinc  = dataCycPhase(:,1)>prctile(dataCycPhase(:,1),percentilethresh(1));
            obj.cellERKlow = dataCycPhase(:,1)<prctile(dataCycPhase(:,1),percentilethresh(2));
            %{
            subplot(1,3,1);histogram(dataCycPhase,'Normalization','probability'); xlim([1.1 1.6]);
            vline(prctile(dataCycPhase(:,1),percentilethresh(1)));vline(prctile(dataCycPhase(:,1),percentilethresh(2)));
            subplot(1,3,2);histogram(dataCycPhase(obj.cellERKinc),'Normalization','probability'); xlim([1.1 1.6]);
            subplot(1,3,3);histogram(dataCycPhase(obj.cellERKlow),'Normalization','probability'); xlim([1.1 1.6]);
            figu(0.5,1);
            %}
        end
        
        function obj = plotIFhist(obj,varargin)
            stain=varargin{1};
            protein=varargin{2};
            DisplayOption=varargin{3};%'Scatter'; %Panel Scatter Box
            Condi=varargin{4};
            resultdir=varargin{5};
            Inhibitiongate=varargin{6};
            option=varargin{7};
            Hoechstval=obj.IFdataSorted(:,3).*obj.IFdataSorted(:,4);
            %figure;histogram(Hoechstval,200); vline([2500000 6000000]);
            
            switch Inhibitiongate
                case 'G1'
                    HoechstGate=2500000<Hoechstval & Hoechstval<6000000; %gate G1 cells
                case 'G2'
                    HoechstGate=1000000<Hoechstval & Hoechstval<7000000; %gate G1 c
                    figure;histogram(Hoechstval,200); vline([1000000 7000000]);
                case 'S'
                    if ismember(Condi,[3 4 9 12 15 16 21 22])
                        HoechstGate=2800000<Hoechstval & Hoechstval<6300000;
                        %                          figure;histogram(Hoechstval,200); vline([2800000 6300000]);
                    elseif ismember(Condi,[5 6 11 17 18 23 24])
                        HoechstGate=2400000<Hoechstval & Hoechstval<6000000;
                        %                         figure;histogram(Hoechstval,200); vline([2400000 6000000]);
                    else
                        HoechstGate=4000000<Hoechstval & Hoechstval<8400000;
                        %                         figure;histogram(Hoechstval,200); vline([4000000 8400000]);
                    end
                case 'All'
                    HoechstGate=1000000<Hoechstval & Hoechstval<12000000;
            end
            
            if stain==1
                Abvals=obj.IFdataSorted(:,protein); %mRNA
                posval=Abvals>=0;
            elseif stain==2
                Abvals=obj.IFdataSorted(:,protein); %CyclinD or p21
                posval=Abvals>=0;
                %posval=Abvals>=1; Abvals(Abvals<1)=1; Abvals=log2(Abvals);
                
                %                 Abvalsx=obj.IFdataSorted(:,6); %CyclinD or p21
                %                 posval=Abvalsx>=1; Abvalsx(Abvalsx<1)=1; Abvalsx=log2(Abvalsx);
            elseif stain==3
                Abvals=obj.IFdataSorted(:,protein)./obj.IFdataSorted(:,5); %pRb/tRb
                posval=Abvals>=0;
            end
            %bstep=(bmax-bmin)/100; bin=bmin:bstep:bmax; %bi=53; %50; 26
            %bmin=6; bmax=14; bstep=(bmax-bmin)/50; bin=bmin:bstep:bmax; ymax=30; bi=20; %50; 26
            switch DisplayOption
                case {'Panel'}
                    if stain==1 && ismember(protein,[6:13])
                        data=Abvals(posval & HoechstGate);
                        obj.mRNAarea=data;
                        %                         figure;
                        %                         subplot(1,2,1);
                        %                         bmin=0; bmax=150; binNum=35;
                        %                         [h]=histogram(data,binNum,'Normalization','probability','BinLimits',[bmin,bmax]);
                        %                         xlim([bmin bmax]);
                        %                         ylabel('mRNA area');
                        %
                        %                         subplot(1,2,2);
                        %                         boxplot(data,'notch','on','symbol','');
                        %                         set(gca,'XTickLabel',{'Cyclin D2 mRNA'});
                        %                         ylim([-5 bmax]);
                        %                         ylabel('mRNA area');
                    elseif stain==1 && protein==9
                        data=Abvals(posval & HoechstGate);
                        obj.mRNAintensity=data;
                        figure;
                        subplot(1,2,1);
                        bmin=100; bmax=3500; binNum=35;
                        [h]=histogram(data,binNum,'Normalization','probability','BinLimits',[bmin,bmax]);
                        xlim([bmin bmax]);
                        xlabel('mRNA intensity');
                        ylabel('cell fraction');
                        
                        subplot(1,2,2);
                        boxplot(data,'notch','on','symbol','');
                        set(gca,'XTickLabel',{'Cyclin D2 mRNA'});
                        ylim([bmin bmax]);
                        ylabel('mRNA intensity');
                    elseif stain==1 && protein==11
                        data=Abvals(posval & HoechstGate);
                        obj.mRNAcount=data;
                        figure;
                        subplot(1,2,1);
                        bmin=0; bmax=30; binNum=35;
                        [h]=histogram(data,binNum,'Normalization','probability','BinLimits',[bmin,bmax]);
                        xlim([bmin bmax]);
                        ylabel('mRNA count');
                        
                        subplot(1,2,2);
                        boxplot(data,'notch','on','symbol','');
                        set(gca,'XTickLabel',{'Cyclin D2 mRNA'});
                        ylim([-5 bmax]);
                        ylabel('mRNA count');
                    elseif stain==2
                        data=Abvals(posval & HoechstGate);
                        obj.Protein=data;
                        %                         figure;
                        %                         subplot(1,2,1);
                        %                         bmin=9; bmax=15; binNum=35;
                        %                         [h]=histogram(data,binNum,'Normalization','probability','BinLimits',[bmin,bmax]);
                        %                         xlim([bmin bmax]);
                        %                         ylabel('cell fraction');
                        %                         xlabel('protein intensity');
                        %
                        %                         subplot(1,2,2);
                        %                         boxplot(data,'notch','on','symbol','');
                        %                         set(gca,'XTickLabel',{'MYC'});
                        %                         ylim([bmin bmax]);
                        %                         ylabel('MYC Protein log2');
                    elseif stain==3
                        %                         data=Abvals(posval & HoechstGate);
                        %                         obj.Protein=data;
                        %                         bmin=0; bmax=6; binNum=50;
                        %                         ymax=1.12; thresh=1.44;
                        %                         [h]=histogram(data,binNum,'Normalization','probability','BinLimits',[bmin,bmax]);
                        %                         hypoRb=sum(h.Values(h.BinEdges>=0 & h.BinEdges<=thresh));
                        %                         hyperRb=sum(h.Values(h.BinEdges>thresh & h.BinEdges<h.BinLimits(2)));
                        %                         xlim([bmin bmax]);%ylim([0 ymax]);
                        %
                        %                         vline(thresh); title(['hypo=',num2str(round(hypoRb,2)),' hyper=',num2str(round(hyperRb,2))],'FontSize',9);
                        
                        data=Abvals(posval & HoechstGate);
                        obj.Protein=data;
                        bmin=0; bmax=1.6; binNum=50;
                        ymax=1.12; thresh=0.8;
                        [h]=histogram(data,binNum,'Normalization','probability','BinLimits',[bmin,bmax]);
                        hypoRb=sum(h.Values(h.BinEdges>=0 & h.BinEdges<=thresh));
                        obj.hyperpRb=sum(h.Values(h.BinEdges>thresh & h.BinEdges<h.BinLimits(2)));
                        xlim([bmin bmax]);%ylim([0 ymax]);
                        
                        vline(thresh); title(['hypo=',num2str(round(hypoRb,2)),' hyper=',num2str(round(obj.hyperpRb,2))],'FontSize',9);
                    end
                case 'Scatter'
                    figure;
                    AbvalTotal=Abvals(posval & HoechstGate);
                    %Abvalsx=Hoechstval(posval & HoechstGate);
                    Abvalsx=Abvalsx(posval & HoechstGate);
                    dscatter(Abvalsx,AbvalTotal);
                    ylabel('Cyclin D2 mRNA count'); xlabel('Hoechst');
                    xlim([1500000 11000000]); ylim([0 100]);
                case 'Box'
                    figure;
                    %                     subplot(1,3,1);
                    data=Abvals(posval & HoechstGate);
                    obj.Protein=data;
                    %temp2=Abvals(posval & HoechstGate & obj.cellERKinc);
                    %temp3=Abvals(posval & HoechstGate & obj.cellERKlow);
                    %data=nan(length(temp1),3);
                    
                    %data(1:length(temp1),1)=temp1;
                    %data(1:length(temp2),2)=temp2;
                    %data(1:length(temp3),3)=temp3;
                    boxplot(data,'notch','on','symbol','');
                    %set(gca,'XTickLabel',{'Total';'ERK high';'ERK low'});
                    %p(1)=ranksum(temp1,temp2);p(2)=ranksum(temp1,temp3);p(3)=ranksum(temp2,temp3);
                    %                     if protein==6 %CyclinD
                    %                         ylim([8 14]); ylabel('CyclinD (log2)');
                    %                     else protein==7 %p21
                    %                         ylim([6 14]); ylabel('p21 (log2)');
                    %                     end
            end
        end
        
        function add_plot_info(obj)
            if ~isempty(obj.xlabels);
                if length(obj.xlabels)==1
                    xlabel(obj.xlabels{1});
                else
                    xlabel(obj.xlabels{obj.ifor});
                end
            end
            if ~isempty(obj.ylabels);
                if length(obj.ylabels)==1
                    ylabel(obj.ylabels{1});
                else
                    ylabel(obj.ylabels{obj.ifor});
                end
            end
            if ~isempty(obj.titles);
                if length(obj.titles)==1
                    title(obj.titles{1});
                else
                    title(obj.titles{obj.ifor});
                end
            end
        end
        
        function boxplot_median_erk(obj,varargin)
            if nargin==4
                G2 = varargin{1};
                S = varargin{1} + varargin{2};
                G1 = varargin{1} + varargin{2} + varargin{3};
            else
                sprintf('enter G2, S and G1 duration')
            end
            
            cycleFrame = [5000-round(G1/obj.timeInterval),5000-round(S/obj.timeInterval)-1;...    % G1 frame
                5000-round(S/obj.timeInterval),5000-round(G2/obj.timeInterval)-1;...      % S frame
                5000-round((G2+2)/obj.timeInterval),5000-11];                                  % G2 frame
            
            for i = 1:3
                data(:,i) = nanmedian(obj.dataAligned(:,cycleFrame(i,1):cycleFrame(i,2),2),2);
            end
            boxplot(data,'notch','on','symbol','');
            set(gca,'XTickLabel',{'G1';'S';'G2'}); ylim([1.15 1.8])
            
            figure;
            anovavalues = [data(:,1);data(:,2);data(:,3)];
            anovalabels = [repmat('G1',size(data(:,1)));repmat('S ',size(data(:,2)));repmat('G2',size(data(:,3)))];
            
            [p,s,stats] = anova1(anovavalues,anovalabels,'off');
            [c,m,h,nms] = multcompare(stats);
        end
        
        function boxplot_median_erk_cyc_qui_separated(obj,varargin)
            if nargin==4
                G2 = varargin{1};
                S = varargin{1} + varargin{2};
                G1 = varargin{1} + varargin{2} + varargin{3};
            else
                sprintf('enter G2, S and G1 duration')
            end
            
            
            dataCyc = obj.dataAligned(obj.cellState,:,2);
            dataQui = obj.dataAligned(~obj.cellState,:,2);
            
            cycleFrame = [5000-round(G1/obj.timeInterval),5000-round(S/obj.timeInterval)-1;...    % G1 frame
                5000-round(S/obj.timeInterval),5000-round(G2/obj.timeInterval)-1;...      % S frame
                5000-round((G2+2)/obj.timeInterval),5000-11];                                  % G2 frame
            
            for i = 1:3
                dataCycPhase(:,i) = nanmedian(dataCyc(:,cycleFrame(i,1):cycleFrame(i,2)),2);
                dataQuiPhase(:,i) = nanmedian(dataQui(:,cycleFrame(i,1):cycleFrame(i,2)),2);
            end
            figure;
            subplot(1,2,1)
            boxplot(dataCycPhase,'notch','on','symbol','');
            set(gca,'XTickLabel',{'G1';'S';'G2'});title('Cycling');ylim([1.25 1.75])
            subplot(1,2,2)
            boxplot(dataQuiPhase,'notch','on','symbol','');
            set(gca,'XTickLabel',{'G1';'S';'G2'});title('Quiescent');ylim([1.25 1.75])
            
            figure;
            anovavalues = [dataCycPhase(:,1);dataCycPhase(:,2);dataCycPhase(:,3)];
            anovalabels = [repmat('G1',size(dataCycPhase(:,1)));repmat('S ',size(dataCycPhase(:,2)));repmat('G2',size(dataCycPhase(:,3)))];
            [p,s,stats] = anova1(anovavalues,anovalabels,'off');
            [c,m,h,nms] = multcompare(stats);
            title('Cycling')
            
            figure;
            anovavalues = [dataQuiPhase(:,1);dataQuiPhase(:,2);dataQuiPhase(:,3)];
            anovalabels = [repmat('G1',size(dataQuiPhase(:,1)));repmat('S ',size(dataQuiPhase(:,2)));repmat('G2',size(dataQuiPhase(:,3)))];
            [p,s,stats] = anova1(anovavalues,anovalabels,'off');
            [c,m,h,nms] = multcompare(stats);
            title('Quiescent')
        end
        
        function bargraph_mean_erk_cyc_qui_separated(obj,varargin)
            if nargin==4
                G2 = varargin{1};
                S = varargin{1} + varargin{2};
                G1 = varargin{1} + varargin{2} + varargin{3};
            else
                sprintf('enter G2, S and G1 duration')
            end
            
            
            dataCyc = obj.dataAligned(obj.cellState,:,2);
            dataQui = obj.dataAligned(~obj.cellState,:,2);
            
            cycleFrame = [5000-round(G1/obj.timeInterval),5000-round(S/obj.timeInterval)-1;...    % G1 frame
                5000-round(S/obj.timeInterval),5000-round(G2/obj.timeInterval)-1;...      % S frame
                5000-round((G2+2)/obj.timeInterval),5000-11];                                  % G2 frame
            
            for i = 1:3
                dataCycPhase(:,i) = nanmedian(dataCyc(:,cycleFrame(i,1):cycleFrame(i,2)),2);
                mediandataCycPhase(i) = nanmedian(dataCycPhase(:,i));
                stddataCycPhase(i) = nanstd(dataCycPhase(:,i))/sqrt(10);
                dataQuiPhase(:,i) = nanmedian(dataQui(:,cycleFrame(i,1):cycleFrame(i,2)),2);
                medianandataQuiPhase(i) = nanmedian(dataQuiPhase(:,i));
                stddataQuiPhase(i) = nanstd(dataQuiPhase(:,i))/sqrt(10);
            end
            barweb([mediandataCycPhase' medianandataQuiPhase'], [stddataCycPhase' stddataQuiPhase'], 1, [{'G1';'S';'G2'}], ['ERK activity'], [' '], ['ERK activity'], [], [], [{'CDKinc';'CDKlow'}], 1, 'axis')
            ylim([1.4 1.6]);
            
            figure;
            anovavalues = [dataCycPhase(:,1);dataCycPhase(:,2);dataCycPhase(:,3);dataQuiPhase(:,1);dataQuiPhase(:,2);dataQuiPhase(:,3)];
            anovalabels = [repmat('CDK inc G1',size(dataCycPhase(:,1)));repmat('CDK inc S ',size(dataCycPhase(:,2)));repmat('CDK inc G2',size(dataCycPhase(:,3)));repmat('CDK low G1',size(dataQuiPhase(:,1)));repmat('CDK low S ',size(dataQuiPhase(:,2)));repmat('CDK low G2',size(dataQuiPhase(:,3)))];
            [p,s,stats] = anova1(anovavalues,anovalabels,'off');
            [c,m,h,nms] = multcompare(stats);
            title('ERK activity')
        end
        
        function obj = plot_heatmap_alignByMitosis(obj)
            numFramesToPlot = obj.timeToPlotBeforeAfterMitosis./obj.timeInterval;
            %             timeToPlot = linspace(-obj.timeToPlotBeforeAfterMitosis,obj.timeToPlotBeforeAfterMitosis,numFramesToPlot*2);
            frameToPlot = 5000-numFramesToPlot+1:5000+numFramesToPlot;
            for ifor = 1:size(obj.dataAligned,3)
                subplot(1,size(obj.dataAligned,3),ifor)
                imagesc([obj.dataAligned(:,frameToPlot,ifor)],obj.colorRange{ifor});
                set(gca,'XTick',[1,ceil(length(frameToPlot)/2),length(frameToPlot)]);
                set(gca,'XTickLabel',[-obj.timeToPlotBeforeAfterMitosis,0,obj.timeToPlotBeforeAfterMitosis]);
                obj.ifor = ifor;add_plot_info(obj)
                colorbar
            end
        end
        
        
         %% assorted
        function obj = judgecdk2_detect_rising_pts(obj,varargin)
            % some codes are arbitrary
            % This script will fit a linear line to CDK2 from different time point.
            % It will determine CDK2Inc/CDK2Low based on the diffence between data and the best-fitted line.
            % It will give you a rising point. At the end you will get
            % several plots, and you may have to suspend it for the last one
            % It will do some cleaning as well
            ptsToFit = 6/obj.timeInterval;
            slope = 0.04;
            lineToFit = [1:ptsToFit].*slope;    % This line will be fitted by changing starting point.
            startThres = 1.5;                   % start point must be lower than this
            riseThres = 1.5;                    % cdkInc cells must have a point above this
            errorIncLowThres = 15;              % Threshold of absolute error between data and the best-fitted line
            cdklowMaxThres = 1;                 % used for cleaning. If Cdk2Low, it shouldn't go above this
            cdkDropThres = -0.5;
            cdkShouldRiseAtLeast = 12;
            
            divFrame = obj.divFrameAfterAlignment(:,2);
            
            smoothCdk2 = tsmovavg(obj.dataAligned(:,:,1),'s',ceil(1/obj.timeInterval),2);
            %             smoothCdk2 = obj.dataAlignedSmooth(:,:,1);
            
            % first remove traces if it's second division is somehow within ptsToFit or
            % if it doesn't have ptsToFit points after the first mitosis
            secDivisionTooSoon = divFrame < 5000+ptsToFit+2+ceil(6./obj.timeInterval);
            smoothCdk2 = smoothCdk2(~secDivisionTooSoon,:);
            divFrame = divFrame(~secDivisionTooSoon,:);
            badCell = find(~secDivisionTooSoon==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            atLeastLengthNan = ~logical(sum(isnan(smoothCdk2(:,5000:5000+ptsToFit)),2));
            smoothCdk2 = smoothCdk2(atLeastLengthNan,:);
            divFrame = divFrame(atLeastLengthNan,:);
            badCell = find(atLeastLengthNan==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            
            % fitting
            for j = 1:size(smoothCdk2,1)
                wholedata = smoothCdk2(j,5000:end);% traces after first division
                firstNan = find(isnan(wholedata),1,'first');
                firstNan = min(firstNan, divFrame(j)-5000-ceil(6./obj.timeInterval));% assume rising point should be at least 6hours before next division
                fit = zeros(firstNan-1-ptsToFit,1);
                % find the best fit and where it starts
                for starti = 1:firstNan-1-ptsToFit
                    fit(starti) = sum(abs(wholedata(starti:starti+ptsToFit-1) - (lineToFit + wholedata(starti))));
                    if wholedata(starti)>startThres   % start point shouldn't be too large
                        fit(starti) = Inf;
                    end
                    if max(wholedata) < riseThres     % if Cdk2Inc it should go above this
                        fit(starti) = Inf;
                    end
                end
                [fitDif(j), fitId(j)] = min(fit);     % the error and the rising point of the best fit
                
                % Adjust the rising point
                % after finding the best fit, find local minima on the left. If there is, update.
                localMinima = find([0;0;diff(sign(diff(wholedata')))]>0)-1;
                localMinima = localMinima(localMinima <= fitId(j));
                [~, id] = min(abs(localMinima - fitId(j)));
                if isempty(localMinima)
                    correctedFitId(j) = fitId(j);
                else
                    correctedFitId(j) = localMinima(id);
                end
                
            end
            correctedFitId = correctedFitId-1;% relative to 5000. 5000+id is the absolute point
            fitId = fitId-1;
            
            % Determine Cdk2Inc and Cdk2low cells
            cellInc = fitDif<errorIncLowThres;
            
            % some cleaning. if it is Cdk2low, then maximum value after several frames
            % shouldn't be above 1 or so.
            maxAfterFirstDivision = max(obj.dataAlignedSmooth(:,5000+4/obj.timeInterval:end,1),[],2);
            withinCdklowMaxThres = maxAfterFirstDivision < cdklowMaxThres;
            badLogical = ~cellInc' & ~withinCdklowMaxThres;
            badCell = find(badLogical==1);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            fitDif = fitDif(~badLogical);
            fitId = fitId(~badLogical);
            correctedFitId = correctedFitId(~badLogical);
            cellInc = cellInc(~badLogical);
            
            %if rise, then it shouldn't go down too early
            smoothCdk2_2 = tsmovavg(obj.dataAlignedSmooth(:,:,1),'s',ceil(1/obj.timeInterval),2);
            cdk2Drop = circshift(smoothCdk2_2',-ceil(0.6/obj.timeInterval))' - smoothCdk2_2 < cdkDropThres; % recognize as drop when change in 30 min exceeds 0.5
            cdk2Drop = cdk2Drop(:,5000:end);
            
            for i = 1:size(cdk2Drop,1)
                trace = cdk2Drop(i,:);
                traceShouldRise = trace(correctedFitId(i)+1:correctedFitId(i)+1+ceil(cdkShouldRiseAtLeast/obj.timeInterval));
                badCell(i) = sum(traceShouldRise);
            end
            badCell = find(badCell>0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            fitDif = fitDif(~badLogical);
            fitId = fitId(~badLogical);
            correctedFitId = correctedFitId(~badLogical);
            cellInc = cellInc(~badLogical);
            
            obj.cellState = cellInc;
            
            risingPoint = correctedFitId;
            risingPoint(~cellInc) = NaN;
            obj.risingPoint = risingPoint;
            
            % Visualization 1
            subplot(1,3,1)
            plot(smoothCdk2(obj.cellState,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','b');hold on
            title('rising')
            subplot(1,3,2)
            plot(smoothCdk2(~obj.cellState,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','r');
            title('quiescent')
            subplot(1,3,3)
            plot(smoothCdk2(risingPoint<2/obj.timeInterval,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','b');
            title('rising within 2 hours')
            
            %             % Visualization 2
            %             figure;
            %             p = panel;
            %             p.margin = [1 1 1 1];
            %             p.pack(3,5)
            %             dataInc = obj.dataAligned(obj.cellState,:,1);
            %             idInc = correctedFitId(obj.cellState)-1;
            %             frame = [1:20:size(dataInc,1)];
            
            %             panelm = 1;
            %             for i = 1:length(frame)-1
            %                 paneln = rem(i,5);paneln(paneln==0) = 5;
            %                 p(panelm,paneln).select();
            %                 plot(dataInc(frame(i):frame(i+1),5000:5000+ceil(30/obj.timeInterval),1)');hold on
            %                 temp = dataInc(frame(i):frame(i+1),5000+idInc(frame(i):frame(i+1))-1);
            %                 temp = diag(diag(temp));temp(temp==0)=NaN;   % just to make colors better
            %                 plot(idInc(frame(i):frame(i+1)),temp,'o','markersize',8);
            %                 plot(idInc(frame(i):frame(i+1)),temp,'o','markersize',4);
            %                 xlim([0,ceil(30/obj.timeInterval)]);
            %                 if paneln ==5
            %                     panelm = panelm+1;
            %                 end
            %                 if panelm == 4
            %                     figure;
            %                     p = panel;
            %                     p.margin = [1 1 1 1];
            %                     p.pack(3,5)
            %                     panelm = 1;
            %                 end
            %             end
            
            % Visualization 3
            %             figure;
            %             data = obj.dataAligned;
            %             risingPoint = obj.risingPoint;
            %             cdkdata = data(:,5000:5200,1);
            %             erkdata = data(:,5000:5200,2);
            %             p = panel;
            %             p.margin = [1 1 1 1];
            %             p.position
            %             p.pack(2,5)
            %             for i = 1:size(cdkdata,1)
            %                 pn = rem(i,5);pn(pn==0) = 5;
            %                 p(1,pn).select()
            %                 plot(cdkdata(i,:));
            %                 hold on;plot([risingPoint(i),risingPoint(i)],[0 10],'r');ylim([0,2.5]);xlim([0 200])
            %                 p(2,pn).select()
            %                 plot(erkdata(i,:));
            %                 hold on;plot([risingPoint(i),risingPoint(i)],[0 10],'r');ylim([1.1,2]);xlim([0 200])
            %                 if rem(i,5) ==0
            %                     waitforbuttonpress
            %                     p = panel;
            %                     p.margin = [1 1 1 1];
            %                     p.pack(2,5)
            %                 end
            %             end
        end
        
        function obj = judgecdk2_detect_rising_pts_2(obj,varargin) %for serum starvation
            % some codes are arbitrary
            % This script will fit a linear line to CDK2 from different time point.
            % It will determine CDK2Inc/CDK2Low based on the diffence between data and the best-fitted line.
            % It will give you a rising point. At the end you will get
            % several plots, and you may have to suspend it for the last one
            % It will do some cleaning as well
            channel=varargin{1};
            ptsToFit = 20/obj.timeInterval;
            slope = 0.04;
            lineToFit = [1:ptsToFit].*slope;    % This line will be fitted by changing starting point.
            startThres = 1.5;                   % start point must be lower than this
            riseThres = 1.5;                    % cdkInc cells must have a point above this
            errorIncLowThres = 15;              % Threshold of absolute error between data and the best-fitted line
            cdklowMaxThres = 1;                 % used for cleaning. If Cdk2Low, it shouldn't go above this
            cdkDropThres = -0.5;
            cdkShouldRiseAtLeast = 12;
            
            smoothCdk2 = tsmovavg(obj.dataSorted(:,:,channel),'s',ceil(1/obj.timeInterval),2);
            %             smoothCdk2 = obj.dataAlignedSmooth(:,:,1);
                     
            lastFrame=size(smoothCdk2,2);
            % fitting
            for j = 1:size(smoothCdk2,1)
                wholedata = smoothCdk2(j,1:end);% traces after first division
                %firstNan = find(isnan(wholedata),1,'last');
               fit = zeros(ptsToFit,1);
                % find the best fit and where it starts
                for starti = 1:lastFrame-ptsToFit
                    fit(starti) = sum(abs(wholedata(starti:starti+ptsToFit-1) - (lineToFit + wholedata(starti))));
                    if wholedata(starti)>startThres   % start point shouldn't be too large
                        fit(starti) = Inf;
                    end
                    if max(wholedata) < riseThres     % if Cdk2Inc it should go above this
                        fit(starti) = Inf;
                    end
                end
                [fitDif(j), fitId(j)] = min(fit);     % the error and the rising point of the best fit
                
                % Adjust the rising point
                % after finding the best fit, find local minima on the left. If there is, update.
                localMinima = find([0;0;diff(sign(diff(wholedata')))]>0)-1;
                localMinima = localMinima(localMinima <= fitId(j));
                [~, id] = min(abs(localMinima - fitId(j)));
                if isempty(localMinima)
                    correctedFitId(j) = fitId(j);
                else
                    correctedFitId(j) = localMinima(id);
                end
                
            end
            correctedFitId = correctedFitId-1;% relative to 5000. 5000+id is the absolute point
            fitId = fitId-1;
            
            % Determine Cdk2Inc and Cdk2low cells
            cellInc = fitDif<errorIncLowThres;
            
            % some cleaning. if it is Cdk2low, then maximum value after several frames
            % shouldn't be above 1 or so.
            maxAfterFirstDivision = max(obj.dataSorted(:,1+20/obj.timeInterval:end,1),[],2);
            withinCdklowMaxThres = maxAfterFirstDivision < cdklowMaxThres;
            badLogical = ~cellInc' & ~withinCdklowMaxThres;
            badCell = find(badLogical==1);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            fitDif = fitDif(~badLogical);
            fitId = fitId(~badLogical);
            correctedFitId = correctedFitId(~badLogical);
            cellInc = cellInc(~badLogical);
                   
            obj.cellState = cellInc;
            
            risingPoint = correctedFitId;
            risingPoint(~cellInc) = NaN;
            obj.risingPoint = risingPoint;
            
            % Visualization 1
            subplot(1,3,1)
            plot(smoothCdk2(obj.cellState,:,1)','b');hold on
            title('rising')
            subplot(1,3,2)
            plot(smoothCdk2(~obj.cellState,:,1)','r');
            title('quiescent')
            subplot(1,3,3)
            plot(smoothCdk2(risingPoint<10/obj.timeInterval,:,1)','b');
            title('rising within 10 hours')
            
                        % Visualization 2
                        figure;
                        p = panel;
                        p.margin = [1 1 1 1];
                        p.pack(3,5)
                        dataInc = obj.dataSorted(obj.cellState,:,channel);
                        idInc = correctedFitId(obj.cellState)-1;
                        frame = [1:20:size(dataInc,1)];
            
                        panelm = 1;
                        for i = 1:length(frame)-1
                            paneln = rem(i,5);paneln(paneln==0) = 5;
                            p(panelm,paneln).select();
                            plot(dataInc(frame(i):frame(i+1),1:1+ceil(30/obj.timeInterval),1)');hold on
                            temp = dataInc(frame(i):frame(i+1),1+idInc(frame(i):frame(i+1))-1);
                            temp = diag(diag(temp));temp(temp==0)=NaN;   % just to make colors better
                            plot(idInc(frame(i):frame(i+1)),temp,'o','markersize',8);
                            plot(idInc(frame(i):frame(i+1)),temp,'o','markersize',4);
                            xlim([0,ceil(30/obj.timeInterval)]);
                            if paneln ==5
                                panelm = panelm+1;
                            end
                            if panelm == 4
                                figure;
                                p = panel;
                                p.margin = [1 1 1 1];
                                p.pack(3,5)
                                panelm = 1;
                            end
                        end
            
            % Visualization 3
            %             figure;
            %             data = obj.dataAligned;
            %             risingPoint = obj.risingPoint;
            %             cdkdata = data(:,5000:5200,1);
            %             erkdata = data(:,5000:5200,2);
            %             p = panel;
            %             p.margin = [1 1 1 1];
            %             p.position
            %             p.pack(2,5)
            %             for i = 1:size(cdkdata,1)
            %                 pn = rem(i,5);pn(pn==0) = 5;
            %                 p(1,pn).select()
            %                 plot(cdkdata(i,:));
            %                 hold on;plot([risingPoint(i),risingPoint(i)],[0 10],'r');ylim([0,2.5]);xlim([0 200])
            %                 p(2,pn).select()
            %                 plot(erkdata(i,:));
            %                 hold on;plot([risingPoint(i),risingPoint(i)],[0 10],'r');ylim([1.1,2]);xlim([0 200])
            %                 if rem(i,5) ==0
            %                     waitforbuttonpress
            %                     p = panel;
            %                     p.margin = [1 1 1 1];
            %                     p.pack(2,5)
            %                 end
            %             end
        end
        
        %         function obj = judgecdk2_detect_rising_pts2(obj,varargin)
        %             slopewindow=10;
        %             for i=1:length(obj.dataAligned(:,:,1),1)
        %                 signal_total= obj.dataAligned(i,:,1);
        %                             firstframe(i)=samplestats(i,1);
        %                 lastframe(i)=find(~isnan(signal_total),1,'last')
        %                 %signal=signal(firstframe(i):lastframe(i));
        %                 signal_total_smooth(i,:)=smooth(signal_total,5);
        %                 signal=signal_total(5000:(i):lastframe(i));
        %                 signal_smooth=smooth(signal)'; %incoming signal always raw (unsmoothened)
        %                                 numframes=lastframe(i)-firstframe(i)+1;
        %                             sigstore(i,firstframe(i):lastframe(i))=signal_smooth;
        %                 %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 %     badtraces(i)=min(signal(1:slopewindow))>1; %will remove permanently high signals
        %
        %                 badtraces(i)=min(signal_total(5000:5000+slopewindow))>1; %will remove permanently high signals
        %
        %                  badCell = find(~secDivisionTooSoon==0);
        %             cumGoodTrace = cumsum(obj.goodTrace);
        %             for i = 1:length(badCell)
        %                 obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
        %             end
        %             end
        %
        %         end
        
        function obj = plot_compare_risingpoint(obj,varargin)
            individual=1;
            FASTINC = 3;         % in hour
            SLOWINC = [3,40];  % in hour, give a range of rising point. You can change this
            if numel(varargin)>0
                SLOWINC = varargin{1};
            end
            
            frameToPlot = 5000+[-20/obj.timeInterval:25/obj.timeInterval];
            fastIncTime = FASTINC/obj.timeInterval;
            slowIncTime = [SLOWINC(1)/obj.timeInterval, SLOWINC(2)/obj.timeInterval];
            %             fastCdk = obj.dataAlignedSmooth(obj.risingPoint<fastIncTime,:,:);
            %             slowCdk = obj.dataAlignedSmooth(obj.risingPoint>slowIncTime(1) & obj.risingPoint<slowIncTime(2),:,:);
            %             lowCdk = obj.dataAlignedSmooth(~obj.cellState,:,:);
            
            cdkInc=(obj.cellState' & obj.risingPoint<fastIncTime);
            cdkLow=(isnan(obj.risingPoint));% | (obj.risingPoint>slowIncTime(3) & obj.risingPoint<slowIncTime(4)));% ~obj.cellDelay); % isnan(obj.risingPoint)'
            cdkDelay=(obj.risingPoint>slowIncTime(1) & obj.risingPoint<slowIncTime(2));% | (~obj.cellState' & ~obj.cellLow'); %obj.cellDelay' & ~obj.cellLow' &
            
            fastCdk = obj.dataAlignedSmooth(cdkInc,:,:); %cellDelay
            slowCdk = obj.dataAlignedSmooth(cdkDelay,:,:);
            lowCdk = obj.dataAlignedSmooth(cdkLow,:,:);
            dataByStateSmooth{1} = lowCdk;
            dataByStateSmooth{2} = fastCdk;
            dataByStateSmooth{3} = slowCdk;
            
            %             fastCdk = obj.dataAligned(obj.risingPoint<fastIncTime,:,:);
            %             slowIncTime = [SLOWINC(1)/obj.timeInterval, SLOWINC(2)/obj.timeInterval];
            %             slowCdk = obj.dataAligned(obj.risingPoint>slowIncTime(1) & obj.risingPoint<slowIncTime(2),:,:);
            %             lowCdk = obj.dataAligned(~obj.cellState,:,:);
            
            fastCdk = obj.dataAligned(cdkInc,:,:);
            slowCdk = obj.dataAligned(cdkDelay,:,:);
            lowCdk = obj.dataAligned(cdkLow,:,:);
            dataByState{1} = lowCdk;
            dataByState{2} = fastCdk;
            dataByState{3} = slowCdk;
            
            obj.cellState=cdkInc;
            obj.cellDelay=cdkDelay;
            obj.cellLow=cdkLow;
            
            p = panel;p.pack(2,3);p.margin = [10,10,10,10];
            time = (frameToPlot-5000).*obj.timeInterval;
            colName = {'cdkLow','cdkHigh','cdkDelay'};
            color={[1 0 0], [0 0 1], [0 1 0]};
            for i = 1:3
                y1 = dataByStateSmooth{i}(:,frameToPlot,1);
%                 y2 = dataByState{i}(:,frameToPlot,2);
                if individual
                    p(1,i).select()
                    title(colName{i});
                    plot(time,y1,'Color',color{i}); %hold on;
                    xlim([time(1),time(end)]); ylim([0 3]);
                    h1(i) = gca;
                    p(2,i).select()
                    
%                     numberOfWell=23;
%                     Mean=nanmean(y2);
%                     SEM=nanstd(y2)/sqrt(numberOfWell);
%                     shadedErrorBar(time,Mean,SEM,{'Color',color{i}}); %hold on;
%                     plot(time,y2,'Color',color{i}); hold on;
%                     xlim([time(1),time(end)]); ylim([1.35 2]);
%                     h2(i) = gca;
                else
                    p(1,i).select()
                    title(colName{i});
                    shadedErrorBar(time,y1,{@nanmean,@nanstd},{'r-','markerfacecolor','r'});
                    xlim([time(1),time(end)]); ylim(obj.colorRange{1});
                    h1(i) = gca;
                    p(2,i).select()
                    shadedErrorBar(time,y2,{@nanmean,@nanstd},{'r-','markerfacecolor','r'});
                    xlim([time(1),time(end)]); ylim(obj.colorRange{2});
                    h2(i) = gca;
                end
            end
            linkaxes(h1,'xy');%linkaxes(h2,'xy');
            %h3 = [h1,h2];
            %linkaxes(h3,'x')
        end
        
        function obj = judgeCDK2Inc(obj,varargin)
            tempGoodTrace = ~isnan(obj.dataAligned(:,5000,1));
            tempGoodTrace = tempGoodTrace & ~isnan(obj.dataAligned(:,5000+ceil(16/obj.timeInterval),1));
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            
            smoothCdk2 = tsmovavg(obj.dataAligned(:,:,1),'s',ceil(6/obj.timeInterval),2);
            frame6to16 = 5000+ceil(6/obj.timeInterval):5000+ceil(16/obj.timeInterval);
            for j = 1:size(smoothCdk2,1)
                p = polyfit(1:length(frame6to16),smoothCdk2(j,frame6to16),1);
                slope(j,:) = p(1);
            end
            
            thres = 0.005;
            if length(varargin)~=0;thres = varargin{1};end
            cellInc = slope>thres;
            obj.cellState = cellInc;
            %             plot(obj.dataAligned(cellInc,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','b');hold on
            %             plot(obj.dataAligned(~cellInc,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','r');
            plot(smoothCdk2(cellInc,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','r');hold on
            plot(smoothCdk2(~cellInc,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','b');
        end
        
        function obj = judgeCDK2IncByTimeAfterMitosis(obj,varargin)
            figure;
            hoursAfterMitosisToMeasure=varargin{1};
            gating=varargin{3};
            channel=varargin{4};
            cutofftime=4;
            tempGoodTrace = ~isnan(obj.dataAligned(:,5000,channel));
            tempGoodTrace = tempGoodTrace & ~isnan(obj.dataAligned(:,5000+ceil(cutofftime/obj.timeInterval),channel));
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            
            smoothCdk2 = tsmovavg(obj.dataAligned(:,:,channel),'s',ceil(6/obj.timeInterval),2);
            judgeFrame = 5000+ceil(hoursAfterMitosisToMeasure/obj.timeInterval);
            judgeFrame_1 = 5000+ceil((hoursAfterMitosisToMeasure+2)/obj.timeInterval);
            judgeFrame_2 = 5000+ceil((hoursAfterMitosisToMeasure+5)/obj.timeInterval);
            judgeFrame_3 = 5000+ceil((hoursAfterMitosisToMeasure+2)/obj.timeInterval);
            lastframe=5000+ceil(cutofftime/obj.timeInterval);
            %             for j = 1:(obj.dataAligned,1)
            %                 if isnan(obj.dataAligned(j,judgeFrame,1))
            %                     cellInc(idx(j),1) = obj.dataAligned(idx(j),lastframe,1)>1;
            %                 else
            %                     cellInc(j,1) = obj.dataAligned(j,judgeFrame,1)>1;
            %                 end
            %             end
            %             cellInc(:,1) = obj.dataAligned(:,judgeFrame,1)>1;
            if varargin{2}
                switch gating
                    case 'G1'
                        cellInc(:,1) = obj.dataAlignedSmooth(:,judgeFrame,channel)>1 & obj.dataAlignedSmooth(:,judgeFrame_1,channel)>1.4;% & obj.dataAlignedSmooth(:,judgeFrame_2,1)>1.3;
                    case 'G2'
%                         cellInc(:,1) = obj.dataAlignedSmooth(:,judgeFrame,1)>1 & obj.dataAlignedSmooth(:,judgeFrame_1,1)>0.6;
                end
                idx = isnan(obj.dataAlignedSmooth(:,judgeFrame,channel));
                cellInc(idx,1) = obj.dataAlignedSmooth(idx,lastframe,channel)>0.8;
                temp_Delay(:,1) = obj.dataAlignedSmooth(:,judgeFrame_2,channel)>1 | obj.dataAlignedSmooth(:,judgeFrame_3,channel)>0.8;
                cellDelay(:,1)=(temp_Delay==1 & cellInc==0);
                cellLow=(temp_Delay==0 & cellInc==0);
                obj.cellState = cellInc;
                obj.cellLow = cellLow;
                obj.cellDelay = cellDelay;
                plot(obj.dataAlignedSmooth(cellInc,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),channel)','b');hold on
                plot(obj.dataAlignedSmooth(cellLow,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),channel)','r');
                plot(obj.dataAlignedSmooth(cellDelay,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),channel)','k');
            else
                cellInc(:,1) = obj.dataAligned(:,judgeFrame,channel)>1;
                idx = isnan(obj.dataAligned(:,judgeFrame,channel));
                cellInc(idx,1) = obj.dataAligned(idx,lastframe,channel)>0.7;
                temp_Delay(:,1) = obj.dataAligned(:,judgeFrame_2,channel)>1 | obj.dataAligned(:,judgeFrame_3,channel)>0.6;
                cellDelay(:,1)=(temp_Delay==1 & cellInc==0);
                cellLow=(temp_Delay==0 & cellInc==0);
                obj.cellState = cellInc;
                obj.cellLow = cellLow;
                obj.cellDelay = cellDelay;
                plot(obj.dataAligned(cellInc,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),channel)','b');hold on
                plot(obj.dataAligned(cellLow,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),channel)','r');
                plot(obj.dataAligned(cellDelay,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),channel)','k');
            end
        end
        
        function obj = judgeCDK2IncByTimeAfterMitosis_2(obj,varargin) % serum starvation
            figure;
            hoursAfterMitosisToMeasure=varargin{1};
            gating=varargin{3};
            cutofftime=4;
            tempGoodTrace = ~isnan(obj.dataAligned(:,5000,1));
            tempGoodTrace = tempGoodTrace & ~isnan(obj.dataAligned(:,5000+ceil(cutofftime/obj.timeInterval),1));
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            
            smoothCdk2 = tsmovavg(obj.dataAligned(:,:,1),'s',ceil(6/obj.timeInterval),2);
            judgeFrame = 5000+ceil(hoursAfterMitosisToMeasure/obj.timeInterval);
            judgeFrame_1 = 5000+ceil((hoursAfterMitosisToMeasure+2)/obj.timeInterval);
            judgeFrame_2 = 5000+ceil((hoursAfterMitosisToMeasure+5)/obj.timeInterval);
            judgeFrame_3 = 5000+ceil((hoursAfterMitosisToMeasure+2)/obj.timeInterval);
            lastframe=5000+ceil(cutofftime/obj.timeInterval);
            %             for j = 1:(obj.dataAligned,1)
            %                 if isnan(obj.dataAligned(j,judgeFrame,1))
            %                     cellInc(idx(j),1) = obj.dataAligned(idx(j),lastframe,1)>1;
            %                 else
            %                     cellInc(j,1) = obj.dataAligned(j,judgeFrame,1)>1;
            %                 end
            %             end
            %             cellInc(:,1) = obj.dataAligned(:,judgeFrame,1)>1;
            if varargin{2}
                switch gating
                    case 'G1'
                        cellInc(:,1) = obj.dataAlignedSmooth(:,judgeFrame,1)>1 & obj.dataAlignedSmooth(:,judgeFrame_1,1)>1.4;% & obj.dataAlignedSmooth(:,judgeFrame_2,1)>1.3;
                    case 'G2'
%                         cellInc(:,1) = obj.dataAlignedSmooth(:,judgeFrame,1)>1 & obj.dataAlignedSmooth(:,judgeFrame_1,1)>0.6;
                end
                idx = isnan(obj.dataAlignedSmooth(:,judgeFrame,1));
                cellInc(idx,1) = obj.dataAlignedSmooth(idx,lastframe,1)>0.8;
                temp_Delay(:,1) = obj.dataAlignedSmooth(:,judgeFrame_2,1)>1 | obj.dataAlignedSmooth(:,judgeFrame_3,1)>0.8;
                cellDelay(:,1)=(temp_Delay==1 & cellInc==0);
                cellLow=(temp_Delay==0 & cellInc==0);
                obj.cellState = cellInc;
                obj.cellLow = cellLow;
                obj.cellDelay = cellDelay;
                plot(obj.dataAlignedSmooth(cellInc,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','b');hold on
                plot(obj.dataAlignedSmooth(cellLow,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','r');
                plot(obj.dataAlignedSmooth(cellDelay,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','k');
            else
                cellInc(:,1) = obj.dataAligned(:,judgeFrame,1)>1;
                idx = isnan(obj.dataAligned(:,judgeFrame,1));
                cellInc(idx,1) = obj.dataAligned(idx,lastframe,1)>0.7;
                temp_Delay(:,1) = obj.dataAligned(:,judgeFrame_2,1)>1 | obj.dataAligned(:,judgeFrame_3,1)>0.6;
                cellDelay(:,1)=(temp_Delay==1 & cellInc==0);
                cellLow=(temp_Delay==0 & cellInc==0);
                obj.cellState = cellInc;
                obj.cellLow = cellLow;
                obj.cellDelay = cellDelay;
                plot(obj.dataAligned(cellInc,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','b');hold on
                plot(obj.dataAligned(cellLow,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','r');
                plot(obj.dataAligned(cellDelay,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','k');
            end
        end
        
        function obj = judgeCDK2Inc_meyerlab(obj,varargin)
            
            minlengthdaughter = 15;
            earlytime = 14;
            
            daughterDivFrame = obj.divFrameAfterAlignment(:,2);
            daughterDivFrame(isnan(daughterDivFrame))= size(obj.dataAligned,2);
            tempGoodTrace = daughterDivFrame > 5000+earlytime;
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            daughterDivFrame = daughterDivFrame(tempGoodTrace);
            
            tracesCdk2 = obj.dataAligned(:,:,1);
            numtraces=size(tracesCdk2,1);
            earlyval=ones(numtraces,1)*NaN;
            lateval=ones(numtraces,1)*NaN;
            maxval=ones(numtraces,1)*NaN;
            minval=ones(numtraces,1)*NaN;
            
            for i=1:numtraces
                earlyval(i)=tracesCdk2(i,5000+earlytime);
                lateval(i)=tracesCdk2(i,5000+minlengthdaughter-1);
                maxval(i)=nanmax(tracesCdk2(i,5000+earlytime:daughterDivFrame(i))); %max from minlength-end
                minval(i)=nanmin(tracesCdk2(i,5000+earlytime:daughterDivFrame(i)));
            end
            earlycutoff=0.6; %default 0.55
            latecutoff=earlycutoff+0.01*(minlengthdaughter-earlytime);
            Cdk2inc=earlyval>earlycutoff & lateval>latecutoff & minval>earlycutoff;
            Cdk2low=earlyval<earlycutoff & maxval<earlycutoff;  % THIS IS NOT EXCLUSIVE
            
            % throw away something not reliable
            tempGoodTrace = (Cdk2inc|Cdk2low);
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            Cdk2inc = Cdk2inc(tempGoodTrace);
            Cdk2low = Cdk2low(tempGoodTrace);
            
            obj.cellState = Cdk2inc;
            %%
            cellInc = Cdk2inc;
            smoothCdk2 = tsmovavg(obj.dataAligned(:,:,1),'s',ceil(6/obj.timeInterval),2);
            
            plot(smoothCdk2(cellInc,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','b');hold on
            plot(smoothCdk2(~cellInc,5000-ceil(24/obj.timeInterval):5000+ceil(24/obj.timeInterval),1)','r');
        end
        
        %% data cleaning
        function obj = clean_manually_remove_data(obj,varargin)
            channel = 1;
            figure;
            if numel(varargin)>0; channel = varargin{1}; align = varargin{2};end
            frameToPlot = 5000-ceil(16/obj.timeInterval):5000+ceil(16/obj.timeInterval);
            if align
                dataToPlot = obj.dataAligned(:,frameToPlot,channel);
            else
                dataToPlot = obj.dataSorted(:,:,channel);
            end
            plot(dataToPlot');
            [x,y] = getpts(gcf);
            xD = repmat(1:size(dataToPlot,2),size(dataToPlot,1),1);
            rowStore = [];
            for i = 1:length(x)
                dist = abs(xD - x(i)) + abs(dataToPlot-y(i));
                [row,col] = find(dist==nanmin(dist(:)));
                rowStore = [rowStore;row];
            end
            tempGoodTrace = ones(size(dataToPlot,1),1);
            tempGoodTrace(rowStore) = 0;
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
            if align
                dataToPlot = obj.dataAligned(:,frameToPlot,channel);
            else
                dataToPlot = obj.dataSorted(:,:,channel);
            end
            plot(dataToPlot');
        end
        
        function obj = clean_remove_suddendropCDK2(obj)
            %remove dirty Cdk2 from storeStoreData
            dataSorted = obj.dataSorted;
            tempGoodTrace = sum(dataSorted(:,:,1)>6,2)==0 & sum(dataSorted(:,:,1)<-0.5,2)==0; % ratio should be -0.5-6
            badTrace = max(diff(tsmovavg(dataSorted(:,:,1),'s',5,2),[],2)>1,[],2);% delete the sudden change. it was 2.
            tempGoodTrace = tempGoodTrace&~badTrace;
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        function obj = clean_leastnum_frames(obj,varargin)
            dataSorted = obj.dataSorted;
            tempGoodTrace = sum(~isnan(dataSorted(:,:,1)),2)>=varargin{1};
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        function obj = clean_leastnum_frames_2(obj,varargin)
            dataAligned = obj.dataAligned;
            tempGoodTrace = sum(~isnan(dataAligned(:,1:5000,1)),2)>=varargin{1};
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        
        function obj = clean_leastnum_frames_3(obj,varargin)
            dataAligned = obj.dataAligned;
            tempGoodTrace = sum(~isnan(dataAligned(:,5001:10000,1)),2)>=varargin{1};
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        function obj = clean_leastnum_frames_4(obj,varargin)
            dataAligned = obj.dataAligned;
            tempGoodTrace = sum(~isnan(dataAligned(:,5001:10000,1)),2)<=varargin{1};
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
         function obj = clean_leastnum_frames_5(obj,varargin)
            dataSorted = obj.dataSorted;
            tempGoodTrace = sum(~isnan(dataSorted(:,varargin{1}:varargin{2},varargin{4})),2)>=varargin{3};
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        
        function obj = clean_remove_short_after_mitosis(obj,varargin)
            dataAligned = obj.dataAligned;
            tempGoodTrace = sum(~isnan(dataAligned(:,5000:5200,1)),2)>=varargin{1};
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        function obj = clean_remove_short_before_mitosis(obj,varargin)
            dataAligned = obj.dataAligned;
            tempGoodTrace = sum(~isnan(dataAligned(:,5000-200:5000,1)),2)>=varargin{1};
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        function obj = clean_remove_mean_outside_n_std(obj,varargin)
            % varargin{1}: which channel to use
            % varargin{2}: n. remove outside the mean+-n*std
            dataSorted = obj.dataSorted;
            nStd = varargin{2};
            traces = dataSorted(:,:,varargin{1});
            meanTrace = nanmean(traces,2);
            dataMean = nanmean(meanTrace);
            dataStd = nanstd(meanTrace);
            
            tempGoodTrace = meanTrace<dataMean+nStd*dataStd & meanTrace>dataMean-nStd*dataStd;
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        function obj = clean_remove_CDK2_not_dropped(obj,varargin)
            thres = 0;
            if numel(varargin)>0;
                thres = varargin{1};
            end
            cdk2Drop = obj.dataAligned(:,5000-ceil(2/obj.timeInterval):5000-ceil(0.5/obj.timeInterval),1);
            for j = 1:size(cdk2Drop,1)
                p = polyfit(1:size(cdk2Drop,2),cdk2Drop(j,:),1);
                slope(j,:) = p(1);
            end
            tempGoodTrace = slope<thres;
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        function obj = clean_remove_CDK2_above_thres_at_mitosis(obj,varargin)
            thres = 0;
            if numel(varargin)>0;
                thres = varargin{1};
            end
            cdk2atMitosis = obj.dataAligned(:,5000,1);
            
            tempGoodTrace = cdk2atMitosis < thres;
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        function obj = find_NO_divnum_frame(obj,varargin)
            for i=1:length(obj.rawDivFrameSorted);
                if isnan(obj.rawDivFrameSorted(i,:))
                    tempGoodTrace(i,1)=1;
                end
            end
            %             tempGoodTrace = rawDivFrame>=varargin{1} & rawDivFrame<=varargin{2};
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        function obj = clean_divnum_frame(obj,varargin)
            for i=1:length(obj.rawDivFrameSorted);
                if sum(~isnan(obj.rawDivFrameSorted(i,:)))==1
                    if obj.rawDivFrameSorted(i,1)>=varargin{1} & obj.rawDivFrameSorted(i,1)<=varargin{2};
                        tempGoodTrace(i,1)=1;
                    else
                        tempGoodTrace(i,1)=0;
                    end
                end
            end
            %             tempGoodTrace = rawDivFrame>=varargin{1} & rawDivFrame<=varargin{2};
            badCell = find(tempGoodTrace==0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
            end
        end
        
        function clean_remove_dataaligned_indicate_timeandrange(obj,varargin)
            % first argument is time point, second argument is range to
            % remove e.g. function(4950, [0, 1.1]) will remove traces in
            % this range at 4950th. First check dataAligned looks like by
            % passing no argument
            channel = 1;
            if numel(varargin)>0
                channel = varargin{1};
            end
            if numel(varargin)<2
                figure;
                if obj.nodivision
                    plot(1:300,obj.dataAligned(:,1:300,channel)')
                else
                    plot(4800:5200,obj.dataAligned(:,4800:5200,channel)')
                end
            end
            if numel(varargin)==3
                timePoint = varargin{2};
                rangeToRemove = varargin{3};
                data = obj.dataAligned(:,timePoint,channel);
                removeTrace = data>rangeToRemove(1) & data<rangeToRemove(2);
                badCell = find(removeTrace==1);
                cumGoodTrace = cumsum(obj.goodTrace);
                for i = 1:length(badCell)
                    obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
                end
            end
        end
        
        function clean_remove_dataaligned_indicate_timeandrange_2(obj,varargin)
            % first argument is time point, second argument is range to
            % remove e.g. function(4950, [0, 1.1]) will remove traces in
            % this range at 4950th. First check dataAligned looks like by
            % passing no argument
            channel = 1;
            if numel(varargin)>0
                channel = varargin{1};
            end
            if numel(varargin)<2
                figure;
                if obj.nodivision
                    plot(1:300,obj.dataAligned(:,1:300,channel)')
                else
                    plot(4800:5200,obj.dataAligned(:,4800:5200,channel)')
                end
            end
            if numel(varargin)==3
                timePoint = varargin{2};
                rangeToRemove = varargin{3};
                data = obj.dataAligned(:,timePoint,channel);
                removeTrace = data==rangeToRemove;
                badCell = find(removeTrace==1);
                cumGoodTrace = cumsum(obj.goodTrace);
                for i = 1:length(badCell)
                    obj.goodTrace(find(cumGoodTrace==badCell(i))) = 0;
                end
            end
        end
        %         function obj = clean_manual_test(obj,varargin)
        %             channel = 1;
        %             if numel(varargin)>0;channel = varargin{1};end
        %             dataToPlot = obj.dataSorted(:,:,channel);
        %             figObj = figure('ButtonDownFcn',@(obj,evt) 0);
        %             plot(dataToPlot');
        %             fig_obj = datacursormode(gcf);
        %             set(fig_obj,'DisplayStyle','datatip','SnapToDataVertex','on','Enable','on')
        %             waitfor(gcf);
        %             c_info = getCursorInfo(fig_obj);
        %             % Make selected line wider
        %             set(c_info.Target,'LineWidth',5)
        %             set(gca,'buttondownFcn',@figObj.test)
        %         end
        % %         function test(obj,srchand,eventdata)
        % %             5
        % %         end
        % %
        function obj = clean_remove_by_lower_upper_limit(obj,varargin)
            % first argument is channel, second argument is lower and upper
            % limit (e.g. [0.05 Inf])
            
            channel = varargin{1};
            lowerLimit = varargin{2}(1);
            upperLimit = varargin{2}(2);
            
            data = obj.dataSorted(:,:,channel);
            badData = data<lowerLimit | data>upperLimit;
            badCell = find(nansum(badData,2)>0);
            cumGoodTrace = cumsum(obj.goodTrace);
            for i = 1:length(badCell)
                obj.goodTrace(cumGoodTrace==badCell(i)) = 0;
            end
        end
    end
end

