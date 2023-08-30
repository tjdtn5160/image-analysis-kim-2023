classdef calculate_next_cell_cycle_entry_prob < sort_or_align_by_mitosis
    
    properties
        inhibitionTime % in hours
        timeWindow = 3;
        
        divideCellsNumStore
        nondivideCellsNumStore
        startTimeStore
        timeAndProbability
    end
    
    properties (Dependent=true)
        inhibitionFrame
    end
    
    methods
        function obj = calculate_next_cell_cycle_entry_prob(rawData,rawDivFrame)
            obj@sort_or_align_by_mitosis(rawData,rawDivFrame);
        end
        
        function obj = set.inhibitionTime(obj,varargin)
             obj.inhibitionTime = varargin{1};
        end
        function inhibitionFrame = get.inhibitionFrame(obj)
             inhibitionFrame = round(obj.inhibitionTime/obj.timeInterval);
        end
        
        function calc_prob(obj,varargin)
            % to calculate probability, you first need to judge CDK2 to get
            % cellState
            maximumTime = (size(obj.dataSorted,2)-1)*obj.timeInterval;
            
            firstDivFrame = obj.divFrameSorted(:,1);
            firstDivTime = (firstDivFrame-1)*obj.timeInterval;
            divTimeRelativeToInhibition = firstDivTime - obj.inhibitionTime;
            
            % here is just to determine where to start the calculation
            start = 0;int = 0;
            while start==0
                start = (0 - int*obj.timeWindow)<=min(divTimeRelativeToInhibition);
                int = int+1;
            end
            startTime = 0 - (int-1)*obj.timeWindow;
            
            
            divideCellsNumStore =[];
            nondivideCellsNumStore = [];
            startTimeStore = []; % when to start window
            while startTime < maximumTime - 20 % ad hoc. assume we judge CDK2inc by 20h after first mitosis
                    cellsInWindow = divTimeRelativeToInhibition>=startTime & divTimeRelativeToInhibition < startTime + obj.timeWindow;
                    divideCellsInWindow  = obj.cellState(find(cellsInWindow));
                    divideCellsNum = sum(divideCellsInWindow);
                    nondivideCellsNum = sum(~divideCellsInWindow);
                    
                    divideCellsNumStore = [divideCellsNumStore,divideCellsNum];
                    nondivideCellsNumStore = [nondivideCellsNumStore,nondivideCellsNum];
                    startTimeStore = [startTimeStore,startTime];
                    startTime = startTime + obj.timeWindow;
            end
            obj.divideCellsNumStore = divideCellsNumStore;
            obj.nondivideCellsNumStore = nondivideCellsNumStore;
            obj.startTimeStore = startTimeStore;

        end

        function plot_total_cell_divide(obj,varargin)
            startTime = obj.startTimeStore;
            endTime = obj.startTimeStore + obj.timeWindow;
            barh([obj.divideCellsNumStore',obj.nondivideCellsNumStore'],'stacked');
%             prob = obj.divideCellsNumStore./(obj.divideCellsNumStore+obj.nondivideCellsNumStore)
            % make labels
            for i = 1:length(startTime)
%               'MEK inhibition started'
              str{i} = sprintf('%dh:%dh',startTime(i),endTime(i))
            end
            set(gca,'YTickLabel',str,'YTick',1:length(startTime));
            xlabel('cell number');
            legend('cycling','quiescent')
            set(gca,'View',[0 -90])
            ylabel('time MEK inh started relative to 1st mitosis (h)');
        end

        function plot_probability_cell_divide(obj,varargin)
            startTime = obj.startTimeStore;
            endTime = obj.startTimeStore + obj.timeWindow;
            prob = obj.divideCellsNumStore./(obj.divideCellsNumStore+obj.nondivideCellsNumStore);
            barh(prob);
%             prob = obj.divideCellsNumStore./(obj.divideCellsNumStore+obj.nondivideCellsNumStore)
            % make labels
            for i = 1:length(startTime)
%               'MEK inhibition started'
              str{i} = sprintf('%dh:%dh',startTime(i),endTime(i))
            end
            set(gca,'YTickLabel',str,'YTick',1:length(startTime));
            xlabel('probability to enter next cell cycle');
            set(gca,'View',[0 -90])
            ylabel('time MEK inh started relative to 1st mitosis (h)');
        end
        
        function obj = plot_line_probability_cell_divide(obj,varargin)
            figure;
            condition=varargin{1};
            resultdir=varargin{2};
            startTime = obj.startTimeStore;
            endTime = obj.startTimeStore + obj.timeWindow;
            prob = obj.divideCellsNumStore./(obj.divideCellsNumStore+obj.nondivideCellsNumStore);
            
            x = -fliplr(startTime);
            y = fliplr(1-prob);
            
            plot(x,y,'-o')
%             xlim([-13,7]);
            ylim([0 1])
            xlabel('Time of treatment -time of mitosis (hr)');
            ylabel('Probability of quiescence');
            timeAndProbability = [x;y];
            obj.timeAndProbability = timeAndProbability;
%             save([resultdir condition '.mat'],'y');
            end
        
     end
         
end

