classdef detect_erk_pulse < sort_or_align_by_mitosis_demo2
    
    properties
        ERKCHAN
        erkDataClean
        pulse
    end
    
    methods
        function obj = detect_erk_pulse(rawData,rawDivFrame,ERKCHAN)
            obj@sort_or_align_by_mitosis_demo2(rawData,rawDivFrame);
            obj.ERKCHAN = ERKCHAN;
        end
        
        function obj = erkclean_filtering_savitzky(obj,varargin)
            % press ENTER to renew figure
            data = obj.dataSorted(:,:,obj.ERKCHAN);
            dataToPlot = randi(size(data,1),size(data,1),1);
            polynomial = varargin{1};
            frameSize = varargin{2};            
            dataFir = sgolayfilt(data',polynomial,frameSize);
            dataFir = dataFir';
            
            obj.erkDataClean = dataFir;
            p = panel();
            p.pack(5,2);
            p.de.margin = 0.05;
            w = 1;
            for i = 1:length(dataToPlot)
                h(1) = p(w,1).select();
                plot(data(dataToPlot(i,:),:));
                h(2) = p(w,2).select();
                plot(dataFir(dataToPlot(i,:),:));      
                linkaxes(h,'xy')
                if rem(i,5)==0
                pause
                w = 0;
                p = panel();
                p.pack(5,2);
                p.de.margin = 0.05;
                end
                w = w+1;
            end
        end
        
        function erkclean_remove_trend(obj,varargin)
            
            
        end
        
        function obj = detect_by_absolute_threshold(obj,varargin)
            thres = prctile(obj.dataSorted(:,:,obj.ERKCHAN),75);
            if numel(varargin)>0;thres = varargin{1};end
            obj.pulse = obj.dataSorted(:,:,obj.ERKCHAN)>thres;
        end
        
        function obj = detect_by_fitting_sine(obj,varargin)
%             data = obj.dataSorted(:,:,obj.ERKCHAN);
%             for i = 1:size(data,1)
%                 ppx = fft(data(i,:));
%                 absppx = abs(ppx);
%                 ppx(absppx<2.5) = 0;
% %                 ppx(1:5,:) = 0;
%                 dataCleaned(i,:) = ifft(ppx);
%             end
%             
%             p = panel();
%             p.pack(5,2);
%             p.de.margin = 0.05;
%             w = 1;
%             for i = 1:size(dataCleaned,1)
%                 h(1) = p(w,1).select();
%                 plot(data(i,:));
%                 h(2) = p(w,2).select();
%                 plot(dataCleaned(i,:));
%                 linkaxes(h,'xy')
%                 if rem(i,5)==0
%                 pause
%                 w = 0;
%                 p = panel();
%                 p.pack(5,2);
%                 p.de.margin = 0.05;
%                 end
%                 w = w+1;
%             end
        end
        
     function obj = plot_erk_pulse(obj,varargin)
            % press ENTER to renew figure
            data = obj.dataSorted(:,:,obj.ERKCHAN);
            if ~isempty(obj.erkDataClean);
                data = obj.erkDataClean;
            end
            dataToPlot = randi(size(data,1),size(data,1),1);
            
            pulseNan = obj.pulse;
            pulseNan = single(pulseNan);
            pulseNan(~pulseNan)=NaN;
            
            p = panel();
            p.pack(5,1);
            p.de.margin = 0.25;
            w = 1;
            for i = 1:length(dataToPlot)
                h(1) = p(w,1).select();
                plot(obj.dataSorted(dataToPlot(i),:,obj.ERKCHAN));
%                 h(2) = p(w,1).select();
                hold on
                plot(obj.dataSorted(dataToPlot(i),:,obj.ERKCHAN).*pulseNan(dataToPlot(i,:),:),'or');    
                ylim(obj.colorRange{obj.ERKCHAN})
                linkaxes(h,'xy')
                if rem(i,5)==0
                pause
                w = 0;
                p = panel();
                p.pack(5,1);
                p.de.margin = 0.25;
                end
                w = w+1;
            end
     end
     
     function detect_erk_pulse_freq_domain(obj,varargin)
        Fs = 1/0.2;
        t = 0:1/Fs:39.8;
        % t = t(30:80+zeroPat);
        y1 = store(i,20:80);
        y1 = y1';         
        
        [pxx,f] = pwelch(y1,40,10,30,Fs);
        Y1(i,:) = pxx';
        % semilogy(f,abs(Y1));
        % loglog(f,abs(Y1(i,:)),'color',col(condi,:),'linewidth',0.05);hold on
        Y1(i,:) = smooth(Y1(i,:),3);
     end
        
     function obj = detect_findpeaks(obj,varargin)
         data = obj.dataSorted(:,:,obj.ERKCHAN);
         pulseNan = zeros(size(data));
         for i = 1:size(data,1)
             if numel(varargin)==0
                [pks,loc] = findpeaks(data(i,:));
             elseif numel(varargin)==1
                [pks,loc] = findpeaks(data(i,:),'MINPEAKHEIGHT',varargin{1});
             elseif numel(varargin)==2
                [pks,loc] = findpeaks(data(i,:),'MINPEAKHEIGHT',varargin{1},'MINPEAKDISTANCE',varargin{2});
             end
             pulseNan(i,loc) = 1;
         end
            dataToPlot = randi(size(data,1),size(data,1),1);
            
            obj.pulse = pulseNan;
            
            pulseNan = single(pulseNan);
            pulseNan(~pulseNan)=NaN;
            
%             p = panel();
%             p.pack(5,1);
%             p.de.margin = 0.25;
%             w = 1;
%             rawData = obj.dataSorted(:,:,obj.ERKCHAN);
%             for i = 1:length(dataToPlot)
%                 h(1) = p(w,1).select();
%                 plot(rawData(dataToPlot(i,:),:));
% %                 h(2) = p(w,1).select();
%                 hold on
%                 plot(rawData(dataToPlot(i,:),:).*pulseNan(dataToPlot(i,:),:),'or');    
%                 ylim(obj.colorRange{obj.ERKCHAN})
%                 linkaxes(h,'xy')
%                 if rem(i,5)==0
%                 pause
%                 w = 0;
%                 p = panel();
%                 p.pack(5,1);
%                 p.de.margin = 0.25;
%                 end
%                 w = w+1;
%             end
     end
     
        
    end
    
end

