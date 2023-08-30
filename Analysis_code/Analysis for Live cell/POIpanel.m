function POIpanel(traceIDs,firstframe,lastframe,sigstore,badtraces,ylims,varargin)
numvar=length(varargin); %number of different POIs to mark
for i=1:numel(traceIDs)
    figure(ceil(i/24)); set(gcf,'color','w');
    subaxis(4,6,mod(i-1,24)+1,'ML',0.05,'MR',0.02,'MT',0.03,'MB',0.05,'SH',0.03);
    frames=firstframe(traceIDs(i)):lastframe(traceIDs(i));
    plot(1:length(frames),sigstore(traceIDs(i),frames));
    ylim(ylims);
    hold on;
    if badtraces(traceIDs(i))==1
        midx=round(length(frames)/2);
        plot(midx,0.5,'rx','markersize',24);
        continue;
    end
    for vc=1:numvar
        plotPOI(varargin{vc},i,traceIDs,sigstore,frames);
    end
%     if nargin==8
%         plot(1:length(frames),altstore(traceIDs(i),frames),'k--');
%     end
end
end

function plotPOI(POI,IDnumber,traceIDs,sigstore,frames)
    if ~isnan(POI(traceIDs(IDnumber)))
        plot(POI(traceIDs(IDnumber)),sigstore(traceIDs(IDnumber),frames(POI(traceIDs(IDnumber)))),'go','markerfacecolor','g','markersize',6);
    end
end