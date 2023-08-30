function [risetime,badtraces]=getGemininfeatures_1(sampletraces,samplestats,risethresh)
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
risetime=ones(samplesize,1)*NaN;
risetimedb=risetime;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;
altstore=sigstore;
prebuffer=6;
for i=1:samplesize
    signal=sampletraces(i,:);
    sigstore(i,:)=signal;
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    numframes=lastframe(i)-firstframe(i)+1;
    %signal=signal(firstframe(i):lastframe(i));
    signal=smooth(signal(firstframe(i):lastframe(i)))'; %incoming signal always raw (unsmoothened)
    %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% general requirements
    postbuffer=10;
    futurebuffer=10;
    earlycutoff=20;
    gate=zeros(1,numframes);
    for j=1:numframes-futurebuffer
        presentheight=signal(j)<risethresh; %default 0.05
        %presentheight=max(signal(1:j))<risethresh;
        %minfutureheight=min(signal(j+futurebuffer:end))>signal(j+futurebuffer-1);
        minfutureheight=min(signal(j+futurebuffer:end))>signal(j);
        bufferfutureheight=signal(j+postbuffer)>signal(j)+2*risethresh; %50
        %bufferfutureheight=signal(j+postbuffer)>signal(j)+125;
        gate(j)=presentheight & minfutureheight & bufferfutureheight;
    end
    gate(1:earlycutoff)=0;
    %sig_revslope=10*getslope_reverse(signal,1:10);
    %sig_fwdslope=10*getslope_forward_avg(signal,1:5);
    sig_fwdslope=10*getslope_forward_avg(signal,1:postbuffer);
    sig_time=1:length(signal);
    filterscore=sig_fwdslope-signal+sig_time+200;
    filterscore=filterscore.*gate;
    altstore(i,firstframe(i):lastframe(i))=sig_fwdslope;
    %tempsearch=find(sig_fwdslope>0.05 & abs(signal)<0.03,1,'last');
    %tempsearch=find(filterscore>0.05,1,'first');
    filtermax=max(filterscore);
    tempsearch=find(filterscore==filtermax,1,'first');
    %if isempty(tempsearch) || signal(end)<0.05 %0.05
    if ~isempty(tempsearch) && filtermax>0
        risetimedb(i)=tempsearch;
        risetime(i)=tempsearch+firstframe(i)-1; %return absolute POI rather than relative to mitosis
    %elseif signal(end)<0.1
    elseif max(signal(prebuffer:end))<risethresh
        risetime(i)=NaN;
    else
        badtraces(i)=1;
    end
end
% keyboard;
%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
traceIDs=1:96;
ylims=[0 200];
POIpanel(traceIDs,firstframe,lastframe,sigstore,badtraces,ylims,risetimedb);
%}