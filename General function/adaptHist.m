function [out] = adaptHist(data,xLabel,nbins)
if nargin<2
    xLabel='';
end

if nargin<3
    nbins=50; %default for number of bins
end

figure;
if sum(mod(data,1))==0 && length(unique(data)) < 500
    histogram(data,'Normalization','probability',...
        'BinMethod','integers'); %binmethod sets binning.
    out=1;
else
    histogram(data,nbins,'Normalization','probability');
    out=2;
end
xlabel(xLabel);
ylabel('Fraction');
end

