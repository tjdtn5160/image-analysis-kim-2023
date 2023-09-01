function res = evalBySlidingBins(xIn,yIn,bins,binwidth,func)
%evalBySlidingBins: runs the designated function on bins of your data. 
%First assign a vector of your bin sizes, and a value for bin width. Inputs are x values, y values, bins, bin width, function to run. 
%binsbrdu=0:0.2:20; binwidthnbrdu=0.4;
%brduavg_off=evalBySlidingBins(time_since_APC_off,antibody_log2_new,binsbrdu,binwidthnbrdu,'nanmedian');

%syntax: res = evalBySlidingBins(xIn,yIn,bins,func)  Evaluates the function
%specified by func on the data (xIn,yIn) by binning the x values. It bins
%the x values using sliding bins centered at the locations specified by
%bins and with width binwidth. The function func should take a list and
%return a value.

x2=xIn(:); y2=yIn(:);
tempInd= isnan(x2) | isnan(y2);
x3=x2(~tempInd);
y3=y2(~tempInd);
mat=[x3 y3];
smat=sortrows(mat,1);
b1=1; b2=1;
for i=1:length(bins)
    b1=findmin(smat(:,1),bins(i),binwidth,b1);
    b2=findmax(smat(:,1),bins(i),binwidth,b2);
    list=smat(b1:b2,2);
    res(i)=feval(func,list);
end

end

%*************************************************************************

function ind=findmin(list,middle,binsize,prev)
target=middle-binsize/2;
ind=prev;
while (list(ind)<target)&(ind<length(list))
    ind=ind+1;
end
end

%*************************************************************************

function ind=findmax(list,middle,binsize,prev)
target=middle+binsize/2;
ind=prev;
while (list(ind)<target)&(ind<length(list))
    ind=ind+1;
end
end

