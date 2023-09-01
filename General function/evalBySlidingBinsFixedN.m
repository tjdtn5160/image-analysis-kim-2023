function [xOut,yOut]=evalBySlidingBinsFixedN(xIn,yIn,n,func,step)
%evalBySlidingBins: runs the designated function on bins of your data. 
%First assign a vector of your bin sizes, and a value for bin width. Inputs are x values, y values, bins, bin width, function to run. 
%binsbrdu=0:0.2:20; binwidthnbrdu=0.4;
%brduavg_off=evalBySlidingBins(time_since_APC_off,antibody_log2_new,binsbrdu,binwidthnbrdu,'nanmedian');

%Evaluates the function func on a sliding bin of width n (n being a number
%of data points) over the paired data (xIn,yIn)

if nargin<5
    step=1;
end
temp=[xIn(:) yIn(:)];
temp(isnan(xIn(:)+yIn(:)),:)=[];
temp=sortrows(temp,1);
xIn=temp(:,1); yIn=temp(:,2);
for i=1:((length(xIn)-n)/step +1);
    start=((i-1)*step+1);
    range=start:(start+n-1);
    xOut(i)=median(xIn(range));
    yOut(i)=feval(func,yIn(range));
end

