function [ trace1_smooth ] = nansmooth(trace1,smooth_factor)
%Smooths single vector, ignoring NaNs
if nargin<2
    smooth_factor=5;
end
trace1_smooth=trace1;
nan_index=find(~isnan(trace1));
smooth_trace_temp=smooth(trace1(nan_index),smooth_factor);
trace1_smooth(nan_index)=smooth_trace_temp;
end