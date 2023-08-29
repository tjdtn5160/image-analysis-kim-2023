function [risetime_cdk2 , cdk2_low , badtraces_cdk2] = getCdk2rise_yilin(traces , tracestats)
%% Initialization
n_traces = size(traces , 1);
risetime_cdk2 = NaN(n_traces , 1);
% risetime_cdk2_cusum = risetime_cdk2;
badtraces_cdk2 = zeros(n_traces , 1);
cdk2_low = zeros(n_traces , 1);

%% Parameters
% peak_lower = 0.35;
% trough = 0.8;
% peak_higher = 1.2;
binary_cutoff = 0.8;

% Minimum time points after rise
% min_after_rise = 5;

% cusum_climit = 4; % Default in cusum is 5

% Average number of frames per block
average_frames_per_block = 10;

% Threshold for detecting slope change
threshold = 0.02;

% Pooled change in slope
pooled_slope_change = [];
    

%% Loops through and finds rise points
for trace_id = 1 : n_traces
    %% Extraction
    sig = traces(trace_id , :);
    trace_start = tracestats(trace_id , 1);
    trace_end = tracestats(trace_id , 2);
    
    % Take only real signal (all subsequent time points are relative to start of trace)
    sig = sig(trace_start : trace_end);
    n_frames = numel(sig);
    sig_smooth = smooth(sig , 5);
    
    
    %% Gating
    % Traces with small change
%     if range(sig_smooth) < range_cutoff
%         cdk2_low(trace_id) = 1;
%         badtraces_cdk2(trace_id) = 1;
%         continue;
%     end
    

    % If trace starts high
    if (max(sig(1 : 5)) > binary_cutoff)
        badtraces_cdk2(trace_id) = 1;
        continue;
    end
    
    % Whether there is NaN in signal
    if sum(isnan(sig)) > 0
        badtraces_cdk2(trace_id) = 1;
        continue;
    end
    
    % If traces doesn't go high
    if (max(sig) < binary_cutoff)
        cdk2_low(trace_id) = 1;
        continue;
    end
    
%     disp(trace_id);
     

    %% Findchangepts
    

    % Find points of change
    pts = findchangepts(sig_smooth , 'Statistic' , 'linear' , 'maxnumchanges' , round(numel(sig) / average_frames_per_block));
    pts = [1 ; pts ; n_frames];
    
    % Smoothed signal value at points of change
    pts_sig_smooth = sig_smooth(pts);
    
    % Change in smoothed signal value between start and end of each block
    pts_sig_smooth_change = diff(pts_sig_smooth);
    
    % Number of frames in each block
    frames_per_block = diff(pts);
    
%     % Approximate slope of each block (not by linear fitting to save time)
%     slope_block = pts_sig_smooth_change ./ frames_per_block;
    
    % Slope of each block by fitting
    slope_block = NaN(numel(pts) - 1 , 1);
    for i = 1 : numel(slope_block)
        xx = 1 : (frames_per_block(i) + 1);
        yy = sig_smooth(pts(i) : pts(i + 1))';
        temp_fit = polyfit(xx , yy , 1);
        slope_block(i) = temp_fit(1);
    end
    
    slope_block(slope_block < -0.01) = -0.01;
    
    % Pts below cutoff, first stretch
    pts_low = (pts_sig_smooth < binary_cutoff);
    pts_low = pts_low(1 : find(pts_low == 0 , 1 , 'first'));   
    
    slope_block_low = slope_block(pts_low);
    change_slope = diff(slope_block_low);
    
    % Find maximum slope change
%     max_change_slope = find(change_slope == max(change_slope)) + find(pts_low , 1 , 'first');
%     rise_pts = pts(max_change_slope);
    
    % Find first change point with major slope change
    first_change_slope = find(change_slope > threshold , 1 , 'first') + find(pts_low , 1 , 'first');
    rise_pts = pts(first_change_slope);
    
    pooled_slope_change = [pooled_slope_change ; change_slope];

    try
        risetime_cdk2(trace_id) = rise_pts(1) + trace_start - 1;
%         risetime_cdk2_cusum(trace_id) = rise_cusum_real + trace_start - 1;
    catch
    end
    
end