function histogram_3D_v2(yvalues,bmin,bmax,num_bins,offset_increment,color1,color2,xlims,mthd)
%% Plots multiple histograms on the same axis to give a 3D-like image

% Input is a cell array of conditions, where each object in the cell is a
% vector of data points from each condition.
% Yaxis of the histogram is % of population 
% num_bins is how many bins the histogram will have
% offset_increment is how much each histogram is offset from the previous
% histogram (units are % population)
% color1 and color2 should be a 3-element vector represnting RGB values
% xlims should be a 2-element vector representing the lower xlimit and the
% upper xlimit respectively (eg. [4 14]). Default is the min/max values of
% the data set. 
% mthd designated which plotting method you want to use (ie 'histogram' or
% 'surface' or 'mesh'). Default is 'histogram'
%
% Written by Steve Cappell
% 2014_08_28

%% Default Parameters
if nargin<9 | isempty(mthd)
    mthd='histogram';
end
if nargin<8 
    xlims=[];
    mthd='histogram';
end

if nargin<7
    color2=[0 0 0];
    mthd='histogram';
end

if nargin<6
    color1=[1 0 0];
    color2=[0 0 1];
    mthd='histogram';
end

if nargin<5
    offset_increment=0.5;
    color1=[1 0 0];
    color2=[0 0 1];
    mthd='histogram';
end

if nargin<4
    num_bins=100;
    offset_increment=0.5;
    color1=[1 0 0];
    color2=[0 0 1];
    mthd='histogram';
end

if isempty(num_bins)
    num_bins=100;
end
if isempty(color1) & isempty(color2)
    color1=[1 0 0];
    color2=[0 0 1];
end
if isempty(offset_increment)
    offset_increment=0.5;
end
if isempty(xlims)
    xlims=[];
end
%%
number_of_conditions=length(yvalues);

%%% Create a Color Transition Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
color_col_1=linspace(color1(1),color2(1),number_of_conditions);
color_col_2=linspace(color1(2),color2(2),number_of_conditions);
color_col_3=linspace(color1(3),color2(3),number_of_conditions);
color_transition=[color_col_1', color_col_2',color_col_3'];

%%% Make Histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate max and min of values to determine bins
for condition=1:number_of_conditions; 
    min_value_each_condition(condition)=prctile(yvalues{condition},2);
    max_value_each_condition(condition)=prctile(yvalues{condition},99.5);
end
if isempty(bmin)
min_value=min(min_value_each_condition);
max_value=max(max_value_each_condition);
else
    min_value=bmin;
    max_value=bmax;
end
bins=linspace(min_value,max_value,num_bins);

%%% Compute histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for condition=1:number_of_conditions;
    histogram_of_data=hist(yvalues{condition},bins);
    percent_total=100*histogram_of_data/sum(histogram_of_data);
    smooth_percent_total(condition,:)=smooth(percent_total,2); %smooths histogram slightly
end

%%% Plot histograms as a filled polygon, offset slightly %%%%%%%%%%%%%%%%%%
%Extends the values to make a connected polygon from the original line
if isempty(xlims)
    xlim1=1;
    xlim2=length(bins);
else
    xlim1=find(bins>xlims(1),1,'first');
    xlim2=find(bins>xlims(2),1,'first');
    xlim2=xlim2-1;
end

if ~isempty(xlims) & xlims(2)>max(bins)
    error('histogram_3D:xlims_out_of_bounds', 'Upper X limit exceeds data size, try smaller value');
end
if ~isempty(xlims) & xlims(1)<min(bins)
    error('histogram_3D:xlims_out_of_bounds', 'Lower X limit is too small, try larger value');
end
    
new_xvalues=[bins(xlim1), bins(xlim1:xlim2), bins(xlim2), bins(xlim1)]; %Adds extra values to the Xaxis to make the line a closed shape
column_of_zeros=zeros(number_of_conditions,1);  %Create a vector of zeros to append to the yvalues to make the line a close shape
new_yvalues=[column_of_zeros, smooth_percent_total(:,xlim1:xlim2), column_of_zeros, column_of_zeros]; %Addes zeros to the beginning and end of the yvalues to make them a closed shape
zaxis=zeros(size(new_xvalues,1),size(new_xvalues,2));
%% Plot the histogram as a filled object with a vertical offset and changing color
xaxis_flip=fliplr(new_xvalues);
yaxis_flip=fliplr(new_yvalues);

switch mthd
    case 'histogram'
        for k=1:number_of_conditions;
            hold on
            inv_k=abs(k-(number_of_conditions+1));
            %     fill(new_xvalues,new_yvalues(inv_k,:)+((inv_k-1)*offset_increment),color_transition(inv_k,:),'EdgeColor',[0 0 0]);
            fill3(zaxis+((k-1)*offset_increment),xaxis_flip,yaxis_flip(inv_k,:),color_transition(inv_k,:),'EdgeColor',[0 0 0]);

            ylabel('Bins');zlabel('% population');xlim([0 ((number_of_conditions-1)*offset_increment)]);
            set(gca,'TickDir','out','xtick',[]);
            az=105;el=22;
            view(az, el);
        end
        hold off
        
    case 'surface'

        addition=zeros(1,size(new_yvalues,2));
        new_yvalues_2=[addition; new_yvalues;addition]; %add zeros as first and last timepoint to close off surface plot
        y=1:number_of_conditions+2; %extend yaxis so it matches dimensions of new_yvalues_2

        surf(new_xvalues,y,new_yvalues_2,'FaceColor',color1,'Edgecolor','none');
        xlabel('bins');ylabel('Conditions');zlabel('% population');ylim([0 number_of_conditions+2]);
        set(gca,'TickDir','out');
        camlight right;lighting phong
        az=10;el=26;
        view(az, el);
    
    case 'mesh'

        addition=zeros(1,size(new_yvalues,2));
        new_yvalues_2=[addition; new_yvalues;addition]; %add zeros as first and last timepoint to close off surface plot
        y=1:number_of_conditions+2; %extend yaxis so it matches dimensions of new_yvalues_2

        surf(new_xvalues,y,new_yvalues_2);
        xlabel('bins');ylabel('Conditions');zlabel('% population');ylim([0 number_of_conditions+2]);
        set(gca,'TickDir','out');
        az=10;el=26;
        view(az, el);
end
