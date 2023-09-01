%%%
%resize the figure window.
%example 1:simply maximize the window
%figu;
%example 2:specifity the ratio
%figu(0.5,0.5);

function [] = figu(varargin)
    y = 1;x = 1;
    if size(varargin,2)>0
       y = varargin{1};
       x = varargin{2};
    end
    param = get(0,'ScreenSize');
    left = param(1);
    bottom = param(2);
    screenWidth = param(3);
    screenHeight = param(4);
    width = round(x*screenWidth);
    height = round(y*screenHeight);
    set(gcf,'position',[left, bottom, width, height]);
%     width = round(x*1600);
%     height = round(y*770);
%     set(gcf,'position',[0 47 width height]);
end