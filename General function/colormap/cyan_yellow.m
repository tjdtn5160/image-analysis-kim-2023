function h = cyan_yellow(m)
%CYAN_YELLOW    Cyan-yellow color map
%   CYAN_YELLOW(M) returns an M-by-3 matrix containing a "cyan-yellow" colormap.
%   CYAN_YELLOW, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(hot)
%
%   See also HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, 
%   COLORMAP, RGBPLOT.

%   written by Sean Collins

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n = floor(m/2);
n2=m-n-1;

r = [zeros(n,1); (0:n2)'/n2];
g = [0.6824*(n:-1:1)'/n; 0.949*(0:n2)'/n2];  %Trying to optimize for cyan and yellow in Illustrator cmyk
b = [0.9373*(n:-1:0)'/n; zeros(n2,1)];

h = [r g b];
