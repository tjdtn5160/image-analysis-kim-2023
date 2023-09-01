function h=parula_gray(m)
% Makes a variant of the jet colormap with black for the first row

if nargin<1
    m=255;
end
h=parula(m);
h=[0.85 0.85 0.85; h];

