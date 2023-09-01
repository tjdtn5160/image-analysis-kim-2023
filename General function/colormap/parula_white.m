function h=parula_black(m)
% Makes a variant of the jet colormap with black for the first row

if nargin<1
    m=255;
end
h=parula(m);
h=[1 1 1; h];

