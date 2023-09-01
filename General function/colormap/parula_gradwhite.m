function h=parula_gradwhite(m)
% Makes a variant of the jet colormap with black for the first row

if nargin<1
    m=255;
end
h=parula(m);
a=makeColorMap([0.95 0.95 0.95],[0.2081    0.1663    0.5292],56);
h=[a(1:55,:); h];
end

