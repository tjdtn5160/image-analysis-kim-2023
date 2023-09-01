function h=parula_gradwhiteNum(num_colors)

if nargin<1
num_colors=24;
end
if num_colors<2
    num_colors=2;
end

h=parula(200);
a=makeColorMap([0.95 0.95 0.95],[0.2081    0.1663    0.5292],num_colors); 
h=[a(1:num_colors-1,:); h];
end