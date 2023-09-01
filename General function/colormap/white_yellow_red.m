function h=white_yellow_red(m)
if nargin < 1, m = 255; end

a=makeColorMap([1 1 1],[1 1 0],125);
b=makeColorMap([1 1 0],[1 0 0],125);
h=[0 0 0;a;b];
end