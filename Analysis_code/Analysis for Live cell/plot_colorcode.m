function h = plot_colorcode(x , y , c , linewidth)
% Plot color-coded lines
z = zeros(size(x));

h = surface([x ; x] , [y ; y] , [z ; z] , [c ; c] , 'facecol' , 'no' , 'edgecol' , 'interp' , 'linew' , linewidth);
end
