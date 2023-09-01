function ImageNew=bgsub_fast(ImageOld,Radius,fra,ratio)
% fra: fraction; if zero, will be generated automatically
% ratio: ratio to shrink the image, 0-1

%% shrink the image
ImageOld_s=imresize(ImageOld,ratio);
Radius_s=ceil(Radius*ratio);
%% get background image
Imagebg_s=getbg_simple(ImageOld_s,Radius_s,fra);
Imagebg=imresize(Imagebg_s,size(ImageOld));

%% flatten the image
ImageNew=ImageOld-Imagebg;