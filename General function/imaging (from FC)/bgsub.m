function ImageNew=bgsub(ImageOld,Radius,fra)
% fra: fraction; if zero, will be generated automatically

%% flatten the image
ImageNew=ImageOld-getbg_simple(ImageOld,Radius,fra);