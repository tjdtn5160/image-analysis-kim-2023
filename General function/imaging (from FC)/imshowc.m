function imshowc(im1,im2,im3,contrast_inc,contrast_dec,brightness_inc,brightness_dec,option)

if isempty(contrast_dec)
    contrast_dec=1;
elseif isempty(contrast_inc)
    contrast_inc=0;
elseif isempty(brightness_dec)
    brightness_dec=1;
elseif isempty(brightness_inc)
    brightness_inc=0;
end


I(:,:,1)=im2double(im1)*1;
I(:,:,2)=zeros(size(im1));
I(:,:,3)=zeros(size(im1));

if (~isempty(im2))
    I(:,:,2)=im2double(im2)*1;
end

if (~isempty(im3))
    I(:,:,3)=im2double(im3)*1;
end

if (option==1)
    for i=1:3
        I(:,:,i) =imadjust(mat2gray(I(:,:,i)),[contrast_inc contrast_dec],[brightness_inc brightness_dec]);
    end
end

imshow(I);
