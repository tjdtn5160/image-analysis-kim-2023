function im1=imageSubtractBackground(im,p_tile,blocksize)
%creates a new image background subtracted image. It uses the prctile
%percentile pixel intensity in a 32x32 neighborhood as the background
%intensity. prctile should be specified as a decimal (e.g. 0.10)

if nargin<2
    p_tile=10;    %default - use the 10th percentile
end
if nargin<3
    blocksize=16;
end

if ndims(im)==3
    [r c s]=size(im);
else
    [r c]=size(im);
end
if ndims(im)==2 || s==1     %single image
    im1=processSingleImage(im,p_tile,blocksize);
end
if ndims(im)==3 && s>1      %image stack
    im1=im*0;
    for j=1:s
        im1(:,:,j)=processSingleImage(im(:,:,j),p_tile,blocksize);
    end    
end
end

%------------------------------------------------------------------------

function im1=processSingleImage(im,p_tile,blocksize)
%this subfunction does the work (one for one image)

%tic
im=squeeze(im); %remove the singleton dimension if one exists
im=im-double(repmat(prctile(double(im),15,2),[1 size(im,2)]));  %remove any horizontal stripes of background
temp=blkproc(im,[blocksize blocksize],@(x) prctile(x(:),p_tile));
bg=imresize(temp,size(im),'bilinear');
%toc
%size(im)
%bg=ordfilt2(im,round(prctile*32*32),ones(32,32),'symmetric');
im1=im-bg;

% fprintf('Median bg = %.0f\n',median(double(bg(:))));
% showImagesWithLinkedAxes({im,bg,im1});
% figure; imagesc(bg); colorbar;

end