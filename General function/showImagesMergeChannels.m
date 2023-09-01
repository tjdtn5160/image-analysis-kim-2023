function im=showImagesMergeChannels(output,r,g,b,adjustFlag)

if nargin<5
    adjustFlag=1;
end
if nargin<4
    b=zeros(size(r));
end
if nargin<3
    g=zeros(size(r));
end

im=0*repmat(r,[1 1 3]);

if adjustFlag>0
    im(:,:,1)=mat2gray(r);
    im(:,:,2)=mat2gray(g);
    im(:,:,3)=mat2gray(b);
else
    im(:,:,1)=r;
    im(:,:,2)=g;
    im(:,:,3)=b;
end

if output
    imshow(im);
end