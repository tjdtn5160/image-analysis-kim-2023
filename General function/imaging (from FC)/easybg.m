function bg=easybg(img,mask,dmb)
% img: original image
% mask: mask of original image
% dm: direction of background

%%
sz=size(img);

%% easy background
signalsum=sum(img.*(~mask),dmb);
masksum=sum(~mask,dmb);
signalavg=signalsum./masksum;

%% fill gap if any
if sum(isnan(signalavg))>0
    signalavg=getbg(signalavg,~isnan(signalavg),round(sz(2)/10));
end
    
%% smoothen
signalsmo=imfilter(signalavg,fspecial('gaussian',round(sz(2)/10),sz(2)/20),'replicate');

%% output
if dmb==1
    bg=ones(sz(1),1)*signalsmo;
elseif dmb==2
    bg=signalsmo*ones(1,sz(2));
end