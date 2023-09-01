function cir_r=getringmask(szo)

%% set up mask
szi=round(szo*sqrt(0.5));szd=szo-szi;
cir_o=strel('disk',szo,0);cir_o=getnhood(cir_o);
cir_i=strel('disk',szi,0);cir_i=getnhood(cir_i);
cir_r=cir_o;
cir_r(1+szd:end-szd,1+szd:end-szd)=~cir_i;
cir_r(1+szo,:)=[];cir_r(:,1+szo)=[];
cir_r=cir_r/sum(cir_r(:));