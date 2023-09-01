function [th0_1a,th1_1a,th0_2a,th1_2a]=getmembraneth_weak(tempMIX0,sz1b,sz1s)

%%
tempMIX1=bgsub(tempMIX0,sz1b,0.05);
tempMIX1a=mat2gray(tempMIX1,prctile(double(tempMIX1(:)),[5,90]));
tempMIX2=bgsub(tempMIX0,sz1s,0.05);
tempMIX2a=mat2gray(tempMIX2,prctile(double(tempMIX2(:)),[5,90]));

%%
tempseries_2a=tempMIX2a((tempMIX2a>0)&(tempMIX2a<1));
[~,th0_2a,bg_2a]=ThreshImage(tempseries_2a);
qq_2a=bg_2a-median(tempseries_2a(tempseries_2a<=bg_2a));
th1_2a=bg_2a+4*qq_2a;

%%
tempseries_1a=tempMIX1a((tempMIX1a>0)&(tempMIX1a<1));
[~,th0_1a,bg_1a]=ThreshImage(tempseries_1a);
qq_1a=bg_1a-median(tempseries_1a(tempseries_1a<=bg_1a));
th1_1a=bg_1a+4*qq_1a;