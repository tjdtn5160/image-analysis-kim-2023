function tempMIX4=generatemembranemask_weak(tempMIX0,sz1b,sz1s,th0_1a,th1_1a,th0_2a,th1_2a)

%%
tempMIX1=bgsub(tempMIX0,sz1b,0.05);
tempMIX1a=mat2gray(tempMIX1,prctile(double(tempMIX1(:)),[5,90]));
tempMIX2=bgsub(tempMIX0,sz1s,0.05);
tempMIX2a=mat2gray(tempMIX2,prctile(double(tempMIX2(:)),[5,90]));

%%
th_2a=1/mean([1/th0_2a,1/th1_2a]);
msk1=im2bw(tempMIX2a,th_2a);
th_1a=1/mean([1/th0_1a,1/th1_1a]);
msk2=im2bw(tempMIX1a,th_1a);
msk_no_hole=imfill(msk1|msk2,'holes');
msk3=im2bw(tempMIX2a,min([th0_2a,th1_2a])) & msk_no_hole;
msk4=im2bw(tempMIX1a,min([th0_1a,th1_1a])) & msk_no_hole;
tempMIX4=msk1|msk2|msk3|msk4;