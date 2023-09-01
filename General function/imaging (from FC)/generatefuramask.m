function tempMIX4=generatefuramask(tempMIX0,sz1b,sz1s,th0_1a,th1_1a,th0_2a,th1_2a)

%%
tempMIX1=bgsub(tempMIX0,sz1b,0.05);
tempMIX1a=imadjust(mat2gray(tempMIX1));
tempMIX2=bgsub(tempMIX0,sz1s,0.05);
tempMIX2a=imadjust(mat2gray(tempMIX2));

%%
th_1a=2/(2*1/th0_1a+1/th1_1a);
th_2a=(2*th0_2a+th1_2a)/2;
msk1=im2bw(tempMIX1a,th_1a);
msk2=im2bw(tempMIX2a,th_2a);
tempMIX4=msk1|msk2;