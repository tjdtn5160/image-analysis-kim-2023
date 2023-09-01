function [alligned_cdk,alligned_geminin,TgridInFrames]=allignMyTraces_timeOfMitosis(cdk,geminin,timeofMitosis)
%Alligns traces
%default input


         %TgridInFrames = -60:110;     %change
         TgridInFrames = -100:200;
         
        %% Defines size of Index
        x=1:size(cdk,2);
        x2=1:size(geminin,2);


        %% Plots alligned traces
        alligned_cdk2=NaN(size(cdk,1), length(TgridInFrames));
        alligned_geminin_1=NaN(size(geminin,1), length(TgridInFrames));

        %%
        timeShiftPerCellFromZero=NaN(size(cdk,1),1);
        for i = 1:size(cdk, 1)
            ind=timeofMitosis(i);    %alligns traces to geminin turning on
            if isempty(ind)
                ind=280;
            end
            if isnan(ind)
                ind=280;
            end
            timeShiftPerCellFromZero(i)=ind;
            xx=x-ind;
            xxmin=max(min(TgridInFrames),min(xx));
            xxmax=min(max(TgridInFrames),max(xx));
            inds1=find(xx==xxmin);
            inds2=find(xx==xxmax);
            inda1=find(TgridInFrames==xxmin);
            inda2=find(TgridInFrames==xxmax);
%             smoothtrj=conv(cdk(i,:),normpdf(-3:3,0,2),'same');    %change
%             smoothtrj=smooth(cdk(i,:),1);
            smoothtrj=cdk(i,:);

            alligned_cdk2(i,inda1:inda2)=smoothtrj(inds1:inds2);
   
            xx2=x2-ind;
            xx2min=max(min(TgridInFrames),min(xx2));
            xx2max=min(max(TgridInFrames),max(xx2));
            inds21=find(xx2==xx2min);
            inds22=find(xx2==xx2max);
            inda21=find(TgridInFrames==xx2min);
            inda22=find(TgridInFrames==xx2max);
%             smoothtrj2=conv(geminin(i,:),normpdf(-3:3,0,2),'same');     %Change
%             smoothtrj2=smooth(geminin(i,:),1);
            smoothtrj2=geminin(i,:);
            
            alligned_geminin_1(i,inda21:inda22)=smoothtrj2(inds21:inds22);
        end

         alligned_cdk=alligned_cdk2;
         alligned_geminin=alligned_geminin_1;


    
    end


