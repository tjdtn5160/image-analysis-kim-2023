load('Data4');
condi = 1;
conditions = {'EGF20Serum5';'EGF10Serum2.5';'EGF5Serum1.25';'MitogenWithdraw';'DMSO';'MEKi';'CDK12i10uM';'CDK2i';'CDK4i';'MEKiCDK12i';'MEKi CDK2i';'MEKi CDK4i'};
data = ret{condi}.store.storeData;
divFrame = ret{condi}.store.allParentDivFrame;
D = sort_or_align_by_mitosis_demo2(data,divFrame)
if condi==8|condi==7|condi==10|condi==11
    D.clean_leastnum_frames(100); % this will remove frames shorter than 320.
else
    D.clean_leastnum_frames(200); % this will remove frames shorter than 320.
    D.clean_remove_suddendropCDK2
    D.clean_remove_CDK2_not_dropped
end
D.clean_remove_by_lower_upper_limit(2,[0.95,2.6])
D.clean_remove_mean_outside_n_std(2,3); % 1st arg means we use second channel, 2nd arg means discard outside mean+-2.5*std
D.clean_remove_CDK2_above_thres_at_mitosis(1.4)
D.clean_remove_dataaligned_indicate_timeandrange
D.clean_remove_dataaligned_indicate_timeandrange(1,4950,[0,1.1])
D.clean_remove_dataaligned_indicate_timeandrange(1,5050,[3,4])
% D.judgeCDK2Inc(0.002)
D.judgecdk2_detect_rising_pts

obj = D;
frameToPlot = 5000+[-12/obj.timeInterval:24/obj.timeInterval];
fastIncTime = 4/obj.timeInterval;
fastCdk = obj.dataAligned(obj.risingPoint<fastIncTime,:,:);
obj = D;
slowIncTime = [12/obj.timeInterval, 18/obj.timeInterval];
slowCdk = obj.dataAligned(obj.risingPoint>slowIncTime(1) & obj.risingPoint<slowIncTime(2),:,:);
lowCdk = obj.dataAligned(~obj.cellState,:,:);
dataByState{1} = lowCdk;
dataByState{2} = fastCdk;
dataByState{3} = slowCdk;

p = panel;p.pack(2,3);p.margin = [10,10,10,10];
time = (frameToPlot-5000).*obj.timeInterval;
for i = 1:3
p(1,i).select()
y = dataByState{i}(:,frameToPlot,1);
shadedErrorBar(time,y,{@nanmean,@nanstd},{'r-','markerfacecolor','r'});
xlim([time(1),time(end)]);ylim(obj.colorRange{1});h1(i) = gca;
p(2,i).select()
y = dataByState{i}(:,frameToPlot,2);
shadedErrorBar(time,y,{@nanmean,@nanstd},{'r-','markerfacecolor','r'});
xlim([time(1),time(end)]);ylim(obj.colorRange{2});h2(i) = gca;
end
linkaxes(h1,'xy');linkaxes(h2,'xy');
h3 = [h1,h2];
linkaxes(h3,'x')
