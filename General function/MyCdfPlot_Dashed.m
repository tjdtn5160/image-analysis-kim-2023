function [Y,n_cumsum] =  MyCdfPlot_Dashed(X,nbins,colorIndex,dashed)

if(length(X)>2) 

Y = [min(X):((max(X)-min(X))/nbins):max(X)];
n_elements=histc(X,Y);
n_elements=n_elements./length(X);
n_cumsum = cumsum(n_elements);



%i = find(n_cumsum>0.9,1,'first')
%Y(i)
%n_cumsum(i)

ColOrd=[1 0 0 ;0 1 0; 0 0 1;0 1 1;1 0 1;1 1 0;...
    0.5 0 0;0 0.5 0; 0 0 0.5;0 0.5 0.5;0.5 0 0.5;0.5 0.5 0;...
    .25 0 0;0 .25 0; 0 0 .25;0 .25 .25;.25 0 .25;.25 .25 0;...
    0.1 0 0;0 0.1 0; 0 0 0.1;0 0.1 0.1;0.1 0 0.1;0.1 0.1 0];
%thisColVec=ColOrd(mod(colorIndex,12)+1,:);
thisColVec=ColOrd(colorIndex,:);
if (dashed==1) 
    plot(Y,n_cumsum,'--','Color',thisColVec,'LineWidth',2);
    %plot(Y,n_cumsum,'--','Color',[0.25 0 0],'LineWidth',1);
else 
    plot(Y,n_cumsum,'-','Color',thisColVec,'LineWidth',2);
    %plot(Y,n_cumsum,'-','Color',[1 0 0],'LineWidth',2);
end

axis([min(X),max(X),0,1])
%hold on
%line([Y(i),Y(i)],[0,n_cumsum(i)],'LineWidth',4,'Color',[1 0 0]);
%line([0.5,5000],[1,10000],'LineWidth',4,'Color',[1 0 0]);



else 
    disp('Error empty length given');
end