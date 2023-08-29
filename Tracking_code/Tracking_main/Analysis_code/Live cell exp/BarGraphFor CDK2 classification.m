root='J:\2. Data\IXmicro (New Axon)\1. Project\1. ERK and Cell cycle\2014-05-16 (MCF10A; ERK-EV-nls, hDHB-mCherry; 12min, 48hr)\';
Resultdir=[root 'Results\'];
load([Resultdir 'Cdk2 classification based on ERK.mat']);
label={'CDKinc'; 'CDKlow'}; 

figure;
bar(DMSO','BarWidth',1); ylim([0 100]); xlim([0.5 5.5]);
legend(label); ylabel('Cell Fraction'); title('DMSO');
xaxis={'Total' 'G1 10%' 'G1 90%' 'G2 10%' 'G2 90%'}; set(gca,'XTick',1:5,'XTickLabel',xaxis);
% export_fig([Resultdir 'CDK classification (DMSO)'],'-eps','-transparent','-nocrop');

figure;
bar(Mitogen','BarWidth',1); ylim([0 100]); xlim([0.5 5.5]);
legend(label); ylabel('Cell Fraction'); title('Mitogen withdraw');
xaxis={'Total' 'G1 10%' 'G1 90%' 'G2 10%' 'G2 90%'}; set(gca,'XTick',1:5,'XTickLabel',xaxis);
% export_fig([Resultdir 'CDK classification (Mitogen)'],'-eps','-transparent','-nocrop');

figure;
bar(MEKi','BarWidth',1); ylim([0 100]); xlim([0.5 5.5]);
legend(label); ylabel('Cell Fraction'); title('MEK inhibition');
xaxis={'Total' 'G1 10%' 'G1 90%' 'G2 10%' 'G2 90%'}; set(gca,'XTick',1:5,'XTickLabel',xaxis);
% export_fig([Resultdir 'CDK classification (MEKi)'],'-eps','-transparent','-nocrop');

figure;
bar(CDK4i','BarWidth',1); ylim([0 100]); xlim([0.5 5.5]);
legend(label); ylabel('Cell Fraction'); title('CDK4 inhibition');
xaxis={'Total' 'G1 10%' 'G1 90%' 'G2 10%' 'G2 90%'}; set(gca,'XTick',1:5,'XTickLabel',xaxis);
% export_fig([Resultdir 'CDK classification (CDK4i)'],'-eps','-transparent','-nocrop');