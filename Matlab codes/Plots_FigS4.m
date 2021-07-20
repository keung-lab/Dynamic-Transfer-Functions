%%Plots_FigureS4
%Plots different gating strategies and resulting histograms

file_names={'20200121/dark/export_JL2020-01-21_dark_E2.0001.mqd.csv';'20200121/P2 light/export_JL2020-01-21_light_A2.0001.mqd.csv'};%Dark condition, light condition'
cd('Flow cytometry data')
for i=2
 

Y=csvread(string(file_names(i)),1,0);
mCherry=Y(:,6);
FSC_A=Y(:,1);
FSC_H=Y(:,2);
SSC_A=Y(:,3);

log_FSC=log10(FSC_A);
log_FSC_H=log10(FSC_H);
log_SSC=log10(SSC_A);
radius=.7;
%med_FSC=median(FSC_A);
%med_SSC=median(SSC_A);
med_FSC=median(log_FSC);
med_SSC=median(log_SSC);
k=log_SSC./log_FSC;
med_k=median(k);
%Gate=FSC_A>med_FSC-radius & FSC_A<med_FSC+radius;%Gate for FSC-A<=radius au
Gate=(log_FSC-med_FSC).^2+(k-med_k).^2<=radius^2;
Gate2=(log_FSC_H)./log_FSC>median(log_FSC_H./log_FSC)-.1 & (log_FSC_H)./log_FSC<median(log_FSC_H./log_FSC)+.1;
Both_Gates=(Gate);
%Gating strategy-large ratius of 0.7
figure
subplot(2,2,1)
scatter(FSC_A,SSC_A,1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,2)
scatter(FSC_A(Gate),SSC_A(Gate),1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,3)
scatter(FSC_A(Gate),FSC_H(Gate),1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,4)
scatter(FSC_A(Both_Gates),FSC_H(Both_Gates),1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
%saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/Revision figures/Gating_r_.7.png')
%Histogram for large radius=0.7
figure
edges=logspace(-2,4);
histogram(mCherry,edges,'Normalization','probability')
set(gca,'Xscale','log')
%saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/Revision figures/Gating_histogram_r_.7.png')

%Gating strategy-small ratius of 0.1
mCherry=Y(:,6);
FSC_A=Y(:,1);
FSC_H=Y(:,2);
SSC_A=Y(:,3);
radius=.1;
%med_FSC=median(FSC_A);
%med_SSC=median(SSC_A);
med_FSC=median(log_FSC);
med_SSC=median(log_SSC);
k=log_SSC./log_FSC;
med_k=median(k);
%Gate=FSC_A>med_FSC-radius & FSC_A<med_FSC+radius;%Gate for FSC-A<=radius au
Gate=(log_FSC-med_FSC).^2+(k-med_k).^2<=radius^2;
Gate2=(log_FSC_H)./log_FSC>median(log_FSC_H./log_FSC)-.1 & (log_FSC_H)./log_FSC<median(log_FSC_H./log_FSC)+.1;
Both_Gates=(Gate);
figure
subplot(2,2,1)
scatter(FSC_A,SSC_A,1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,2)
scatter(FSC_A(Gate),SSC_A(Gate),1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,3)
scatter(FSC_A(Gate),FSC_H(Gate),1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,4)
scatter(FSC_A(Both_Gates),FSC_H(Both_Gates),1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
%saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/Revision figures/Gating_r_.1.png')

%Histogram for small radius=0.1
figure
edges=logspace(-2,4);
histogram(mCherry,edges,'Normalization','probability')
set(gca,'Xscale','log')
%saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/Revision figures/Gating_histogram_r_.1.png')


%Gating strategy-large radius of 0.7, normalized to FSC-A
mCherry=Y(:,6);
FSC_A=Y(:,1);
FSC_H=Y(:,2);
SSC_A=Y(:,3);
mCherry=mCherry./FSC_A; %size correction used for Figures 2,3, 5,6
radius=.7;
%med_FSC=median(FSC_A);
%med_SSC=median(SSC_A);
med_FSC=median(log_FSC);
med_SSC=median(log_SSC);
k=log_SSC./log_FSC;
med_k=median(k);
%Gate=FSC_A>med_FSC-radius & FSC_A<med_FSC+radius;%Gate for FSC-A<=radius au
Gate=(log_FSC-med_FSC).^2+(k-med_k).^2<=radius^2;
Gate2=(log_FSC_H)./log_FSC>median(log_FSC_H./log_FSC)-.1 & (log_FSC_H)./log_FSC<median(log_FSC_H./log_FSC)+.1;
Both_Gates=(Gate);
figure
subplot(2,2,1)
scatter(FSC_A,SSC_A,1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,2)
scatter(FSC_A(Gate),SSC_A(Gate),1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,3)
scatter(FSC_A(Gate),FSC_H(Gate),1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')
subplot(2,2,4)
scatter(FSC_A(Both_Gates),FSC_H(Both_Gates),1,'filled')
xlim([1,200])
ylim([1,200])
set(gca,'Xscale','log')
set(gca,'Yscale','log')

%Histogram for large radius=0.7, mCherry normalized by FSC-A
figure
edges=logspace(-2,4);
histogram(mCherry,edges,'Normalization','probability')
set(gca,'Xscale','log')
%saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/Revision figures/Gating_histogram_areacorrected.png')

end