%CROverallAnalysis. Written on 20200219 by Jessica Lee
%Modified 20201029
%Modified 20210512 to be used for revision experiments done on
%20210425-20210501
%Requires additional m files: CRAnalysisFunction.m, fwer_hombonf.m, and
%ProfileClusterAnalysis.m
%Required Excel files: CRAnalysis.
close all
%Run CRAnalysisFunction to get foldChanges and CVs. 
experiment=input('What revision experiment? 1-AM, 2-PWM, 3-Truncated ZF');
if experiment==1
    fileName='Excel Files/CRAnalysis_AM';
[CRs_freq30,data_freq30,~,fc2_freq30,fc3_freq30,CV_freq30,darkCV_freq30,fc2_dark_freq30,file_freq30]=CRAnalysisFunction_revisions(fileName,'Cond2'); %Table for condition 2:1.2e10
close all
[CRs_freq1500,data_freq1500,~,fc2_freq1500,fc3_freq1500,CV_freq1500,darkCV_freq1500,fc2_dark_freq1500,file_freq1500]=CRAnalysisFunction_revisions(fileName,'Cond1'); %Table for condition 1: 6e9
close all
[CRs_freq10,data_freq10,~,fc2_freq10,fc3_freq10,CV_freq10,darkCV_freq10,fc2_dark_freq10,file_freq10]=CRAnalysisFunction_revisions(fileName,'Cond3'); %Table for condition 3:6e10
close all
[~,names]=xlsread('Excel Files/CRAnalysis_AM','List of CRs'); %Gives list of all CRs in column 2
elseif experiment==2
        fileName='Excel Files/CRAnalysis_PWM';
[CRs_freq30,data_freq30,~,fc2_freq30,fc3_freq30,CV_freq30,darkCV_freq30,fc2_dark_freq30,file_freq30]=CRAnalysisFunction_revisions(fileName,'Cond2'); %Table for condition 2
close all
[CRs_freq1500,data_freq1500,~,fc2_freq1500,fc3_freq1500,CV_freq1500,darkCV_freq1500,fc2_dark_freq1500,file_freq1500]=CRAnalysisFunction_revisions(fileName,'Cond1'); %Table for condition 1
close all
[CRs_freq10,data_freq10,~,fc2_freq10,fc3_freq10,CV_freq10,darkCV_freq10,fc2_dark_freq10,file_freq10]=CRAnalysisFunction_revisions(fileName,'Cond3'); %Table for condition 3
close all
[~,names]=xlsread('Excel Files/CRAnalysis_AM','List of CRs'); %Gives list of all CRs in column 2
elseif experiment==3
    fileName='Excel Files/Truncated_FM';
    [CRs_freq30,data_freq30,~,fc2_freq30,fc3_freq30,CV_freq30,darkCV_freq30,fc2_dark_freq30,file_freq30]=CRAnalysisFunction_revisions(fileName,'Cond2'); %Table for condition 2
close all
[CRs_freq1500,data_freq1500,~,fc2_freq1500,fc3_freq1500,CV_freq1500,darkCV_freq1500,fc2_dark_freq1500,file_freq1500]=CRAnalysisFunction_revisions(fileName,'Cond1'); %Table for condition 1
close all
[CRs_freq10,data_freq10,~,fc2_freq10,fc3_freq10,CV_freq10,darkCV_freq10,fc2_dark_freq10,file_freq10]=CRAnalysisFunction_revisions(fileName,'Cond3'); %Table for condition 3
close all
[~,names]=xlsread('Excel Files/Truncated_FM','List of CRs');
end

%% 
global number_of_clusters face_color face_color_CV
%fc2_dark_freq30=fc2_dark_freq30';
%fc2_dark_freq1500=fc2_dark_freq1500';
%fc2_dark_freq10=fc2_dark_freq10';
%Reorder data so they correspond
namesOfCRs=names(:,1);
numberOfCRs=length(namesOfCRs);

%averaged over biological replicates
FC1_cond2=zeros(numberOfCRs,1);
FC1_cond1=zeros(numberOfCRs,1);
FC1_cond3=zeros(numberOfCRs,1);

lightCV_cond2=zeros(numberOfCRs,1);
lightCV_cond1=zeros(numberOfCRs,1);
lightCV_cond3=zeros(numberOfCRs,1);
darkCV_cond2=zeros(numberOfCRs,1); %averaged over biological replicates
darkCV_cond1=zeros(numberOfCRs,1); %averaged over biological replicates
darkCV_cond3=zeros(numberOfCRs,1);

FC2_dark_cond1=zeros(numberOfCRs,1);
FC2_light_cond1=zeros(numberOfCRs,1);
FC2_dark_cond2=zeros(numberOfCRs,1);
FC2_light_cond2=zeros(numberOfCRs,1);
FC2_dark_cond3=zeros(numberOfCRs,1);
FC2_light_cond3=zeros(numberOfCRs,1);

FC3_cond1=zeros(numberOfCRs,1);
FC3_cond2=zeros(numberOfCRs,1);
FC3_cond3=zeros(numberOfCRs,1);



%not averaged
fc3_cond1=zeros(numberOfCRs,80);
fc3_cond2=zeros(numberOfCRs,80);
fc3_cond3=zeros(numberOfCRs,80);

fc2_cond1=zeros(numberOfCRs,80);
fc2_cond2=zeros(numberOfCRs,80);
fc2_cond3=zeros(numberOfCRs,80);

fc2_dark_cond1=zeros(numberOfCRs,80);
fc2_dark_cond2=zeros(numberOfCRs,80);
fc2_dark_cond3=zeros(numberOfCRs,80);

CV_cond1=zeros(numberOfCRs,80);
CV_cond2=zeros(numberOfCRs,80);
CV_cond3=zeros(numberOfCRs,80);

CV_dark_cond1=zeros(numberOfCRs,80); %not averaged
CV_dark_cond2=zeros(numberOfCRs,80); %not averaged
CV_dark_cond3=zeros(numberOfCRs,80); %not averaged
files_light_cond1=num2cell(NaN(numberOfCRs,80));
files_light_cond2=num2cell(NaN(numberOfCRs,80));
files_light_cond3=num2cell(NaN(numberOfCRs,80));
for i=1:numberOfCRs
    tempCR=namesOfCRs(i);
    cond2_ind=strcmp(CRs_freq30,tempCR);
    existvar=sum(cond2_ind);
    if existvar==0
        FC1_cond2(i)=NaN;
        FC2_dark_cond2(i)=NaN;
        FC2_light_cond2(i)=NaN;
        FC3_cond2(i)=NaN;
        lightCV_cond2(i)=NaN;
        darkCV_cond2(i)=NaN;
     files_light_cond2(i,:)=num2cell(NaN(1,80));
       fc3_cond2(i,:)=NaN(1,80);
       fc2_cond2(i,:)=NaN(1,80);
       CV_cond2(i,:)=NaN(1,80);
       fc2_dark_cond2(i,:)=NaN(1,80);
       CV_dark_cond2(i,:)=NaN(1,80);
    else
   FC1_cond2(i)=data_freq30(cond2_ind,1);
   FC2_dark_cond2(i)=data_freq30(cond2_ind,2);
   FC2_light_cond2(i)=data_freq30(cond2_ind,3);
   files_light_cond2(i,:)=file_freq30(:,cond2_ind);
   FC3_cond2(i)=data_freq30(cond2_ind,4);
   darkCV_cond2(i)=data_freq30(cond2_ind,5);
   lightCV_cond2(i)=data_freq30(cond2_ind,6);
   fc3_cond2(i,:)=fc3_freq30(:,cond2_ind);
  fc2_cond2(i,:)=fc2_freq30(:,cond2_ind);
  CV_cond2(i,:)=CV_freq30(:,cond2_ind);
  fc2_dark_cond2(i,:)=fc2_dark_freq30(:,cond2_ind);
  CV_dark_cond2(i,:)=darkCV_freq30(:,cond2_ind);
    end
    cond1_ind=strcmp(CRs_freq1500,tempCR);
    existvar=sum(cond1_ind);
    if existvar==0
        FC1_cond1(i)=NaN;
        FC2_dark_cond1(i)=NaN;
        FC2_light_cond1(i)=NaN;
        lightCV_cond1(i)=NaN;
        FC3_cond1(i)=NaN;
        darkCV_cond1(i)=NaN;
        fc3_cond1(i,:)=NaN(1,80);
        fc2_cond1(i,:)=NaN(1,80);
        CV_cond1(i,:)=NaN(1,80);
        fc2_dark_cond1(i,:)=NaN(1,80);
        CV_dark_cond1(i,:)=NaN(1,80);
        files_light_cond1(i,:)=num2cell(NaN(1,80));
    else
        FC1_cond1(i)=data_freq1500(cond1_ind,1);
        FC2_dark_cond1(i)=data_freq1500(cond1_ind,2);
        FC2_light_cond1(i)=data_freq1500(cond1_ind,3);
        FC3_cond1(i)=data_freq1500(cond1_ind,4);
        darkCV_cond1(i)=data_freq1500(cond1_ind,5);
        lightCV_cond1(i)=data_freq1500(cond1_ind,6);
        fc3_cond1(i,:)=fc3_freq1500(:,cond1_ind);
        fc2_cond1(i,:)=fc2_freq1500(:,cond1_ind);
        CV_cond1(i,:)=CV_freq1500(:,cond1_ind);
        fc2_dark_cond1(i,:)=fc2_dark_freq1500(:,cond1_ind);
        CV_dark_cond1(i,:)=darkCV_freq1500(:,cond1_ind);
        files_light_cond1(i,:)=file_freq1500(:,cond1_ind);
    end
        cond3_ind=strcmp(CRs_freq10,tempCR);
    existvar=sum(cond3_ind);
    if existvar==0
        FC1_cond3(i)=NaN;
        FC2_dark_cond3(i)=NaN;
        FC2_light_cond3(i)=NaN;
        lightCV_cond3(i)=NaN;
        FC3_cond3(i)=NaN;
        darkCV_cond3(i)=NaN;
        fc3_cond3(i,:)=NaN(1,80);
        fc2_cond3(i,:)=NaN(1,80);
        CV_cond3(i,:)=NaN(1,80);
        fc2_dark_cond3(i,:)=NaN(1,80);
        CV_dark_cond3(i,:)=NaN(1,80);
        files_light_cond3(i,:)=num2cell(NaN(1,80));
    else
        FC1_cond3(i)=data_freq10(cond3_ind,1);
        FC2_dark_cond3(i)=data_freq10(cond3_ind,2);
        FC2_light_cond3(i)=data_freq10(cond3_ind,3);
        FC3_cond3(i)=data_freq10(cond3_ind,4);
        darkCV_cond3(i)=data_freq10(cond3_ind,5);
        lightCV_cond3(i)=data_freq10(cond3_ind,6);
        fc3_cond3(i,:)=fc3_freq10(:,cond3_ind);
        fc2_cond3(i,:)=fc2_freq10(:,cond3_ind);
        CV_cond3(i,:)=CV_freq10(:,cond3_ind);
        fc2_dark_cond3(i,:)=fc2_dark_freq10(:,cond3_ind);
        CV_dark_cond3(i,:)=darkCV_freq10(:,cond3_ind);
        files_light_cond3(i,:)=file_freq10(:,cond3_ind);
    end
end
%% 
face_color=[0,0,0.29;
    0.13,.25,.6;
    .3,.56,.35;
    .99,0.68,.38;
    .74,.12,.18;
    .61,.4,.65;
    0,0,0;
    0,.8,1;
    0,1,0;];
face_color_CV=[0,0,0.29;
    0.13,.25,.6;
    .3,.56,.35;
    .74,.12,.18;
    .61,.4,.65;
    0,0,0;
    0,.8,1;
    0,1,0;];
namesOfCRs=names(:,1);
Save_condition=input('Save figure? 0=no, 1=yes');
AUC=[5*6e10/1500*14,5*6e10/30*14,6*6e10/10*14];
Freq_1=[1/1500, 1/30, 1/10];
[n,Albert]=xlsread('Excel Files/CRAnalysis','GO_complex');
%Synergy
Synergy_cond1=((FC2_light_cond1-1)-(FC2_dark_cond1-1)-(FC2_light_cond1(1)-1))./(FC2_light_cond1-1);
Synergy_cond2=((FC2_light_cond2-1)-(FC2_dark_cond2-1)-(FC2_light_cond2(1)-1))./(FC2_light_cond2-1);
Synergy_cond3=((FC2_light_cond3-1)-(FC2_dark_cond3-1)-(FC2_light_cond3(1)-1))./(FC2_light_cond3-1);
%calculation of slopes
AUC_cond1=5*6e10/1500*14*ones(numberOfCRs,1);
AUC_cond2=5*6e10/30*14*ones(numberOfCRs,1);
AUC_cond3=5*6e10/10*14*ones(numberOfCRs,1);
slope1_FC1_cond1_cond2=(FC1_cond2-FC1_cond1)./(AUC_cond2-AUC_cond1);
slope1_FC1_cond2_cond3=(FC1_cond3-FC1_cond2)./(AUC_cond3-AUC_cond2);
slope2_FC2_dark_cond1=(FC2_light_cond1-FC2_dark_cond1)./(AUC_cond1);
slope2_FC2_cond1_cond2=(FC2_light_cond2-FC2_light_cond1)./(AUC_cond2-AUC_cond1);
slope2_FC2_cond2_cond3=(FC2_light_cond3-FC2_light_cond2)./(AUC_cond3-AUC_cond2);
slope2_overall=[slope2_FC2_dark_cond1,slope2_FC2_cond1_cond2,slope2_FC2_cond2_cond3];
slope3_FC3_cond1_cond2=(FC3_cond2-FC3_cond1)./(AUC_cond2-AUC_cond1);
slope3_FC3_cond2_cond3=(FC3_cond3-FC3_cond2)./(AUC_cond3-AUC_cond2);
slope3_overall=[slope3_FC3_cond1_cond2, slope3_FC3_cond2_cond3];
slope5=(Synergy_cond2-Synergy_cond1)./(AUC_cond2-AUC_cond1);
slope6=(darkCV_cond2-darkCV_cond2)./(AUC_cond2-AUC_cond1);
slope7=(FC2_light_cond2-FC2_light_cond1)./(AUC_cond2-AUC_cond1);
FC3_overall=[FC3_cond1 FC3_cond2 FC3_cond3];%all conditions based on average
%Remove empty (NaN) conditions
nanIndices=any(isnan(FC3_overall),2);
FC3_overall=FC3_overall(~nanIndices,:);
namesOfCRs=namesOfCRs(~nanIndices);
fc3_cond1=fc3_cond1(:);
fc3_cond2=fc3_cond2(:);
fc3_cond3=fc3_cond3(:);
files_light_cond1=files_light_cond1(:);
files_light_cond2=files_light_cond2(:);
files_light_cond3=files_light_cond3(:);
fc2_dark_cond1=fc2_dark_cond1(:);
fc2_dark_cond2=fc2_dark_cond2(:);
fc2_dark_cond3=fc2_dark_cond3(:);
fc2_cond1=fc2_cond1(:);
fc2_cond2=fc2_cond2(:);
fc2_cond3=fc2_cond3(:);
CV_cond1=CV_cond1(:);
CV_cond2=CV_cond2(:);
CV_cond3=CV_cond3(:);
CV_dark_cond1=CV_dark_cond1(:);
CV_dark_cond2=CV_dark_cond2(:);
CV_dark_cond3=CV_dark_cond3(:);
fc3_overall=[fc3_cond1 fc3_cond2 fc3_cond3];%not averaged
fc2_dark_overall=[fc2_dark_cond1 fc2_dark_cond2 fc2_dark_cond3];%not averaged
fc2_overall=[fc2_cond1 fc2_cond2 fc2_cond3];%not averaged
files_overall=[files_light_cond1 files_light_cond2 files_light_cond3];
CV_overall=[CV_cond1 CV_cond2 CV_cond3];
CV_dark_overall=[CV_dark_cond1 CV_dark_cond2 CV_dark_cond3];

%Set anything less than 0, to a really small value
fc3_overall(fc3_overall<0)=.001;
fc2_dark_overall(fc2_dark_overall<0)=.001;
fc2_overall(fc2_overall<0)=.001;
FC3_overall(FC3_overall<0)=.001;
CV_overall(CV_overall<0)=.000001;

%Make vector for fc3_overall with CR names
CR_80=[];
CRs=names(:,1);
for i=1:80
  CR_80=[CR_80; CRs];  
end

VP16_JY30_fc3=fc3_overall(strcmp(CR_80,'JY30'),:);
VP16_JY145_fc3=fc3_overall(strcmp(CR_80,'JY145'),:);
VP16_JY145_fc2_dark=fc2_dark_overall(strcmp(CR_80,'JY145'),:);
VP16_JY145_CV=CV_overall(strcmp(CR_80,'JY145'),:);
VP16_JY145_fc2=fc2_overall(strcmp(CR_80,'JY145'),:);
VP16_JY145_CV_dark=CV_dark_overall(strcmp(CR_80,'JY145'),:);

nanIndices=any(isnan(fc3_overall),2);
fc3_overall(nanIndices,:)=[];
CV_overall(nanIndices,:)=[];
CR_80(nanIndices,:)=[];
files_overall(nanIndices,:)=[];
fc2_dark_overall(nanIndices,:)=[];
fc2_overall(nanIndices,:)=[];
CV_dark_overall(nanIndices,:)=[];

%All CRs without masking (masking will be done for clustering)
CR_80_not_masked=CR_80;
fc3_overall_notmasked=fc3_overall;
CV_overall_notmasked=CV_overall;
files_overall_notmasked=files_overall;
CV_dark_overall_notmasked=CV_dark_overall;
fc2_dark_overall_notmasked=fc2_dark_overall;
fc2_overall_notmasked=fc2_overall;

%remove VP16 from clustering (will be added back later)
fc3_overall(strncmp(CR_80,'VP16 ',5),:)=[];
files_overall(strncmp(CR_80,'VP16 ',5),:)=[];
fc2_dark_overall(strncmp(CR_80,'VP16 ',5),:)=[];
fc2_overall(strncmp(CR_80,'VP16 ',5),:)=[];
CV_overall(strncmp(CR_80,'VP16 ',5),:)=[];
CV_dark_overall(strncmp(CR_80,'VP16 ',5),:)=[];
CR_80(strncmp(CR_80,'VP16 ',5))=[];


%Apply filters to remove low variability and low values of FC
mask=genevarfilter(fc3_overall);
fc3_overall=fc3_overall(mask,:);
fc2_dark_overall=fc2_dark_overall(mask,:);
fc2_overall=fc2_overall(mask,:);
CV_overall=CV_overall(mask,:);
CV_dark_overall=CV_dark_overall(mask,:);
CR_80=CR_80(mask,:);
files_overall=files_overall(mask,:);


[mask,fc3_overall,CR_80]=genelowvalfilter(fc3_overall,CR_80,'absval',2);
fc2_dark_overall=fc2_dark_overall(mask,:);
fc2_overall=fc2_overall(mask,:);
CV_overall=CV_overall(mask,:);
CV_dark_overall=CV_dark_overall(mask,:);
files_overall=files_overall(mask,:);

%Add VP16 back into the values for clustering (just one though)
CR_80=['VP16 (JY145)';CR_80];
CV_overall=[nanmean(VP16_JY145_CV);CV_overall];
fc3_overall=[nanmean(VP16_JY145_fc3);fc3_overall];
fc2_dark_overall=[nanmean(VP16_JY145_fc2_dark);fc2_dark_overall];
fc2_overall=[nanmean(VP16_JY145_fc2);fc2_overall];
CV_dark_overall=[nanmean(VP16_JY145_CV_dark);CV_dark_overall];
files_overall=[{'hi','hi','hi'};files_overall];

%Take the log of everything
CV_overall=log10(CV_overall);
CV_dark_overall=log10(CV_dark_overall);
CV_overall_notmasked=log10(CV_overall_notmasked);
CV_dark_overall_notmasked=log10(CV_dark_overall_notmasked);
fc3_overall=log10(fc3_overall);
fc2_dark_overall=log10(fc2_dark_overall);
fc2_overall=log10(fc2_overall);
fc3_overall_notmasked=log10(fc3_overall_notmasked);
FC3_overall=log10(FC3_overall);

%Adjust value based on dark fc2, ie set the fc2 equal for each condition.
%(This does not affect fc3)
mean_fc2_dark=nanmean(fc2_dark_overall,2);
mean_fc2_dark_notmasked=nanmean(fc2_dark_overall_notmasked,2);
fc2_overall=fc2_overall./fc2_dark_overall.*mean_fc2_dark;
fc2_dark_overall=fc2_dark_overall./fc2_dark_overall.*mean_fc2_dark;
fc2_dark_overall_notmasked=fc2_dark_overall_notmasked./fc2_dark_overall_notmasked.*mean_fc2_dark_notmasked;
% mean_CV_dark=nanmean(CV_dark_overall,2);
% CV_overall=CV_overall./CV_dark_overall.*mean_CV_dark;
% CV_dark_overall=CV_dark_overall./CV_dark_overall.*mean_CV_dark;

axes_positions=[.05,.7567, .45,.1683;.5203,.7565,.45,.1683;
                .05,.456,.45,.1683;  .5203,.456,.45,.1683;
                .05,.1554,.45,.1683; .5203,.1554,.45,.1683];

%Cluster profile of all replicates for FC3
number_of_clusters=5;
[cdis_FC3,~,Distance_matrix_fc3,ctrs]=ProfileClusterAnalysis(CRs,Freq_1,(fc3_overall),number_of_clusters);%based on each replicate
[~,cdis_FC3_not_masked] = pdist2(ctrs,(fc3_overall_notmasked),'correlation','Smallest',1);
%Values based off of k-means but reordered and slightly changed for more
%specific filtering
ctrs=[0.756341421875678,-0.3,-0.338099818683020;
    -0.736510820962514,0.102405013200597,0.65;
    -0.351484685217366,0.769986528990965,-0.418501843773599;
    -0.775512363928304,0.550835601762814,0.3;
    0.0178768566235528,-0.685625005390339,0.667748148766786];
[~,cdis_FC3_not_masked] = pdist2(ctrs,(fc3_overall_notmasked),'correlation','Smallest',1);
[~,cdis_FC3]=pdist2(ctrs,fc3_overall,'correlation','Smallest',1);
for c=1:number_of_clusters
    subplot(3,2,c);
    %set(gca,'position',positions_profile_plots(c,:))
    %plot((ind_variable),(dep_variable(cidx==c,:))','color',face_color(c,:),'linewidth',3);
    plot(1:3,(ctrs(c,:))','color',face_color(c,:),'linewidth',3);
    axis tight
    set(gca,'position',axes_positions(c,:))
    %ylim([-.8,.8])
    set(gca,'yticklabel','')
end
%saveas(gcf,'LogPlots/Profile_FC3.png')

%Cluster profiles of all replicates for CV
[cdis_CV,~,~,ctrsCV]=ProfileClusterAnalysis(CRs,Freq_1,(CV_overall),number_of_clusters);
[~,cdis_CV_not_masked] = pdist2(ctrs,(CV_overall_notmasked),'correlation','Smallest',1);
%Adjust profiles to make more distinct filtering behaviors (1-saturation,
%2-linear, 3-band-pass, 4-band-stop, 5-low-pass
ctrsCV=[0.788353724500856,-0.394924286856965,-0.393429437643890;
    -0.516116430390331,-0.213976559640784,0.730092990031115;
    -0.202378178954878,0.748614469598275,-0.546236290643397;
    0.322669678434008,-0.772813223608021,0.450143545174013;
    0.498197985466675,0.287914018610248,-0.786112004076924];
[~,cdis_CV]=pdist2(ctrsCV,CV_overall,'correlation','Smallest',1);
subplot(3,2,4)

for c=1:number_of_clusters
    subplot(3,2,c);
    %set(gca,'position',positions_profile_plots(c,:))
    %plot((ind_variable),(dep_variable(cidx==c,:))','color',face_color(c,:),'linewidth',3);
    plot(1:3,(ctrsCV(c,:))','color',face_color_CV(c,:),'linewidth',3);
    axis tight
    set(gca,'position',axes_positions(c,:))
    %ylim([-.8,.8])
    set(gca,'yticklabel','')
end

fc2_comp=[fc2_overall(:,1);fc2_overall(:,2);fc2_overall(:,3)];
fc3_comp=[fc3_overall(:,1);fc3_overall(:,2);fc3_overall(:,3)];
CV_comp=[CV_overall(:,1);CV_overall(:,2);CV_overall(:,3)];
CV_dark_comp=[CV_dark_overall(:,1);CV_dark_overall(:,2);CV_dark_overall(:,3)];
cdis_FC3_comp=[cdis_FC3';cdis_FC3';cdis_FC3'];
CR_comp=[CR_80;CR_80;CR_80];
D='Dark';
A='Cond1';
B='Cond2';
C='Cond3';
cond_comp=[repmat({A},length(cdis_FC3),1); repmat({B},length(cdis_FC3),1); repmat({C},length(cdis_FC3),1)];

%% Find CRs of interest with split profiles/MIs
%load in MI
Save_condition=input('Save figure? 0=no, 1=yes');
[num,text2]=xlsread(fileName,'MI');
CR_compiled=text2(2:end,1);
MI_compiled=num(1:end,1);
VP16_MI=nanmean(MI_compiled(strcmp(CR_compiled,'JY145')));
MI_compiled=MI_compiled./VP16_MI;
files_compiled=text2(2:end,4);
test2=zeros(length(CRs),1);
test1=zeros(length(CRs),1);
for i=1:length(CRs)
   Indices=strcmp(CR_80_not_masked,CRs(i));
   Indices2=strcmp(CR_compiled,CRs(i));
   if sum(Indices)==0
   else
   CR_mode=mode(cdis_FC3_not_masked(Indices));
   if sum(cdis_FC3_not_masked(Indices)==CR_mode)/sum(Indices)<.5
       test1(i)=1;
   else
   end
   MI_CV=nanstd(MI_compiled(Indices2))/nanmean(MI_compiled(Indices2));
   if MI_CV>.8
       test2(i)=1;
   else
   end
   
   end
   
end
%Plot average MI for each cluster. Based on individual replicates
cdis_reordered=NaN(length(MI_compiled),1);
cdis_CV_reordered=NaN(length(MI_compiled),1);
fc3_reordered=NaN(length(MI_compiled),3);
CV_reordered=NaN(length(MI_compiled),3);
CV_dark_reordered=NaN(length(MI_compiled),3);
fc2_dark_reordered=NaN(length(MI_compiled),3);
fc2_reordered=NaN(length(MI_compiled),3);
MI_CRs=unique(CR_compiled);
k=1;
for j=1:length(CR_compiled)
    Indices1=strcmp(CR_80_not_masked,CR_compiled(j));
    Indices=strcmp(files_overall_notmasked(:,1),files_compiled(j));
    Indices_overall=Indices1+Indices==2;
    Indices3=strcmp(CR_80,CR_compiled(j));
    Indices4=strcmp(files_overall(:,1),files_compiled(j));
    Indices_overall2=Indices3+Indices4==2;
    if sum(Indices_overall2)==0 
        if sum(Indices_overall)==0
        cdis_reordered(k)=NaN;
        cdis_CV_reordered(k)=NaN;
        else
        cdis_reordered(k)=0;
        cdis_CV_reordered(k)=0;
        end
    else
      cdis_reordered(k)=cdis_FC3(Indices_overall2);  
        cdis_CV_reordered(k,:)=cdis_CV(Indices_overall2);
    end
    if sum(Indices_overall)==0
        fc3_reordered(k,:)=NaN(1,3);
        CV_reordered(k,:)=NaN(1,3);
        fc2_dark_reordered(k,:)=NaN(1,3);
        fc2_reordered(k,:)=NaN(1,3);
        CV_dark_reordered(k,:)=NaN(1,3);
    else
    fc3_reordered(k,:)=fc3_overall_notmasked(Indices_overall,:);
    CV_reordered(k,:)=CV_overall_notmasked(Indices_overall,:);
    fc2_dark_reordered(k,:)=fc2_dark_overall_notmasked(Indices_overall,:);
    fc2_reordered(k,:)=fc2_overall_notmasked(Indices_overall,:);
    CV_dark_reordered(k,:)=CV_dark_overall_notmasked(Indices_overall,:);
    end
 k=k+1;  
end
%Remove any MIs that don't have FC3 (due to outlier or less than 500 cells)
indices=(~any(isnan(fc3_reordered),2));
CR_compiled=CR_compiled(indices);
MI_compiled=MI_compiled(indices);
fc3_reordered=fc3_reordered(indices,:);
CV_reordered=CV_reordered(indices,:);
CV_dark_reordered=CV_dark_reordered(indices,:);
cdis_reordered=cdis_reordered(indices);
cdis_CV_reordered=cdis_CV_reordered(indices);
fc2_dark_reordered=fc2_dark_reordered(indices,:);



 %% Graph of MI average for each CR
 if experiment==3
     [num,text2]=xlsread('Excel Files/Supplemental Table 3','Truncated ZF');
 elseif experiment==2
     [num,text2]=xlsread('Excel Files/Supplemental Table 3','PWM');
 elseif experiment==1
     [num,text2]=xlsread('Excel Files/Supplemental Table 3','AM');
 end
 CR=text2(3:end,1);
 %CR=CR_compiled;
 
 MI=num(1:end,1);
 MI_VP16=nanmean(MI(strcmp(CR,'JY30')));
 MI(strcmp(CR,'JY145'))=[];
 CR(strcmp(CR,'JY145'))=[];
 MI(strcmp(CR,'JY28'))=[];
 CR(strcmp(CR,'JY28'))=[];
 CR(strcmp(CR,'JY30'))={'VP16 only'};
 MI=MI./MI_VP16;
 %MI=MI_compiled;
 
 if experiment==3
     [num,text2]=xlsread('Excel Files/Supplemental Table 3','FM');
     CR_FM=text2(3:end,1);
     MI_FM=num(1:end,3);
 MI_VP16_FM=nanmean(MI_FM(strcmp(CR_FM,'VP16 (JY30)')));
 MI_FM(strcmp(CR_FM,'VP16 (JY145)'))=[];
 CR_FM(strcmp(CR_FM,'VP16 (JY145)'))=[];
 MI_FM(strcmp(CR_FM,'JY28'))=[];
 CR_FM(strcmp(CR_FM,'JY28'))=[];
 CR_FM(strcmp(CR_FM,'VP16 (JY30)'))={'VP16 only'};
 MI_FM=MI_FM./MI_VP16_FM;
 index=find(ismember(CR_FM,CR));%only keep those that are the same as truncated ZF screen
 CR_FM=CR_FM(index);
 MI_FM=MI_FM(index);
 Truncated_temp=[ones(length(CR),1);zeros(length(CR_FM),1)];
 Truncated=string(length(Truncated_temp));
 Truncated(Truncated_temp==1)='Truncated';
 Truncated(Truncated_temp==0)='Not truncated';
 t=table(CR,MI,'VariableNames',{'Gene','MI'});
 %t=table([CR;CR_FM],[MI;MI_FM],Truncated','VariableNames',{'Gene','MI','Truncated'});
 
 t2=table(CR_FM,MI_FM,'VariableNames',{'Gene','MI'});
t = sortrows(t,'Gene','descend');
t2 = sortrows(t2,'Gene','descend');
 stats_truncated=grpstats(t,{'Gene'},{'mean','std','meanci','sem'});
 stats_not_truncated=grpstats(t2,{'Gene'},{'mean','std','meanci','sem'});
stats_truncated = sortrows(stats_truncated,'Gene','descend');
stats_not_truncated = sortrows(stats_not_truncated,'Gene','descend');
 figure
 hold on
 %plot mean values
 %not truncated does not have Y11
 y_position=[stats_truncated.mean_MI';stats_not_truncated.mean_MI'];
 x_position=1:length(stats_truncated.mean_MI);
 b=bar(y_position');
 b(1).FaceColor=[.7,.7,.7];
 b(2).FaceColor=[1,1,1];
 %Plot individual points and get Welch t-test p-values
 %Add significance *
  unique_CRs=stats_truncated.Gene;
 p_values=NaN(length(unique_CRs),1);
 for i=1:length(unique_CRs)
     X=t2.MI(strcmp(t2.Gene,unique_CRs(i)));%truncated
     Y=t.MI(strcmp(t.Gene,unique_CRs(i)));%not truncated
     [h,p_values(i)]=ttest2(X,Y);
     plot(x_position(i).*ones(length(Y),1)-.15,Y','.','color',[.3,.3,.3])
     plot(x_position(i).*ones(length(X),1)+.15,X','.','color',[.3,.3,.3])
     if p_values(i)<0.05
         plot(x_position(i),max(stats_not_truncated.mean_MI(i),stats_truncated.mean_MI(i))+.2,'*k')
     else
     end
 end
 
 set(gca,'xtick',1:length(stats_truncated.mean_MI))
set(gca,'XTickLabel',stats_truncated.Gene)
set(gca,'XTickLabelRotation',90)
ylabel('MI_F_M')
legend({'Truncated ZF','Full ZF'}) 
ylim([0,1.5])
%  [~,~,stats]=anovan(t.MI,{t.Gene,t.Truncated},'model','interaction','varnames',{'CR','Truncated'})
%  multcompare(stats)
 else
   t=table(CR,MI,'VariableNames',{'Gene','MI'});  
 
 
 t = sortrows(t,'MI','ascend');
 [~,~,stats]=anova1(t.MI,t.Gene);
 multcompare(stats)
 
 Stats_MI=grpstats(t,{'Gene'},{'mean','std','meanci','sem'});
 Stats_MI=Stats_MI(Stats_MI.GroupCount>=3,:);
 Stats_MI = sortrows(Stats_MI,'mean_MI','ascend');
 figure
 hold on
 set(gcf,'position',[440,225,1000,450]);
 bar(1:length(Stats_MI.Gene(Stats_MI.mean_MI<=0.5)), Stats_MI.mean_MI(Stats_MI.mean_MI<=.5),'facecolor',[1,.85,.9])
 bar(length(Stats_MI.Gene(Stats_MI.mean_MI<=0.5))+1:sum(Stats_MI.mean_MI<1), Stats_MI.mean_MI(Stats_MI.mean_MI>.5 & Stats_MI.mean_MI<1),'facecolor',[.7,.7,.7])
 bar(length(Stats_MI.Gene(Stats_MI.mean_MI<1))+1:(length(Stats_MI.mean_MI)), (Stats_MI.mean_MI(Stats_MI.mean_MI>=1)),'facecolor',[.8,.8,.99])
    errorbar(1:length(Stats_MI.Gene),Stats_MI.mean_MI,Stats_MI.std_MI,'linestyle','none','color',[.75,.75,.75],'Capsize',0)

 for i=1:length(Stats_MI.Gene)
     plot(i,MI(strcmp(CR,Stats_MI.Gene(i))),'.','color',[0.3,.3,.3],'Markersize',10)
    % boxplot(MI(strcmp(CR,Stats_MI.Gene(i))),'positions',i,'plotstyle','compact')
 end
 %errorbar(1:length(Stats_MI.Gene),Stats_MI.mean_MI,Stats_MI.meanci_MI(:,2)-Stats_MI.mean_MI,'linestyle','none','color',[.5,.5,.5],'Capsize',0)
 set(gca,'XTick',1:length(Stats_MI.Gene))
    set(gca,'XTickLabel',Stats_MI.Gene)
    set(gca,'XTickLabelRotation',90)
    ylim([0,1.8])
 end
    if experiment==1
    ylabel('MI_A_M')
    saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/CR screen/LogPlots/MI_AM.png')
    elseif experiment==2
        ylabel('MI_P_W_M')
    saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/CR screen/LogPlots/MI_PWM.png')
    elseif experiment==3
        ylabel('MI_F_M')
        saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/CR screen/LogPlots/MI_truncatedFM.png')
    end
    
%  %% Determining link between GO/GA and MI
%  %remove VP16
%  MI_compiled_noVP16=MI_compiled(~strncmp(CR_compiled,'VP16',4));
%  %Average MI for each GO
%  Average_MI_GO=zeros(length(Complex_matrix(1,:)),1);
%  MI_GO=NaN(length(MI_compiled_noVP16),length(Complex_matrix(1,:)));
%  figure
%  hold on
%  for i=1:length(Complex_matrix(1,:))
%      Average_MI_GO(i)=nanmean(MI_compiled_noVP16(logical(Complex_matrix(:,i))));
%      MI_GO(logical(Complex_matrix(:,i)),i)=MI_compiled_noVP16(logical(Complex_matrix(:,i)));
%      err(i)=nanstd(MI_compiled_noVP16(logical(Complex_matrix(:,i))))./sqrt(length(MI_compiled_noVP16(logical(Complex_matrix(:,i)))));
% 
%  end
%  [values,I]=sort(Average_MI_GO,'descend');
%  err=err(I);
%  labels=[];
%  labels=unique_complex(I);
%  labels(isnan(values))=[];
%  err(isnan(values))=[];
%  values(isnan(values))=[];
%       bar(1:length(values(values>=1.13)),values(values>=1.13),'FaceColor',[.9,.7,.7]);
%       bar(length(values(values>=1.13))+1:length(values(values>=.5)),values(values<=1.13 & values>=.5),'FaceColor',[.9,.9,.9]);
%       bar(length(values(values>=.5))+1:length(values),values(values<.5),'FaceColor',[.5,.5,.5]);
%      errorbar(1:length(values),values,err,'linestyle','none','color','k','linewidth',2,'capsize',0);
%      set(gca,'XTick',1:length(values))
%     set(gca,'XTickLabel',labels)
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     ylim([0,1.5])
%     saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/CR screen/LogPlots/averageMI_GO.png')
%  figure
% [p_value_GO_MI_anova,~,stats]=anova1(MI_GO)
% multcompare_CV_MI=multcompare(stats)
%  %Fisher Exact test between GO and MI above 1
%  enrichment=NaN(length(Complex_matrix(1,:)),1);
% p_MI_GO=NaN(length(Complex_matrix(1,:)),1);
% for i=1:length(Complex_matrix(1,:))
%     
%     temp=crosstab(Complex_matrix(:,i),MI_increased);
%  if sum(Complex_matrix(:,i))==0
%  else
%     [h,p_MI_GO(i),stats]=fishertest(temp,'tail','both');
%     enrichment(i)=temp(2,2)./sum(temp(2,:))./(sum(temp(:,2))./sum(temp,'all'));
%  end
%     
% end
% [p_MI_GO_corrected,c_alpha,h]=fwer_holmbonf(p_MI_GO,0.05); %Holm-Bonferroni correction
% 
%     p_MI_GO_transformed=-log10(p_MI_GO_corrected);
%  figure
%  hold on
%  [values,I]=sort(p_MI_GO_transformed,'descend');
%     bar((1:10),values(1:10),'facecolor',[.9,.7,.7]);
%     labels=unique_complex(I);
%     set(gca,'XTick',1:10)
%     set(gca,'XTickLabel',labels(1:10))
%     set(gca,'XTickLabelRotation',45)
%     ylim([0,2])
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_GO_pvalues.png')
%  end
%   sig=p_MI_GO_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(enrichment,'descend');
%     labels=unique_complex(I);
%     sigtemp=sig(I);
%     sigtemp=sigtemp(1:10);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:10,values(1:10),'FaceColor',[.9,.7,.7]);
%     set(gca,'XTick',1:10)
%     set(gca,'XTickLabel',labels(1:10))
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,7])
%     scatter(x,values(sigtemp)+max(14)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_GO_enrichment.png')
%  end
%  
%   %Fisher Exact test between GO and MI below 0.5
%   p_MI_decrease_GO=NaN(length(Complex_matrix(1,:)),1);
%   enrichment=NaN(length(Complex_matrix(1,:)),1);
% for i=1:length(Complex_matrix(1,:))
%     
%     temp=crosstab(Complex_matrix(:,i),MI_decreased);
%  if sum(Complex_matrix(:,i))==0
%  else
%     [h,p_MI_decrease_GO(i),stats]=fishertest(temp,'tail','both');
%     enrichment(i)=temp(2,2)./sum(temp(2,:))./(sum(temp(:,2))./sum(temp,'all'));
%  end
%     
% end
%     [p_MI_decrease_GO_corrected,c_alpha,h]=fwer_holmbonf(p_MI_decrease_GO,0.05); %Holm-Bonferroni correction
% 
%     p_MI_decrease_GO_transformed=-log10(p_MI_decrease_GO_corrected);
%  figure
%  hold on
%  [values,I]=sort(p_MI_decrease_GO_transformed,'descend');
%     bar((1:10),values(1:10),'facecolor',[.5,.5,.5]);
%     labels=unique_complex(I);
%     set(gca,'XTick',1:10)
%     set(gca,'XTickLabel',labels(1:10))
%     set(gca,'XTickLabelRotation',45)
%     ylim([0,3])
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_decrease_GO_pvalues.png')
%  end
%  
%   sig=p_MI_decrease_GO_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(enrichment,'descend');
%     labels=unique_complex(I);
%     sigtemp=sig(I);
%     sigtemp=sigtemp(1:10);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:10,values(1:10),'FaceColor',[.5,.5,.5]);
%     set(gca,'XTick',1:10)
%     set(gca,'XTickLabel',labels(1:10))
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,7])
%     scatter(x,values(sigtemp)+max(5)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_decrease_GO_enrichment.png')
%  end
%   %Average MI for each GA
%  Average_MI_GA=zeros(length(Activity_matrix(1,:)),1);
%  MI_GA=NaN(length(MI_compiled_noVP16),length(Activity_matrix(1,:)));
%  figure
%  hold on
%  for i=1:length(Activity_matrix(1,:))
%      Average_MI_GA(i)=nanmean(MI_compiled_noVP16(logical(Activity_matrix(:,i))));
%      MI_GA(logical(Activity_matrix(:,i)),i)=MI_compiled_noVP16(logical(Activity_matrix(:,i)));
%      err(i)=nanstd(MI_compiled_noVP16(logical(Activity_matrix(:,i))))./sqrt(length(MI_compiled_noVP16(logical(Activity_matrix(:,i)))));
% 
%  end
%     [values,I]=sort(Average_MI_GA);
%     
%      bar(1:length(Average_MI_GA),values,'FaceColor',[.5,.5,.5]);
%      errorbar(1:length(Average_MI_GA),values,err(I),'linestyle','none','color','k','linewidth',2,'capsize',0);
%      set(gca,'XTick',1:length(Average_MI_GA))
%     set(gca,'XTickLabel',unique_activity(I))
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     ylim([0,1.5])
%     %saveas(gcf,'LogPlots/AverageMI_GA.png')
%  figure
% [p_value_GA_MI_anova,~,stats]=anova1(MI_GA)
% multcompare_CV_MI=multcompare(stats)
% 
%  %Fisher Exact test between GA and MI above 1
%  enrichment=NaN(length(Activity_matrix(1,:)),1);
% p_MI_GA=NaN(length(Activity_matrix(1,:)),1);
% for i=1:length(Activity_matrix(1,:))
%     
%     temp=crosstab(Activity_matrix(:,i),MI_increased);
%  if sum(Activity_matrix(:,i))==0
%  else
%     [h,p_MI_GA(i),stats]=fishertest(temp,'tail','both');
%     enrichment(i)=temp(2,2)./sum(temp(2,:))./(sum(temp(:,2))./sum(temp,'all'));
%  end
%     
% end
%     [p_MI_GA_corrected,c_alpha,h]=fwer_holmbonf(p_MI_GA,0.05); %Holm-Bonferroni correction
% 
%     p_MI_GA_transformed=-log10(p_MI_GA_corrected);
%  figure
%  hold on
%  p_MI_GA_transformed(strcmp(unique_activity,''))=[];
%  labels=unique_activity(~strcmp(unique_activity,''))
%  [values,I]=sort(p_MI_GA_transformed,'descend');
%     bar(1:length(values),values,'FaceColor',[.9,.7,.7]);
%     set(gca,'XTick',1:length(Activity_matrix(1,:)))
%     set(gca,'XTickLabel',labels(I))
%      set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     ylim([0,1.5])
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_GA_pvalues.png')
%  end
%  
%    sig=p_MI_GA_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(enrichment,'descend');
%     labels=unique_activity(I);
%     sigtemp=sig(I);
%     %sigtemp=sigtemp(1:10);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:length(values),values,'FaceColor',[.9,.7,.7]);
%     set(gca,'XTick',1:length(values))
%     set(gca,'XTickLabel',labels)
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,5])
%     scatter(x,values(sigtemp)+max(5)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_GA_enrichment.png')
%  end
%  
%   %Fisher Exact test between GA and MI below 0.5
%  enrichment=NaN(length(Activity_matrix(1,:)),1);
% p_MI_decrease_GA=NaN(length(Activity_matrix(1,:)),1);
% for i=1:length(Activity_matrix(1,:))
%     
%     temp=crosstab(Activity_matrix(:,i),MI_decreased);
%  if sum(Activity_matrix(:,i))==0
%  else
%     [h,p_MI_decrease_GA(i),stats]=fishertest(temp,'tail','both');
%     enrichment(i)=temp(2,2)./sum(temp(2,:))./(sum(temp(:,2))./sum(temp,'all'));
%  end
%     
% end
%     [p_MI_decrease_GA_corrected,c_alpha,h]=fwer_holmbonf(p_MI_GA,0.05); %Holm-Bonferroni correction
% 
%     p_MI_decrease_GA_transformed=-log10(p_MI_decrease_GA_corrected);
%  figure
%  hold on
%  p_MI_GA_transformed(strcmp(unique_activity,''))=[];
%  labels=unique_activity(~strcmp(unique_activity,''))
%  [values,I]=sort(p_MI_decrease_GA_transformed,'descend');
%     bar(1:length(values),values,'FaceColor',[.5,.5,.5]);
%     set(gca,'XTick',1:length(Activity_matrix(1,:)))
%     set(gca,'XTickLabel',labels(I))
%      set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     ylim([0,1.5])
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_decrease_GA_pvalues.png')
%  end
%  
%    sig=p_MI_decrease_GA_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(enrichment,'descend');
%     labels=unique_activity(I);
%     sigtemp=sig(I);
%     %sigtemp=sigtemp(1:10);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:length(values),values,'FaceColor',[.5,.5,.5]);
%     set(gca,'XTick',1:length(values))
%     set(gca,'XTickLabel',labels)
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,5])
%     scatter(x,values(sigtemp)+max(5)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_decrease_GA_enrichment.png')
%  end
%  %% Connection between MI and process
%  Average_MI_Process=zeros(length(Process_matrix(1,:)),1);
%  MI_Process=NaN(length(MI_compiled_noVP16),length(Process_matrix(1,:)));
%  figure
%  hold on
%  for i=1:length(Process_matrix(1,:))
%      Average_MI_Process(i)=nanmean(MI_compiled_noVP16(logical(Process_matrix(:,i))));
%      MI_Process(logical(Process_matrix(:,i)),i)=MI_compiled_noVP16(logical(Process_matrix(:,i)));
%      err(i)=nanstd(MI_compiled_noVP16(logical(Process_matrix(:,i))))./sqrt(length(MI_compiled_noVP16(logical(Process_matrix(:,i)))));
% 
%  end
%  [values,I]=sort(Average_MI_Process,'descend');
%  err=err(I);
%  labels=[];
%  labels=unique_Process(I);
%   labels(isnan(values))=[];
%  err(isnan(values))=[];
%  values(isnan(values))=[];
%       bar(1:length(values(values>=1.13)),values(values>=1.13),'FaceColor',[.8,.8,.99]);
%       bar(length(values(values>=1.13))+1:length(values(values>=.5)),values(values<=1.13 & values>=.5),'FaceColor',[.9,.9,.9]);
%       bar(length(values(values>=.5))+1:length(values),values(values<.5),'FaceColor',[1,.85,.9]);
%      errorbar(1:length(values),values,err,'linestyle','none','color','k','linewidth',1.5,'capsize',0);
%      set(gca,'XTick',1:length(values))
%     set(gca,'XTickLabel',labels)
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     ylim([0,1.5])
%     saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/CR screen/LogPlots/averageMI_Process.png')
%  figure
% [p_value_Process_MI_anova,~,stats]=anova1(MI_Process)
% multcompare_Process_MI=multcompare(stats)
%  
% %Fisher Exact test between GO and High MI 
%  enrichment=NaN(length(Process_matrix(1,:)),1);
% p_MI_Process=NaN(length(Process_matrix(1,:)),1);
% for i=1:length(Process_matrix(1,:))
%     
%     temp=crosstab(Process_matrix(:,i),MI_increased);
%  if sum(Process_matrix(:,i))==0
%  else
%     [h,p_MI_Process(i),stats]=fishertest(temp,'tail','both');
%     enrichment(i)=temp(2,2)./sum(temp(2,:))./(sum(temp(:,2))./sum(temp,'all'));
%  end
%     
% end
% [p_MI_Process_corrected,c_alpha,h]=fwer_holmbonf(p_MI_Process,0.05); %Holm-Bonferroni correction
% 
%     p_MI_Process_transformed=-log10(p_MI_Process_corrected);
%  figure
%  hold on
%  [values,I]=sort(p_MI_Process_transformed,'descend');
%     bar((1:10),values(1:10),'facecolor',[.9,.7,.7]);
%     labels=unique_Process(I);
%     set(gca,'XTick',1:10)
%     set(gca,'XTickLabel',labels(1:10))
%     set(gca,'XTickLabelRotation',45)
%     ylim([0,2.5])
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_Process_pvalues.png')
%  end
%     sig=p_MI_Process_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(enrichment,'descend');
%     labels=unique_Process(I);
%     sigtemp=sig(I);
%     sigtemp=sigtemp(1:16);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:16,values(1:16),'FaceColor',[.9,.7,.7]);
%     set(gca,'XTick',1:16)
%     set(gca,'XTickLabel',labels(1:16))
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,15])
%     scatter(x,values(sigtemp)+max(5)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_Process_enrichment.png')
%  end
%   %Fisher Exact test between GO and MI below 0.5
%   p_MI_decrease_Process=NaN(length(Process_matrix(1,:)),1);
%   enrichment=NaN(length(Process_matrix(1,:)),1);
%   diminishment=NaN(length(Function_matrix(1,:)),1);
% for i=1:length(Process_matrix(1,:))
%     temp2=crosstab(~Process_matrix(:,i),MI_decreased);
%     temp=crosstab(Process_matrix(:,i),MI_decreased);
%  if sum(Process_matrix(:,i))==0
%  else
%     [h,p_MI_decrease_Process(i),stats]=fishertest(temp,'tail','both');
%     [h,p_MI_decrease_Process2(i),stats]=fishertest(temp2,'tail','both');
%     enrichment(i)=temp(2,2)./sum(temp(2,:))./(sum(temp(:,2))./sum(temp,'all'));
%     diminishment(i)=temp2(2,2)./sum(temp2(2,:))./(sum(temp2(:,2))./sum(temp2,'all'));
%  end
%     
% end
%     [p_MI_decrease_Process_corrected,c_alpha,h]=fwer_holmbonf(p_MI_decrease_Process,0.05); %Holm-Bonferroni correction
%     [p_MI_decrease_Process2_corrected,c_alpha,h]=fwer_holmbonf(p_MI_decrease_Process2,0.05); %Holm-Bonferroni correction
%     p_MI_decrease_Process_transformed=-log10(p_MI_decrease_Process_corrected);
%  figure
%  hold on
%  [values,I]=sort(p_MI_decrease_Process_transformed,'descend');
%     bar((1:10),values(1:10),'facecolor',[.5,.5,.5]);
%     labels=unique_Process(I);
%     set(gca,'XTick',1:10)
%     set(gca,'XTickLabel',labels(1:10))
%     set(gca,'XTickLabelRotation',45)
%     ylim([0,2.5])
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_decrease_Process_pvalues.png')
%  end
%   sig=p_MI_decrease_Process_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(enrichment,'descend');
%     labels=unique_Process(I);
%     sigtemp=sig(I);
%     sigtemp=sigtemp(1:16);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:16,values(1:16),'FaceColor',[.5,.5,.5]);
%     set(gca,'XTick',1:16)
%     set(gca,'XTickLabel',labels(1:16))
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,12])
%     scatter(x,values(sigtemp)+max(5)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_decrease_Process_enrichment.png')
%  end
%   sig=p_MI_decrease_Process2_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(diminishment,'descend');
%     labels=unique_Process(I);
%     sigtemp=sig(I);
%     sigtemp=sigtemp(1:16);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:16,values(1:16),'FaceColor',[.5,.5,.5]);
%     set(gca,'XTick',1:16)
%     set(gca,'XTickLabel',labels(1:16))
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,12])
%     scatter(x,values(sigtemp)+max(5)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_decrease_Process_diminishment.png')
%  end
%  %% Connection between MI and function
%  Average_MI_Function=zeros(length(unique_Function),1);
%  MI_Function=NaN(length(MI_compiled_noVP16),length(unique_Function));
%  err=[];
%  figure
%  hold on
%  for i=1:length(unique_Function)
%      Average_MI_Function(i)=nanmean(MI_compiled_noVP16(logical(Function_matrix(:,i))));
%      MI_Function(logical(Function_matrix(:,i)),i)=MI_compiled_noVP16(logical(Function_matrix(:,i)));
%      err(i)=nanstd(MI_compiled_noVP16(logical(Function_matrix(:,i))))./sqrt(length(MI_compiled_noVP16(logical(Function_matrix(:,i)))));
% 
%  end
%  [values,I]=sort(Average_MI_Function,'descend');
%  err=err(I);
%  labels=[];
%  labels=unique_Function(I);
%   labels(isnan(values))=[];
%  err(isnan(values))=[];
%  values(isnan(values))=[];
%       bar(1:length(values(values>=1)),values(values>=1),'FaceColor',[.9,.7,.7]);
%       bar(length(values(values>=1))+1:length(values(values>=.5)),values(values<=1 & values>=.5),'FaceColor',[.9,.9,.9]);
%       bar(length(values(values>=.5))+1:length(values),values(values<.5),'FaceColor',[.5,.5,.5]);
%      errorbar(1:length(values),values,err,'linestyle','none','color','k','linewidth',2,'capsize',0);
%      set(gca,'XTick',1:length(values))
%     set(gca,'XTickLabel',labels)
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     ylim([0,1.5])
%     %saveas(gcf,'LogPlots/averageMI_Function.png')
%  figure
% [p_value_Function_MI_anova,~,stats]=anova1(MI_Function)
% multcompare_Function_MI=multcompare(stats)
% 
%  %Fisher Exact test between GO and high MI
%  
% p_MI_Function=NaN(length(Function_matrix(1,:)),1);
% enrichment=NaN(length(Function_matrix(1,:)),1);
% diminishment=NaN(length(Function_matrix(1,:)),1);
% for i=1:length(Function_matrix(1,:))
%     temp=crosstab(Function_matrix(:,i),MI_increased);
%     temp2=crosstab(~Function_matrix(:,i),MI_increased);
%  if sum(Function_matrix(:,i))==0
%  else
%     [h,p_MI_Function(i),stats]=fishertest(temp,'tail','both');
%     [h,p_MI_Function2(i),stats]=fishertest(temp,'tail','both');
%     enrichment(i)=temp(2,2)./sum(temp(2,:))./(sum(temp(:,2))./sum(temp,'all'));
%     diminishment(i)=temp2(2,2)./sum(temp2(2,:))./(sum(temp2(:,2))./sum(temp2,'all'));
%  end
%     
% end
% [p_MI_Function_corrected,c_alpha,h]=fwer_holmbonf(p_MI_Function,0.05); %Holm-Bonferroni correction
% [p_MI_Function2_corrected,c_alpha,h]=fwer_holmbonf(p_MI_Function2,0.05);
%     p_MI_Function_transformed=-log10(p_MI_Function_corrected);
%  figure
%  hold on
%  [values,I]=sort(p_MI_Function_transformed,'descend');
%     bar((1:10),values(1:10),'facecolor',[.9,.7,.7]);
%     labels=unique_Function(I);
%     set(gca,'XTick',1:10)
%     set(gca,'XTickLabel',labels(1:10))
%     set(gca,'XTickLabelRotation',45)
%     ylim([0,2])
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_Function_pvalues.png')
%  end
%  sig=p_MI_Function_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(enrichment,'descend');
%     labels=unique_Function(I);
%     sigtemp=sig(I);
%     sigtemp=sigtemp(1:16);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:16,values(1:16),'FaceColor',[.9,.7,.7]);
%     set(gca,'XTick',1:16)
%     set(gca,'XTickLabel',labels(1:16))
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,15])
%     scatter(x,values(sigtemp)+max(5)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_Function_enrichment.png')
%  end
%  sig=p_MI_Function2_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(diminishment,'descend');
%     labels=unique_Function(I);
%     sigtemp=sig(I);
%     sigtemp=sigtemp(1:16);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:16,values(1:16),'FaceColor',[.9,.7,.7]);
%     set(gca,'XTick',1:16)
%     set(gca,'XTickLabel',labels(1:16))
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,12])
%     scatter(x,values(sigtemp)+max(5)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_Function_diminishment.png')
%  end
%   %Fisher Exact test between Function and MI below 0.5
%   p_MI_decrease_Function=NaN(length(Function_matrix(1,:)),1);
%   enrichment=NaN(length(Function_matrix(1,:)),1);
%   diminishment=NaN(length(Function_matrix(1,:)),1);
% for i=1:length(Function_matrix(1,:))
%     temp2=crosstab(~Function_matrix(:,i),MI_decreased);
%     temp=crosstab(Function_matrix(:,i),MI_decreased);
%  if sum(Function_matrix(:,i))==0
%  else
%     [h,p_MI_decrease_Function(i),stats]=fishertest(temp,'tail','both');
%     [h,p_MI_decrease_Function2(i),stats]=fishertest(temp,'tail','both');
%     enrichment(i)=temp(2,2)./sum(temp(2,:))./(sum(temp(:,2))./sum(temp,'all'));
%     diminishment(i)=temp2(2,2)./sum(temp2(2,:))./(sum(temp2(:,2))./sum(temp2,'all'));
%  end
%     
% end
%     [p_MI_decrease_Function_corrected,c_alpha,h]=fwer_holmbonf(p_MI_decrease_Function,0.05); %Holm-Bonferroni correction
% [p_MI_decreased_Function2_corrected,c_alpha,h]=fwer_holmbonf(p_MI_decrease_Function2,0.05); %Holm-Bonferroni correction
%     p_MI_decrease_Function_transformed=-log10(p_MI_decrease_Function_corrected);
%  figure
%  hold on
%  [values,I]=sort(p_MI_decrease_Function_transformed,'descend');
%     bar((1:10),values(1:10),'facecolor',[.5,.5,.5]);
%     labels=unique_Function(I);
%     set(gca,'XTick',1:10)
%     set(gca,'XTickLabel',labels(1:10))
%     set(gca,'XTickLabelRotation',45)
%     ylim([0,2])
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_decrease_Function_pvalues.png')
%  end
%  sig=p_MI_decrease_Function_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(enrichment,'descend');
%     labels=unique_Function(I);
%     sigtemp=sig(I);
%     sigtemp=sigtemp(1:16);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:16,values(1:16),'FaceColor',[.5,.5,.5]);
%     set(gca,'XTick',1:16)
%     set(gca,'XTickLabel',labels(1:16))
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,12])
%     scatter(x,values(sigtemp)+max(5)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_decrease_Function_enrichment.png')
%  end
%  sig=p_MI_decreased_Function2_corrected<=0.05;
%  figure
%  set(gcf,'position',[440,225,560,570]);
%  
%     [values,I]=sort(diminishment,'descend');
%     labels=unique_Function(I);
%     sigtemp=sig(I);
%     sigtemp=sigtemp(1:16);
%     x=find(sigtemp==1);
%     
%      hold on
%     bar(1:16,values(1:16),'FaceColor',[.5,.5,.5]);
%     set(gca,'XTick',1:16)
%     set(gca,'XTickLabel',labels(1:16))
%     set(gca,'XTickLabelRotation',45)
%     set(gca,'fontsize',8)
%     
%     ylim([0,12])
%     scatter(x,values(sigtemp)+max(5)./25,'*k')
%  
%  if Save_condition==1
%      %saveas(gcf,'LogPlots/MI_decrease_Function_diminishment.png')
%  end
 %% Example CRs for each cluster
 indices=any(~isnan(fc3_reordered),2);
 fc3_filtered=fc3_reordered(indices,:);
 CR_filtered=CR_compiled(indices,:);
 CV_dark_filtered=CV_dark_reordered(indices,:);
 CV_filtered=CV_reordered(indices,:);
 cdis_fc3_filtered=cdis_reordered(indices,:);
 t=table(CR_filtered,fc3_filtered(:,1),fc3_filtered(:,2),fc3_filtered(:,3),CV_dark_filtered(:,2),CV_filtered(:,1),CV_filtered(:,2),CV_filtered(:,3),'VariableNames',{'Gene','FC1','FC2','FC3','CV0','CV1','CV2','CV3'});
Stats_fc3=grpstats(t,{'Gene'},{'mean','std','meanci','sem'});
max_CI=max(abs([(Stats_fc3.meanci_FC1(:,1)-Stats_fc3.meanci_FC1(:,2))./Stats_fc3.mean_FC1, (Stats_fc3.meanci_FC2(:,1)-Stats_fc3.meanci_FC2(:,2))./Stats_fc3.mean_FC2, (Stats_fc3.meanci_FC3(:,1)-Stats_fc3.meanci_FC3(:,2))./Stats_fc3.mean_FC3])');
table2=Stats_fc3;

[distance,cdis_fc3_average]=pdist2(ctrs,[table2.mean_FC1,table2.mean_FC2,table2.mean_FC3],'correlation','Smallest',1);
slope1=abs(table2.mean_FC1-table2.mean_FC2);
slope2=abs(table2.mean_FC2-table2.mean_FC3);
fc3_range=max([table2.mean_FC1 table2.mean_FC2 table2.mean_FC3]')-min([table2.mean_FC1 table2.mean_FC2 table2.mean_FC3]');
table2=[table(cdis_fc3_average', slope1, slope2,fc3_range','VariableNames',{'Cluster','Slope1','Slope2','RangeofFC3'}) table2]; 
figure
hold on
gene=string(zeros(height(table2),1));
%table2=table2(distance'<.3,:);
%table2=table2(table2.RangeofFC3<1,:);
%table2=table2(max_CI<2,:);
k=1;
for j=1:5
    temp_table=table2(table2.Cluster==j,:);
    
   
if j==1

    ngrps=0;
elseif j==2
%     Indices=find(min(temp_table.mean_FC2));
%     temp_table=temp_table(Indices,:);
ngrps=[10];
elseif j==3
%     Indices=find(max(temp_table.mean_FC2));
%     temp_table=temp_table(Indices,:);
ngrps=[1];
elseif j==4
%     Indices=find(min(abs(temp_table.meanci_FC3(:,2)-temp_table.mean_FC3)));
ngrps=[1];
elseif j==5
% temp_table=Stats_fc3(strcmp(Stats_fc3.Gene,'hda3'),:);
ngrps=1;
end
% % temp_table=temp_table(Indices,:);
 ngrps=height(temp_table);
    if ngrps==0
       
    else
 subplot(3,2,j)
 hold on
for i=1:ngrps
errorbar(Freq_1,[temp_table.mean_FC1(i) temp_table.mean_FC2(i) temp_table.mean_FC3(i)],[temp_table.sem_FC1(i) temp_table.sem_FC2(i) temp_table.sem_FC3(i)],'Linewidth',2,'Color',face_color(j,:),'capsize',0)
ylim([-2,4])
gene(k)=temp_table.Gene(i);
k=k+1;
end

    end
end
if experiment==1
genes={'ies4','rxt2','ecm5','ecm5','gcn5'};
elseif experiment==2
    genes={'ies4','ecm5','gcn5','hda3','hda3'};
elseif experiment==3
  genes={'ies4','rxt2','ecm5','ecm5','gcn5'};  
end

% test=[fc3_comp(strcmp(CR_comp, 'cdc73'),:); fc3_comp(strcmp(CR_comp,'arp8')); fc3_comp(strcmp(CR_comp,'ecm5'));fc3_comp(strcmp(CR_comp,'taf9'));fc3_comp(strcmp(CR_comp,'hda3'))];
% testCR=[CR_comp(strcmp(CR_comp, 'cdc73'),:); CR_comp(strcmp(CR_comp,'arp8')); CR_comp(strcmp(CR_comp,'ecm5')); CR_comp(strcmp(CR_comp,'taf9'));CR_comp(strcmp(CR_comp,'hda3'))];
% testCond=[cond_comp(strcmp(CR_comp, 'cdc73'),:); cond_comp(strcmp(CR_comp,'arp8')); cond_comp(strcmp(CR_comp,'ecm5')); cond_comp(strcmp(CR_comp,'taf9'));cond_comp(strcmp(CR_comp,'hda3'))];
% %[p,h,stats]=anova1(test(strcmp(testCR,'hda3')),testCond(strcmp(testCR,'hda3')));
% %results=multcompare(stats)
figure
set(gcf,'position',[440,225,560,570]);
hold on

for i=1:5
indices=strcmp(table2.Gene,genes(i));
subplot(3,2,table2.Cluster(genes(i)))
hold on
% bar(3*(i-1)+1:3*i,[table2.mean_FC1(indices) table2.mean_FC2(indices) table2.mean_FC3(indices)],'facecolor',face_color(table2.Cluster(genes(i)),:))
%     errorbar(3*(i-1)+1:3*i,[table2.mean_FC1(indices) table2.mean_FC2(indices) table2.mean_FC3(indices),],...
%     [table2.meanci_FC1(indices,2)-table2.mean_FC1(indices) table2.meanci_FC2(indices,2)-table2.mean_FC2(indices) table2.meanci_FC3(indices,2)-table2.mean_FC3(indices)],'Linewidth',1,'Color','k','Linestyle','none','capsize',2)
values=fc3_overall_notmasked(strcmp(CR_80_not_masked,genes(i)),:);
% xvalues=(3*(i-1)+1:3*i);

errorbar(Freq_1,[table2.mean_FC1(indices) table2.mean_FC2(indices) table2.mean_FC3(indices)],...
    [table2.sem_FC1(indices) table2.sem_FC2(indices) table2.sem_FC3(indices)],'Linewidth',2,'Color',face_color(table2.Cluster(genes(i)),:),'capsize',0)
plot(Freq_1,values','.','color',[0.5,.5,.5],'Markersize',10)
%xlabel('Frequency (sec^-^1)')
%ylabel('log_1_0(Fold change)')
%ylim([-1,4])
legend(genes(i),'location','southeast')
end
set(gca,'Xtick',1:5)
set(gca,'linewidth',1,'ticklength',[0.005 0.005])
%saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/CR screen/LogPlots/FC_exampleCR.png')
xlabel('Frequency (sec^-^1)')
ylabel('log_1_0(Fold change)')



figure
limits=[0,3;0,3;0,5;-1.5,3;-1,1];
set(gcf,'position',[440,225,560,570]);
for i=1:5
indices=strcmp(table2.Gene,genes(i));
subplot(3,2,table2.Cluster(genes(i)))
hold on
bar(1:3,[table2.mean_FC1(indices) table2.mean_FC2(indices) table2.mean_FC3(indices)],'facecolor',face_color(table2.Cluster(genes(i)),:))
    errorbar(1:3,[table2.mean_FC1(indices) table2.mean_FC2(indices) table2.mean_FC3(indices),],...
    [table2.sem_FC1(indices) table2.sem_FC2(indices) table2.sem_FC3(indices)],'Linewidth',2,'Color','k','Linestyle','none','Capsize',0)
values=fc3_overall_notmasked(strcmp(CR_80_not_masked,genes(i)),:);
% xvalues=(3*(i-1)+1:3*i);

%errorbar(Freq_1,[table2.mean_FC1(indices) table2.mean_FC2(indices) table2.mean_FC3(indices)],...
    %[table2.meanci_FC1(indices,2)-table2.mean_FC1(indices) table2.meanci_FC2(indices,2)-table2.mean_FC2(indices) table2.meanci_FC3(indices,2)-table2.mean_FC3(indices)],'Linewidth',2,'Color',face_color(table2.Cluster(genes(i)),:))
plot(1:3,values','.','color',[0.5,.5,.5],'Markersize',10)
legend(genes(i),'location','north')
ylim(limits(i,:));
set(gca,'Xtick',1:3)
set(gca,'XTickLabel',{''})
set(gca,'linewidth',1.25,'ticklength',[0.02 0.02])
end

%     set(gca,'fontsize',8)

    %saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/CR screen/LogPlots/FC_CRexample_bar.png')
figure
indices=strcmp(table2.Gene,'JY145');
bar(1:3,[table2.mean_FC1(indices) table2.mean_FC2(indices) table2.mean_FC3(indices)],'facecolor',face_color(table2.Cluster('JY145'),:))
hold on   
errorbar(1:3,[table2.mean_FC1(indices) table2.mean_FC2(indices) table2.mean_FC3(indices),],...
    [table2.sem_FC1(indices) table2.sem_FC2(indices) table2.sem_FC3(indices)],'Linewidth',2,'Color','k','Linestyle','none','Capsize',0)
values=fc3_overall_notmasked(strcmp(CR_80_not_masked,'JY145'),:);
% xvalues=(3*(i-1)+1:3*i);

%errorbar(Freq_1,[table2.mean_FC1(indices) table2.mean_FC2(indices) table2.mean_FC3(indices)],...
    %[table2.meanci_FC1(indices,2)-table2.mean_FC1(indices) table2.meanci_FC2(indices,2)-table2.mean_FC2(indices) table2.meanci_FC3(indices,2)-table2.mean_FC3(indices)],'Linewidth',2,'Color',face_color(table2.Cluster(genes(i)),:))
plot(1:3,values','.','color',[0.5,.5,.5],'Markersize',10)
legend('VP16 only','location','north')
%ylim(limits(i,:));
set(gca,'Xtick',1:3)
set(gca,'XTickLabel',{''})
set(gca,'linewidth',1.25,'ticklength',[0.02 0.02])
figure
hold on
values=[];

test=NaN.*ones(4,5);
for i=1:number_of_clusters
    values=MI_compiled(strcmp(CR_compiled,genes(i)));
    bar(table2.Cluster(genes(i)),nanmean(values),'facecolor',face_color(table2.Cluster(genes(i)),:))
    plot(table2.Cluster(genes(i)),values,'.','color',[.5,.5,.5],'Markersize',15)
    errorbar(table2.Cluster(genes(i)),nanmean(values),std(values)./sqrt(length(values)),'color','k','linewidth',2,'capsize',0)
    test(1:length(values),table2.Cluster(genes(i)))=values;
    
end
set(gca,'XTickLabel',{''})
set(gca,'Xtick',1:number_of_clusters)
ylabel('MI_F_M')
set(gca,'linewidth',1.5,'ticklength',[0.02 0.02])
    %saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/CR screen/LogPlots/FC_CRexample_MI.png')

figure
[p,tbl,stats]=anova1(test);
multcompare(stats)
%CV for examples
figure
set(gcf,'position',[440,225,560,570]);
for i=1:5
indices=strcmp(table2.Gene,genes(i));
subplot(3,2,table2.Cluster(genes(i)))
hold on
bar(1:4,[table2.mean_CV0(indices) table2.mean_CV1(indices) table2.mean_CV2(indices) table2.mean_CV3(indices)],'facecolor',face_color(table2.Cluster(genes(i)),:))
    errorbar(1:4,[table2.mean_CV0(indices) table2.mean_CV1(indices) table2.mean_CV2(indices) table2.mean_CV3(indices)],...
    [table2.sem_CV0(indices) table2.sem_CV1(indices) table2.sem_CV2(indices) table2.sem_CV3(indices)],'Linewidth',2,'Color','k','Linestyle','none','capsize',0)
values=[CV_dark_filtered(strcmp(CR_filtered,genes(i)),2) CV_filtered(strcmp(CR_filtered,genes(i)),:)];
% xvalues=(3*(i-1)+1:3*i);

%errorbar(Freq_1,[table2.mean_FC1(indices) table2.mean_FC2(indices) table2.mean_FC3(indices)],...
    %[table2.meanci_FC1(indices,2)-table2.mean_FC1(indices) table2.meanci_FC2(indices,2)-table2.mean_FC2(indices) table2.meanci_FC3(indices,2)-table2.mean_FC3(indices)],'Linewidth',2,'Color',face_color(table2.Cluster(genes(i)),:))
plot(1:4,values','.','color',[0.5,.5,.5],'Markersize',10)
legend(genes(i),'location','north')
ylim([-.5,.5])
set(gca,'Xtick',1:4)
set(gca,'XTickLabel',{''})
set(gca,'linewidth',1.25,'ticklength',[0.02 0.02])
end

%     set(gca,'fontsize',8)

    %saveas(gcf,'/Volumes/GoogleDrive/My Drive/Shared folder/YeastLightExperiments/CR screen/LogPlots/FC_CRexample_CVbar.png')

test=[CV_comp(strcmp(CR_comp, 'cdc73'),:); CV_comp(strcmp(CR_comp,'arp8')); CV_comp(strcmp(CR_comp,'ecm5'));CV_comp(strcmp(CR_comp,'taf9'));CV_comp(strcmp(CR_comp,'hda3'))];
testCR=[CR_comp(strcmp(CR_comp, 'cdc73'),:); CR_comp(strcmp(CR_comp,'arp8')); CR_comp(strcmp(CR_comp,'ecm5')); CR_comp(strcmp(CR_comp,'taf9'));CR_comp(strcmp(CR_comp,'hda3'))];
testCond=[cond_comp(strcmp(CR_comp, 'cdc73'),:); cond_comp(strcmp(CR_comp,'arp8')); cond_comp(strcmp(CR_comp,'ecm5')); cond_comp(strcmp(CR_comp,'taf9'));cond_comp(strcmp(CR_comp,'hda3'))];
[p,h,stats]=anovan(test,{testCR,testCond},'varnames',{'CR','Cond'},'model','interaction');
results=multcompare(stats,'Dimension',[1 2])
%% Histograms for example CRs for FC clusters
genes={'cdc73','arp8','ecm5','taf9','hda3'};
[~,d]=xlsread('Excel Files/CRMIsamples.xlsx');
CR_MIsamples=d(3:end,1);
Cond0_folder=d(3:end,3);
Cond0_file=d(3:end,4);
Cond1_folder=d(3:end,5);
Cond1_file=d(3:end,6);
Cond2_folder=d(3:end,7);
Cond2_file=d(3:end,8);
Cond3_folder=d(3:end,9);
Cond3_file=d(3:end,10);
%
replicate=[2,3,2,3,3];
figure
set(gcf,'position',[440,225,560,570]);
for l=1:5
    subplot(3,2,table2.Cluster(genes(l)))
    hold on
    xlim([0.01,3000])
indices=find(strcmp(CR_MIsamples,genes(l))==1);
indices2=find(strcmp(CR_80_not_masked,genes(l))==1);
mCherry=[];
for j=replicate(l)
name=char(strcat('Flow cytometry data/',Cond0_folder(indices(j)),Cond0_file(indices(j))));
Y1=csvread(name,1,0);
mCherry2=Y1(:,6);
FSC_A=Y1(:,1);
FSC_H=Y1(:,2);
SSC_A=Y1(:,3);

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
Both_Gates=(Gate+Gate2==2);
mCherry1=mCherry2(Both_Gates)./FSC_A(Both_Gates);%correct for size
%mCherry1=mCherry(Both_Gates);
mCherry=[mCherry;mCherry1];
end
edges=logspace(-2,4,50);
%edges=[0,edges];
histogram(mCherry,edges,'normalization','probability','facecolor',face_color(1,:),'facealpha',1)

mCherry=[];
for j=replicate(l)
name=char(strcat('Flow cytometry data/',Cond1_folder(indices(j)),Cond1_file(indices(j))));
Y2=csvread(name,1,0);
mCherry2=Y2(:,6);
FSC_A=Y2(:,1);
FSC_H=Y2(:,2);
SSC_A=Y2(:,3);
log_FSC=log10(FSC_A);
log_FSC_H=log10(FSC_H);
log_SSC=log10(SSC_A);
%med_FSC=median(FSC_A);
%med_SSC=median(SSC_A);
med_FSC=median(log_FSC);
med_SSC=median(log_SSC);
k=log_SSC./log_FSC;
med_k=median(k);
%Gate=FSC_A>med_FSC-radius & FSC_A<med_FSC+radius;%Gate for FSC-A<=radius au
Gate=(log_FSC-med_FSC).^2+(k-med_k).^2<=radius^2;
Gate2=(log_FSC_H)./log_FSC>median(log_FSC_H./log_FSC)-.1 & (log_FSC_H)./log_FSC<median(log_FSC_H./log_FSC)+.1;
Both_Gates=(Gate+Gate2==2);
mCherry1=mCherry2(Both_Gates)./FSC_A(Both_Gates);%correct for size
%mCherry1=mCherry(Both_Gates);

mCherry=[mCherry; mCherry1];
end

histogram(mCherry,edges,'normalization','probability','facecolor',face_color(2,:),'facealpha',1)

mCherry=[];
for j=replicate(l)
name=char(strcat('Flow cytometry data/',Cond2_folder(indices(j)),Cond2_file(indices(j))));
Y3=csvread(name,1,0);
mCherry2=Y3(:,6);
FSC_A=Y3(:,1);
FSC_H=Y3(:,2);
SSC_A=Y3(:,3);
log_FSC=log10(FSC_A);
log_FSC_H=log10(FSC_H);
log_SSC=log10(SSC_A);
%med_FSC=median(FSC_A);
%med_SSC=median(SSC_A);
med_FSC=median(log_FSC);
med_SSC=median(log_SSC);
k=log_SSC./log_FSC;
med_k=median(k);
%Gate=FSC_A>med_FSC-radius & FSC_A<med_FSC+radius;%Gate for FSC-A<=radius au
Gate=(log_FSC-med_FSC).^2+(k-med_k).^2<=radius^2;
Gate2=(log_FSC_H)./log_FSC>median(log_FSC_H./log_FSC)-.1 & (log_FSC_H)./log_FSC<median(log_FSC_H./log_FSC)+.1;
Both_Gates=(Gate+Gate2==2);
mCherry1=mCherry2(Both_Gates)./FSC_A(Both_Gates);%correct for size
%mCherry1=mCherry(Both_Gates);
mCherry=[mCherry1;mCherry];
end
histogram(mCherry,edges,'normalization','probability','facecolor',face_color(3,:),'facealpha',1)

mCherry=[];
for j=replicate(l)
name=char(strcat('Flow cytometry data/',Cond3_folder(indices(j)),Cond3_file(indices(j))));
Y4=csvread(name,1,0);
mCherry2=Y4(:,6);
FSC_A=Y4(:,1);
FSC_H=Y4(:,2);
SSC_A=Y4(:,3);
log_FSC=log10(FSC_A);
log_FSC_H=log10(FSC_H);
log_SSC=log10(SSC_A);
%med_FSC=median(FSC_A);
%med_SSC=median(SSC_A);
med_FSC=median(log_FSC);
med_SSC=median(log_SSC);
k=log_SSC./log_FSC;
med_k=median(k);
%Gate=FSC_A>med_FSC-radius & FSC_A<med_FSC+radius;%Gate for FSC-A<=radius au
Gate=(log_FSC-med_FSC).^2+(k-med_k).^2<=radius^2;
Gate2=(log_FSC_H)./log_FSC>median(log_FSC_H./log_FSC)-.1 & (log_FSC_H)./log_FSC<median(log_FSC_H./log_FSC)+.1;
Both_Gates=(Gate+Gate2==2);
mCherry1=mCherry2(Both_Gates)./FSC_A(Both_Gates);%correct for size
%mCherry1=mCherry(Both_Gates);
mCherry=[mCherry;mCherry1];
end
histogram(mCherry,edges,'normalization','probability','facecolor',face_color(4,:),'facealpha',1)

xlim([0,1000])
ylim([0,.3])
set(gca,'XScale','log')
end
%saveas(gcf,'LogPlots/histograms_CRexamples.png')
%%  %% Example CRs for each CV cluster
indices=~any(isnan(CV_reordered),2);
 CV_filtered=CV_reordered(indices,:);
 CR_filtered=CR_compiled(indices,:);
 cdis_CV_filtered=cdis_CV_reordered(indices,:);
 t=table(CR_filtered,CV_filtered(:,1),CV_filtered(:,2),CV_filtered(:,3),'VariableNames',{'Gene','CV1','CV2','CV3'});
Stats_CV=grpstats(t,{'Gene'},{'mean','sem','meanci'});
max_CI=max(abs([(Stats_CV.meanci_CV1(:,1)-Stats_CV.meanci_CV1(:,2))./Stats_CV.mean_CV1, (Stats_CV.meanci_CV2(:,1)-Stats_CV.meanci_CV2(:,2))./Stats_CV.mean_CV2, (Stats_CV.meanci_CV3(:,1)-Stats_CV.meanci_CV3(:,2))./Stats_CV.mean_CV3])');
table2=Stats_CV;

[~,cdis_CV_average]=pdist2(ctrsCV,[table2.mean_CV1,table2.mean_CV2,table2.mean_CV3],'correlation','Smallest',1);
slope1=abs(table2.mean_CV1-table2.mean_CV2);
slope2=abs(table2.mean_CV2-table2.mean_CV3);
CV_range=max([table2.mean_CV1 table2.mean_CV2 table2.mean_CV3]')-min([table2.mean_CV1 table2.mean_CV2 table2.mean_CV3]');
table2=[table(cdis_CV_average', slope1, slope2,CV_range','VariableNames',{'Cluster','Slope1','Slope2','RangeofCV'}) table2]; 
figure
hold on
gene=string(zeros(height(table2),1));
%table2=table2(table2.RangeofFC3<4,:);
% k=1;
% for j=1:5
%     temp_table=table2(table2.Cluster==j,:);
%     
%    
% if j==1
% 
%     ngrps=[6];
% elseif j==2
% %     Indices=find(min(temp_table.mean_FC2));
% %     temp_table=temp_table(Indices,:);
% ngrps=[4];
% elseif j==3
% %     Indices=find(max(temp_table.mean_FC2));
% %     temp_table=temp_table(Indices,:);
% ngrps=[5];
% elseif j==4
% %     Indices=find(min(abs(temp_table.meanci_FC3(:,2)-temp_table.mean_FC3)));
% ngrps=[6];
% elseif j==5
% ngrps=[8];
% end
% % % temp_table=temp_table(Indices,:);
% ngrps=height(temp_table);
%     if ngrps==0
%        
%     else
%         
% for i=ngrps
% errorbar(Freq_1,[temp_table.mean_CV1(i) temp_table.mean_CV2(i) temp_table.mean_CV3(i)],[temp_table.meanci_CV1(i,2)-temp_table.mean_CV1(i) temp_table.meanci_CV2(i,2)-temp_table.mean_CV2(i) temp_table.meanci_CV3(i,2)-temp_table.mean_CV3(i)],'Linewidth',3,'Color',face_color(j,:))
% ylim([1,3])
% gene(k)=temp_table.Gene(i);
% k=k+1;
% end
% 
%     end
% end
% gene(strcmp(gene,'0'))=[];
% legend(gene,'location','northeastoutside')
% xlabel('Frequency (sec^-^1)')
% ylabel('log_1_0(CV)')
genes={'bur2','cdc73','swc5','hir3','taf1'};
figure
test=[CV_comp(strcmp(CR_comp,'rxt3'),:); CV_comp(strcmp(CR_comp, 'cdc73'),:); CV_comp(strcmp(CR_comp,'swc5')); CV_comp(strcmp(CR_comp,'hir3'))];
testCR=[CR_comp(strcmp(CR_comp,'rxt3'),:); CR_comp(strcmp(CR_comp, 'cdc73'),:); CR_comp(strcmp(CR_comp,'swc5')); CR_comp(strcmp(CR_comp,'hir3'))];
testCond=[cond_comp(strcmp(CR_comp,'rxt3'),:); cond_comp(strcmp(CR_comp, 'cdc73'),:); cond_comp(strcmp(CR_comp,'swc5')); cond_comp(strcmp(CR_comp,'hir3'))];
[p,h,stats]=anovan(test,{testCR,testCond},'varnames',{'CR','Cond'},'model','interaction');
results=multcompare(stats,'Dimension',[1 2])
figure
set(gcf,'position',[440,225,560,570]);
hold on

for i=1:5
indices=strcmp(table2.Gene,genes(i));
subplot(3,2,table2.Cluster(genes(i)))
hold on
 bar(1:3,[table2.mean_CV1(indices) table2.mean_CV2(indices) table2.mean_CV3(indices)],'facecolor',face_color_CV(table2.Cluster(genes(i)),:))
 errorbar(1:3,[table2.mean_CV1(indices) table2.mean_CV2(indices) table2.mean_CV3(indices),],...
     [table2.sem_CV1(indices) table2.sem_CV2(indices) table2.sem_CV3(indices)],'Linewidth',2,'Color','k','Linestyle','none','capsize',0)
values=CV_reordered(strcmp(CR_compiled,genes(i)),:);
% xvalues=(3*(i-1)+1:3*i);

%errorbar(Freq_1,[table2.mean_CV1(indices) table2.mean_CV2(indices) table2.mean_CV3(indices)],...
    %[table2.meanci_CV1(indices,2)-table2.mean_CV1(indices) table2.meanci_CV2(indices,2)-table2.mean_CV2(indices) table2.meanci_CV3(indices,2)-table2.mean_CV3(indices)],'Linewidth',2,'Color',face_color(table2.Cluster(genes(i)),:))
plot(1:3,values','.','color',[0.5,.5,.5],'Markersize',10)
%xlabel('Frequency (sec^-^1)')
%ylabel('log_1_0(CV)')
ylim([-.5,.5])
legend(genes(i),'location','north')
set(gca,'xtick',1:3)
set(gca,'XTickLabel',{''})
end

if Save_condition==1
    %saveas(gcf,'LogPlots/CV_CRexample_bar.png')
end
%Average MI for each example
figure
hold on
values=[];
for i=1:number_of_clusters
    values=MI_compiled(strcmp(CR_compiled,genes(i)));
    bar(table2.Cluster(genes(i)),nanmean(values),'facecolor',face_color(table2.Cluster(genes(i)),:))
    plot(table2.Cluster(genes(i)),values,'.','color',[.5,.5,.5],'Markersize',10)
    errorbar(table2.Cluster(genes(i)),nanmean(values),std(values)./sqrt(length(values)),'color','k','capsize',0,'linewidth',2)
end
ylabel('MI_F_M')
set(gca,'XTickLabel',{''})
if Save_condition==1
    %saveas(gcf,'LogPlots/CV_CRexample_MI.png')
end
%% %% Histograms for example CRs for CV clusters
genes={'bur2','cdc73','swc5','hir3','taf1'};
[~,d]=xlsread('Excel Files/CRMIsamples.xlsx');
CR_MIsamples=d(3:end,1);
Cond0_folder=d(3:end,3);
Cond0_file=d(3:end,4);
Cond1_folder=d(3:end,5);
Cond1_file=d(3:end,6);
Cond2_folder=d(3:end,7);
Cond2_file=d(3:end,8);
Cond3_folder=d(3:end,9);
Cond3_file=d(3:end,10);
%
replicate=[2,3,1,3,3];
figure
set(gcf,'position',[440,225,560,570]);
for l=1:5
    subplot(3,2,table2.Cluster(genes(l)))
    hold on
    xlim([0.01,3000])
indices=find(strcmp(CR_MIsamples,genes(l))==1);
indices2=find(strcmp(CR_80_not_masked,genes(l))==1);
mCherry=[];
for j=replicate(l)
name=char(strcat('Flow cytometry data/',Cond0_folder(indices(j)),Cond0_file(indices(j))));
Y1=csvread(name,1,0);
mCherry2=Y1(:,6);
FSC_A=Y1(:,1);
FSC_H=Y1(:,2);
SSC_A=Y1(:,3);

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
Both_Gates=(Gate+Gate2==2);
mCherry1=mCherry2(Both_Gates)./FSC_A(Both_Gates);%correct for size
%mCherry1=mCherry(Both_Gates);
mCherry=[mCherry;mCherry1];
end
edges=logspace(-2,4,50);
%edges=[0,edges];
histogram(mCherry,edges,'normalization','probability','facecolor',face_color(1,:))

mCherry=[];
for j=replicate(l)
name=char(strcat('Flow cytometry data/',Cond1_folder(indices(j)),Cond1_file(indices(j))));
Y2=csvread(name,1,0);
mCherry2=Y2(:,6);
FSC_A=Y2(:,1);
FSC_H=Y2(:,2);
SSC_A=Y2(:,3);
log_FSC=log10(FSC_A);
log_FSC_H=log10(FSC_H);
log_SSC=log10(SSC_A);
%med_FSC=median(FSC_A);
%med_SSC=median(SSC_A);
med_FSC=median(log_FSC);
med_SSC=median(log_SSC);
k=log_SSC./log_FSC;
med_k=median(k);
%Gate=FSC_A>med_FSC-radius & FSC_A<med_FSC+radius;%Gate for FSC-A<=radius au
Gate=(log_FSC-med_FSC).^2+(k-med_k).^2<=radius^2;
Gate2=(log_FSC_H)./log_FSC>median(log_FSC_H./log_FSC)-.1 & (log_FSC_H)./log_FSC<median(log_FSC_H./log_FSC)+.1;
Both_Gates=(Gate+Gate2==2);
mCherry1=mCherry2(Both_Gates)./FSC_A(Both_Gates);%correct for size
%mCherry1=mCherry(Both_Gates);

mCherry=[mCherry; mCherry1];
end

histogram(mCherry,edges,'normalization','probability','facecolor',face_color(2,:))

mCherry=[];
for j=replicate(l)
name=char(strcat('Flow cytometry data/',Cond2_folder(indices(j)),Cond2_file(indices(j))));
Y3=csvread(name,1,0);
mCherry2=Y3(:,6);
FSC_A=Y3(:,1);
FSC_H=Y3(:,2);
SSC_A=Y3(:,3);
log_FSC=log10(FSC_A);
log_FSC_H=log10(FSC_H);
log_SSC=log10(SSC_A);
%med_FSC=median(FSC_A);
%med_SSC=median(SSC_A);
med_FSC=median(log_FSC);
med_SSC=median(log_SSC);
k=log_SSC./log_FSC;
med_k=median(k);
%Gate=FSC_A>med_FSC-radius & FSC_A<med_FSC+radius;%Gate for FSC-A<=radius au
Gate=(log_FSC-med_FSC).^2+(k-med_k).^2<=radius^2;
Gate2=(log_FSC_H)./log_FSC>median(log_FSC_H./log_FSC)-.1 & (log_FSC_H)./log_FSC<median(log_FSC_H./log_FSC)+.1;
Both_Gates=(Gate+Gate2==2);
mCherry1=mCherry2(Both_Gates)./FSC_A(Both_Gates);%correct for size
%mCherry1=mCherry(Both_Gates);
mCherry=[mCherry1;mCherry];
end
histogram(mCherry,edges,'normalization','probability','facecolor',face_color(3,:))

mCherry=[];
for j=replicate(l)
name=char(strcat('Flow cytometry data/',Cond3_folder(indices(j)),Cond3_file(indices(j))));
Y4=csvread(name,1,0);
mCherry2=Y4(:,6);
FSC_A=Y4(:,1);
FSC_H=Y4(:,2);
SSC_A=Y4(:,3);
log_FSC=log10(FSC_A);
log_FSC_H=log10(FSC_H);
log_SSC=log10(SSC_A);
%med_FSC=median(FSC_A);
%med_SSC=median(SSC_A);
med_FSC=median(log_FSC);
med_SSC=median(log_SSC);
k=log_SSC./log_FSC;
med_k=median(k);
%Gate=FSC_A>med_FSC-radius & FSC_A<med_FSC+radius;%Gate for FSC-A<=radius au
Gate=(log_FSC-med_FSC).^2+(k-med_k).^2<=radius^2;
Gate2=(log_FSC_H)./log_FSC>median(log_FSC_H./log_FSC)-.1 & (log_FSC_H)./log_FSC<median(log_FSC_H./log_FSC)+.1;
Both_Gates=(Gate+Gate2==2);
mCherry1=mCherry2(Both_Gates)./FSC_A(Both_Gates);%correct for size
%mCherry1=mCherry(Both_Gates);
mCherry=[mCherry;mCherry1];
end
histogram(mCherry,edges,'normalization','probability','facecolor',face_color(4,:))

xlim([0,1000])
ylim([0,.3])
set(gca,'XScale','log')
end
%saveas(gcf,'LogPlots/CV_histograms_CRexamples.png')

%% 

function p_value_corrected=EnrichmentPlot(group,condition,type)
global number_of_clusters face_color face_color_CV
%p_GO_cluster_corrected=EnrichmentPlot(Complex_matrix,cdis_reordered_noVP16,'GO');
if strcmp(type,'Cluster')==1
    number_iterations=number_of_clusters;
    color=face_color;
elseif strcmp(type,'ClusterCV')==1
    number_iterations=number_of_clusters;
    color=face_color_CV;
elseif ('GO_cluster')
    number_iterations=length(group(:,1));
    color=face_color;
    number_plots=number_of_clusters;
else
    number_iterations=length(unique(group));
end

p_value=NaN(number_iterations,1);
for i=1:number_iterations
    temp=crosstab(group==i,condition);
    [h,p_value(i),stats]=fishertest(temp,'tail','both');
    enrichment(i)=temp(2,2)./sum(temp(2,:))./(sum(temp(:,2))./sum(temp,'all'));
end
    [p_value_corrected,c_alpha,h]=fwer_holmbonf(p_value,0.05); %Holm-Bonferroni correction
    sig=p_value_corrected<=0.05;
    x=1:length(group);
 figure
 hold on
 for i=1:number_iterations
    bar(i,enrichment(i),'FaceColor',color(i,:));
    %bar(i,p_MI_cluster_transformed(i),'FaceColor',face_color(i,:));
 end
    scatter(x(sig),enrichment(sig)+.1,'*k')
    ylim([0,3])
end

