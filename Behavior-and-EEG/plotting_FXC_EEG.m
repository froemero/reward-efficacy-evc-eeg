% Sanity check FXCEEG data

%%

%% SPECIFY YOUR WORKING DIRECTORY! %%
PATH = './';
% load behavioral Data
load(sprintf('%sBehavior/FXCallSubDataTable.mat', PATH))
% load EEG Data

load(sprintf('%sExport/CDAT.mat', PATH))

%%

load(sprintf('%sExport/chanlocs.mat', PATH));
nchans = 35;%size(LDAT,1);
for n = 1:nchans
    
    expression = [chanlocs(n).labels '=' sprintf('%d',n) ';'];
    eval(expression)
end




%% 2) get unique subject names for later processing and preselection of partially existing data
[~,idx] = unique(FXCallSubDataTable.SubID,'first');
vps = FXCallSubDataTable.SubID(sort(idx));
outSubs = [1, 21, 35, 38, 39, 43, 46, 53, 63];

dropoutVP = ismember(FXCallSubDataTable.SubID, outSubs); 

%% load & define stuff for plotting


PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;

%[ALLEEG EEG CURRENTSET ALLCOM] = eeglab; % for topoplot function


%% Grand average all conditions

cTIME    = linspace(-200, 1500, size(CDAT,2));
XMIN = -200;
XMAX =  1500;

CERP=nanmean(CDAT,3);


%% prepare filters
% 1) FXC by REWARD for c
E1R0 = find(~dropoutVP & FXCallSubDataTable.EffLvl==1 & FXCallSubDataTable.RewLvl==1 & FXCallSubDataTable.IsMiss~=1 & FXCallSubDataTable.RT>200);
E0R0 = find(~dropoutVP & FXCallSubDataTable.EffLvl==0 & FXCallSubDataTable.RewLvl==1 & FXCallSubDataTable.IsMiss~=1 & FXCallSubDataTable.RT>200);
E1R1 = find(~dropoutVP & FXCallSubDataTable.EffLvl==1 & FXCallSubDataTable.RewLvl==2 & FXCallSubDataTable.IsMiss~=1 & FXCallSubDataTable.RT>200);
E0R1 = find(~dropoutVP & FXCallSubDataTable.EffLvl==0 & FXCallSubDataTable.RewLvl==2 & FXCallSubDataTable.IsMiss~=1 & FXCallSubDataTable.RT>200);

%% 1) Cue
YMIN=-6;
YMAX=2;

XMIN = -200;
XMAX =  1500;
figure();
hold on;
title('Cz Cue Efficacy by reward', 'fontsize', 14)
 f = fill([1000 1000 1500 1500],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
 f.EdgeColor = 'none';
 f.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(cTIME,nanmean(CDAT(Cz,:,E1R1),3),'-','linewidth',2, 'color', [247/255, 147/255, 30/255]);
Leg2=plot(cTIME,nanmean(CDAT(Cz,:,E1R0),3),'-','linewidth',2, 'color', [0/255, 114/255, 178/255]); % rgb(83.5%,36.9%,0%)
Leg3=plot(cTIME,nanmean(CDAT(Cz,:,E0R1),3),'-','linewidth',2, 'color', [252/255, 238/255, 33/255]); %rgb(80%,47.5%,65.5%)
Leg4=plot(cTIME,nanmean(CDAT(Cz,:,E0R0),3),'-','linewidth',2, 'color', [41/255, 171/255, 226/255]); %rgb(80%,47.5%,65.5%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('Time [ms]')
ylabel('Amplitude [�V]')
set(gca, 'fontsize', 12);
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
t=gridxy(0, 'Linestyle',':','LineWidth',1);
t2=gridxy(1500, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2, Leg3, Leg4],'High efficacy, High reward','High efficacy, Low reward',  'Low efficacy, High reward', 'Low efficacy, Low reward');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 12, 'position',[.30,.20,lp(3:4)])%
set(lh,'box','off')
set(gcf,'units','centimeters','position',[0 0 12 10])
%%
saveas(gcf, 'Figures/CNV_time_course.pdf')
%% P3b
YMIN=-6;
YMAX=8;

XMIN = -200;
XMAX =  1500;
figure();
hold on;
title('Pz Cue Efficacy by reward', 'fontsize', 14)
 f = fill([250 250 550 550],[YMIN YMAX YMAX YMIN], [.8 .8 .8]);
 f.EdgeColor = 'none';
 f.EdgeColor = 'none';
%  f2 = fill([700 700 1000 1000],[YMIN YMAX YMAX YMIN], [.7 .7 .7]);
%  f2.EdgeColor = 'none';
%  f2.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
Leg1=plot(cTIME,nanmean(CDAT(Pz,:,E1R1),3),'-','linewidth',2, 'color', [247/255, 147/255, 30/255]);
Leg2=plot(cTIME,nanmean(CDAT(Pz,:,E1R0),3),'-','linewidth',2, 'color', [0/255, 114/255, 178/255]); % rgb(83.5%,36.9%,0%)
Leg3=plot(cTIME,nanmean(CDAT(Pz,:,E0R1),3),'-','linewidth',2, 'color', [252/255, 238/255, 33/255]); %rgb(80%,47.5%,65.5%)
Leg4=plot(cTIME,nanmean(CDAT(Pz,:,E0R0),3),'-','linewidth',2, 'color', [41/255, 171/255, 226/255]); %rgb(80%,47.5%,65.5%)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabel('Time [ms]')
ylabel('Amplitude [�V]')
set(gca, 'fontsize', 12);
set(gca,'TickDir','out')
xlabh = get(gca,'XLabel');
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
t=gridxy(0, 'Linestyle',':','LineWidth',1);
t2=gridxy(1500, 'Linestyle',':','LineWidth',1);
h = legend([Leg1, Leg2, Leg3, Leg4],'high efficacy high reward','high efficacy low reward',  'low efficacy high reward', 'low efficacy low reward');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,  'fontsize', 12)%, 'position',[.65,.85,lp(3:4)]
set(lh,'box','off')
set(gcf,'units','centimeters','position',[0 0 12 10])
saveas(gcf, 'Figures/Cue_P3b_time_course.pdf')





addpath(genpath('EEGfunctions/'));



%% plot LMM topos


%% this is CNV
load(sprintf('%sExport/plotes2.mat', PATH), 'Emat2') % estimates

load(sprintf('%sExport/plotts2.mat', PATH), 'Tmat2') % t-values


%% 

PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-1.5 1.5];
figure();
subplot(1,2,1)
title('Reward Effect')
topoplot(Emat2.Reward,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Fz, FCz, Cz],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
 
 export_fig('Figures/CNV_reward_topo','-pdf','-painters');
%%
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-1.5 1.5];
figure();
subplot(1,2,1)
title('Efficacy Effect')
topoplot(Emat2.Efficacy,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Fz,FCz, Cz],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
 
 
 export_fig('Figures/CNV_efficacy_topo','-pdf','-painters');

 %% Cue P3b
load(sprintf('%sExport/plotes5.mat', PATH), 'Emat5') % estimates

load(sprintf('%sExport/plotts5.mat', PATH), 'Tmat5') % t-values


%% 

PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-1.5 1.5];
figure();
subplot(1,2,1)
title('Reward Effect')
topoplot(Emat5.Reward,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz, P4, P3],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
 
 export_fig('Figures/P3b_reward_topo','-pdf','-painters');

%%
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-1.5 1.5];
figure();
subplot(1,2,1)
title('Efficacy Effect')
topoplot(Emat5.Efficacy,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[Pz, P4, P3],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
 
 export_fig('Figures/P3b_efficacy_topo','-pdf','-painters');
 