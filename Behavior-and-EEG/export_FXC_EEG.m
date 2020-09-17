% Export time window average or time series to access in R for analyses  

%% SPECIFY YOUR WORKING DIRECTORY! %%
PATH = './'; %edit to your needs
EXPORTPATH ='Export/'; % Export folder (segmented data will be saved there)
srate=512; % sampling rate
bl=200; %length of prestimulus interval in ms

%% Cue

load(strcat(PATH,EXPORTPATH,'CDAT.mat'), 'CDAT'); % load whatever your data is
%%
BFS = CDAT; % set BFS (big fucking structure) to whatever your data is
nchans = size(BFS,1);

%% export mean CNV --> 1000 - 1500
mint=1000; %export start in ms
maxt=1500; %export end in ms
ERPNAME='CNV'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates

%%
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% export mean Cue P3b
mint=250; %export start in ms
maxt=550; %export end in ms
ERPNAME='P3b'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates

%%
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');



%% Feedback

load(strcat(PATH,EXPORTPATH,'FDAT.mat'), 'FDAT'); % load whatever your data is
%%

BFS = FDAT; % set BFS to whatever your data is
nchans = size(BFS,1);

%% peak to peak FRN
mint = 200;
maxt = 400;
round((mint+bl)*srate/1000)
round((maxt+bl)*srate/1000)
round((150)*srate/1000)
P2PFCZ = p2p(FDAT, 205, 307, 10, 77);

save(sprintf('%s%sP2PFCZ.mat', PATH,EXPORTPATH), 'P2PFCZ', '-v6');


%% Response

load(strcat(PATH,EXPORTPATH,'RDAT.mat'), 'RDAT'); % load whatever your data is
%%

BFS = RDAT; % set BFS to whatever your data is
nchans = size(BFS,1);

%% ERN 
mint=0; %export start in ms
maxt=100; %export end in ms
ERPNAME='ERN'; % name of the component you want to export

% the rest is automatic
EXPERP=nanmean(BFS(:,round((mint+bl)*srate/1000):round((maxt+bl)*srate/1000),:),2); % adapt for multiple sampling rates

%%
EXPERP=reshape(EXPERP, [size(BFS,1), size(BFS,3)]);
EXPERP=EXPERP';
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% export BASE
load(strcat(PATH,EXPORTPATH,'CBASE.mat'), 'CBASE')
ERPNAME='Baseline';
EXPERP = CBASE';
mint=200; %export start in ms
maxt=0;
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');

%% export just BASE as is
save(sprintf('%s%sBASE.mat', PATH, EXPORTPATH), 'BASE', '-v6');

%% export response BASE
load(strcat(PATH,EXPORTPATH,'RBASE.mat'), 'RBASE')
ERPNAME='rBaseline';
EXPERP = RBASE';
mint=200; %export start in ms
maxt=0;
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');


%% export Feedback BASE
load(strcat(PATH,EXPORTPATH,'FBASE.mat'), 'FBASE')
ERPNAME='fBaseline';
EXPERP = FBASE';
mint=200; %export start in ms
maxt=0;
expression= [ERPNAME,sprintf('%d%d',mint,maxt) '=EXPERP;']; % have one name for all and use eval to create multiple outputs depending on what is set
eval(expression);
save(sprintf('%s%s%s%d%d.mat', PATH, EXPORTPATH,ERPNAME,mint,maxt), sprintf('%s%d%d',ERPNAME,mint,maxt), '-v6');
