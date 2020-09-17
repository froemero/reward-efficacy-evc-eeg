%% Making BFS for FXC EEG

% 1) Cue-locked -200: 1500 ms
%                 xmin: -0.1992
%                 xmax: 1.4980
%                 DP = 870
% 2) Stim-locked: -200:800 ms?
%                 xmin: -0.1992
%                 xmax: 0.7988
%                 DP: 512

% 3) R-locked: -200: 800 ms
% 4) FB locked: -200: 800 ms
% 5) Cue-locked but whole Interval 
%% Info for segmentation analyses:

%% SPECIFY YOUR WORKING DIRECTORY! %%
PATH = './'; 
addpath('EEGfunctions/'); 

%% about the data:
%srate = ugly 512;
% divide by 1.9531 

% **Cue**
% 91: efficacy0 reward1
% 92: efficacy0 reward2
% 93: efficacy1 reward1
% 94: efficacy1 reward2
%
% **stroop stimuli**
% 51: incongruent
% 52: congruent
% 53: neutral
%
% **response**
% 200: wrong
% 210: correct
%
% **feedback**
% 10: $0
% 20: $0.10
% 110: $1.00

alltrigs= {[91:94]; [51:53];[200, 210]; [10, 20, 110]};%
allIL= [-0.2, 1.5; -0.2, 0.8; -0.2, 0.8; -0.2, 0.8];%

%% for response-locked segementation:
load(sprintf('%sBehavior/FXCallSubDataTable.mat', PATH))

%% approach: Get all condition triggers
% just NAN the artifact trials and paste all subject files together.
%% run only when needed, adds EEGLab

% addpath('N:/Software/eeglab13_5_4b')
%% preparation
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%%
SOURCEFILES = dir(strcat(PATH, 'EEG-continuous/*.set')); %all cleaned files
SUBJECTS = 1:numel(SOURCEFILES);
%%
allsegs= { 'c';'s';'r'; 'f'}; % change back! 
contmat=zeros(50,3, 4); % check number of trials found per participant (before & after AR)
%% loop through different segmentations to save memory
for nsegs = 1:length(allsegs)
    %%
    if strcmp(allsegs{nsegs}, 'c')
        DAT=nan(35,870,height(FXCallSubDataTable)); % preallocate memory and preset datamatrix (larger for cue than others)
        BASE = nan(35,height(FXCallSubDataTable));
    elseif strcmp(allsegs{nsegs}, 'r')
        
        DAT=nan(35,512,height(FXCallSubDataTable(~isnan(FXCallSubDataTable.Resp),:))); % preallocate memory and preset datamatrix (larger for cue than others)
        BASE = nan(35,height(FXCallSubDataTable(~isnan(FXCallSubDataTable.Resp),:)));
    elseif strcmp(allsegs{nsegs}, 'l')
        DAT=nan(35,1382,height(FXCallSubDataTable)); %     
    else
        DAT=nan(35,512,height(FXCallSubDataTable));
        BASE = nan(35,height(FXCallSubDataTable));
    end
    
    % get relevant markers
    Trigs =alltrigs{nsegs};
    Mks = cellstr(num2str(Trigs(:)));
    if strcmp(allsegs{nsegs}, 'r')
        [subIDs,uniqueRows] = unique(FXCallSubDataTable.SubID(~isnan(FXCallSubDataTable.Resp))); % get first line for each sub to start putting data into
        
    else
        [subIDs,uniqueRows] = unique(FXCallSubDataTable.SubID); % get first line for each sub to start putting data into
    end
    %% go through all subjects
    for s =1:numel(SOURCEFILES)
        %%
        vpn=SOURCEFILES(s).name(1:3); % read out of dataset name
        s_id = str2double(SOURCEFILES(s).name(1:3)); % reads in subject ID and converts it to double
        fprintf('processing participant number %d\n', s_id)
        contmat(s,1,nsegs)=s_id;
        pcnt=uniqueRows(s); % get first line for each sub to start putting data into
        lns = sum(FXCallSubDataTable.SubID==s_id); % how many trials
        %%
        if lns>0
           
            
            % load dataset
            EEG = pop_loadset('filename',sprintf('%s_ICAPruned_Reref.set',vpn),'filepath',sprintf('%sEEG-continuous',PATH));
            
             % load fixed trigger file for r-locked trials
            if strcmp(allsegs{nsegs}, 'r')
                triggersFixed =readtable(sprintf('%s/TriggersCleaned/%s_triggersFixed.csv',PATH, vpn));
                % change triggercodes to fixed trigger codes
                for eventnum = 1:length(EEG.event)
                    
                    EEG.event(eventnum).type = triggersFixed.eventCodeNew(eventnum);
                    
                end
                clear triggersFixed;
                
            end
            %% get event-locked data
            fprintf('%s-locked segmentation of participant number %d\n', allsegs{nsegs}, s_id)
            IL =allIL(nsegs,:); 

            [EEG, indices]=pop_epoch(EEG, Mks,IL);
            
            if strcmp(allsegs{nsegs}, 'r')
                ndif = length(FXCallSubDataTable.RT( FXCallSubDataTable.SubID==s_id & ~isnan(FXCallSubDataTable.Resp))) - size(EEG.data,3);
            else
                ndif = 600- size(EEG.data,3);
            end
            
            % export baseline
            if (strcmp(allsegs{nsegs}, 'c')|| strcmp(allsegs{nsegs}, 'f')|| strcmp(allsegs{nsegs}, 's'))
                BASE(:,pcnt + ndif:pcnt+ndif+size(EEG.data,3)-1)=nanmean(EEG.data(:,1:103,1:size(EEG.data,3)),2);
                
            end
            %% export baseline and baseline correction
            if (strcmp(allsegs{nsegs}, 'r')|| strcmp(allsegs{nsegs}, 's'))% used to be c
                if strcmp(allsegs{nsegs}, 'r')
                    BASE(:,pcnt + ndif:pcnt+ndif+size(EEG.data,3)-1)=nanmean(EEG.data(:,1:103,1:size(EEG.data,3)),2);
                end
                EEG = pop_rmbase( EEG, [-200 0]);
            elseif strcmp(allsegs{nsegs}, 'l')
                      EEG = pop_rmbase( EEG, [-200 0]);
            else
                EEG = pop_rmbase( EEG, [-200 0]);
            end
            
            DP = size(EEG.data,2);
            
            % select all trials based on object onset markers
            % no further segmentation required --> matching with logfiles
            contmat(s,2,nsegs)=EEG.trials;
            %%
            % Artifact detection (artifact trials are naned out later)
            % (Irej: index of rejected trials)
            %% I don't do the artifact rejection on the EOG channels
            [EEGn, Irej] = pop_eegthresh(EEG,1,[1:33] ,-150,150, IL(1),IL(2)-0.002,1,0); 
            
            %% %   >> [rej rejE] = rejtrend( signal, winsize, maxslope, minR, step);
            
            [rej, rejE] = rejtrend(EEG.data(1:33, :,:),DP,50,0.3,1); 
            
            frej= find(rej==1);
            Out=unique(sort([Irej,frej]));
            if isempty(Out)
                remout=0;
            else
                remout = length(Out);
            end
            fprintf('%d trials removed alltogether\n', remout);
            
            %% add subject data to data matrix after naning out artifact segments
            EEG.data(:,:,Out) = nan; % invalidate data with artifacts
            DAT(:,:,pcnt + ndif:pcnt+ndif+size(EEG.data,3)-1)=EEG.data(:,:,1:size(EEG.data,3));
            contmat(s,3,nsegs)=(EEG.trials-size(Out,2));
            
        else
            fprintf('Subject not in table! Run analyzeFXC first!\n')
        end
    end
    
    
    fprintf('%d participants in dataset... saving %s locked data\n', s, allsegs{nsegs})
    if strcmp(allsegs{nsegs}, 'c')
        CDAT = DAT;
        save(sprintf('%sExport/CDAT.mat', PATH), 'CDAT', '-v7.3');
        CBASE = BASE;
        save(sprintf('%sExport/CBASE.mat', PATH), 'CBASE', '-v7.3');
        clear BASE;
        clear CBASE;
        clear CDAT;
    elseif strcmp(allsegs{nsegs}, 's')
        SDAT = DAT;
        SBASE = BASE;
        save(sprintf('%sExport/SDAT.mat', PATH), 'SDAT', '-v7.3');
        save(sprintf('%sExport/SBASE.mat', PATH), 'SBASE', '-v7.3');
        clear SDAT;
        clear BASE;
        clear SBASE;
    elseif strcmp(allsegs{nsegs}, 'r')
        RDAT=DAT;
        save(sprintf('%sExport/RDAT.mat', PATH), 'RDAT', '-v7.3');
        RBASE = BASE;
        save(sprintf('%sExport/RBASE.mat', PATH), 'RBASE', '-v7.3');
        clear RBASE;
        clear BASE;
        clear RDAT;
    elseif strcmp(allsegs{nsegs}, 'l')
        LDAT= DAT;
        save(sprintf('%sExport/LDAT.mat', PATH), 'LDAT', '-v7.3');
    else
        FDAT=DAT;
        save(sprintf('%sExport/FDAT.mat', PATH), 'FDAT', '-v7.3');
        FBASE = BASE;
        save(sprintf('%sExport/FBASE.mat', PATH), 'FBASE', '-v7.3');
        clear FBASE;
        clear BASE;
        clear FDAT;
    end
    
    
end

chanlocs = EEG.chanlocs;
save(sprintf('%sExport/chanlocs.mat', PATH), 'chanlocs', '-v7.3');

save(sprintf('%sExport/contmat.mat', PATH), 'contmat', '-v7.3');